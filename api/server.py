import json
import shutil
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import Depends, FastAPI, File, HTTPException, UploadFile, status
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.security import OAuth2PasswordRequestForm
from pydantic import BaseModel
from sqlalchemy.orm import Session

from . import config
from .orchestrator import Engine
from .db import get_db, get_user, create_user, increment_question_count, delete_user, init_db
from .auth import hash_password, verify_password, create_token, current_user, User
from .history import save_history, load_history, list_sessions, delete_session, rename_session

# MCP + LLM engine — created once at startup, shared across all requests
engine: Optional[Engine] = None


# ----------------------------
# App startup / shutdown
# ----------------------------

@asynccontextmanager
async def lifespan(app: FastAPI):
    init_db()
    global engine
    engine = Engine()
    await engine.connect()
    print(f"[api] ready — connected servers: {engine.connected_servers or '(none)'}")
    yield
    await engine.close()


app = FastAPI(title="Drug Discovery Platform API", lifespan=lifespan)


# ----------------------------
# Per-user directory helpers
# ----------------------------

def _user_dir(username: str) -> Path:
    """Create the user's root folder with all subdirs. Called once at register/login."""
    base = config.USER_DATA_DIR / username
    for sub in ("history", "files", "job_tmp"):
        (base / sub).mkdir(parents=True, exist_ok=True)
    return base

def _seed_case_studies(username: str) -> None:
    """Copy the bundled case studies into a new user's history/ so they see examples.
    Titled via _meta.json (not the message content) so the sidebar shows the
    filename regardless of what the first message in each case study says.

    The source JSONs store tool file references as relative 'files/...' paths
    (the same convention live sessions use before path rewriting), which only
    resolve correctly under the seeded-to user's own files/ dir — so the
    relative prefix is rewritten to that user's absolute path at copy time
    rather than baked into the shared template."""
    case_studies_dir = config.CODE_DIR / "case_studies"
    history_dir       = config.USER_DATA_DIR / username / "history"
    files_dir         = config.USER_DATA_DIR / username / "files"
    user_files_prefix = str(files_dir)
    for src in case_studies_dir.glob("*.json"):
        text = src.read_text().replace('"files/', f'"{user_files_prefix}/')
        (history_dir / src.name).write_text(text)
        rename_session(username, src.stem, src.stem)

    # The actual assets those 'files/...' references point to
    for src in (case_studies_dir / "files").glob("*"):
        if src.is_file():
            shutil.copy2(src, files_dir / src.name)

# Pure path helpers — the dirs already exist (created at register/login)
def _files_dir(username: str) -> Path:
    return config.USER_DATA_DIR / username / "files"

def _job_tmp_dir(username: str) -> Path:
    return config.USER_DATA_DIR / username / "job_tmp"

def _job_ready_file(username: str) -> Path:
    return _job_tmp_dir(username) / "job_ready.json"

def _safe_file(username: str, file_path: str) -> Path:
    """Resolve a path inside the user's files/ dir, rejecting traversal escapes."""
    base   = _files_dir(username).resolve()
    target = (base / file_path).resolve()
    if target != base and base not in target.parents:
        raise HTTPException(status_code=403, detail="Access denied")
    if not target.exists():
        raise HTTPException(status_code=404, detail="File not found")
    return target

def _load_system_prompt(username: str) -> str:
    try:
        text = config.SYSTEM_FILE.read_text()
    except Exception:
        text = ""
    files_dir = str(config.USER_DATA_DIR / username / "files")
    return (
        text
        + f"\n\nCurrent user: {username}"
        + f"\nFile dir: {files_dir}"
        + f"\nAlways pass '{files_dir}' as file_dir when calling any tool that writes files."
    )


# ----------------------------
# Health
# ----------------------------

@app.get("/health")
async def health() -> Dict[str, Any]:
    return {
        "status"           : "ok",
        "connected_servers": engine.connected_servers if engine else [],
        "default_model"    : config.DEFAULT_MODEL,
    }


# ----------------------------
# Auth
# ----------------------------

class RegisterRequest(BaseModel):
    username       : str
    password       : str
    free_questions : int = config.FREE_QUESTIONS_DEFAULT

class TokenResponse(BaseModel):
    access_token : str
    token_type   : str = "bearer"

class UserInfo(BaseModel):
    username       : str
    free_questions : int
    n_questions    : int
    is_active      : bool


@app.post("/auth/register", status_code=status.HTTP_201_CREATED)
async def register(req: RegisterRequest, db: Session = Depends(get_db)) -> Dict[str, str]:
    if get_user(db, req.username):
        raise HTTPException(status_code=400, detail="Username already taken")
    create_user(db, req.username, hash_password(req.password), req.free_questions, req.password)
    _user_dir(req.username)   # create history/ files/ job_tmp/ on the filesystem
    _seed_case_studies(req.username)   # give new users example sessions to look at
    return {"message": f"User '{req.username}' created"}


@app.post("/auth/login", response_model=TokenResponse)
async def login(
    form: OAuth2PasswordRequestForm = Depends(),
    db  : Session = Depends(get_db),
) -> TokenResponse:
    user = get_user(db, form.username)
    if not user or not verify_password(form.password, user.hashed_password):
        raise HTTPException(status_code=401, detail="Invalid credentials")
    _user_dir(user.username)   # ensure dirs exist for existing accounts
    return TokenResponse(access_token=create_token(user.username))


@app.get("/auth/me", response_model=UserInfo)
async def me(user: User = Depends(current_user)) -> UserInfo:
    return UserInfo(
        username       = user.username,
        free_questions = user.free_questions,
        n_questions    = user.n_questions,
        is_active      = user.is_active,
    )


class KeysRequest(BaseModel):
    openai_key     : Optional[str] = None
    openrouter_key : Optional[str] = None
    ollama_key     : Optional[str] = None

class KeysResponse(BaseModel):
    openai_key     : str
    openrouter_key : str
    ollama_key     : str

@app.get("/auth/keys", response_model=KeysResponse)
async def get_keys(user: User = Depends(current_user)) -> KeysResponse:
    def mask(k: str) -> str:
        return k[:8] + "…" if len(k) > 8 else k
    return KeysResponse(
        openai_key     = mask(user.openai_key     or ""),
        openrouter_key = mask(user.openrouter_key or ""),
        ollama_key     = mask(user.ollama_key     or ""),
    )

@app.patch("/auth/keys")
async def update_keys(
    req  : KeysRequest,
    user : User    = Depends(current_user),
    db   : Session = Depends(get_db),
) -> Dict[str, str]:
    if req.openai_key     is not None: user.openai_key     = req.openai_key
    if req.openrouter_key is not None: user.openrouter_key = req.openrouter_key
    if req.ollama_key     is not None: user.ollama_key     = req.ollama_key
    db.commit()
    return {"message": "Keys updated"}


@app.delete("/auth/account")
async def delete_account(
    user : User    = Depends(current_user),
    db   : Session = Depends(get_db),
) -> Dict[str, str]:
    """Delete the account (frees the username) but keep files — the user's folder
    is renamed to <username>_deleted instead of removed."""
    src = config.USER_DATA_DIR / user.username
    if src.exists():
        dst = config.USER_DATA_DIR / f"{user.username}_deleted"
        i   = 1
        while dst.exists():                 # avoid clobbering an earlier deletion
            dst = config.USER_DATA_DIR / f"{user.username}_deleted_{i}"
            i  += 1
        src.rename(dst)
    delete_user(db, user)
    return {"message": "Account deleted"}


# ----------------------------
# Models
# ----------------------------

@app.get("/models")
async def list_models() -> Dict[str, Any]:
    try:
        lines  = (config.CONFIG_DIR / "llm_models").read_text().splitlines()
        models = [l.strip() for l in lines if l.strip() and "(err)" not in l]
    except Exception:
        models = [config.DEFAULT_MODEL]
    return {"models": models, "default": config.DEFAULT_MODEL}


# ----------------------------
# Notifications — job completion alerts
# ----------------------------

@app.get("/notifications")
async def get_notifications(user: User = Depends(current_user)) -> Dict[str, Any]:
    """Return pending job-ready notifications, then clear them."""
    path = _job_ready_file(user.username)
    if not path.exists():
        return {"notifications": []}
    try:
        items = json.loads(path.read_text())
    except Exception:
        items = []
    path.unlink(missing_ok=True)
    return {"notifications": items}


# ----------------------------
# User files — list, upload, download, delete
# ----------------------------

@app.get("/files/list")
async def list_files(user: User = Depends(current_user)) -> Dict[str, Any]:
    """Recursively list all files in the user's files/ dir with metadata."""
    base  = _files_dir(user.username)
    files = []
    for f in sorted(base.rglob("*"), key=lambda p: p.stat().st_mtime, reverse=True):
        if f.is_file():
            stat = f.stat()
            files.append({
                "path"    : str(f.relative_to(base)),
                "name"    : f.name,
                "size"    : stat.st_size,
                "modified": stat.st_mtime,
            })
    return {"files": files}


@app.post("/files/upload", status_code=status.HTTP_201_CREATED)
async def upload_file(
    file : UploadFile = File(...),
    user : User       = Depends(current_user),
) -> Dict[str, str]:
    """Upload a file into the user's files/ dir."""
    dest = _files_dir(user.username) / (file.filename or "upload")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(await file.read())
    return {"path": str(dest.relative_to(_files_dir(user.username)))}


@app.delete("/files/{file_path:path}")
async def delete_file(
    file_path : str,
    user      : User = Depends(current_user),
) -> Dict[str, str]:
    _safe_file(user.username, file_path).unlink()
    return {"message": "Deleted"}


@app.get("/files/{file_path:path}")
async def serve_file(
    file_path : str,
    user      : User = Depends(current_user),
) -> FileResponse:
    """Serve a file from the current user's files/ directory."""
    return FileResponse(_safe_file(user.username, file_path))


# ----------------------------
# Sessions
# ----------------------------

@app.get("/sessions")
async def get_sessions(user: User = Depends(current_user)) -> Dict[str, Any]:
    """List all conversations for the current user, newest first."""
    sessions = list_sessions(user.username)
    for s in sessions:
        s["running"] = engine.turns.is_running(s["id"])
    return {"sessions": sessions}


# ----------------------------
# History
# ----------------------------

@app.get("/history/{session_id}")
async def get_history(
    session_id : str,
    user       : User = Depends(current_user),
) -> Dict[str, Any]:
    """Return all messages for a session — only if it belongs to the current user."""
    msgs = load_history(user.username, session_id)
    if not msgs:
        raise HTTPException(status_code=404, detail="Session not found")
    # Strip the system prompt before sending to the client
    visible = [m for m in msgs if m.get("role") not in ("system", "feedback")]
    return {"session_id": session_id, "messages": visible}


class RenameRequest(BaseModel):
    title : str

@app.patch("/history/{session_id}")
async def rename_history(
    session_id : str,
    req        : RenameRequest,
    user       : User = Depends(current_user),
) -> Dict[str, str]:
    if not rename_session(user.username, session_id, req.title):
        raise HTTPException(status_code=404, detail="Session not found")
    return {"message": "Renamed"}

@app.delete("/history/{session_id}")
async def remove_history(
    session_id : str,
    user       : User = Depends(current_user),
) -> Dict[str, str]:
    if not delete_session(user.username, session_id):
        raise HTTPException(status_code=404, detail="Session not found")
    return {"message": "Deleted"}


# ----------------------------
# Feedback
# ----------------------------

FEEDBACK_OPTIONS = {
    "I did not get what I wanted",
    "The agent job was good but needed a lot of instructions",
    "The agent job was good with reasonable instructions",
    "The agent worked better than I expected",
}

class FeedbackRequest(BaseModel):
    rating : str

@app.post("/feedback/{session_id}", status_code=status.HTTP_201_CREATED)
async def submit_feedback(
    session_id : str,
    req        : FeedbackRequest,
    user       : User = Depends(current_user),
) -> Dict[str, str]:
    if req.rating not in FEEDBACK_OPTIONS:
        raise HTTPException(status_code=400, detail="Invalid rating value")
    msgs = load_history(user.username, session_id)
    if not msgs:
        raise HTTPException(status_code=404, detail="Session not found")
    msgs.append({
        "role"     : "feedback",
        "rating"   : req.rating,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    })
    save_history(user.username, session_id, msgs)
    return {"message": "Feedback recorded"}


# ----------------------------
# Chat — SSE, decoupled from the HTTP request (see turns.py / orchestrator.py):
# starting a turn and watching it are two different things. /chat starts (or
# 409s if one's already running) and streams; /chat/stream/{id} only attaches,
# for reconnecting without submitting a new message; the request disconnecting
# never stops the turn itself, only that one watcher.
# ----------------------------

class ChatMessage(BaseModel):
    role    : str
    content : Optional[str] = None

class ChatRequest(BaseModel):
    messages   : List[ChatMessage]
    model      : Optional[str] = None
    session_id : str

@app.post("/chat")
async def chat(
    req  : ChatRequest,
    user : User    = Depends(current_user),
    db   : Session = Depends(get_db),
) -> StreamingResponse:
    """Start a chat turn and stream it as SSE. Each frame: data: {json} with
    type token/tool/done/error. 409 if a turn is already running for this
    session_id — reattach via GET /chat/stream/{session_id} instead."""

    # Quota check — users with their own API key bypass the free limit.
    # Only charged for a genuinely new turn — never on a race where this
    # fires again for a session that's already running (see 409 below).
    has_quota   = True
    has_own_key = bool(user.openai_key or user.openrouter_key or user.ollama_key)
    if engine.turns.get(req.session_id) is None:
        has_quota = increment_question_count(db, user)
    if not has_quota and not has_own_key:
        raise HTTPException(status_code=402,
                            detail="Free question quota reached. Please supply your own API key.")

    model_full = req.model or config.DEFAULT_MODEL

    # Build message list: system prompt + history from file + latest user message
    system_msg = {"role": "system", "content": _load_system_prompt(user.username)}
    history    = [m for m in load_history(user.username, req.session_id)
                  if m.get("role") not in ("system", "feedback")]
    incoming   = [m.model_dump(exclude_none=True) for m in req.messages]
    messages : List[Dict[str, Any]] = (
        [system_msg]
        + (history if history else incoming[:-1])
        + [incoming[-1]]
    )

    turn, is_new = engine.start_turn(user.username, req.session_id, messages, model_full)
    if not is_new:
        raise HTTPException(status_code=409, detail="A turn is already running for this session.")

    async def event_stream():
        async for event in turn.subscribe():
            yield f"data: {json.dumps(event)}\n\n"

    return StreamingResponse(event_stream(), media_type="text/event-stream")


@app.get("/chat/stream/{session_id}")
async def chat_stream(session_id: str, user: User = Depends(current_user)) -> StreamingResponse:
    """Reattach to an already-running turn — after a session switch,
    logout/login, dropped connection, or another tab — without submitting a
    new message. 404 means nothing is currently running for this session
    (either it finished already, or it was never started)."""
    turn = engine.turns.get(session_id)
    if turn is None or turn.username != user.username:
        raise HTTPException(status_code=404, detail="No turn is currently running for this session.")

    async def event_stream():
        async for event in turn.subscribe():
            yield f"data: {json.dumps(event)}\n\n"

    return StreamingResponse(event_stream(), media_type="text/event-stream")


@app.post("/chat/{session_id}/cancel")
async def cancel_chat(session_id: str, user: User = Depends(current_user)) -> Dict[str, bool]:
    turn = engine.turns.get(session_id)
    if turn is not None and turn.username != user.username:
        raise HTTPException(status_code=404, detail="No turn is currently running for this session.")
    return {"ok": engine.turns.cancel(session_id)}
