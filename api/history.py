import json
from pathlib import Path
from typing import Any, Dict, List

from . import config


# ----------------------------
# Path helpers (dirs are created at register/login, not here)
# ----------------------------

def _history_dir(username: str) -> Path:
    return config.USER_DATA_DIR / username / "history"

def _history_file(username: str, session_id: str) -> Path:
    return _history_dir(username) / f"{session_id}.json"

def _meta_file(username: str) -> Path:
    return _history_dir(username) / "_meta.json"


# ----------------------------
# Per-session metadata (custom title + model), keyed by session_id
# ----------------------------

def _load_meta(username: str) -> Dict[str, Dict[str, str]]:
    path = _meta_file(username)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}

def _save_meta(username: str, meta: Dict[str, Dict[str, str]]) -> None:
    path = _meta_file(username)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(meta, indent=2))

def _set_meta(username: str, session_id: str, **fields: str) -> None:
    """Merge fields into one session's metadata entry."""
    meta = _load_meta(username)
    meta.setdefault(session_id, {}).update(fields)
    _save_meta(username, meta)


# ----------------------------
# Read / write
# ----------------------------

def save_history(username: str, session_id: str, messages: List[Dict[str, Any]],
                 model: str = "") -> None:
    path = _history_file(username, session_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(messages, indent=2, default=str))
    if model:
        _set_meta(username, session_id, model=model)

def load_history(username: str, session_id: str) -> List[Dict[str, Any]]:
    path = _history_file(username, session_id)
    if not path.exists():
        return []
    try:
        return json.loads(path.read_text())
    except Exception:
        return []

def list_sessions(username: str) -> List[Dict[str, Any]]:
    """All sessions for a user, newest-first. Title/model come from _meta.json,
    falling back to the first user message for the title."""
    meta     = _load_meta(username)
    sessions = []
    for f in sorted(_history_dir(username).glob("*.json"),
                    key=lambda p: p.stat().st_mtime, reverse=True):
        if f.name == "_meta.json":
            continue
        try:
            msgs    = json.loads(f.read_text())
            default = next(
                (m["content"][:60] for m in msgs
                 if m.get("role") == "user" and m.get("content")),
                "Untitled",
            )
            info = meta.get(f.stem, {})
            sessions.append({
                "id"      : f.stem,
                "title"   : info.get("title") or default,
                "model"   : info.get("model", ""),
                "modified": f.stat().st_mtime,
            })
        except Exception:
            pass
    return sessions

def delete_session(username: str, session_id: str) -> bool:
    path = _history_file(username, session_id)
    if not path.exists():
        return False
    path.unlink()
    meta = _load_meta(username)
    if meta.pop(session_id, None) is not None:
        _save_meta(username, meta)
    return True

def rename_session(username: str, session_id: str, title: str) -> bool:
    if not _history_file(username, session_id).exists():
        return False
    _set_meta(username, session_id, title=title[:60])
    return True
