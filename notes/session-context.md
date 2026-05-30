# Session Context — BCAI (code repo)

Read this at the start of any BCAI coding session. Full workspace-level context (version diffs, broader decisions) is at `research/bcai/notes/session-context.md`.

---

## Three-layer architecture

```
Layer 1 — biochem          atomic science functions; no orchestration, no MCP
Layer 2 — pipeline.py      sync orchestration; composes biochem calls; no MCP; testable standalone
Layer 3 — server.py        thin async MCP wrapper; job registry + thread; delegates to pipeline
```

**Exception:** structure sub-agent has no `pipeline.py` because `bc.run_af3` / `bc.run_esm3`
already are complete pipelines in `biochem` — `server.py` wraps them directly.

### File tree

```
bcai/main/server.py
bcai/subagents/
  jobs.py                           shared job registry + check_status()
  screening/
    pipeline.py                     prepare_library(), run_virtual_screening()
    server.py                       MCP wrappers + list_libraries, get_top_hits
  structure/
    server.py                       MCP wrappers for bc.run_af3, bc.run_esm3
  literature/
    pipeline.py                     ← TODO
    server.py                       ← drafted, PAUSED
bcai/data/libraries/                prepared compound libraries
```

---

## Standard pattern for a new sub-agent

### pipeline.py
```python
def my_pipeline(arg1, arg2, file_dir="files", progress_cb=None):
    """Sync. All logic here. No MCP, no threads (internal threads OK for monitoring)."""
    def _p(msg):
        if progress_cb: progress_cb(msg)
    _p("step 1...")
    result = bc.some_function(arg1, file_dir=file_dir)
    _p("step 2...")
    ...
    return {...}
```

### server.py
```python
from subagents.jobs import _new_job, _update_job, check_status
from subagents.myagent import pipeline

@mcp.tool()
def my_tool(arg1, arg2, file_dir: str = "files") -> dict:
    """One-line description for the LLM."""
    file_dir = os.getenv("FILE_DIR", file_dir)
    job_id = _new_job("my_tool")
    def _run():
        try:
            result = pipeline.my_pipeline(arg1, arg2, file_dir=file_dir,
                         progress_cb=lambda msg: _update_job(job_id, progress=msg))
            _update_job(job_id, status="done", progress="done", result=result)
        except Exception as e:
            _update_job(job_id, status="failed", error=str(e))
    threading.Thread(target=_run, daemon=True).start()
    return {"job_id": job_id, "status": "running", "next_steps": ["check_status(job_id)"]}

mcp.tool()(check_status)   # standalone server use
```

---

## _ASYNC_TOOLS

In `main/server.py`:
```python
_ASYNC_TOOLS = {"run_af3", "run_esm3"}
for fn in TOOLS:
    if fn.__name__ in _ASYNC_TOOLS:
        continue
    mcp.tool()(fn)
```
When a `biochem` function is promoted to async (wrapped in a sub-agent), add its name here so the sync version is not registered.

---

## file_dir / env var pattern

In every function that writes files:
```python
def fn(..., file_dir: str = 'files'):
    file_dir = os.getenv("FILE_DIR", file_dir)
```
Priority: env var (MCP server) → explicit arg → default `'files'`.
Apply at top of each `pipeline.py` function and each `server.py` tool that accepts `file_dir`.

---

## Conda envs

| Env | Use |
|---|---|
| `bcai-dev` | Active dev — editable `biochem` install |
| `bcai` | Fallback / production snapshot |
| `biochem` | Standalone `biochem` package work |

---

## Literature sub-agent (done 2026-05-25)

`subagents/literature/pipeline.py` + `subagents/literature/server.py` — complete, wired into mother.

- `search_literature(query, source="both", max_results, file_dir)` — PubMed + Semantic Scholar, saves JSON
- `extract_information(papers_file, question, file_dir)` — general LLM extraction, saves CSV with `answer / evidence / confidence`
- LLM model from `LLM_MODEL_LITERATURE` env var; `_make_client()` is sync mirror of `client.py`'s `make_client`
- `check_status` removed from all sub-agent files — lives only in `jobs.py`, registered once in `main/server.py`

## Async job notification (done 2026-05-27, revised 2026-05-27)

When a tool returns a `job_id`, `client.py` auto-notifies the UI when the job completes — no manual `check_status` needed.

- `client.py`: after each query, scans tool messages for `job_id` → stores in `_pending_jobs` / `_seen_job_ids`. In the `while True` loop, polls `check_status` via MCP for each pending job. When done, writes `session/<user>/job_ready.json`. `_seen_job_ids` prevents re-registering completed job_ids from history.
- `BioChemAIgent_chat_UI.py`: `job_ready_changed()` checked on every rerun. When flag appears → sets `job_notification` in session state, deletes the flag file, reruns. A dismissable `st.success` banner is shown until the user clicks ✕.
- Flag file: `session/<user>/job_ready.json` — contains `job_id`, `type`, `message`.
- All new code marked with `# claude ---` / `# ---` blocks.

### Known limitations — requires replacing the Streamlit UI

Streamlit's architecture makes the following impossible or impractical to fix properly.
All three point to the same solution: replace the Streamlit frontend with a proper web stack (e.g. FastAPI + WebSocket + React/Vue, or any JS-based UI).

1. **Async job notification requires a page interaction**
   Streamlit is single-threaded; `time.sleep()` blocks all UI interaction. There is no server-push mechanism.
   The job-done banner appears on the user's next natural interaction, not the instant the job finishes.

2. **No secure user accounts — passwords are visible in plain text**
   The current `USER_FILE` CSV stores passwords unhashed and readable by anyone with file access.
   Proper auth (hashed passwords, session tokens, HTTPS) cannot be layered on top of Streamlit safely.

3. **No streaming — full response must be ready before it is shown**
   LLM responses are generated token-by-token but Streamlit can only display the complete message.
   Streaming requires a persistent connection (WebSocket or SSE) between the browser and the backend,
   which Streamlit does not expose.

---

## Docker / biochem (done 2026-05-27)

- `biochem` published to `https://github.com/imsb-uke/biochem.git` (fresh repo, clean history)
- `bcai/environment.yml`: replaced `-e ../biochem` with `git+https://github.com/imsb-uke/biochem.git`
- `biochem/protein_ligand_basics.py`: `from admet_ai import ADMETModel` moved inside `admet_predict()` (lazy import) to fix `pkg_resources` crash at Docker startup
- Docker build: `cd bcai && docker build -t bcai-image .`
- Docker run: `docker run --rm -it -p 8501:8501 -v "$PWD":/workspace bcai-image bash`

---

## Code quality backlog

- **`main/client.py` — replace `os.system()` logging with Python `logging` module**
  All `os.system(f"echo ... >> client.log")` calls → `logging.info(...)`.
  Reasons: `os.system` deprecated in favour of `subprocess`; shell injection risk; not cross-platform.
