# API backend

FastAPI backend with streaming SSE chat, MCP tool orchestration, user auth, and per-user session history.

## Layout

- `config.py` — paths, env loading, model + MCP server defaults.
- `orchestrator.py` — `Engine`: connects to MCP servers at startup, builds OpenAI tool specs, runs streaming round-1 (tool calls) + streamed final answer. Delta-based tool-call argument collection preserves nested arguments (e.g. `render_structures.style_rules`) correctly.
- `server.py` — FastAPI app: auth (JWT + SQLite), `/chat` (SSE stream), `/history`, `/sessions`, `/files`, `/notifications`.
- `history.py` — per-user JSON session history. System messages are stripped on load and excluded on save to prevent accumulation.
- `auth.py` / `db.py` — JWT auth, SQLite user table, per-user API key storage.

## Run (Docker — recommended)

```bash
cd code/
docker compose up -d          # starts api + web
# browser: http://localhost:3000
```

Code changes to `api/`, `mcp/`, `config/`, `docs/` take effect after:
```bash
docker compose restart api
```

Changes to the `biochem` package (installed at build time from GitHub) require:
```bash
docker compose build --no-cache api && docker compose up -d
```

## Key design decisions

- **Round-1 streaming:** the first LLM round (which may call tools) is streamed. Tool-call argument deltas are collected by index and concatenated. This means any text the model generates *before* a tool call streams live as a visible bubble, rather than being silently swallowed.
- **System message hygiene:** system prompt is injected fresh each request; never stored in history. Prevents unbounded context growth across sessions.
- **Tool schema:** `render_structures` uses `TypedDict` for `style_rules`/`surface_rules` so FastMCP generates a proper schema with named properties (`select`, `style`). Without this, gpt-5 collapses nested dicts to `[{}, {}]`.

## SSE event types

| type | payload | meaning |
|------|---------|---------|
| `token` | `text` | one streamed text chunk |
| `tool` | `name`, `content` | tool was called and returned |
| `done` | `content` | full answer complete |
| `error` | `message` | something failed |
| `session` | `session_id` | session ID for this turn |
