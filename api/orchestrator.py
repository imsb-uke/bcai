import asyncio
import os
import re
import json
from contextlib import AsyncExitStack
from pathlib import Path
from typing import Any, AsyncGenerator, Dict, List, Tuple

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from openai import AsyncOpenAI
from ollama import Client as OllamaClient

from . import config
from .history import save_history
from .turns import RunningTurn, TurnRegistry


# ----------------------------
# Small pure helpers
# ----------------------------
def slug(s: str) -> str:
    """Sanitize a string to match ^[a-zA-Z0-9_-]+$."""
    return re.sub(r"[^a-zA-Z0-9_-]", "_", s)

def str2dict(s: str) -> dict:
    try:
        return json.loads(s)
    except Exception:
        return {}

def load_servers_from_config(config_file_path: str) -> List[Tuple[str, Dict[str, Any]]]:
    """Read one MCP config file -> list of (unique_server_key, server_def)."""
    with open(config_file_path, "r") as f:
        cfg = json.load(f)
    base = os.path.splitext(os.path.basename(config_file_path))[0]
    return [
        (f"{slug(base)}__{slug(key)}", server_def)
        for key, server_def in cfg.get("mcpServers", {}).items()
    ]

def message_to_dict(msg) -> dict:
    """Convert an OpenAI ChatCompletionMessage to a plain JSON-serialisable dict."""
    d = {"role": getattr(msg, "role", None), "content": getattr(msg, "content", None)}
    for attr in ("refusal", "name", "function_call"):
        v = getattr(msg, attr, None)
        if v is not None:
            d[attr] = v
    # tool_calls are Pydantic objects — convert to plain dicts so they survive
    # JSON round-trips in the DB and can be replayed correctly on the next turn
    tcs = getattr(msg, "tool_calls", None)
    if tcs:
        d["tool_calls"] = [
            tc.model_dump() if hasattr(tc, "model_dump") else vars(tc)
            for tc in tcs
        ]
    return d


# ----------------------------
# Engine
# ----------------------------
class Engine:
    """Holds live MCP sessions and the tool routing table for the app lifetime."""

    def __init__(self) -> None:
        self._exit_stack     = AsyncExitStack()
        self.sessions        : Dict[str, ClientSession]       = {}  # server_key -> session
        self.route_map       : Dict[str, Tuple[str, str]]     = {}  # public tool name -> (server_key, bare name)
        self.tool_schemas    : Dict[str, Dict[str, Any]]      = {}  # public tool name -> inputSchema
        self.connected_servers : List[str]                    = []
        self.turns            = TurnRegistry()

    # ----------------------------
    # Lifecycle
    # ----------------------------
    async def connect(self) -> None:
        """Connect to every MCP server listed in config. Failed servers are skipped."""
        for cfg_path in config.MCP_CONFIG_FILES:
            try:
                servers = load_servers_from_config(cfg_path)
            except Exception as e:
                print(f"[engine] cannot read config {cfg_path}: {e}")
                continue
            for server_key, server_def in servers:
                try:
                    params = StdioServerParameters(
                        command = server_def.get("command"),
                        args    = server_def.get("args", []),
                        env     = {**os.environ, **server_def.get("env", {})},
                        cwd     = config.MCP_SERVER_CWD,   # server.py runs from mcp/
                    )
                    stdio, write = await self._exit_stack.enter_async_context(stdio_client(params))
                    session = await self._exit_stack.enter_async_context(ClientSession(stdio, write))
                    await session.initialize()
                    self.sessions[server_key] = session
                    self.connected_servers.append(server_key)
                    print(f"[engine] connected to '{server_key}'")
                except Exception as e:
                    print(f"[engine] FAILED to connect '{server_key}': {e}")

    async def close(self) -> None:
        await self._exit_stack.aclose()

    # ----------------------------
    # Tool surface
    # ----------------------------
    async def build_tools(self) -> List[Dict[str, Any]]:
        """Merge tools from all connected servers into OpenAI tool specs."""
        self.route_map.clear()
        self.tool_schemas.clear()
        tools: List[Dict[str, Any]] = []
        for server_key, session in self.sessions.items():
            for tool in (await session.list_tools()).tools:
                # Public name: <server_key>__<tool_name>, both sanitized
                public = f"{slug(server_key)}__{slug(tool.name)}"
                if public in config.DISABLED_TOOLS:
                    continue
                self.route_map[public] = (server_key, tool.name)
                self.tool_schemas[public] = tool.inputSchema or {}
                tools.append({
                    "type": "function",
                    "function": {
                        "name"       : public,
                        "description": f"[{server_key}] {tool.description or ''}".strip(),
                        "parameters" : tool.inputSchema or {"type": "object", "properties": {}},
                    },
                })
        return tools

    @staticmethod
    def _rewrite_paths(obj: Any, base: Path) -> Any:
        """Recursively rewrite any relative path starting with 'files' to absolute."""
        if isinstance(obj, str):
            s = obj.strip("/")
            if s and not s.startswith("/") and (s == "files" or s.startswith("files/")):
                return str(base / s)
            return obj
        if isinstance(obj, list):
            return [Engine._rewrite_paths(i, base) for i in obj]
        if isinstance(obj, dict):
            return {k: Engine._rewrite_paths(v, base) for k, v in obj.items()}
        return obj

    async def _call_tool(self, public_name: str, args: Dict[str, Any],
                         username: str = "") -> str:
        """Route a tool call to its MCP server and return text output."""
        # Rewrite all relative paths (files/...) to absolute user paths so the
        # agent never needs to know the server's directory structure.
        if username:
            args = self._rewrite_paths(args, config.USER_DATA_DIR / username)

        # Force the real authenticated username onto any tool that accepts one —
        # never trust whatever (if anything) the LLM filled into that field, since
        # async job tools key their completion notification off it.
        if username and "username" in (self.tool_schemas.get(public_name, {}).get("properties") or {}):
            args = {**args, "username": username}

        if public_name not in self.route_map:
            return f"Error: unknown tool '{public_name}'."
        server_key, bare = self.route_map[public_name]
        session = self.sessions.get(server_key)
        if session is None:
            return f"Error: server '{server_key}' not connected."
        result = await session.call_tool(bare, arguments=args or {})
        parts: List[str] = []
        for part in getattr(result, "content", []) or []:
            if hasattr(part, "text") and isinstance(part.text, str):
                parts.append(part.text)
            else:
                parts.append(str(part))
        output = "\n".join(p for p in parts if p) or "(no content)"

        return output

    # ----------------------------
    # Chat orchestration
    # ----------------------------
    async def stream_query(
        self, messages: List[Dict[str, Any]], model_full: str, username: str = "",
    ) -> AsyncGenerator[Dict[str, Any], None]:
        """Run one chat turn and yield events upstream:
          {"type": "tool",  "name": <tool>, "content": <result>}  -- a tool was called
          {"type": "token", "text": <chunk>}                      -- one piece of the streamed answer
          {"type": "done",  "content": <full>}                    -- answer complete
          {"type": "error", "message": <msg>}                     -- something failed

        Loops LLM <-> tool hops (up to config.MAX_TOOL_HOPS) so the model can
        chain multiple tool calls before answering, instead of stopping after
        a single round. On the hop past the limit, tools are withheld to force
        a final answer; if the model still returns tool calls at that point
        (ignoring tool_choice="none"), the turn ends with an error rather than
        looping forever.
        """
        llm_source, model = model_full.split("|", 1)
        client, client_type = self._make_client(llm_source)
        tools = await self.build_tools()

        try:
            hops = 0
            while True:
                hops += 1
                force_final = hops > config.MAX_TOOL_HOPS

                if client_type == "openai":
                    # Streaming round. If the model returns tool calls we collect
                    # the argument deltas by index; if it returns plain text we
                    # stream tokens directly and finish.
                    full_content            = ""
                    collected_tool_calls   : Dict[int, Dict[str, Any]] = {}
                    finish_reason           = None
                    native_finish_reason    = None
                    usage                   = None

                    stream = await client.chat.completions.create(
                        model=model, messages=messages,
                        tools=None if force_final else (tools or None),
                        tool_choice="none" if force_final else ("auto" if tools else "none"),
                        stream=True,
                    )
                    async for chunk in stream:
                        if getattr(chunk, "usage", None):
                            usage = chunk.usage
                        if not chunk.choices:
                            continue
                        choice = chunk.choices[0]
                        if choice.finish_reason:
                            finish_reason = choice.finish_reason
                            native_finish_reason = getattr(choice, "native_finish_reason", None) \
                                or (choice.model_extra or {}).get("native_finish_reason")
                        delta = choice.delta
                        if delta.content:
                            full_content += delta.content
                            yield {"type": "token", "text": delta.content}
                        if delta.tool_calls:
                            for tc_delta in delta.tool_calls:
                                i = tc_delta.index
                                if i not in collected_tool_calls:
                                    collected_tool_calls[i] = {
                                        "id": "", "type": "function",
                                        "function": {"name": "", "arguments": ""},
                                    }
                                if tc_delta.id:
                                    collected_tool_calls[i]["id"] += tc_delta.id
                                if tc_delta.function:
                                    if tc_delta.function.name:
                                        collected_tool_calls[i]["function"]["name"] += tc_delta.function.name
                                    if tc_delta.function.arguments:
                                        collected_tool_calls[i]["function"]["arguments"] += tc_delta.function.arguments

                    tool_calls_list = [collected_tool_calls[i] for i in sorted(collected_tool_calls)]

                    if not tool_calls_list:
                        if full_content:
                            # Normal case: the model answered in text, no more tools needed.
                            yield {"type": "done", "content": full_content}
                            return
                        # The completion produced neither visible text nor a tool call —
                        # e.g. OpenRouter's finish_reason:"error" after burning the
                        # generation budget on reasoning tokens. HTTP 200 does not mean
                        # the completion succeeded, so this must not be saved as a normal
                        # empty answer (see server.py's history-save guard).
                        print(f"[engine] EMPTY COMPLETION model={model} hop={hops} "
                              f"finish_reason={finish_reason} native_finish_reason={native_finish_reason} "
                              f"usage={usage}")
                        yield {
                            "type": "error",
                            "message": (
                                "The LLM did not generate any text response. "
                                "Please review the agent's tool calls above and continue the chat."
                            ),
                        }
                        return

                    # Persist the assistant message built from stream deltas
                    messages.append({
                        "role"      : "assistant",
                        "content"   : full_content or None,
                        "tool_calls": tool_calls_list,
                    })

                    if force_final:
                        # Reaching here means the model returned tool calls even
                        # though tools were withheld this hop (tools=None) — i.e.
                        # it hallucinated a call rather than answering. Stop
                        # rather than looping past the hop budget.
                        yield {
                            "type": "error",
                            "message": f"Exceeded max_tool_hops ({config.MAX_TOOL_HOPS}).",
                        }
                        return

                    # Execute each tool call, then loop to the next hop
                    for tc in tool_calls_list:
                        public_name = tc["function"]["name"]
                        args        = str2dict(tc["function"]["arguments"] or "{}")
                        try:
                            output = await self._call_tool(public_name, args, username)
                        except Exception as e:
                            output = f"Routing/Execution error: {e}"
                        yield {"type": "tool", "name": public_name, "content": output, "args": args}
                        messages.append({
                            "role"        : "tool",
                            "name"        : public_name,
                            "content"     : output,
                            "tool_call_id": tc["id"],
                        })

                else:
                    # Ollama — sync client, no streaming
                    resp          = client.chat(
                        model=model, messages=messages,
                        tools=None if force_final else tools, stream=False,
                    )
                    assistant_msg = resp.message
                    messages.append(message_to_dict(assistant_msg))

                    tool_calls = getattr(assistant_msg, "tool_calls", None)

                    if not tool_calls:
                        content = assistant_msg.content or ""
                        yield {"type": "token", "text": content}
                        yield {"type": "done",  "content": content}
                        return

                    if force_final:
                        yield {
                            "type": "error",
                            "message": f"Exceeded max_tool_hops ({config.MAX_TOOL_HOPS}).",
                        }
                        return

                    for tc in tool_calls:
                        public_name = tc.function.name
                        args        = str2dict(tc.function.arguments or "{}")
                        try:
                            output = await self._call_tool(public_name, args, username)
                        except Exception as e:
                            output = f"Routing/Execution error: {e}"
                        yield {"type": "tool", "name": public_name, "content": output}
                        messages.append({"role": "tool", "name": public_name, "content": output})
                    # loop continues to the next hop

        except Exception as e:
            print(f"[engine] EXCEPTION in stream_query: {type(e).__name__}: {e}")
            yield {"type": "error", "message": "No response from the LLM. Please try again."}

    # ----------------------------
    # Background turns — decoupled from the HTTP request (see turns.py):
    # a turn runs as an independent asyncio.Task, not as part of whatever
    # request happens to be watching it, so a client that disconnects
    # (session switch, logout/login, closed tab, lost connectivity) never
    # stops or loses an in-flight turn. Any number of requests can attach
    # to the same turn's live output at once via RunningTurn.subscribe().
    # ----------------------------

    def start_turn(
        self, username: str, session_id: str,
        messages: List[Dict[str, Any]], model_full: str,
    ) -> Tuple[RunningTurn, bool]:
        """
        Starts a new background turn for (username, session_id), unless one is
        already running for this session_id — in which case the existing turn
        is returned unchanged and `messages` is ignored (callers must check the
        returned `is_new` flag: False means nothing new was submitted).
        """
        existing = self.turns.get(session_id)
        if existing is not None:
            return existing, False

        async def run(turn: RunningTurn) -> None:
            await self._run_turn(turn, username, session_id, messages, model_full)

        return self.turns.start(session_id, username, run), True

    async def _run_turn(
        self, turn: RunningTurn, username: str, session_id: str,
        messages: List[Dict[str, Any]], model_full: str,
    ) -> None:
        """
        One chat turn: run stream_query, publishing each event to every
        current/future subscriber, then save history once the turn concludes
        — on a normal finish, an error event, or cancellation alike (the
        `finally` below), so a checkpoint is never lost regardless of how the
        turn ended. Cancellation (see TurnRegistry.cancel) is the only thing
        that stops a turn early — a disconnected watcher never does.
        """
        full_answer = ""
        had_error   = False
        try:
            async for event in self.stream_query(messages, model_full, username):
                if event["type"] == "token":
                    full_answer += event["text"]
                elif event["type"] == "error":
                    had_error = True
                turn.publish(event)
        except asyncio.CancelledError:
            turn.publish({"type": "error", "message": "Cancelled."})
            raise
        finally:
            # Mirrors the previous save-once-at-the-end behavior from
            # server.py's /chat handler, just no longer tied to that request.
            history = [m for m in messages if m.get("role") != "system"]
            if not (had_error and not full_answer):
                history = history + [{"role": "assistant", "content": full_answer}]
            save_history(username, session_id, history, model=model_full)

    # ----------------------------
    # LLM client factory
    # ----------------------------
    @staticmethod
    def _make_client(llm_source: str):
        """Define client"""
        if llm_source == "openai":
            return AsyncOpenAI(), "openai"
        if llm_source == "ollama_local":
            # Local Ollama; install models via: ollama run <model_name>
            return AsyncOpenAI(base_url="http://localhost:11434/v1", api_key="ollama"), "openai"
        if llm_source == "openrouter":
            return AsyncOpenAI(base_url="https://openrouter.ai/api/v1",
                               api_key=os.getenv("OPENROUTER_API_KEY")), "openai"
        if llm_source == "ollama_cloud":
            return OllamaClient(host="https://ollama.com",
                                headers={"Authorization": "Bearer " + os.getenv("OLLAMA_API_KEY", "")}), "ollama"
        raise ValueError(f"Unknown LLM source '{llm_source}'")
