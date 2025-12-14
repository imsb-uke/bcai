import os
import json
import time
import asyncio
from datetime import datetime
from dotenv import load_dotenv
from contextlib import AsyncExitStack
from typing import Any, Dict, List, Tuple, Optional

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from openai import AsyncOpenAI
from ollama import Client

import threading

load_dotenv("../env")

import src.utils_client as utils
import src.user as user
from src.history import save_history, load_history, message_to_dict

os.system(f"echo 'Run client.py {datetime.now()}' > client.log")

# ----------------------------
# Configurations and initiations
# ----------------------------
# Set working directory
client = None
client_type = None

SYSTEM_FILE_DIR = os.getenv("SYSTEM_FILE_DIR")
WORKING_DIR = "work_space"
SESSION_PARENT_DIR = "session"

# Server config files
SERVER_CONFIG_FILES = [
    "../config/ChEMBL_MCP_Server_docker.json",
    "../config/PDB-MCP-Server_docker.json",
    "../config/BioChemAIgent_MCP_Server.json",
]

# List LLM models
llm_model_set = []
with open("../llm_models", 'r') as file:
    for l in file:
        l = l.strip()
        if l:
            llm_model_set.append(l)

# Multi-server state
exit_stack = AsyncExitStack()
sessions: Dict[str, ClientSession] = {}           # server_key -> ClientSession
writers: Dict[str, Any] = {}                      # server_key -> write handle (unused but kept for completeness)
# Tool routing map: public tool name -> (server_key, original_tool_name)
route_map: Dict[str, Tuple[str, str]] = {}

# ----------------------------
# Server utilities
# ----------------------------
def make_client(llm_source):
    """Define client"""
    if llm_source == 'openai':             
        # OpenAI client/model
        client = AsyncOpenAI()
        client_type = 'openai'
    elif llm_source == 'ollama_local':         
        # Ollama local client/model; Install each model via $ ollama run model_name
        client = AsyncOpenAI(base_url="http://localhost:11434/v1", api_key="ollama")
        client_type = 'openai'
    elif llm_source == 'ollama_cloud':     
        # Ollama cloud client/model
        client = Client(host="https://ollama.com", headers={'Authorization': 'Bearer ' + os.getenv("OLLAMA_API_KEY")})
        client_type = 'ollama'
    elif llm_source == 'openrouter':
        client = AsyncOpenAI(base_url="https://openrouter.ai/api/v1", api_key=os.getenv("OPENROUTER_API_KEY"))
        client_type = 'openai'
    return client, client_type

def load_system_message():
    """Load the system message for system"""
    try:
        with open (SYSTEM_FILE_DIR, 'r') as file:
            text = file.read()
    except:
        text = ""
    return text
    
def route(public_tool_name: str) -> Tuple[str, str]:
    """
    Look up (server_key, original_tool_name) from a sanitized public tool name.
    """
    if public_tool_name not in route_map:
        raise ValueError(f"Unknown tool '{public_tool_name}'.")
    return route_map[public_tool_name]
    
async def connect_to_all_servers(config_files: List[str], list_tools: bool = True) -> None:
    """
    Connect to all servers across all config files. Populates global 'sessions' and 'writers'.
    """
    global sessions, writers, exit_stack

    for cfg_path in config_files:
        for server_key, server_def in utils.load_servers_from_config(cfg_path):
            # StdioServerParameters expects command + args + env
            params = StdioServerParameters(
                command=server_def.get("command"),
                args=server_def.get("args", []),
                env=server_def.get("env", {}),
            )
            stdio_transport = await exit_stack.enter_async_context(stdio_client(params))
            stdio, write = stdio_transport

            session = await exit_stack.enter_async_context(ClientSession(stdio, write))
            await session.initialize()

            sessions[server_key] = session
            writers[server_key] = write

            if list_tools:
                tools_result = await session.list_tools()
                print(f"\nConnected to '{server_key}' with tools:")
                for tool in tools_result.tools:
                    print(f"  - {tool.name}: {tool.description}")

async def build_tools() -> List[Dict[str, Any]]:
    """
    Merge tools from all connected MCP servers, exposing them to OpenAI with **sanitized** names.
    Fills the global 'route_map' so we can route calls back to the correct MCP server/tool.
    """
    global route_map
    route_map.clear()

    merged_tools: List[Dict[str, Any]] = []
    for server_key, session in sessions.items():
        tools_result = await session.list_tools()
        for tool in tools_result.tools:
            # Public name: <server_key>__<tool_name>, both sanitized
            public_name = f"{utils.slug(server_key)}__{utils.slug(tool.name)}"
            route_map[public_name] = (server_key, tool.name)

            merged_tools.append({
                "type": "function",
                "function": {
                    "name": public_name,  # MUST match ^[a-zA-Z0-9_-]+$
                    "description": f"[{server_key}] {tool.description or ''}".strip(),
                    "parameters": tool.inputSchema or {"type": "object", "properties": {}},
                },
            })
    return merged_tools


async def mcp_call(
    server_key: str,
    bare_tool_name: str,
    arguments: Optional[Dict[str, Any]],
) -> str:
    """
    Call an MCP tool on a specific session and return a best-effort text rendering of the result.
    """
    session = sessions.get(server_key)
    if session is None:
        return f"Error: server '{server_key}' is not connected."

    args = arguments or {}
    result = await session.call_tool(bare_tool_name, arguments=args)

    # return result.content[0].text      # simple version, for advanced -> below
    # Normalize tool result content to a simple string
    text_parts: List[str] = []
    if getattr(result, "content", None):
        for part in result.content:
            # Common MCP part types have .text, .type, etc.
            if hasattr(part, "text") and isinstance(part.text, str):
                text_parts.append(part.text)
            else:
                # Fallback to JSON-ish string
                try:
                    text_parts.append(json.dumps(getattr(part, "__dict__", str(part)), ensure_ascii=False))
                except Exception:
                    text_parts.append(str(part))

    return "\n".join([p for p in text_parts if p]) or "(no content)"

async def cleanup() -> None:
    global exit_stack
    await exit_stack.aclose()

# ----------------------------
# Chat Orchestration
# ----------------------------
async def process_query(history_file: str, model: str) -> str:
    """
    Send the user's query to LLM, allow tool calls, execute them across MCP sessions,
    then get a final answer.
    """
    # Prepare tool surface
    tools = await build_tools()

    # Build message history
    chat_history = load_history(history_file)
    if chat_history[0]['role'] == 'system':
        messages = chat_history
    else:
        system_message = load_system_message()
        messages = [{"role": "system", "content": system_message}] + chat_history
    # messages.append({"role": "user", "content": query})

    # First round (tools allowed)
    if client_type == "openai":
        first = await client.chat.completions.create(
            model=model,
            messages=messages,
            tools=tools,
            tool_choice="auto",
        )
        assistant_msg = first.choices[0].message
    else:
        first = client.chat(
            model=model,
            messages=messages,
            tools=tools,
            stream=False,
        )
        assistant_msg = first.message

    messages.append(message_to_dict(assistant_msg))

    # Handle tool calls
    if getattr(assistant_msg, "tool_calls", None):
        for tc in assistant_msg.tool_calls:
            
            public_name = tc.function.name
            args = tc.function.arguments or {}
            if isinstance(args, str):                 # For OpnenAI, args is string, then needs json.loads(raw_args)
                args = utils.str2dict(args)

            try:
                server_key, bare_tool_name = route(public_name)
                tool_output = await mcp_call(server_key, bare_tool_name, args)
            except Exception as e:
                tool_output = f"Routing/Execution error: {e}"

            tool_message = {
                "role": "tool",
                "name": public_name,
                "content": tool_output,
            }
            if client_type == "openai":
                tool_message["tool_call_id"] = tc.id
            messages.append(tool_message)

        # Final round (no further tool calls)
        if client_type == "openai":
            final = await client.chat.completions.create(
                model=model,
                messages=messages,
                tools=tools,
                tool_choice="none",
            )
            assistant_msg = final.choices[0].message
        else:
            final = client.chat(
                model=model,
                messages=messages,
                tools=tools,
                stream=False,
            )
            assistant_msg = final.message
        messages.append(message_to_dict(assistant_msg))

    save_history(messages, history_file)
    return messages[-1]["content"]

# ----------------------------
# Main Loop
# ----------------------------
async def main() -> None:
    print("\n=========================")
    print("Welcome to BioChemAIgent Chat, a Drug Developer Agent ....")
    print(f"Version: 0.9.0")
    print("=========================\n")
    print("Here is the list of availbe tools\n")
        
    # Connect to all declared servers
    await connect_to_all_servers(SERVER_CONFIG_FILES, list_tools=True)

    global client, client_type
    
    os.makedirs(WORKING_DIR, exist_ok=True)

    while True:
        # List files at work space and choose the first one
        file_list = [f for f in os.listdir(WORKING_DIR) if f.split('.')[0]]
        file = os.path.join(WORKING_DIR, file_list[0]) if file_list else None
        time.sleep(0.05)
        if file and os.path.getsize(file) > 0:
            # print(f'\nRunning {file}')
            os.system(f"echo '{datetime.now()} | Running {file}' >> client.log")

            try:
                # Load session info of the user
                user_name = os.path.splitext(os.path.basename(file))[0]
                session_dir = os.path.join(SESSION_PARENT_DIR, user_name)
                session_file = os.path.join(session_dir, "chat_session.json")
                with open(session_file, 'r') as f:
                    session_info = json.load(f)
    
                if user_name != session_info['user_name']:
                    # print('Username does not match')
                    os.system(f"echo 'Username does not match' >> client.log")
    
                # Set the model and user tokens
                llm_source, model = session_info['model_full'].split("|", 1)
                # print(f'LLM model: {llm_source} {model}')
                os.system(f"echo 'LLM model: {llm_source} {model}' >> client.log")
                os.environ["OPENAI_API_KEY"] = session_info['openai_api_key']
                os.environ["OLLAMA_API_KEY"] = session_info['ollama_api_key']
                os.environ["OPENROUTER_API_KEY"] = session_info['openrouter_api_key']
                client, client_type = make_client(llm_source)
                
                response = ''
                counter = 0
                while response == '' and counter < 5:
                    counter += 1
                    response = await process_query(file, model)
                    if response == '':
                        user_input = ''
                        print('Processing ...')
                        os.system(f"echo 'Processing ...' >> client.log")

                if response == '' and counter < 5:
                    response = 'I could not solve the task in 5 iterations. Please simplify the tesk.'
                    
                    # Move the history file back
                    # os.system(f"mv {file} {session_info['history_file']}")
                    
            except Exception as e:
                response = f"Error: {e}"
                os.system(f"echo {response} >> client.log")
                save_history(load_history(file) + [{"role": "assistant", "content": response}], file)

            # Move the history file back
            os.system(f"mv {file} {session_info['history_file']}")

    # # Cleanup
    # try:
    #     if os.path.exists(HISTORY_FILE):
    #         os.remove(HISTORY_FILE)
    # except Exception:
    #     pass
    # await cleanup()

if __name__ == "__main__":
    asyncio.run(main())
