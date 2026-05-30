import os
import re
import json
import time
import streamlit as st
from typing import Any, Dict, List, Tuple, Optional


def slug(s: str) -> str:
    """Sanitize any string to match ^[a-zA-Z0-9_-]+$ by replacing invalid chars with '_'."""
    return re.sub(r'[^a-zA-Z0-9_-]', '_', s)

def _unique_server_key(config_file_path: str, server_key_in_file: str) -> str:
    """
    Build a unique (and sanitized) server key from the file name and the key inside that file.
    This prevents collisions across multiple files with similarly named servers.
    """
    base = os.path.splitext(os.path.basename(config_file_path))[0]
    return f"{slug(base)}__{slug(server_key_in_file)}"

def load_servers_from_config(config_file_path: str) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Read one config file and return list of (unique_server_key, server_def).
    The config must have a top-level "mcpServers" object mapping keys -> server defs.
    """
    with open(config_file_path, "r") as f:
        cfg = json.load(f)
    servers = cfg.get("mcpServers", {})
    out: List[Tuple[str, Dict[str, Any]]] = []
    for key, server_def in servers.items():
        out.append((_unique_server_key(config_file_path, key), server_def))
    return out

def str2dict(args: str) -> dict:
    try:
        args = json.loads(args)
    except Exception:
        args = {}
    return args

def escape_brackets_outside_code(text):
    # Regex: split into code spans (`...`) and normal text
    parts = re.split(r'(`[^`]*`)', text)

    escaped_parts = []
    for part in parts:
        if part.startswith("`") and part.endswith("`"):
            # It's code → keep unchanged
            escaped_parts.append(part)
        else:
            # It's normal text → escape brackets
            escaped_parts.append(
                part.replace("[", r"\[").replace("]", r"\]")
            )
    return "".join(escaped_parts)

def message_box(text, bg="#f5f5f5", border="#dddddd", height=10):
    st.markdown(
        f"""
        <div style="
            border:1px solid {border};
            background-color:{bg};
            padding:{height}px 18px;
            border-radius:30px;   /* more circular */
            margin-bottom:10px;
            display:inline-block;
        ">
            {text}
        </div>
        """,
        unsafe_allow_html=True
    )

# Not used
def check_file_change(path: str):
    # Get initial modification time
    last_mtime = os.path.getmtime(path)
    time.sleep(0.05)
    # Get next modification time
    current_mtime = os.path.getmtime(path)
    return current_mtime != last_mtime and os.path.getsize(path) > 0

