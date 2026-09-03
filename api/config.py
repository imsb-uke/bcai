import os
from pathlib import Path
from dotenv import load_dotenv

# Directories relative to code/ (Docker workspace root: /workspace)
CODE_DIR   = Path(__file__).resolve().parent.parent
CONFIG_DIR = CODE_DIR / "config"
MCP_DIR    = CODE_DIR / "mcp"
DOCS_DIR   = CODE_DIR / "docs"

# Load API keys and settings from config/env
load_dotenv(os.getenv("ENV_FILE", str(CONFIG_DIR / "env")), override=False)

# System prompt fed to the LLM agent
SYSTEM_FILE = Path(os.getenv("SYSTEM_FILE", str(DOCS_DIR / "system.md")))

# MCP servers to connect on startup (comma-separated paths override the default)
_DEFAULT_CONFIGS = [
    str(CONFIG_DIR / "BioChemAIgent_MCP_Server.json"),
    str(CONFIG_DIR / "PDB-MCP-Server_docker.json"),
    str(CONFIG_DIR / "ChEMBL_MCP_Server_docker.json"),
    str(CONFIG_DIR / "UniProt-MCP-Server_docker.json"),
]
MCP_CONFIG_FILES = [
    p.strip() for p in os.getenv("MCP_CONFIG_FILES", ",".join(_DEFAULT_CONFIGS)).split(",")
    if p.strip()
]

# Working directory for the MCP server subprocess
MCP_SERVER_CWD = str(MCP_DIR)

# Per-user data root (history, files, job_tmp) — absolute so MCP tools resolve correctly
USER_DATA_DIR = Path(os.getenv("USER_DATA_DIR", str(CODE_DIR / "user_data")))

# Default LLM model — format: "<source>|<model>"
DEFAULT_MODEL = os.getenv("LLM_MODEL", "openrouter|openai/gpt-4.1-mini")

# Free questions granted to a new user at registration
FREE_QUESTIONS_DEFAULT = int(os.getenv("FREE_QUESTIONS_DEFAULT", "50"))

# Max LLM<->tool round trips in one chat turn before forcing a final answer
MAX_TOOL_HOPS = int(os.getenv("MAX_TOOL_HOPS", "8"))

# Public tool names never exposed to the LLM (redundant/superseded by another
# tool the agent kept ignoring in favor of this one — see system.md)
DISABLED_TOOLS = {
    "PDB-MCP-Server_docker__pdb-server__download_structure",  # use get_pdb instead
}
