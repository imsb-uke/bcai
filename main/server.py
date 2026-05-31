import os
import sys
from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP

## ======= Import tools and subagents =======

# Import biochem tools
load_dotenv("../env")
from biochem import *
# TOOLS is already defined in biochem/__init__.py

# Make bcai/ importable so subagents can be reached as subagents.*
_BCAI_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _BCAI_ROOT not in sys.path:
    sys.path.insert(0, _BCAI_ROOT)

from subagents.jobs import check_status, cancel_job

from subagents.screening.server import (
    list_libraries,
    prepare_library,
    run_screening,
    get_top_hits,
)

from subagents.structure.server import (
    run_af3,
    run_esm3,
)

from subagents.literature.server import (
    search_literature,
    extract_information,
)

from subagents.workflow.server import (
    save_workflow,
    list_workflows,
    delete_workflow,
    run_workflow,
)

SCREENING_TOOLS   = [list_libraries, prepare_library, run_screening, get_top_hits]
STRUCTURE_TOOLS   = [run_af3, run_esm3]
LITERATURE_TOOLS  = [search_literature, extract_information]
WORKFLOW_TOOLS    = [save_workflow, list_workflows, delete_workflow, run_workflow]

# Tools handled by async sub-agents — skip from the sync mother registration
_ASYNC_TOOLS = {"run_af3", "run_esm3"}

## ======= Register tools and subagents =======
# Create an MCP server
mcp = FastMCP(name="BioChemAIgent")

# Register biochem tools (as the Mother agent), skipping those promoted to async sub-agents
for fn in TOOLS:
    if fn.__name__ in _ASYNC_TOOLS:
        continue
    mcp.tool()(fn)

# Register screening subagent tools
for fn in SCREENING_TOOLS:
    mcp.tool()(fn)

# Register structure subagent tools (async wrappers for run_af3 / run_esm3)
for fn in STRUCTURE_TOOLS:
    mcp.tool()(fn)

# Register literature subagent tools
for fn in LITERATURE_TOOLS:
    mcp.tool()(fn)

# Register workflow subagent tools
for fn in WORKFLOW_TOOLS:
    mcp.tool()(fn)

# Single check_status / cancel_job cover jobs from all sub-agents
mcp.tool()(check_status)
mcp.tool()(cancel_job)

# Run the server
if __name__ == "__main__":
    mcp.run(transport="stdio")