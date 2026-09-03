"""
Workflow Sub-Agent — async MCP server for user-defined multi-step workflows.

Tools exposed:
  save_workflow(name, description, steps)  → saves workflow JSON (sync)
  list_workflows()                          → lists available workflows (sync)
  run_workflow(workflow_name, file_dir)     → job_id (async)

Typical workflow:
  1. save_workflow("protein_screening", "Screen a library against a target", steps=[...])
  2. run_workflow("protein_screening")  →  job_id
  3. check_status(job_id)              →  progress / final result
     (or click Refresh in the UI sidebar to see step-by-step progress)

Run standalone from bcai/subagents/workflow/:
    python server.py
"""

import os
import sys
import threading
from pathlib import Path

from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP

# ── Config ────────────────────────────────────────────────────────────────────

_HERE = os.path.dirname(os.path.abspath(__file__))
load_dotenv(os.path.join(_HERE, "../../env"))

_BCAI_ROOT = os.path.dirname(os.path.dirname(_HERE))
if _BCAI_ROOT not in sys.path:
    sys.path.insert(0, _BCAI_ROOT)

from subagents.jobs import _new_job, _update_job, get_stop_event  # claude
from subagents.workflow import pipeline

# ── MCP server ────────────────────────────────────────────────────────────────

mcp = FastMCP(name="WorkflowAgent")


# ── Tool 1: save_workflow (sync) ──────────────────────────────────────────────

@mcp.tool()
def save_workflow(name: str, description: str, steps: list) -> dict:
    """
    Save a named workflow for later re-use.

    name:        short identifier used as filename (e.g. 'protein_screening').
                 Must be filesystem-safe (no spaces or special characters).
    description: human-readable description of what this workflow does.
    steps:       list of step dicts, each with:
                   id   — unique step identifier (e.g. 'step1')
                   tool — function name (e.g. 'prepare_protein', 'prepare_library')
                   args — dict of arguments. Use ${stepId.key} to reference
                          the output of a prior step, e.g.:
                          "protein_df_file": "${step1.protein_df_file}"
    """
    return pipeline.save_workflow(name=name, description=description, steps=steps)


# ── Tool 2: list_workflows (sync) ─────────────────────────────────────────────

@mcp.tool()
def list_workflows() -> dict:
    """List all saved workflows with name, description, step count, and creation date."""
    return pipeline.list_workflows()


# ── Tool 3: delete_workflow (sync) ────────────────────────────────────────────

@mcp.tool()
def delete_workflow(name: str) -> dict:
    """
    Delete a saved workflow by name.

    name: the workflow name to delete (same identifier used in save_workflow).
    """
    return pipeline.delete_workflow(name=name)


# ── Tool 4: run_workflow (async) ──────────────────────────────────────────────

@mcp.tool()
def run_workflow(
    workflow_name: str,
    username:      str,
    file_dir:      str,
) -> dict:
    """
    Run a saved workflow asynchronously. Returns a job_id immediately.
    Progress is written to the user's session/ folder (derived from file_dir).
    Click Refresh in the UI sidebar to see the current step.
    Use check_status(job_id) for the final result, or cancel_job(job_id) to stop.

    workflow_name: name of the workflow to run (as saved by save_workflow).
    username:      the logged-in user's name (used to write per-user progress file).
    file_dir:      base output directory for generated files.
    """
    # Derive from file_dir: .../user_data/<user>/files/ → .../user_data/<user>/job_tmp/
    session_dir = str(Path(file_dir).parent / "job_tmp")

    # Validate before starting the thread so the error is returned immediately
    if not any(w["name"] == workflow_name for w in pipeline.list_workflows()["workflows"]):
        return {
            "message": (
                f"Workflow '{workflow_name}' not found. "
                f"Use list_workflows() to see available workflows."
            )
        }

    job_id     = _new_job("run_workflow", username=username)
    stop_event = get_stop_event(job_id)  # claude

    def _run():
        try:
            result = pipeline.execute_workflow(
                workflow_name = workflow_name,
                file_dir      = file_dir,
                session_dir   = session_dir,  # claude
                stop_event    = stop_event,
                progress_cb   = lambda msg: _update_job(job_id, progress=msg),
            )
            _update_job(job_id, status="done", progress="done", result=result)
        except Exception as e:
            _update_job(job_id, status="failed", error=str(e))

    threading.Thread(target=_run, daemon=True).start()

    return {
        "message":    f"Workflow '{workflow_name}' started.",
        "job_id":     job_id,
        "status":     "running",
        "next_steps": [
            "check_status(job_id)",
            "click Refresh in the UI sidebar to see step progress",
        ],
    }


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    mcp.run(transport="stdio")
