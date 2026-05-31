"""
Workflow pipeline — synchronous orchestration functions.

Composes tool calls from biochem and sub-agent pipelines into user-defined
multi-step workflows. No MCP, no threads — usable standalone or imported
by the workflow sub-agent server.

Workflow JSON format:
  {
    "name": "my_workflow",
    "description": "...",
    "created_at": "...",
    "steps": [
      {"id": "step1", "tool": "prepare_protein", "args": {...}},
      {"id": "step2", "tool": "prepare_library",  "args": {...}},
      {"id": "step3", "tool": "run_virtual_screening",
       "args": {"protein_df_file": "${step1.protein_df_file}",
                "library_dir":     "${step2.library_dir}"}}
    ]
  }

progress_cb: optional callable(str) for progress reporting.
             The sub-agent passes lambda msg: _update_job(job_id, progress=msg).
             Standalone callers can pass None or their own callback.
stop_event:  optional threading.Event for cooperative cancellation.
"""

import os
import sys
import json
import inspect
import re
from datetime import datetime

_HERE = os.path.dirname(os.path.abspath(__file__))
_BCAI_ROOT = os.path.dirname(os.path.dirname(_HERE))
if _BCAI_ROOT not in sys.path:
    sys.path.insert(0, _BCAI_ROOT)

import biochem as bc
from subagents.screening import pipeline as screening_pipeline
from subagents.literature import pipeline as lit_pipeline

_WORKFLOWS_DIR = os.getenv("WORKFLOWS_DIR", os.path.join(_HERE, "../../workflows"))


# ── Helpers ───────────────────────────────────────────────────────────────────

def _resolve(value, step_outputs: dict):
    """
    Recursively resolve ${stepId.key} placeholders in args using prior step outputs.

    Example: "${step1.protein_df_file}" → the actual path returned by step1.
    Unresolved placeholders (step not found or key missing) are left as-is.
    Works on strings, dicts, and lists recursively.
    """
    if isinstance(value, str):
        def _sub(m):
            step_id, key = m.group(1), m.group(2)
            return str(step_outputs.get(step_id, {}).get(key, m.group(0)))
        return re.sub(r'\$\{(\w+)\.(\w+)\}', _sub, value)
    elif isinstance(value, dict):
        return {k: _resolve(v, step_outputs) for k, v in value.items()}
    elif isinstance(value, list):
        return [_resolve(v, step_outputs) for v in value]
    return value


def _get_callable(tool_name: str):
    """
    Resolve a tool name to a callable.

    Resolution order:
      1. subagents.screening.pipeline  — gets full orchestration (tmp dirs, resume, etc.)
      2. subagents.literature.pipeline — full literature search + extraction
      3. biochem (bc)                  — all atomic biochem functions

    Checking pipeline modules before bc ensures that tools like prepare_library
    resolve to the full orchestration version, not the raw biochem primitive.
    """
    for module in (screening_pipeline, lit_pipeline, bc):
        fn = getattr(module, tool_name, None)
        if fn is not None and callable(fn):
            return fn
    raise ValueError(
        f"Unknown tool in workflow: '{tool_name}'. "
        f"Tool must be a function in biochem, screening pipeline, or literature pipeline."
    )


def _call_step(fn, args: dict, progress_cb, stop_event):
    """
    Call a step function with the resolved args, injecting progress_cb and
    stop_event only if the function's signature accepts them.

    This allows workflow steps to work with both biochem functions (which don't
    have progress_cb / stop_event) and pipeline functions (which do).
    """
    sig       = inspect.signature(fn)
    call_args = dict(args)
    if 'progress_cb' in sig.parameters:
        call_args['progress_cb'] = progress_cb
    if 'stop_event' in sig.parameters:
        call_args['stop_event'] = stop_event
    return fn(**call_args)


def _write_progress(progress_file: str, step_num: int, total: int, step_id: str, message: str):
    """
    Accumulate workflow progress in a JSON file in FILE_DIR.
    Each step is appended; sub-step messages update the current step in place.
    The UI reads this file when the user clicks Refresh in the sidebar.
    """
    # claude ---
    try:
        with open(progress_file) as f:
            data = json.load(f)
    except Exception:
        data = {"total": total, "steps": []}

    steps = data.get("steps", [])

    if step_id == "done":
        if steps and steps[-1]["status"] == "running":
            steps[-1]["status"] = "done"
        steps.append({"step": step_num, "step_id": "done", "message": "workflow complete", "status": "done"})
    elif steps and steps[-1]["step"] == step_num:
        steps[-1]["message"] = message          # sub-step progress: update in place
    else:
        if steps and steps[-1]["status"] == "running":
            steps[-1]["status"] = "done"        # previous step finished
        steps.append({"step": step_num, "step_id": step_id, "message": message, "status": "running"})

    data = {"total": total, "updated_at": datetime.now().isoformat(), "steps": steps}
    with open(progress_file, 'w') as f:
        json.dump(data, f, indent=2)
    # ---


# ── save_workflow ─────────────────────────────────────────────────────────────

def save_workflow(name: str, description: str, steps: list) -> dict:
    """
    Save a named workflow JSON to the workflows directory.

    Args:
        name:        short identifier used as filename (e.g. 'protein_screening').
                     Must be filesystem-safe (no spaces or special characters).
        description: human-readable description of what this workflow does.
        steps:       list of step dicts, each with:
                       id   — unique step identifier (e.g. 'step1')
                       tool — function name (e.g. 'prepare_protein', 'prepare_library')
                       args — dict of arguments. Use ${stepId.key} to reference
                              the output of a prior step, e.g.:
                              "protein_df_file": "${step1.protein_df_file}"

    Returns:
        {"message": str, "path": str, "n_steps": int}
    """
    os.makedirs(_WORKFLOWS_DIR, exist_ok=True)
    workflow = {
        "name":        name,
        "description": description,
        "created_at":  datetime.now().isoformat(),
        "steps":       steps,
    }
    path = os.path.join(_WORKFLOWS_DIR, f"{name}.json")
    with open(path, 'w') as f:
        json.dump(workflow, f, indent=2)
    return {
        "message": f"Workflow '{name}' saved.",
        "path":    path,
        "n_steps": len(steps),
    }


# ── delete_workflow ───────────────────────────────────────────────────────────

def delete_workflow(name: str) -> dict:
    """
    Delete a saved workflow JSON from the workflows directory.

    Args:
        name: workflow name to delete (same identifier used in save_workflow).

    Returns:
        {"message": str}
    """
    path = os.path.join(_WORKFLOWS_DIR, f"{name}.json")
    if not os.path.exists(path):
        return {"message": f"Workflow '{name}' not found."}
    os.remove(path)
    return {"message": f"Workflow '{name}' deleted."}


# ── list_workflows ────────────────────────────────────────────────────────────

def list_workflows() -> dict:
    """
    List all saved workflows in the workflows directory.

    Returns:
        {
            "workflows": [{"name", "description", "n_steps", "created_at"}, ...],
            "n_workflows": int,
        }
    """
    os.makedirs(_WORKFLOWS_DIR, exist_ok=True)
    workflows = []
    for fname in sorted(os.listdir(_WORKFLOWS_DIR)):
        if not fname.endswith('.json'):
            continue
        try:
            with open(os.path.join(_WORKFLOWS_DIR, fname)) as f:
                w = json.load(f)
            workflows.append({
                "name":        w.get("name", fname[:-5]),
                "description": w.get("description", ""),
                "n_steps":     len(w.get("steps", [])),
                "created_at":  w.get("created_at", ""),
            })
        except Exception:
            pass
    return {"workflows": workflows, "n_workflows": len(workflows)}


# ── execute_workflow ──────────────────────────────────────────────────────────

def execute_workflow(
    workflow_name: str,
    file_dir:      str  = None,
    session_dir:   str  = None,  # claude: per-user dir for progress file
    stop_event           = None,
    progress_cb          = None,
) -> dict:
    """
    Execute a saved workflow synchronously, running each step in sequence.

    After each step completes, its output dict is stored and made available
    for ${stepId.key} placeholder resolution in all subsequent steps.
    If a step fails, execution stops immediately with a RuntimeError that
    identifies the failing step.

    Args:
        workflow_name: name of the workflow to run (must match a saved JSON).
        file_dir:      base output directory for generated files.
        session_dir:   per-user session directory; workflow_progress.json is
                       written here so multiple users don't overwrite each other.
        stop_event:    optional threading.Event; checked between steps.
                       When set, execution stops cleanly at the next boundary.
        progress_cb:   optional callable(str) for progress messages.
                       The server passes lambda msg: _update_job(job_id, progress=msg).

    Returns:
        {
            "workflow_name": str,
            "n_steps":       int,
            "step_outputs":  {step_id: result_dict, ...},
        }
    """
    file_dir = file_dir or os.getenv("FILE_DIR", "files")

    def _progress(msg):
        if progress_cb:
            progress_cb(msg)

    # Load workflow JSON
    workflow_path = os.path.join(_WORKFLOWS_DIR, f"{workflow_name}.json")
    if not os.path.exists(workflow_path):
        raise FileNotFoundError(
            f"Workflow '{workflow_name}' not found in '{_WORKFLOWS_DIR}'. "
            f"Use list_workflows() to see available workflows."
        )

    with open(workflow_path) as f:
        workflow = json.load(f)

    steps         = workflow.get("steps", [])
    total         = len(steps)
    # claude: write progress to session dir (per-user); fall back to file_dir
    _progress_dir = session_dir or file_dir
    os.makedirs(_progress_dir, exist_ok=True)
    progress_file = os.path.join(_progress_dir, "workflow_progress.json")
    # Clear any stale progress from a previous run  # claude
    if os.path.exists(progress_file):
        os.remove(progress_file)
    step_outputs  = {}

    os.makedirs(file_dir, exist_ok=True)

    for i, step in enumerate(steps, start=1):
        # Cooperative cancellation — checked between every step
        if stop_event and stop_event.is_set():
            _progress("cancelled")
            break

        step_id   = step["id"]
        tool_name = step["tool"]
        args      = _resolve(step.get("args", {}), step_outputs)

        msg = f"step {i}/{total}: {tool_name}"
        _progress(msg)
        _write_progress(progress_file, i, total, step_id, tool_name)  # claude: UI prepends "step N/M:"

        try:
            fn     = _get_callable(tool_name)
            result = _call_step(
                fn, args,
                progress_cb = lambda m: _write_progress(progress_file, i, total, step_id, m),
                stop_event  = stop_event,
            )
            # Store result for placeholder resolution in later steps
            step_outputs[step_id] = result if isinstance(result, dict) else {"result": result}

        except Exception as e:
            raise RuntimeError(f"Step '{step_id}' ({tool_name}) failed: {e}") from e

    _write_progress(progress_file, total, total, "done", "workflow complete")
    _progress("done")

    return {
        "workflow_name": workflow_name,
        "n_steps":       total,
        "step_outputs":  step_outputs,
    }
