"""
Shared job registry for all BCAI sub-agents.

All async sub-agents (screening, structure, literature, ...) import
_new_job and _update_job from here, so a single check_status() covers
every job type regardless of which sub-agent created it.
"""

import threading
import uuid
from datetime import datetime

_jobs: dict = {}
_jobs_lock = threading.Lock()
_stop_events: dict = {}  # claude: job_id -> threading.Event for cooperative cancellation


def _new_job(job_type: str) -> str:
    job_id = str(uuid.uuid4())[:8]
    with _jobs_lock:
        _jobs[job_id] = {
            "type":       job_type,
            "status":     "running",
            "progress":   "starting...",
            "result":     None,
            "error":      None,
            "started_at": datetime.now().isoformat(),
        }
    _stop_events[job_id] = threading.Event()  # claude
    return job_id


def get_stop_event(job_id: str):  # claude
    """Return the threading.Event for the given job (used by worker threads to check for cancellation)."""
    return _stop_events.get(job_id)


def _update_job(job_id: str, **kwargs):
    with _jobs_lock:
        if _jobs[job_id].get("status") == "cancelled":
            return   # job was cancelled; discard any late updates from the thread
        _jobs[job_id].update(kwargs)
        # claude: clean up stop_event when job reaches a terminal state
        if kwargs.get("status") in ("done", "failed"):
            _stop_events.pop(job_id, None)


def cancel_job(job_id: str) -> dict:
    """
    Cancel a running job. Sets the stop_event so the worker thread exits
    cooperatively at its next checkpoint (e.g. between compounds in screening).
    The current atomic step (e.g. one smina run) finishes before the flag is checked.
    """
    with _jobs_lock:
        job = _jobs.get(job_id)
        if job is None:
            return {"message": f"No job found with id '{job_id}'."}
        if job["status"] != "running":
            return {"message": f"Job '{job_id}' is not running (status: {job['status']})."}
        _jobs[job_id]["status"] = "cancelled"
    # claude ---
    stop_event = _stop_events.pop(job_id, None)  # remove from dict on cancel
    if stop_event:
        stop_event.set()
    # ---
    return {"job_id": job_id, "status": "cancelled",
            "message": "Cancellation requested. The job will stop at its next checkpoint."}


def check_status(job_id: str) -> dict:
    """Check the status of any running or completed job (screening, structure prediction, etc.)."""
    with _jobs_lock:
        job = _jobs.get(job_id)

    if job is None:
        return {"message": f"No job found with id '{job_id}'."}

    out = {
        "job_id":     job_id,
        "type":       job["type"],
        "status":     job["status"],
        "progress":   job["progress"],
        "started_at": job["started_at"],
    }
    if job["status"] == "done":
        out["result"] = job["result"]
    if job["status"] == "failed":
        out["error"] = job["error"]

    return out
