"""
Shared job registry for all BCAI sub-agents.

All async sub-agents (screening, structure, literature, ...) import
_new_job and _update_job from here, so a single check_status() covers
every job type regardless of which sub-agent created it.
"""

import json
import os
import threading
import uuid
from datetime import datetime
from pathlib import Path

_jobs: dict = {}
_jobs_lock = threading.Lock()
_stop_events: dict = {}  # claude: job_id -> threading.Event for cooperative cancellation

_USER_DATA_DIR = os.getenv("USER_DATA_DIR")  # set in Docker env; None in local dev


def _write_notification(username: str, job_id: str, job_type: str, status: str) -> None:
    """Write a job-ready notification so the frontend picks it up within 5 s."""
    if not _USER_DATA_DIR or not username:
        return
    notif_file = Path(_USER_DATA_DIR) / username / "job_tmp" / "job_ready.json"
    try:
        notif_file.parent.mkdir(parents=True, exist_ok=True)
        existing: list = []
        if notif_file.exists():
            try:
                existing = json.loads(notif_file.read_text())
            except Exception:
                pass
        existing.append({
            "job_id":  job_id,
            "type":    job_type,
            "status":  status,
            "message": f"'{job_type}' is ready.",
        })
        notif_file.write_text(json.dumps(existing))
    except Exception:
        pass  # never crash the worker thread over a notification


def _new_job(job_type: str, username: str = "") -> str:
    job_id = str(uuid.uuid4())[:8]
    with _jobs_lock:
        _jobs[job_id] = {
            "type":       job_type,
            "status":     "running",
            "progress":   "starting...",
            "result":     None,
            "error":      None,
            "started_at": datetime.now().isoformat(),
            "username":   username,
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
            username = _jobs[job_id].get("username", "")
            job_type = _jobs[job_id].get("type", job_id)
            status   = kwargs["status"]
    if kwargs.get("status") in ("done", "failed") and username:
        _write_notification(username, job_id, job_type, status)


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
