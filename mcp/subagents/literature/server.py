"""
Literature Sub-Agent — async MCP server for literature search and LLM extraction.

Tools exposed:
  search_literature(query, source, max_results, file_dir)  → job_id  (async)
  extract_information(papers_file, question, file_dir)      → job_id  (async)
  check_status(job_id)                                      → status + result (sync)

Typical workflow:
  1. job1 = search_literature("EGFR inhibitors NSCLC IC50")
  2. check_status(job1)  →  result.papers_file
  3. job2 = extract_information(papers_file, "What compounds were tested and what were their IC50 values?")
  4. check_status(job2)  →  result.csv_file

LLM model for extraction is set via LLM_MODEL_LITERATURE in bcai/env, e.g.:
    LLM_MODEL_LITERATURE = 'openrouter|openai/gpt-4o-mini'

Run standalone from bcai/subagents/literature/:
    python server.py
"""

import os
import sys
import threading

from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP

# ── Config ────────────────────────────────────────────────────────────────────

_HERE = os.path.dirname(os.path.abspath(__file__))
load_dotenv(os.path.join(_HERE, "../../env"))

# Make bcai/ importable (needed for subagents.jobs and standalone run)
_BCAI_ROOT = os.path.dirname(os.path.dirname(_HERE))
if _BCAI_ROOT not in sys.path:
    sys.path.insert(0, _BCAI_ROOT)

from subagents.jobs import _new_job, _update_job, get_stop_event  # claude: added get_stop_event
from subagents.literature import pipeline

# ── MCP server ────────────────────────────────────────────────────────────────

mcp = FastMCP(name="LiteratureAgent")


# ── Tool 1: search_literature (async) ────────────────────────────────────────

@mcp.tool()
def search_literature(
    query:       str,
    file_dir:    str,
    source:      str = "both",
    max_results: int = 20,
    username:    str = "",
) -> dict:
    """
    Search the scientific literature and retrieve papers with abstracts.
    Returns a job_id immediately; use check_status(job_id) to get results.

    query:       free-text search query (e.g. 'EGFR inhibitors IC50 NSCLC').
    source:      'pubmed' | 'semantic_scholar' | 'both' (default: 'both').
    max_results: maximum papers to retrieve per source.
    file_dir:    directory to save the output JSON file.
    """
    job_id     = _new_job("search_literature", username=username)
    stop_event = get_stop_event(job_id)  # claude

    def _run():
        if stop_event and stop_event.is_set():  # claude: exit early if already cancelled
            return
        try:
            result = pipeline.search_literature(
                query       = query,
                source      = source,
                max_results = max_results,
                file_dir    = file_dir,
                progress_cb = lambda msg: _update_job(job_id, progress=msg),
            )
            _update_job(job_id, status="done", progress="done", result=result)
        except Exception as e:
            _update_job(job_id, status="failed", error=str(e))

    threading.Thread(target=_run, daemon=True).start()

    return {
        "message":    f"Literature search started: '{query}' ({source}).",
        "job_id":     job_id,
        "status":     "running",
        "next_steps": ["check_status(job_id)"],
    }


# ── Tool 2: extract_information (async) ──────────────────────────────────────

@mcp.tool()
def extract_information(
    papers_file: str,
    question:    str,
    file_dir:    str,
    username:    str = "",
) -> dict:
    """
    For each paper in papers_file, ask the LLM the given question and save
    a structured CSV of answers. Returns a job_id immediately.

    papers_file: path to the JSON file produced by search_literature().
    question:    the question to answer from each abstract, e.g.:
                 'What compounds were tested and what were their IC50 values?'
                 'What proteins or receptors are associated with this disease?'
                 'What cell lines or animal models were used?'
    file_dir:    used to resolve relative paths.
    """

    if not os.path.exists(papers_file):
        return {"message": f"papers_file not found: '{papers_file}'"}

    if not (os.getenv("OPENAI_API_KEY") or os.getenv("OPENROUTER_API_KEY") or os.getenv("OLLAMA_API_KEY")):
        return {"message": "No LLM API key found. Set OPENAI_API_KEY, OPENROUTER_API_KEY, or OLLAMA_API_KEY in env."}

    job_id     = _new_job("extract_information", username=username)
    stop_event = get_stop_event(job_id)  # claude

    def _run():
        if stop_event and stop_event.is_set():  # claude: exit early if already cancelled
            return
        try:
            result = pipeline.extract_information(
                papers_file = papers_file,
                question    = question,
                file_dir    = file_dir,
                progress_cb = lambda msg: _update_job(job_id, progress=msg),
            )
            _update_job(job_id, status="done", progress="done", result=result)
        except Exception as e:
            _update_job(job_id, status="failed", error=str(e))

    threading.Thread(target=_run, daemon=True).start()

    return {
        "message":    f"Information extraction started for '{os.path.basename(papers_file)}'.",
        "job_id":     job_id,
        "status":     "running",
        "next_steps": ["check_status(job_id)"],
    }


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    mcp.run(transport="stdio")
