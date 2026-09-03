"""
Screening Sub-Agent — async virtual screening MCP server.

Tools exposed:
  list_libraries()                         → available compound libraries (sync)
  prepare_library(source, ...)             → job_id  (async, can take hours)
  run_screening(protein_df, library, ...)  → job_id  (async, can take hours)
  check_status(job_id)                     → status + progress for any job (sync)
  get_top_hits(project_name, ...)          → ranked hits from a finished run (sync)

The two async tools are thin wrappers over pipeline.py, which contains the
actual orchestration logic. Same pattern as the structure sub-agent.

Run standalone from bcai/subagents/screening/:
    python server.py
"""

import os
import sys
import threading

import pandas as pd
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
from subagents.screening import pipeline

# Base directory for all prepared libraries
LIBRARY_DIR = os.getenv("LIBRARY_DIR")

# ── MCP server ────────────────────────────────────────────────────────────────

mcp = FastMCP(name="ScreeningAgent")


# ── Tool 1: list_libraries ────────────────────────────────────────────────────

@mcp.tool()
def list_libraries() -> dict:
    """List all available compound libraries. Shows whether each is prepared (ready for screening) and how many compounds it contains."""

    os.makedirs(LIBRARY_DIR, exist_ok=True)
    entries = []

    for name in sorted(os.listdir(LIBRARY_DIR)):
        lib_dir = os.path.join(LIBRARY_DIR, name)
        if not os.path.isdir(lib_dir):
            continue

        df_ligand_file = os.path.join(lib_dir, "df_ligand.csv")
        prepared       = os.path.exists(df_ligand_file)
        n_compounds    = 0
        if prepared:
            try:
                n_compounds = len(pd.read_csv(df_ligand_file))
            except Exception:
                pass

        entries.append({
            "name":        name,
            "prepared":    prepared,
            "n_compounds": n_compounds,
            "path":        lib_dir,
        })

    if not entries:
        return {"message": "No libraries found. Use prepare_library() to add one.", "libraries": []}

    return {
        "message":   f"{len(entries)} librar{'y' if len(entries) == 1 else 'ies'} available.",
        "libraries": entries,
    }


# ── Tool 2: prepare_library (async) ──────────────────────────────────────────

@mcp.tool()
def prepare_library(
    source:       str,
    library_name: str,
    file_dir:     str,
    n_sample:     int   = 5000,
    ph:           float = 7.4,
    username:     str   = "",
) -> dict:
    """
    Prepare a compound library for screening. Returns a job_id immediately.

    source:       path to a CSV file with 'smiles' and 'id' columns,
                  or a built-in alias: 'enamine_real', 'drugbank', 'zinc'.
    library_name: name to store the prepared library under in data/libraries/.
    n_sample:     randomly sample this many compounds from source (0 = use all).
    ph:           protonation pH for ligand preparation.
    file_dir:     output directory for temporary files.
    """
    job_id     = _new_job("prepare_library", username=username)
    stop_event = get_stop_event(job_id)

    def _run():
        try:
            result = pipeline.prepare_library(
                source       = source,
                library_name = library_name,
                n_sample     = n_sample,
                ph           = ph,
                file_dir     = file_dir,
                job_id       = job_id,
                stop_event   = stop_event,
                progress_cb  = lambda msg: _update_job(job_id, progress=msg),
            )
            _update_job(
                job_id,
                status   = "done",
                progress = f"done: {result['n_prepared']} prepared, {result['n_errors']} errors",
                result   = {**result, "next_steps": ["run_screening(protein_df_file, library_name)"]},
            )
        except Exception as e:
            _update_job(job_id, status="failed", error=str(e))

    threading.Thread(target=_run, daemon=True).start()

    return {
        "message":    f"Library preparation started for '{library_name}'.",
        "job_id":     job_id,
        "status":     "running",
        "next_steps": ["check_status(job_id)"],
    }


# ── Tool 3: run_screening (async) ─────────────────────────────────────────────

@mcp.tool()
def run_screening(
    protein_df_file: str,
    library_name:    str,
    file_dir:        str,
    project_name:    str  = "my_screening",
    docking_method:  str  = "smina",
    exhaustiveness:  str  = "16",
    use_docker:      bool = False,  # claude: set True when running without native binaries (e.g. on Mac)
    username:        str  = "",
) -> dict:
    """
    Run virtual screening of a compound library against a prepared protein.
    Returns a job_id immediately; use check_status() to track progress.

    protein_df_file: path to the protein CSV produced by prepare_protein().
    library_name:    name of a prepared library (see list_libraries()).
    project_name:    name for this screening run (results saved under file_dir/project_name/).
    docking_method:  smina | gnina | vina.
    exhaustiveness:  search exhaustiveness (higher = slower but more thorough).
    """
    job_id     = _new_job("run_screening", username=username)
    stop_event = get_stop_event(job_id)  # claude

    def _run():
        try:
            result = pipeline.run_virtual_screening(
                protein_df_file = protein_df_file,
                library_name    = library_name,
                project_name    = project_name,
                docking_method  = docking_method,
                exhaustiveness  = exhaustiveness,
                use_docker      = use_docker,  # claude
                stop_event      = stop_event,  # claude
                file_dir        = file_dir,
                progress_cb     = lambda msg: _update_job(job_id, progress=msg),
            )
            _update_job(job_id, status="done",
                        progress=f"done: screening complete", result=result)
        except Exception as e:
            _update_job(job_id, status="failed", error=str(e))

    threading.Thread(target=_run, daemon=True).start()

    return {
        "message":    f"Screening started: '{library_name}' × '{project_name}'.",
        "job_id":     job_id,
        "status":     "running",
        "next_steps": ["check_status(job_id)"],
    }


# ── Tool 5: get_top_hits ──────────────────────────────────────────────────────

@mcp.tool()
def get_top_hits(
    project_name:    str,
    file_dir:        str,
    docking_method:  str   = "smina",
    top_n:           int   = 20,
    score_threshold: float = None,
) -> dict:
    """
    Rank and return top hits from a completed screening run.

    project_name:    the project_name used in run_screening().
    docking_method:  must match the method used in run_screening().
    top_n:           number of top hits to return.
    score_threshold: optional upper bound on affinity (kcal/mol); e.g. -7.0
                     means only return compounds with affinity ≤ -7.0.
    """
    result_dir = os.path.join(file_dir, project_name)

    if not os.path.exists(result_dir):
        return {"message": f"Result directory not found: '{result_dir}'"}

    query_csv = os.path.join(file_dir, f"{project_name}_docking_query_{docking_method}.csv")
    if not os.path.exists(query_csv):
        return {"message": f"Query table not found: '{query_csv}'. Was run_screening() completed?"}

    df_query = pd.read_csv(query_csv)
    records  = []

    for _, row in df_query.iterrows():
        pdbqt_path = os.path.join(
            result_dir,
            f"dock_{docking_method}_{row['protein_name']}_{row['ligand_name']}.pdbqt"
        )
        if not os.path.exists(pdbqt_path):
            continue

        best_score = None
        with open(pdbqt_path) as f:
            for line in f:
                if "minimizedAffinity" in line or "VINA RESULT" in line:
                    try:
                        parts      = line.strip().split()
                        best_score = float(parts[1] if "VINA" in line else parts[-1])
                        break
                    except Exception:
                        pass

        if best_score is not None:
            records.append({
                "compound":      row["ligand_name"],
                "best_affinity": best_score,
                "smiles":        row.get("ligand_smiles", ""),
                "pdbqt_file":    pdbqt_path,
            })

    if not records:
        return {"message": "No docking results found. Is the screening finished?"}

    df = pd.DataFrame(records).sort_values("best_affinity")

    if score_threshold is not None:
        df = df[df["best_affinity"] <= score_threshold]

    top      = df.head(top_n)
    hits_csv = os.path.join(result_dir, f"top_{top_n}_hits.csv")
    top.to_csv(hits_csv, index=False)

    return {
        "message":       f"Top {len(top)} hits from '{project_name}' (out of {len(records)} docked).",
        "n_docked":      len(records),
        "n_hits":        len(top),
        "best_affinity": float(top["best_affinity"].iloc[0]) if len(top) else None,
        "hits":          top.to_dict(orient="records"),
        "hits_csv":      hits_csv,
    }


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    mcp.run(transport="stdio")
