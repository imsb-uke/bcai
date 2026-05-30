"""
Structure Prediction Sub-Agent — async MCP server for AF3 and ESM3.

Tools exposed:
  run_af3(protein_sequence, ...)   → job_id  (async, can take hours)
  run_esm3(protein_input, task, …) → job_id  (async, minutes–hours)
  check_status(job_id)             → status + progress for any job

Run from bcai/subagents/structure/:
    python server.py
"""

import os
import sys
import threading

from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP

import biochem as bc

# ── Config ────────────────────────────────────────────────────────────────────

_HERE = os.path.dirname(os.path.abspath(__file__))
load_dotenv(os.path.join(_HERE, "../../env"))

# Make bcai/ importable (needed when running standalone)
_BCAI_ROOT = os.path.dirname(os.path.dirname(_HERE))
if _BCAI_ROOT not in sys.path:
    sys.path.insert(0, _BCAI_ROOT)

from subagents.jobs import _new_job, _update_job, get_stop_event  # claude: added get_stop_event

# ── MCP server ────────────────────────────────────────────────────────────────

mcp = FastMCP(name="StructureAgent")


# ── Tool 1: run_af3 (async) ───────────────────────────────────────────────────

@mcp.tool()
def run_af3(
    protein_sequence:       list,
    n_chain_per_sequence:   list | None = None,
    uniprot_id:             list | None = None,
    ligand_sdf_dir:         str  | None = None,
    n_ligand:               int         = 1,
    project_name:           str         = "my_project",
) -> dict:
    """
    Run AlphaFold 3 structure prediction. Returns a job_id immediately;
    use check_status(job_id) to poll for completion.

    protein_sequence:     list of amino acid sequences (one per unique chain type).
    n_chain_per_sequence: number of copies of each sequence (default: [1, 1, ...]).
    uniprot_id:           optional UniProt IDs for MSA lookup.
    ligand_sdf_dir:       optional path to ligand SDF file (protein–ligand complex).
    n_ligand:             number of ligand copies.
    project_name:         name for this prediction run (used for output files).
    """
    job_id     = _new_job("run_af3")
    stop_event = get_stop_event(job_id)  # claude

    def _run():
        if stop_event and stop_event.is_set():  # claude: exit early if already cancelled
            return
        try:
            _update_job(job_id, progress="submitting AF3 job and waiting for output...")
            result = bc.run_af3(
                protein_sequence     = protein_sequence,
                n_chain_per_sequence = n_chain_per_sequence,
                uniprot_id           = uniprot_id,
                ligand_sdf_dir       = ligand_sdf_dir,
                n_ligand             = n_ligand,
                project_name         = project_name,
            )
            _update_job(job_id, status="done", progress="done", result=result)
        except Exception as e:
            _update_job(job_id, status="failed", error=str(e))

    threading.Thread(target=_run, daemon=True).start()

    return {
        "message":    f"AF3 structure prediction started for '{project_name}'.",
        "job_id":     job_id,
        "status":     "running",
        "next_steps": ["check_status(job_id)"],
    }


# ── Tool 2: run_esm3 (async) ──────────────────────────────────────────────────

@mcp.tool()
def run_esm3(
    protein_input: dict,
    task:          str,
    protein_name:  str = "my_protein",
    model_name:    str = "esm3-large-2024-03",
    file_dir:      str = "files",
) -> dict:
    """
    Run ESM3 structure or sequence prediction. Returns a job_id immediately;
    use check_status(job_id) to poll for completion.

    protein_input: dict with any of: pkl_file, pdb_file, sequence, sasa,
                   function, residue_annotations. Unused keys are ignored.
    task:          'structure' | 'sequence' | 'function' | 'sasa' |
                   'residue_annotations'
    protein_name:  name used for output files.
    model_name:    ESM3 model variant.
    file_dir:      output directory (overridden by FILE_DIR env var if set).
    """
    file_dir   = file_dir or os.getenv("FILE_DIR", "files")  # claude
    job_id     = _new_job("run_esm3")
    stop_event = get_stop_event(job_id)  # claude

    def _run():
        if stop_event and stop_event.is_set():  # claude: exit early if already cancelled
            return
        try:
            _update_job(job_id, progress=f"running ESM3 ({task})...")
            result = bc.run_esm3(
                protein_input = protein_input,
                task          = task,
                protein_name  = protein_name,
                model_name    = model_name,
                file_dir      = file_dir,
            )
            _update_job(job_id, status="done", progress="done", result=result)
        except Exception as e:
            _update_job(job_id, status="failed", error=str(e))

    threading.Thread(target=_run, daemon=True).start()

    return {
        "message":    f"ESM3 {task} prediction started for '{protein_name}'.",
        "job_id":     job_id,
        "status":     "running",
        "next_steps": ["check_status(job_id)"],
    }


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    mcp.run(transport="stdio")
