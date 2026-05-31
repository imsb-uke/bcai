"""
Screening pipeline — synchronous orchestration functions.

Composes biochem's core functions into higher-level screening workflows.
No MCP, no threads from the caller's perspective — usable standalone or
imported by the screening sub-agent server.

progress_cb: optional callable(message: str) for progress reporting.
             The sub-agent passes lambda msg: _update_job(job_id, progress=msg).
             Standalone callers can pass None or their own callback.
"""

import os
import time
import shutil  # claude
import threading

import pandas as pd
import biochem as bc

# claude ---
LIBRARY_DIR = os.getenv("LIBRARY_DIR")
_BUILTIN = {
    "enamine_real": os.path.join(LIBRARY_DIR, "_sources/enamine_real.csv"),
    "drugbank":     os.path.join(LIBRARY_DIR, "_sources/drugbank.csv"),
    "zinc":         os.path.join(LIBRARY_DIR, "_sources/zinc.csv"),
}
# ---


# ── prepare_library ───────────────────────────────────────────────────────────

def prepare_library(
    source:       str,
    library_name: str,
    n_sample:     int   = 5000,
    ph:           float = 7.4,
    file_dir:     str,
    job_id:       str   = None,  # claude: used to name the per-job tmp dir
    stop_event           = None,  # claude: threading.Event; set to cancel between compounds
    progress_cb          = None,
) -> dict:
    """
    Prepare a compound library for virtual screening.

    Reads a CSV with 'smiles' and 'id' columns, protonates and converts each
    compound to .pdbqt, and saves everything under library_dir/.
    Supports resume: compounds whose .pdbqt already exists are skipped.

    Args:
        source:       built-in alias ('drugbank', 'enamine_real', 'zinc') or
                      path to a CSV with 'smiles' and 'id' columns.
        library_name: name for this library; stored under LIBRARY_DIR/library_name/.
        n_sample:     randomly sample this many compounds (0 = use all).
        ph:           protonation pH for ligand preparation.
        progress_cb:  optional callable(str) for progress messages.

    Returns:
        {
            "library_name":   str,
            "library_dir":    str,
            "df_ligand_file": str,   # path to df_ligand.csv
            "n_prepared":     int,
            "n_errors":       int,
        }
    """
    source_csv  = _BUILTIN.get(source, source)  # claude: resolve alias or use as path
    library_dir = os.path.join(LIBRARY_DIR, library_name)  # claude

    if not os.path.exists(source_csv):
        raise FileNotFoundError(
            f"Source not found: '{source_csv}'. "
            f"Valid aliases: {list(_BUILTIN.keys())}."
        )

    def _progress(msg):
        if progress_cb:
            progress_cb(msg)

    _progress("loading source CSV...")
    df = pd.read_csv(source_csv)

    for col in ("smiles", "id"):
        if col not in df.columns:
            raise ValueError(
                f"Source CSV must have a '{col}' column. Found: {list(df.columns)}"
            )

    if n_sample and len(df) > n_sample:
        df = df.sample(n=n_sample, random_state=42).reset_index(drop=True)

    total = len(df)

    ligand_dir = os.path.join(library_dir, "ligands")
    # Put intermediate files in a job-specific tmp dir inside file_dir so they
    # never pollute the main files/ folder. Deleted at the end of this function.
    tmp_dir = os.path.join(file_dir, f"tmp_{job_id}" if job_id else "tmp_lib")
    # ---
    for d in (library_dir, ligand_dir, tmp_dir):
        os.makedirs(d, exist_ok=True)

    df_ligand_rows = []
    errors         = []

    for i, row in df.iterrows():
        if stop_event and stop_event.is_set():  # claude
            break
        smiles = str(row["smiles"])
        name   = str(row["id"])
        _progress(f"preparing ligands: {i + 1}/{total}")

        # Resume support — skip if .pdbqt already exists
        pdbqt_path = os.path.join(ligand_dir, f"{name}.pdbqt")
        if os.path.exists(pdbqt_path):
            df_ligand_rows.append({
                "ligand_name":       name,
                "ligand_smiles":     smiles,
                "ligand_sdf_file":   os.path.join(ligand_dir, f"{name}_ph{str(ph).replace('.', '')}_min.sdf"),
                "ligand_pdb_file":   os.path.join(ligand_dir, f"{name}.pdb"),
                "ligand_pdbqt_file": pdbqt_path,
            })
            continue

        try:
            smiles_clean = bc.cxsmiles2smiles(smiles)
            out = bc.prepare_ligand(
                smiles       = smiles_clean,
                ligand_name  = name,
                project_name = "ligand",
                ph           = ph,
                file_dir     = tmp_dir,
            )

            if out.get("ligand_df_file"):
                df_tmp = pd.read_csv(out["ligand_df_file"])
                for col in ("ligand_sdf_file", "ligand_pdb_file", "ligand_pdbqt_file"):
                    src = df_tmp.at[0, col]
                    if src and os.path.exists(str(src)):
                        dst = os.path.join(ligand_dir, os.path.basename(src))
                        os.rename(src, dst)
                        df_tmp.at[0, col] = dst
                df_ligand_rows.append(df_tmp.iloc[0].to_dict())

            for f in os.listdir(tmp_dir):
                try:    os.remove(os.path.join(tmp_dir, f))
                except Exception: pass

        except Exception as e:
            errors.append({"id": name, "error": str(e)})
            for f in os.listdir(tmp_dir):
                try:    os.remove(os.path.join(tmp_dir, f))
                except Exception: pass

    df_ligand      = pd.DataFrame(df_ligand_rows)
    df_ligand_path = os.path.join(library_dir, "df_ligand.csv")
    df_ligand.to_csv(df_ligand_path, index=False)

    if errors:
        pd.DataFrame(errors).to_csv(os.path.join(library_dir, "errors.csv"), index=False)

    shutil.rmtree(tmp_dir, ignore_errors=True)  # claude: delete job tmp dir

    return {
        "library_name":   library_name,  # claude
        "library_dir":    library_dir,
        "df_ligand_file": df_ligand_path,
        "n_prepared":     len(df_ligand_rows),
        "n_errors":       len(errors),
    }


# ── run_virtual_screening ─────────────────────────────────────────────────────

def run_virtual_screening(
    protein_df_file: str,
    library_name:    str,
    project_name:    str   = "my_screening",
    docking_method:  str   = "smina",
    exhaustiveness:  str   = "16",
    use_docker:      bool  = False,  # claude
    stop_event              = None,  # claude: threading.Event; set to cancel between compounds
    file_dir:        str,
    progress_cb              = None,
) -> dict:
    """
    Run virtual screening: build a query table then dock all library ligands
    against the prepared protein.

    Uses an internal thread to run bc.run_molecular_docking() so that the
    output directory can be polled every 5 s for per-compound progress.
    From the caller's perspective this function is still synchronous (blocking).

    Args:
        protein_df_file: path to protein CSV produced by bc.prepare_protein().
        library_name:    name of a prepared library (under LIBRARY_DIR/).
        project_name:    name for this run; results saved under file_dir/project_name/.
        docking_method:  smina | gnina | vina.
        exhaustiveness:  docking exhaustiveness (higher = more thorough, slower).
        file_dir:        base output directory.
        progress_cb:     optional callable(str) for progress messages.

    Returns the bc.run_molecular_docking() result dict, augmented with
    project_name and result_dir.
    """
    library_dir = os.path.join(LIBRARY_DIR, library_name)  # claude

    def _progress(msg):
        if progress_cb:
            progress_cb(msg)

    df_ligand_file = os.path.join(library_dir, "df_ligand.csv")
    if not os.path.exists(df_ligand_file):
        raise FileNotFoundError(
            f"Library not prepared at '{library_dir}'. "
            f"Run prepare_library() first."
        )

    _progress("building query table...")
    query = bc.make_query_table(
        protein_df_file = protein_df_file,
        ligand_df_file  = df_ligand_file,
        project_name    = project_name,
        docking_method  = docking_method,
        file_dir        = file_dir,
    )

    total      = len(pd.read_csv(query["save_dir"]))
    result_dir = os.path.join(file_dir, project_name)
    _progress(f"docking: 0/{total}")

    # Run docking in a sub-thread so we can monitor output for per-compound progress
    docking_result = [None]
    docking_error  = [None]

    def _dock():
        try:
            docking_result[0] = bc.run_molecular_docking(
                query_table_dir = query["save_dir"],
                docking_method  = docking_method,
                exhaustiveness  = exhaustiveness,
                project_name    = project_name,
                use_docker      = use_docker,  # claude
                stop_event      = stop_event,  # claude
                file_dir        = file_dir,
            )
        except Exception as e:
            docking_error[0] = str(e)

    dock_thread = threading.Thread(target=_dock, daemon=True)
    dock_thread.start()

    while dock_thread.is_alive():
        if os.path.exists(result_dir):
            n_done = len([f for f in os.listdir(result_dir) if f.endswith(".pdbqt")])
            _progress(f"docking: {n_done}/{total}")
        time.sleep(5)

    dock_thread.join()

    if docking_error[0]:
        raise Exception(docking_error[0])

    return {
        **docking_result[0],
        "project_name": project_name,
        "result_dir":   result_dir,
        "next_steps":   ["get_top_hits(project_name)"],
    }
