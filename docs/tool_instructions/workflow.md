# Workflow Tool Instructions

A **workflow** is a saved sequence of tool calls that can be re-run on demand.

- Call `save_workflow` to persist the workflow.
- Call `run_workflow` to execute it — runs all steps automatically with no per-step confirmation.
- Call `list_workflows` to show the user their previously saved workflows (name, description, step count, creation date). 

- The user can click Refresh in the sidebar to monitor progress, or call `check_status(job_id)` at any time.

### save_workflow
Saves a named, reusable multi-step workflow to disk. Steps are defined as a list of tool calls; outputs of one step can be passed as inputs to the next using `${stepId.key}` placeholders. Call this once when the user describes a workflow they want to repeat.

Parameters:
- `name`: filesystem-safe identifier (e.g. `'protein_screening'`)
- `description`: human-readable description of what the workflow does
- `steps`: list of step dicts, each with:
  - `id` — unique identifier (e.g. `'step1'`)
  - `tool` — function name (e.g. `'prepare_protein'`, `'prepare_library'`, `'run_virtual_screening'`)
  - `args` — dict of arguments; use `${stepId.key}` to wire a prior step's output as input

Example step wiring:
```json
{"id": "step3", "tool": "run_virtual_screening",
 "args": {"protein_df_file": "${step1.protein_df_file}",
          "library_dir":     "${step2.library_dir}"}}
```

### list_workflows
Lists all saved workflows with their name, description, step count, and creation date.

### delete_workflow
Deletes a saved workflow by name.

Parameters:
- `name`: workflow name to delete (same identifier used in `save_workflow`)

### run_workflow
Runs a saved workflow asynchronously, executing all steps in sequence without requiring user confirmation. Returns a `job_id` immediately. Progress is written per-user to `session/<username>/workflow_progress.json`; the user can click Refresh in the sidebar to see the current step. Use `check_status(job_id)` for the final result, or `cancel_job(job_id)` to stop between steps.

Parameters:
- `workflow_name`: name of the workflow to run (as saved by `save_workflow`)
- `username`: the logged-in user's name — use the value from `Current user:` in your context
- `file_dir`: base output directory for generated files

## Key rule

Workflow step `"tool"` values are **pipeline function names**, not MCP tool names. The mapping is not always 1-to-1 — see the table below.

| MCP tool name | Pipeline function name (use this in workflows) |
|---|---|
| `prepare_library` | `prepare_library` ✓ (same) |
| `run_screening` | `run_virtual_screening` ← different |
| `search_literature` | `search_literature` ✓ (same) |
| `extract_information` | `extract_information` ✓ (same) |
| `get_pdb` | `get_pdb` ✓ (same) |
| `prepare_protein` | `prepare_protein` ✓ (same) |

> **Note:** the table above and the function list below cover the docking/screening/literature
> pipeline specifically — they are not an exhaustive whitelist. Step resolution falls through to
> *any* function in the `biochem` package if it isn't found in the screening or literature
> pipelines first (see `_get_callable` in `pipeline.py`). So functions like
> `extract_pdb_components`, `get_protein_ligand_interaction`, `interaction_plot`, and `obabel`
> are already valid workflow steps too — just use the biochem function name directly as `"tool"`.

---

## Available workflow step functions

### `get_pdb`

Download a protein PDB file from RCSB.

```
Parameters:
  pdb_id     str  — RCSB PDB identifier (e.g. "6N1K")

Output keys:
  saved_file  str  — path to the downloaded .pdb file
```

---

### `prepare_protein`

Prepare a protein structure for docking (protonation, PDBQT conversion, binding box).

```
Parameters:
  pdb_file        str  — path to the input PDB file
  protein_name    str  — short identifier (e.g. "PAH")
  project_name    str  — name for this run (used in output filenames)
  docking_method  str  — "smina" | "gnina" | "vina" | "diffdock"
  partial_charge  str  — "gasteiger" (default) | "mmff94"
  ph              float|None — protonation pH; None = skip protonation step

Output keys:
  protein_df_file  str  — path to the protein query CSV (pass to run_virtual_screening)
```

---

### `prepare_library`

Protonate and convert a compound library to PDBQT format.

```
Parameters:
  source        str    — built-in alias: 'drugbank' | 'enamine_real' | 'zinc'
                         or a path to a CSV file with 'smiles' and 'id' columns
  library_name  str    — name for this library (e.g. "drugbank20")
  n_sample      int    — number of compounds to randomly sample; 0 = use all (default: 5000)
  ph            float  — protonation pH (default: 7.4)

Output keys:
  library_name    str  — same as input (use this in run_virtual_screening)
  library_dir     str  — full path to the prepared library directory
  df_ligand_file  str  — path to df_ligand.csv inside the library
```

---

### `run_virtual_screening`

Dock all ligands in a prepared library against a target protein.

```
Parameters:
  protein_df_file  str   — path to protein CSV from prepare_protein (use ${stepN.protein_df_file})
  library_name     str   — name of the prepared library (use ${stepN.library_name})
  project_name     str   — name for this screening run
  docking_method   str   — "smina" | "gnina" | "vina" (default: "smina")
  exhaustiveness   str   — docking exhaustiveness (default: "16")
  file_dir         str   — base output directory (default: "files")

Output keys:
  project_name  str  — same as input (pass to get_top_hits after workflow completes)
  result_dir    str  — path to docking result directory
```

---

### `search_literature`

Search PubMed and/or Semantic Scholar for papers.

```
Parameters:
  query        str  — free-text search query
  source       str  — "pubmed" | "semantic_scholar" | "both" (default: "both")
  max_results  int  — max papers per source (default: 20)
  file_dir     str  — directory to save output JSON (default: "files")

Output keys:
  papers_file  str  — path to the JSON file with titles and abstracts
```

---

### `extract_information`

For each paper in a `papers_file`, ask an LLM a question and save a structured CSV.

```
Parameters:
  papers_file  str  — path to JSON from search_literature (use ${stepN.papers_file})
  question     str  — the question to answer from each abstract
  file_dir     str  — directory to save output CSV (default: "files")

Output keys:
  csv_file  str  — path to the output CSV with answer, evidence, confidence columns
```

---

## What CANNOT be a workflow step

| Tool | Reason |
|---|---|
| `get_top_hits` | Only in the MCP server, not in the pipeline. Call it interactively after `run_workflow` completes. |
| `render_structures` | Visualization — must be called interactively, not automated. |
| `run_af3`, `run_esm3` | No pipeline layer; these are direct biochem calls that the workflow executor cannot inject progress into usefully. Call them as separate async jobs. |

---

## Placeholder syntax

Use `${stepId.key}` to pass the output of one step as the input of the next.

- `${step1.saved_file}` — the `saved_file` key from the result of the step with id `"step1"`
- `${step2.protein_df_file}` — the `protein_df_file` key from step2's result
- Unresolved placeholders (wrong step ID or missing key) are passed through as-is, which will cause the step to fail with a clear error.

---

## Full example

Workflow: download a protein, prepare it, prepare a compound library, run screening.

```json
{
  "name": "pah_screening",
  "description": "Prepare PAH (6N1K), prepare a 100-compound library, screen with smina.",
  "steps": [
    {
      "id": "step1",
      "tool": "get_pdb",
      "args": {"pdb_id": "6N1K"}
    },
    {
      "id": "step2",
      "tool": "prepare_protein",
      "args": {
        "pdb_file":       "${step1.saved_file}",
        "protein_name":   "PAH",
        "project_name":   "pah_screening",
        "docking_method": "smina"
      }
    },
    {
      "id": "step3",
      "tool": "prepare_library",
      "args": {
        "source":       "drugbank",
        "library_name": "drugbank100",
        "n_sample":     100,
        "ph":           7.4
      }
    },
    {
      "id": "step4",
      "tool": "run_virtual_screening",
      "args": {
        "protein_df_file": "${step2.protein_df_file}",
        "library_name":    "${step3.library_name}",
        "project_name":    "pah_screening",
        "docking_method":  "smina",
        "exhaustiveness":  "16"
      }
    }
  ]
}
```

After `run_workflow` completes, call `get_top_hits(project_name="pah_screening")` interactively to retrieve ranked hits.
