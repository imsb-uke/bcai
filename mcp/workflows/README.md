# workflows/

Saved user-defined workflows for BCAI. Each `.json` file is a named,
reusable sequence of tool calls that can be re-run via `run_workflow`.
This directory is not tracked by git — workflows are local to each deployment.

## Workflow JSON format

```json
{
  "name": "protein_screening",
  "description": "Screen a compound library against a target protein",
  "created_at": "2026-05-31T10:00:00",
  "steps": [
    {
      "id":   "step1",
      "tool": "prepare_protein",
      "args": {"pdb_id": "1abc", "chain_id": "A"}
    },
    {
      "id":   "step2",
      "tool": "prepare_library",
      "args": {"source_csv": "data/chembl_subset.csv", "library_dir": "data/libraries/my_lib"}
    },
    {
      "id":   "step3",
      "tool": "run_virtual_screening",
      "args": {
        "protein_df_file": "${step1.protein_df_file}",
        "library_dir":     "${step2.library_dir}",
        "project_name":    "my_screening"
      }
    }
  ]
}
```

## Placeholder syntax

Use `${stepId.key}` in `args` to reference the output of a prior step.
The key must match a field in that step's result dict.

## Available tools

Any function in `biochem`, `subagents/screening/pipeline.py`, or
`subagents/literature/pipeline.py` can be used as a workflow step.
