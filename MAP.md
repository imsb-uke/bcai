# BioChemAIgent — Project Map

> Last updated: 2026-05-26

---

## 1. Directory Tree

```
bcai/                               ← repo root
│
├── main/                           ← Mother agent (entry point)
│   ├── server.py                   ★ MCP server — registers all 34 tools
│   ├── client.py                   Async LLM client (OpenAI / OpenRouter / Ollama)
│   ├── BioChemAIgent_chat.py       CLI chat interface
│   ├── BioChemAIgent_chat_UI.py    Streamlit chat interface
│   ├── users.csv                   User registry (username / password / quota)
│   ├── client.log                  Runtime log
│   └── src/
│       ├── user.py                 Login, registration, question quota
│       ├── history.py              Conversation save/load
│       └── utils_client.py         Formatting helpers
│
├── subagents/                      ← Multi-agent async wrappers
│   ├── jobs.py                     ★ Shared job registry (_new_job, _update_job, check_status)
│   │
│   ├── screening/                  Virtual screening sub-agent
│   │   ├── pipeline.py             Sync: prepare_library(), run_virtual_screening()
│   │   └── server.py               Async MCP wrappers (list_libraries, prepare_library,
│   │                                                    run_screening, get_top_hits)
│   ├── structure/                  Structure prediction sub-agent
│   │   └── server.py               Async MCP wrappers for bc.run_af3 / bc.run_esm3
│   │                               (no pipeline.py — biochem functions are already pipelines)
│   └── literature/                 Literature search sub-agent
│       ├── pipeline.py             Sync: search_literature(), extract_information()
│       └── server.py               Async MCP wrappers
│
├── data/
│   └── libraries/                  Prepared compound libraries for virtual screening
│       ├── _sources/               Raw input CSVs (id, name, smiles)
│       │   ├── drugbank.csv        2 726 approved drugs
│       │   └── enamine_real.csv    5 000 Enamine REAL compounds
│       ├── drugbank/
│       │   ├── df_ligand.csv       2 726 rows — absolute paths to .pdbqt files
│       │   └── ligands/            8 178 .pdbqt files (protonated, optimised)
│       └── enamine_real/
│           ├── df_ligand.csv       4 976 rows — absolute paths to .pdbqt files
│           └── ligands/            14 931 .pdbqt files
│
├── docs/
│   ├── system.md                   System prompt for the Mother agent
│   ├── drug_discovery_protocol.md  Step-by-step DD protocol surfaced as a tool
│   ├── collaborative_people.md     Contact list surfaced as a tool
│   └── tool_instructions/         Per-tool usage guides (render_structures, run_af3, …)
│
├── config/                         MCP client configuration files (copy into Claude/Cursor)
│   ├── BioChemAIgent_MCP_Server.json
│   ├── ChEMBL_MCP_Server.json
│   ├── ChEMBL_MCP_Server_docker.json
│   ├── PDB-MCP-Server.json
│   └── PDB-MCP-Server_docker.json
│
├── mcp_external/                   Third-party MCP servers (bundled)
│   ├── ChEMBL-MCP-Server/          Node.js; query ChEMBL compound/bioactivity DB
│   │   ├── src/index.ts
│   │   └── build/                  Compiled JS (ready to run)
│   └── PDB-MCP-Server/             Node.js; fetch PDB structures
│       ├── src/index.ts
│       └── build/
│
├── software/                       Notes on external structure prediction software
│   ├── AlphaFold3.md               AF3 API setup guide
│   ├── DiffDock.md                 DiffDock setup guide
│   └── README.md
│
├── notes/
│   └── session-context.md          Running dev notes (architecture, decisions, todo)
│
├── old_versions/
│   └── main_v1_biochem/            v1 server (single-agent, tools inline) — archived
│
├── env                             ⚠ API keys & config (git-ignored)
├── llm_models                      Available LLM model list
├── environment.yml                 Conda env spec (bcai / bcai-dev)
├── Dockerfile                      Container build for bcai server
├── docker-compose.yml              Orchestrates bcai + external MCPs
├── entrypoint.sh                   Docker entrypoint
├── LICENSE
└── README.md

── sibling repo ──────────────────────────────────────────────────────────────
../biochem/                         ← independent package (editable install)
    biochem/
    ├── __init__.py                 TOOLS list + all imports
    ├── get_data.py                 uniprot2seq, get_pdb, read_csv, write_text_file
    ├── resources.py                drug_discovery_protocol, collaborative_people, get_tools_doc
    ├── protein_ligand_basics.py    smiles_to_3d, obabel, admet_predict, protonate_*, …
    ├── molecular_docking.py        prepare_protein/ligand, make_query_table,
    │                               run_molecular_docking, get_protein_ligand_interaction
    ├── af3.py                      run_af3  (AlphaFold3 API)
    ├── esm3.py                     run_esm3 (ESM3 API)
    ├── render_structures.py        render_structures (3D visualisation)
    ├── interaction_plot.py         interaction_plot
    ├── dock_utilities.py           Docking utility functions
    └── generate_graph.py           Graph generation helpers
```

---

## 3. Tool Catalog — All 34 Registered Tools

### 3a. Biochem Tools (sync, registered directly from `biochem.TOOLS`)

| # | Tool | Source file | What it does |
|---|---|---|---|
| 1 | `uniprot2seq` | `get_data.py` | Fetch protein sequence from UniProt ID |
| 2 | `get_pdb` | `get_data.py` | Download PDB structure by PDB ID |
| 3 | `read_csv` | `get_data.py` | Read a CSV file into structured output |
| 4 | `write_text_file` | `get_data.py` | Write text/data to a file |
| 5 | `drug_discovery_protocol` | `resources.py` | Return the DD protocol as a reference doc |
| 6 | `collaborative_people` | `resources.py` | Return contact list for collaboration |
| 7 | `get_tools_doc` | `resources.py` | Return detailed documentation for any tool |
| 8 | `smiles_to_3d` | `protein_ligand_basics.py` | Convert SMILES string to 3D SDF/PDB |
| 9 | `calculate_ligand_info` | `protein_ligand_basics.py` | MW, logP, HBD/HBA, TPSA, … |
| 10 | `smiles_pattern_search` | `protein_ligand_basics.py` | Substructure search on SMILES |
| 11 | `extract_pdb_components` | `protein_ligand_basics.py` | Split PDB into protein/ligand/water |
| 12 | `obabel` | `protein_ligand_basics.py` | Run OpenBabel format conversions |
| 13 | `extract_pdb_info` | `protein_ligand_basics.py` | Extract chain/residue metadata from PDB |
| 14 | `protonate_and_optimize_ligand` | `protein_ligand_basics.py` | Protonate + minimise ligand at pH 7.4 |
| 15 | `protonate_and_optimize_protein` | `protein_ligand_basics.py` | Protonate protein (Propka / Obabel) |
| 16 | `cxsmiles2smiles` | `protein_ligand_basics.py` | Convert ChemDraw CXSMILES → canonical SMILES |
| 17 | `stereoisomers` | `protein_ligand_basics.py` | Enumerate all stereoisomers for a SMILES |
| 18 | `admet_predict` | `protein_ligand_basics.py` | Predict ADMET properties |
| 19 | `prepare_protein` | `molecular_docking.py` | Convert protein PDB → PDBQT for docking |
| 20 | `prepare_ligand` | `molecular_docking.py` | Convert ligand → PDBQT for docking |
| 21 | `make_query_table` | `molecular_docking.py` | Build protein × ligand docking query table |
| 22 | `run_molecular_docking` | `molecular_docking.py` | Run Smina/Vina docking (single pair) |
| 23 | `get_protein_ligand_interaction` | `molecular_docking.py` | Parse docking output → interaction table |
| 24 | `render_structures` | `render_structures.py` | 3D visualisation (py3Dmol / NGLview) |
| 25 | `interaction_plot` | `interaction_plot.py` | 2D interaction diagram |

> `run_af3` and `run_esm3` are in `biochem.TOOLS` but **skipped** in the sync loop
> (`_ASYNC_TOOLS = {"run_af3", "run_esm3"}`) — the async wrappers below replace them.

---

### 3b. Structure Sub-Agent Tools (async — `subagents/structure/server.py`)

| # | Tool | What it does |
|---|---|---|
| 26 | `run_af3` | AlphaFold3 structure prediction — wraps `bc.run_af3`; returns `job_id` immediately |
| 27 | `run_esm3` | ESM3 protein structure / function prediction — wraps `bc.run_esm3`; async |

---

### 3c. Screening Sub-Agent Tools (`subagents/screening/server.py`)

| # | Tool | Sync/Async | What it does |
|---|---|---|---|
| 28 | `list_libraries` | sync | List prepared compound libraries in `data/libraries/` |
| 29 | `prepare_library` | **async** | Load raw CSV → protonate + convert all compounds → `df_ligand.csv` + `.pdbqt` files |
| 30 | `run_screening` | **async** | Build query table then run virtual screening over a full library |
| 31 | `get_top_hits` | sync | Return ranked top hits from a finished screening run |

---

### 3d. Literature Sub-Agent Tools (`subagents/literature/server.py`)

| # | Tool | Sync/Async | What it does |
|---|---|---|---|
| 32 | `search_literature` | **async** | Query PubMed + Semantic Scholar; deduplicate; save JSON of papers + abstracts |
| 33 | `extract_information` | **async** | Per-abstract LLM extraction; save CSV with `answer / evidence / confidence` |

---

### 3e. Shared Infrastructure

| # | Tool | Source | What it does |
|---|---|---|---|
| 34 | `check_status` | `subagents/jobs.py` | Poll any running job (screening, structure, literature) by `job_id` |

---

## 4. Async Job Pattern

All long-running tools return immediately with a `job_id`. The LLM polls with `check_status`.

```
  LLM calls run_screening(...)
        │
        ▼
  server.py: job_id = _new_job("run_screening")
             Thread(target=_run).start()
             return {"job_id": "a3f7b2c1", "status": "running"}
        │
        │   (background thread running pipeline.run_virtual_screening)
        │
  LLM calls check_status("a3f7b2c1")
        │
        ▼
  jobs.py: returns {"status": "running", "progress": "docking compound 412/4976"}
        │
        │   (…later…)
        │
  LLM calls check_status("a3f7b2c1")
        │
        ▼
  jobs.py: returns {"status": "done", "result": {"csv_file": "...", "n_hits": 4976}}
```

**Job lifecycle states:** `running` → `done` | `failed`

**Shared registry:** `_jobs` dict in `subagents/jobs.py` with `threading.Lock`.
All sub-agents import `_new_job` and `_update_job` from there — one `check_status` covers everything.

---

## 5. UI ↔ Client Communication

The Streamlit UI and `client.py` are **two separate processes** running in parallel. They communicate via the filesystem — no sockets, no shared memory.

```
  main/
  ├── work_space/          ← handoff zone: UI drops files here, client picks them up
  │   └── <user>.pkl       history file with latest user message (temporary)
  ├── session/
  │   └── <user>/
  │       └── chat_session.json   model name, API keys, history file path
  └── history/
      └── <user>/
          └── chat_history.pkl   permanent history; UI reads from here
```

### Step-by-step

```
  USER types a question
        │
        ▼
  BioChemAIgent_chat_UI.py
    1. Appends user message → history/<user>/chat_history.pkl
    2. Copies it           → work_space/<user>.pkl
    3. Writes session info → session/<user>/chat_session.json
    4. Sets waiting_for_response = True
    5. Polls history/<user>/chat_history.pkl for mtime change (every 0.1s)
        │
        │  (file appears in work_space/)
        │
        ▼
  client.py  (while True loop, sleep 0.05s)
    1. Detects work_space/<user>.pkl
    2. Reads session/<user>/chat_session.json  → model, API keys
    3. Calls LLM with full history + MCP tools
    4. Saves updated history (with assistant reply) → work_space/<user>.pkl
    5. Moves work_space/<user>.pkl → history/<user>/chat_history.pkl
        │
        │  (mtime of history file changes)
        │
        ▼
  BioChemAIgent_chat_UI.py
    6. Detects mtime change → sets waiting_for_response = False → st.rerun()
    7. Reloads history and displays assistant response
```

### Key points

- **One `client.py` serves all users** — each user's work file is `<username>.pkl`, so multiple users are naturally isolated
- **Session JSON carries the API keys** — the UI passes the user's tokens to the client at query time; the client sets them as env vars before calling the LLM
- **`mv` not `cp` on the way back** — moving the file atomically updates the mtime, which is what the UI polls for
- **Async job notification** — when a tool returns a `job_id`, `client.py` registers it in `_pending_jobs` and polls `check_status` every loop iteration. When done, it writes `session/<user>/job_ready.json`. The UI watches that file and shows a `st.toast()` notification automatically — no manual status check needed.

---

## 6. Data Flow — End-to-End Drug Discovery Workflow

```
 User question
      │
      ▼
 client.py  ──► OpenAI / OpenRouter / Ollama  (async LLM)
                        │  tool calls (MCP)
                        ▼
              main/server.py  (Mother MCP server)
                        │
         ┌──────────────┼──────────────────────────────────┐
         │              │                                   │
         ▼              ▼                                   ▼
  biochem tools    sub-agent tools                  external MCP servers
  (sync, instant)  (async, long-running)            (ChEMBL, PDB — Node.js)
         │              │
         │        ┌─────┴──────────────────────────┐
         │        │ screening    structure  literature
         │        │ pipeline.py  (direct bc) pipeline.py
         │        └────────────────────────────────┘
         │              │
         └──────────────┘
                  │
                  ▼
           data/libraries/    files/    PubMed / S2 / AF3 / ESM3 APIs
```

---

## 6. Compound Libraries

| Library | Source | Compounds (raw) | Prepared (.pdbqt) | Ligand files |
|---|---|---|---|---|
| `drugbank` | DrugBank approved drugs | 2 726 | 2 726 | 8 178 |
| `enamine_real` | Enamine REAL subset (5 HAC bins) | 5 000 | 4 976 | 14 931 |

**File layout per library:**
```
data/libraries/<name>/
    df_ligand.csv       compound table: id, name, smiles, pdbqt_path (absolute)
    ligands/            .pdbqt files (one per stereoisomer / protonation state)

data/libraries/_sources/
    drugbank.csv        raw: id, name, smiles
    enamine_real.csv    raw: merged from 5 HAC-bin CSVs
```

---

## 7. External MCP Servers

These run as **separate processes** alongside the main server. Configure in Claude/Cursor via `config/`.

| Server | Language | Transport | What it provides |
|---|---|---|---|
| **ChEMBL-MCP-Server** | Node.js (TypeScript) | stdio | ChEMBL compound search, bioactivity data, target queries |
| **PDB-MCP-Server** | Node.js (TypeScript) | stdio | PDB structure search and metadata |

Docker images available via `docker-compose.yml` and per-server `Dockerfile`.

```
config/
    BioChemAIgent_MCP_Server.json   ← point Claude at main/server.py
    ChEMBL_MCP_Server.json          ← point Claude at ChEMBL Node.js server
    PDB-MCP-Server.json             ← point Claude at PDB Node.js server
    *_docker.json                   ← Docker variants of the above
```

---

## 8. Front-End Interfaces

| File | Type | Description |
|---|---|---|
| `main/BioChemAIgent_chat.py` | CLI | Terminal chat loop; uses `client.py` |
| `main/BioChemAIgent_chat_UI.py` | Streamlit | Web UI with login, quota tracking, session history |

**User management:** `main/users.csv` — `username, password, free_questions, n_questions`
**Session history:** saved/loaded via `main/src/history.py`

---

## 9. Configuration & Environment

### `bcai/env` (git-ignored — copy and fill in)

```ini
# LLM
LLM_MODEL          = 'openrouter|openai/gpt-4o'
LLM_MODEL_LITERATURE = 'openrouter|openai/gpt-4o-mini'

# API keys
OPENAI_API_KEY     = sk-...
OPENROUTER_API_KEY = sk-or-...
OLLAMA_API_KEY     = ...
ESM3_API_KEY       = ...
AF3_API_KEY        = ...

# Paths
FILE_DIR           = /absolute/path/to/working/files
SYSTEM_FILE_DIR    = /path/to/bcai/docs
USER_FILE          = /path/to/bcai/main/users.csv
```

### LLM provider format

`LLM_MODEL = 'provider|model_name'`

| Provider string | Client used | Notes |
|---|---|---|
| `openai` | `AsyncOpenAI()` | Uses `OPENAI_API_KEY` |
| `ollama_local` | `AsyncOpenAI(base_url="localhost:11434/v1")` | Local Ollama |
| `ollama_cloud` | `ollama.Client(host="https://ollama.com")` | Cloud Ollama |
| `openrouter` | `AsyncOpenAI(base_url="https://openrouter.ai/api/v1")` | Uses `OPENROUTER_API_KEY` |

### Conda environments

| Env | Use |
|---|---|
| `bcai-dev` | Active development — editable install of `biochem` (`pip install -e ../biochem`) |
| `bcai` | Stable snapshot / production |
| `biochem` | Standalone `biochem` package work |

---

## 10. Key File Relationships

```
env  ──────────────────────────────────► all servers and pipelines (dotenv)
biochem (pip install -e) ─────────────► main/server.py, subagents/*/server.py
subagents/jobs.py ────────────────────► subagents/screening/server.py
                                      ► subagents/structure/server.py
                                      ► subagents/literature/server.py
                                      ► main/server.py  (registers check_status)
subagents/screening/pipeline.py ──────► subagents/screening/server.py
subagents/literature/pipeline.py ─────► subagents/literature/server.py
data/libraries/ ──────────────────────► subagents/screening/server.py (list_libraries, run_screening)
docs/system.md ───────────────────────► client.py (SYSTEM_FILE_DIR)
```

---

## 11. Pending / Backlog

| Item | Status | Notes |
|---|---|---|
| AF3 path validation | ⏳ untested | Full agent → `run_af3` → file output path needs e2e test |
| `client.py` logging cleanup | ⚠️ todo | Replace `os.system("echo ... >> client.log")` with Python `logging` module |
| User-defined workflows | 💡 planned | User defines a pipeline once; agent stores and re-runs it |
| ChEMBL in literature sub-agent | 💡 planned | Add ChEMBL as a third source in `search_literature` |
| End-to-end screening test | ⏳ next | Run `run_screening` with a real protein against `drugbank` library |
