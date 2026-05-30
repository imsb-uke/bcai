## Plot with render_structures

`render_structures` gives you full control over models, chains, residues, elements, HETATM, and more. It makes smart defaults (cartoon for PDB/CIF; sticks for SDF/MOL2) but you can override everything with explicit style rules. 

Help
Render multiple molecular files with fine-grained control.

Parameters
----------
files : list of paths (PDB/CIF/SDF/MOL2). Order becomes model indices [0..N-1].
style_rules : list of {"select": <Sel>, "style": <Sty>} applied in order.
    - Example selections:
        {"model": 0}                         # first model
        {"model": [1,2]}                     # multiple models
        {"chain": "A"}                       # chain A (in PDB/CIF models)
        {"resi": [10,11,120]}                # residue numbers
        {"elem": "C"}                        # element carbon
        {"hetflag": True}                    # HETATM only (ligands, ions, etc.) in PDB/CIF
        {"invert": True, "hetflag": True}    # everything *except* HETATM (protein/nucleic)
    - Example styles:
        {"cartoon": {"color":"#2ca02c", "opacity":0.5}}
        {"stick": {"radius":0.25, "colorscheme":"cyanCarbon"}}
        {"sphere": {"radius":1.5}}
        {"line": {}}
        {"cross": {}}
surface_rules : list of {"select": <Sel>, "surface": {"opacity":0.2, "color":"#aaa"}, "type": py3Dmol.MS}
    - type: py3Dmol.VDW (default), py3Dmol.MS (molecular), py3Dmol.SAS (solvent-accessible)
label_rules : list of {"text": "Active site", "position": {"x":..,"y":..,"z":..}, "style": {...}}
    - or {"text":"Residue 197", "select":{"model":0,"chain":"A","resi":197}}
chain_color_map : optional mapping to auto-add per-chain cartoon for polymeric models.
    - Keys can be model index (int) or "all".
background : backgroundall color.
width: width.
height: height.
zoom : auto zoom to all atoms.
file_name : file_name.

How to use
1) Minimal (auto styles)
```
render_structures(
    files=[
        "files/protein.pdb",                # polymeric -> cartoon + sticks for HETATM
        "files/ligand.sdf" # ligand -> sticks
    ],
)
```
2) Full control with explicit rules
```
rules = [
    # Protein A in model 0 as semi-transparent blue cartoon
    {"select": {"model":0, "chain":"A", "hetflag":False},
     "style":  {"cartoon":{"color":"#1f77b4","opacity":0.6}}},

    # Protein B in model 0 as red cartoon + backbone line overlay
    {"select": {"model":0, "chain":"B", "hetflag":False},
     "style":  {"cartoon":{"color":"#d62728","opacity":0.6}, "line": {}}},

    # All ligands (HETATM) in model 0 as thicker sticks
    {"select": {"model":0, "hetflag":True},
     "style":  {"stick":{"radius":0.35}}},

    # External ligand (model 1) cyan sticks
    {"select": {"model":1},
     "style":  {"stick":{"radius":0.3, "colorscheme":"cyanCarbon"}}},

    # Highlight specific residues on chain A
    {"select": {"model":0, "chain":"A", "resi":[123,197,222]},
     "style":  {"stick":{"radius":0.4}}},
]
surfaces = [
    {"select":{"model":0, "chain":"A", "hetflag":False},
     "surface":{"opacity":0.15}, "type": py3Dmol.MS}
]
render_structures(
    files=["files/protein.pdb","files/ligand.sdf"],
    style_rules=rules,
    surface_rules=surfaces,
    background="white",
)
```
3) Quick per-chain coloring (auto-adds cartoon per chain)
```
render_structures(
    files=["files/protein_1.pdb","files/protein_2.pdb","files/ligand_1.sdf"],
    chain_color_map={
        0: {"A":"#3cb44b","B":"#e6194B"},  # model 0 chains
        1: {"A":"#9467bd","_":"#aaaaaa"},  # model 1 chains (fallback "_" if chain ids missing)
        "all": {"_":"#cccccc"}             # optional global fallback
    },
)
# Keys can be model index integer number (e.g. 0, 1, 2) or a string "all".
```
Tips:
- You never have to declare “this is protein/ligand”. Just select what you want (by model, hetflag, chain, resi, elem, etc.) and style it.
- Use {"invert": True, "hetflag": True} to style everything except ligands.
- Put hetflag and the oxygen-sphere rule last (so it overrides the broad HETATM sticks).
- Want a clean reset? Call with style_rules=[] to avoid defaults and only draw what you specify.


