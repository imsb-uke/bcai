import os
import pathlib
import py3Dmol
import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional, Tuple, Union

FILE_DIR = os.getenv("FILE_DIR")

Sel = Dict[str, Any]   # a 3Dmol selection dict, e.g. {"model":0,"chain":"A","resi":[10,11], "elem":"C","hetflag":True}
Sty = Dict[str, Any]   # a 3Dmol style dict,   e.g. {"cartoon":{"color":"#ff0","opacity":0.6}} or {"stick":{"radius":0.25}}

def _infer_format(path: str) -> str:
    ext = pathlib.Path(path).suffix.lower().lstrip(".")
    return {"pdb":"pdb","cif":"cif","mmcif":"cif","mol2":"mol2","sdf":"sdf"}.get(ext, "pdb")

def _read_file(path: str) -> str:
    with open(path, "r") as fh:
        return fh.read()

def _list_chains_from_pdb_text(txt: str) -> List[str]:
    chains = []
    for line in txt.splitlines():
        if line.startswith(("ATOM","HETATM")) and len(line) > 21:
            c = line[21].strip() or "_"
            if c not in chains:
                chains.append(c)
    return chains

def _is_polymeric_format(fmt: str) -> bool:
    # Formats that can show cartoon/secondary structure if the content is protein/NA
    return fmt in ("pdb","cif")

def _extract_atom_coords_from_pdb_text(txt: str) -> List[tuple[float,float,float]]:
    """
    Return list of (x, y, z) for each ATOM/HETATM in file, in the same order
    3Dmol uses internally (sequential through the PDB).
    """
    coords = {}
    for line in txt.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                idx = int(line[6:11]) 
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords[idx] = (x, y, z)
    return coords
    
def _add_interaction_lines_from_csv(
    viewer: py3Dmol.view,
    atom_coords: List[tuple[float,float,float]],
    csv_path: str,
    color: str = "gray",
    radius: float = 0.15,
    dash_length: float = 0.5,
    gap_length: float = 0.3,
    itype_filter: Optional[List[str]] = None    # e.g. ["hbond"]
) -> None:
    """
    Draw dashed cylinders between ligand_idx and protein_idx from a CSV.
    Assumes ligand_idx/protein_idx refer to the ATOM/HETATM order of the PDB.
    """

    def get_coords(idx):
        try:
            idx = int(idx)
            x, y, z = atom_coords[idx]
        except:
            idx = [int(i.strip()) for i in idx[1:-1].split(',')]
            coords_ = np.array([atom_coords[i] for i in idx])
            x, y, z = coords_.mean(axis=0)
        return x, y, z

    if not csv_path or not os.path.exists(csv_path):
        return

    df = pd.read_csv(csv_path)
    
    if itype_filter:
        df = df.loc[df['itype'].isin(itype_filter)]
    
    for lig_idx_raw, prot_idx_raw in zip(df['ligand_idx'], df['protein_idx']):

        x1, y1, z1 = get_coords(lig_idx_raw)
        x2, y2, z2 = get_coords(prot_idx_raw)

        viewer.addCylinder({
            "start": {"x": x1, "y": y1, "z": z1},
            "end":   {"x": x2, "y": y2, "z": z2},
            "color": color,
            "radius": radius,
            "dashed": True,
            "dashLength": dash_length,
            "gapLength": gap_length,
            "fromCap": 1,
            "toCap": 1,
        })

def _highlight_interaction_residues_from_csv(
    viewer: py3Dmol.view,
    model_index: int,
    csv_path: str,
    residue_style: Dict[str, Any],
    itype_filter: Optional[List[str]] = None,  # e.g. ["hbond"]
) -> None:
    """
    Use the interaction CSV to highlight protein residues involved in interactions.

    residue_style is a 3Dmol style dict, e.g.:
    {
        "cartoon": {"color": "white", "opacity": 1.0},
        "stick": {"color": "yellow", "radius": 0.25}
    }
    """
    if not csv_path or not os.path.exists(csv_path):
        return

    df = pd.read_csv(csv_path)
    
    if itype_filter:
        df = df.loc[df['itype'].isin(itype_filter)]

    for chain, resnum in zip(df['prot_chain'], df['prot_res_num']):
        viewer.setStyle(
            {
                "model": model_index,
                "chain": chain,
                "resi": int(resnum),
                "hetflag": False,
            },
            residue_style,
        )
            
def _parse_box_constraint(path: str):
    """Parse simple key=value lines for center_*, size_*."""
    vals = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"): 
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                vals[k.strip()] = float(v.strip())
    cx = vals.get("center_x"); cy = vals.get("center_y"); cz = vals.get("center_z")
    sx = vals.get("size_x");   sy = vals.get("size_y");   sz = vals.get("size_z")
    if None in (cx, cy, cz, sx, sy, sz):
        raise ValueError("box_constraint missing one of center_{x,y,z} or size_{x,y,z}")
    return (cx, cy, cz), (sx, sy, sz)

def _add_wire_box(viewer, center, size, color="#ff0000", opacity=0.6):
    """Add a wireframe box using 12 cylinders."""
    cx, cy, cz = center
    sx, sy, sz = size
    hx, hy, hz = sx/2.0, sy/2.0, sz/2.0

    viewer.addBox({
        "center": {"x": cx, "y": cy, "z": cz},
        "dimensions": {"w": sx, "h": sy, "d": sz},   # width, height, depth
        "color": color,
        "opacity": opacity,
        "wireframe": False    # or False for filled box
            })
    viewer.addBox({
        "center": {"x": cx, "y": cy, "z": cz},
        "dimensions": {"w": sx, "h": sy, "d": sz},   # width, height, depth
        "color": color,
        "opacity": opacity + 0.1,
        "wireframe": True    # or False for filled box
            })
    return viewer


def render_structures(
    files: List[str],
    style_rules: Optional[List[Dict[str, Union[Sel, Sty]]]] = None,
    surface_rules: Optional[List[Dict[str, Any]]] = None,
    label_rules: Optional[List[Dict[str, Any]]] = None,
    chain_color_map: Optional[Dict[Union[int,str], Dict[str,str]]] = None,
    interaction_tables: Optional[List[Optional[str]]] = None,
    interaction_color: str = "red",
    interaction_types: Optional[List[str]] = None,
    box_constraint: str | None = None,
    box_color: str = "blue",
    box_opacity: float = 0.3,
    background: str = "white",
    width: int = 1000,
    height: int = 800,
    zoom: bool = True,
    file_name: str = 'my_protein',
) -> dict:
    """Render multiple molecular files with fine-grained control. Read the help doc get_tools_doc('render_structures')"""
    # Build viewer
    viewer = py3Dmol.view(width=width, height=height)
    viewer.setBackgroundColor(background)

    # Load models
    model_meta = []  # list of dicts: {"path":..., "fmt":..., "text":..., "chains":[...]}
    for p in files:
        fmt = _infer_format(p)
        txt = _read_file(p)
        viewer.addModel(txt, fmt)
        chains = _list_chains_from_pdb_text(txt) if fmt in ("pdb","cif") else []
        coords = _extract_atom_coords_from_pdb_text(txt) if fmt in ("pdb",) else []
        model_meta.append({"path": p, "fmt": fmt, "text": txt, "chains": chains, "coords": coords})

    # Clear any default style
    viewer.setStyle({}, {})

    # Smart defaults (only if user didn't supply style_rules)
    if not style_rules:
        for i, meta in enumerate(model_meta):
            if _is_polymeric_format(meta["fmt"]) and meta["chains"]:
                # Show polymeric models as cartoon; ligands (HETATM) as sticks
                viewer.setStyle({"model": i, "hetflag": False}, {"cartoon": {"color":"spectrum"}})
                viewer.setStyle({"model": i, "hetflag": True},  {"stick": {"radius": 0.25}})
            else:
                # SDF/MOL2 or non-polymeric -> sticks
                viewer.setStyle({"model": i}, {"stick": {"radius": 0.25}})

    # Optional: chain-level coloring helpers
    if chain_color_map:
        # Prepare a list of target models based on keys
        def _targets_for_key(k):
            if k == "all":
                return list(range(len(model_meta)))
            if isinstance(k, int):
                return [k]
            if isinstance(k, str) and k.isdigit():
                return [int(k)]
            return []
        for key, cmap in chain_color_map.items():
            targets = _targets_for_key(key)
            for mi in targets:
                if _is_polymeric_format(model_meta[mi]["fmt"]):
                    chains = model_meta[mi]["chains"] or ["_"]
                    for ch in chains:
                        col = cmap.get(ch) or cmap.get("_")
                        if col:
                            viewer.setStyle(
                                {"model": mi, "chain": ch, "hetflag": False},
                                {"cartoon": {"color": col}}
                            )

    # Apply explicit style rules (in order; later rules can override earlier)
    if style_rules:
        for rule in style_rules:
            sel = dict(rule.get("select", {}))
            sty = dict(rule.get("style", {}))
            # Expand model lists -> apply once per model for clarity
            models = sel.get("model", None)
            if isinstance(models, list):
                for m in models:
                    s2 = dict(sel); s2["model"] = m
                    viewer.setStyle(s2, sty)
            else:
                viewer.setStyle(sel, sty)

    # Surfaces
    if surface_rules:
        for rule in surface_rules:
            sel = dict(rule.get("select", {}))
            srf = dict(rule.get("surface", {}))
            stype = rule.get("type", py3Dmol.VDW)
            viewer.addSurface(stype, srf, sel)

    # Labels
    if label_rules:
        for rule in label_rules:
            text = rule.get("text","")
            if "position" in rule:
                pos = rule["position"]
                viewer.addLabel(text, {"position": pos, **rule.get("style", {})})
            elif "select" in rule:
                # Label center of selection’s bounding box
                viewer.addLabel(text, {"backgroundOpacity": 0.0, **rule.get("style", {})}, rule["select"])
            else:
                viewer.addLabel(text, rule.get("style", {}))


    # Interactions
    if interaction_tables:
        interaction_radius = 0.1
        interaction_dash_length = .5
        interaction_gap_length = .5
        interaction_residue_style = {
            "stick": {
                "radius": 0.25,
                "colorscheme": "element"   # carbon gray, oxygen red, nitrogen blue, etc.
            }
        }
        
        for mi, (meta, csv_path) in enumerate(zip(model_meta, interaction_tables)):
            if not csv_path:
                continue
            _add_interaction_lines_from_csv(
                viewer=viewer,
                atom_coords=meta["coords"],
                csv_path=csv_path,
                color=interaction_color,
                radius=interaction_radius,
                dash_length=interaction_dash_length,
                gap_length=interaction_gap_length,
                itype_filter=interaction_types
            )
            _highlight_interaction_residues_from_csv(
                viewer=viewer,
                model_index=mi,
                csv_path=csv_path,
                residue_style=interaction_residue_style,
                itype_filter=interaction_types,
            )

            
    # Add box from box_constraint file
    if box_constraint:
        center, dims = _parse_box_constraint(box_constraint)
        viewer = _add_wire_box(viewer, center, dims, color=box_color, opacity=box_opacity)
    
    if zoom:
        viewer.zoomTo()

    file_name = os.path.join(FILE_DIR, file_name + '.html')
    html = viewer.write_html()
    with open(file_name, "w") as fh:
        fh.write(html)
    return {
        'message' : f"3D structure is saved at {file_name}.",
        'html_file' : file_name
    }
