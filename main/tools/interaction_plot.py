import pandas as pd
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import ast
import re
import os

FILE_DIR = os.getenv("FILE_DIR")

BASIC_COLOR_NAMES = {
    # basic set; Plotly understands many more, but this catches clear typos
    "black", "white", "red", "green", "blue", "yellow",
    "cyan", "magenta", "orange", "purple", "grey", "gray",
    "silver", "lightgray", "lightgrey", "darkgray", "darkgrey",
    "brown", "pink", "gold", "lime", "navy", "teal",
    "violet", "indigo", "maroon",
}

METAL_SYMBOLS = {
    "Li", "Na", "K", "Rb", "Cs",
    "Be", "Mg", "Ca", "Sr", "Ba",
    "Fe", "Zn", "Mn", "Cu", "Co", "Ni", "Cr",
    "Al", "Cd", "Hg",
}

DEFAULT_ATOM_COLORS = {
    "C": "black",
    "H": "white",
    "O": "red",
    "N": "blue",
    "S": "yellow",
    "P": "orange",
    "F": "lightgreen",
    "Cl": "green",
    "Br": "brown",
    "I": "violet",
}

DEFAULT_EDGE_FEATURE_COLORS = {
    "single": "#555555",
    "double": "#555555",
    "aromatic": "#555555",
    "peptide": "#cccccc",
    "ligand_hbonf_acceptor": "#0033cc",
    "ligand_hbonf_donnor": "#4f7fff",
    "ionic_ligand_negative": "#ff0000",
    "hydrophobic": "#00ff00",
    "ligand_halogen_donnor": "#ff00ff",
    "waterbridge": "#00ffff",
}

PROTEIN_NODE_RADIUS_FACTOR = 1.3
PEPTIDE_RADIUS_FACTOR = 0.6
PL_RADIUS_FACTOR = 0.6
PL_DASHES = 7
PL_DUTY = 0.5

# -------------------------------------------------------------------------
# HELPER FUNCTIONS
# -------------------------------------------------------------------------

# -----------------
# Read and re-indexing
# -----------------

def _reindex(nodes, edges):
    import ast
    import pandas as pd

    # -------------------------
    # Helpers: parse index cells
    # -------------------------
    def _as_list(x):
        """
        Always return a list of ints.
        Accepts: int/float, list, or string like "[14, 15]".
        """
        if pd.isna(x):
            return []
        if isinstance(x, (list, tuple)):
            return [int(v) for v in x]
        s = str(x).strip()
        if s.startswith("[") and s.endswith("]"):
            return [int(v) for v in ast.literal_eval(s)]
        return [int(s)]

    # -------------------------
    # Basic checks
    # -------------------------
    for col in ["index", "molecule"]:
        if col not in nodes.columns:
            raise ValueError(f"Node file must contain a column '{col}'.")
    for col in ["index1", "index2", "type"]:
        if col not in edges.columns:
            raise ValueError(f"Edge file must contain a column '{col}'.")

    nodes = nodes.copy()
    edges = edges.copy()

    # Create/keep old_index
    if "old_index" not in nodes.columns:
        nodes["old_index"] = nodes["index"]

    lig_mask = nodes["molecule"].astype(str).str.lower() == "ligand"
    prot_mask = nodes["molecule"].astype(str).str.lower() == "protein"

    if not lig_mask.any():
        raise ValueError("No ligand nodes found (molecule == 'ligand').")
    if not prot_mask.any():
        raise ValueError("No protein nodes found (molecule == 'protein').")

    lig_nodes = nodes[lig_mask].copy()
    prot_nodes = nodes[prot_mask].copy()

    lig_old_ids = set(lig_nodes["old_index"].astype(int))
    prot_old_ids = set(prot_nodes["old_index"].astype(int))

    # -------------------------
    # Build maps: old_index -> new contiguous index
    # Ligands first, then proteins
    # -------------------------
    lig_map = {}
    current = 1
    for old in sorted(lig_old_ids):
        lig_map[old] = current
        current += 1

    prot_map = {}
    for old in sorted(prot_old_ids):
        prot_map[old] = current
        current += 1

    # -------------------------
    # Map nodes to new_index
    # -------------------------
    def _map_node_row(row):
        old_idx = int(row["old_index"])
        mol = str(row["molecule"]).lower()
        if mol == "ligand":
            return lig_map[old_idx]
        elif mol == "protein":
            return prot_map[old_idx]
        else:
            raise ValueError(f"Unknown molecule type '{row['molecule']}' in row:\n{row}")

    nodes["new_index"] = nodes.apply(_map_node_row, axis=1)

    nodes = nodes.set_index("new_index").sort_index()
    nodes.index.name = "index"

    front_cols = [c for c in ["old_index", "molecule", "symbol", "x", "y", "z"] if c in nodes.columns]
    other_cols = [c for c in nodes.columns if c not in front_cols]
    nodes = nodes[front_cols + other_cols]

    # -------------------------
    # Map edges (support list-valued index1/index2)
    # -------------------------
    def _map_edge_row(row):
        t = str(row["type"]).lower()

        i1_list = _as_list(row["index1"])
        i2_list = _as_list(row["index2"])

        if t == "ligand-ligand":
            bad1 = [i for i in i1_list if i not in lig_map]
            bad2 = [i for i in i2_list if i not in lig_map]
            if bad1 or bad2:
                raise ValueError(f"Ligand-ligand edge uses non-ligand indices: index1={bad1}, index2={bad2}")

            row["index1"] = [lig_map[i] for i in i1_list]
            row["index2"] = [lig_map[i] for i in i2_list]

        elif t == "protein-protein":
            bad1 = [i for i in i1_list if i not in prot_map]
            bad2 = [i for i in i2_list if i not in prot_map]
            if bad1 or bad2:
                raise ValueError(f"Protein-protein edge uses non-protein indices: index1={bad1}, index2={bad2}")

            row["index1"] = [prot_map[i] for i in i1_list]
            row["index2"] = [prot_map[i] for i in i2_list]

        elif t == "ligand-protein":
            # enforce direction: index1 ligand, index2 protein (in OLD index space)
            bad1 = [i for i in i1_list if i not in lig_old_ids]
            bad2 = [i for i in i2_list if i not in prot_old_ids]
            if bad1 or bad2:
                raise ValueError(
                    "ligand-protein edge does not satisfy (index1 ligand, index2 protein): "
                    f"bad_index1={bad1}, bad_index2={bad2}"
                )

            row["index1"] = [lig_map[i] for i in i1_list]
            row["index2"] = [prot_map[i] for i in i2_list]

        else:
            raise ValueError(f"Unknown edge type '{row['type']}' in row:\n{row}")

        return row

    edges = edges.apply(_map_edge_row, axis=1)

    return nodes.copy(), edges.copy()


def _parse_idx(x):
    if isinstance(x, str) and x.strip().startswith('['):
        return ast.literal_eval(x)
    return int(x)
# -----------------
# Geometry
# -----------------
def _create_sphere_mesh(center, radius=1.0, n_lat=8, n_lon=16):
    cx, cy, cz = center
    phis = np.linspace(0, np.pi, n_lat + 1)
    thetas = np.linspace(0, 2 * np.pi, n_lon, endpoint=False)

    xs, ys, zs = [], [], []
    for phi in phis:
        for theta in thetas:
            x = radius * np.sin(phi) * np.cos(theta) + cx
            y = radius * np.sin(phi) * np.sin(theta) + cy
            z = radius * np.cos(phi) + cz
            xs.append(x); ys.append(y); zs.append(z)

    xs = np.array(xs); ys = np.array(ys); zs = np.array(zs)

    i_list, j_list, k_list = [], [], []
    for lat in range(n_lat):
        for lon in range(n_lon):
            i0 = lat * n_lon + lon
            i1 = lat * n_lon + (lon + 1) % n_lon
            i2 = (lat + 1) * n_lon + lon
            i3 = (lat + 1) * n_lon + (lon + 1) % n_lon

            i_list += [i0, i1]
            j_list += [i2, i2]
            k_list += [i1, i3]

    return xs, ys, zs, np.array(i_list), np.array(j_list), np.array(k_list)

def _create_cylinder_mesh(p0, p1, radius=0.2, n_segments=12, n_height=1):
    p0_ = np.asarray(p0, dtype=float)
    p1_ = np.asarray(p1, dtype=float)
    v = p1_ - p0_
    L = np.linalg.norm(v)
    if L == 0:
        return None

    v_hat = v / L

    if abs(v_hat[0]) < 0.9:
        a = np.array([1.0, 0.0, 0.0])
    else:
        a = np.array([0.0, 1.0, 0.0])

    u = np.cross(v_hat, a)
    u /= np.linalg.norm(u)
    w = np.cross(v_hat, u)

    zs_param = np.linspace(0, L, n_height + 1)
    thetas = np.linspace(0, 2 * np.pi, n_segments, endpoint=False)

    xs, ys, zs = [], [], []
    for z_param in zs_param:
        for theta in thetas:
            rcos = radius * np.cos(theta)
            rsin = radius * np.sin(theta)
            point = p0_ + rcos * u + rsin * w + z_param * v_hat
            xs.append(point[0]); ys.append(point[1]); zs.append(point[2])

    xs = np.array(xs); ys = np.array(ys); zs = np.array(zs)

    i_list, j_list, k_list = [], [], []
    for h in range(n_height):
        for s in range(n_segments):
            i0 = h * n_segments + s
            i1 = h * n_segments + (s + 1) % n_segments
            i2 = (h + 1) * n_segments + s
            i3 = (h + 1) * n_segments + (s + 1) % n_segments

            i_list += [i0, i1]
            j_list += [i2, i2]
            k_list += [i1, i3]

    return xs, ys, zs, np.array(i_list), np.array(j_list), np.array(k_list)

def _add_cylinder_trace(fig, p0, p1, color, radius, opacity, name, showlegend=False):
    res = _create_cylinder_mesh(p0, p1, radius=radius, n_segments=14, n_height=1)
    if res is None:
        return
    xs, ys, zs, i, j, k = res
    fig.add_trace(go.Mesh3d(
        x=xs, y=ys, z=zs,
        i=i, j=j, k=k,
        color=color,
        opacity=opacity,
        flatshading=False,
        showscale=False,
        lighting=dict(
            ambient=0.4, diffuse=0.9,
            specular=0.3, roughness=0.6
        ),
        name=name,
        showlegend=showlegend
    ))

def _add_dashed_cylinders(fig, p0, p1, color, radius, opacity, name,
                          num_dashes=7, duty_cycle=0.5):
    p0_ = np.asarray(p0, dtype=float)
    p1_ = np.asarray(p1, dtype=float)
    v = p1_ - p0_
    L = np.linalg.norm(v)
    if L == 0:
        return
    v_hat = v / L

    dash_len = L / num_dashes
    on_len = dash_len * duty_cycle

    for n in range(num_dashes):
        start = n * dash_len
        end = start + on_len
        if start >= L:
            break
        if end > L:
            end = L

        seg_p0 = p0_ + v_hat * start
        seg_p1 = p0_ + v_hat * end
        _add_cylinder_trace(fig, seg_p0, seg_p1, color, radius, opacity, name)

def _create_ring_fill_mesh(points3d):
    pts = np.asarray(points3d, dtype=float)
    n = pts.shape[0]
    if n < 3:
        return None

    center = pts.mean(axis=0)
    pts_c = pts - center

    _, _, Vt = np.linalg.svd(pts_c, full_matrices=False)
    e1 = Vt[0]
    e2 = Vt[1]

    coords2d = pts_c @ np.vstack([e1, e2]).T
    angles = np.arctan2(coords2d[:, 1], coords2d[:, 0])
    order = np.argsort(angles)
    pts_ord = pts[order]

    xs = np.concatenate([pts_ord[:, 0], [center[0]]])
    ys = np.concatenate([pts_ord[:, 1], [center[1]]])
    zs = np.concatenate([pts_ord[:, 2], [center[2]]])
    center_idx = n

    i_list, j_list, k_list = [], [], []
    for k in range(n):
        i_list.append(k)
        j_list.append((k + 1) % n)
        k_list.append(center_idx)

    return xs, ys, zs, np.array(i_list), np.array(j_list), np.array(k_list)

def _double_offset_vector(p0, p1, scale):
    dvec = p1 - p0
    if np.allclose(dvec, 0):
        return np.array([0, 0, 0], dtype=float)
    v = np.array([0., 0., 1.])
    if np.allclose(np.cross(dvec, v), 0.):
        v = np.array([0., 1., 0.])
    normal = np.cross(dvec, v)
    norm = np.linalg.norm(normal)
    if norm == 0:
        return np.array([0., 0., 0.])
    normal /= norm
    return normal * scale
# -----------------
# Get colors and styles
# -----------------
def _is_valid_color_string(c):
    if not isinstance(c, str):
        return False
    s = c.strip()
    if not s:
        return False
    # hex #rgb or #rrggbb
    if re.fullmatch(r"#([0-9a-fA-F]{3}|[0-9a-fA-F]{6})", s):
        return True
    # basic named colors
    if s.lower() in BASIC_COLOR_NAMES:
        return True
    return False

def _get_atom_color(symbol, atom_global_color = None):
    # global override
    if atom_global_color is not None:
        return atom_global_color

    s_raw = str(symbol).strip()
    # metals
    if s_raw in METAL_SYMBOLS:
        return "lightgray"

    # default table
    if s_raw in DEFAULT_ATOM_COLORS:
        return DEFAULT_ATOM_COLORS[s_raw]
    s_cap = s_raw.capitalize()
    if s_cap in DEFAULT_ATOM_COLORS:
        return DEFAULT_ATOM_COLORS[s_cap]
    return "silver"


# -------------------------------------------------------------------------
# HELPERS: color validation & symbol/feature resolution
# -------------------------------------------------------------------------

def _ask_color(prompt, default):
    """Ask for a color; empty = default; invalid → re-ask."""
    while True:
        val = input(prompt).strip()
        if not val:
            return default
        if _is_valid_color_string(val):
            return val
        print("  Invalid color. Use a named color (e.g. 'red') or hex like '#ff0000'. Try again.")

def _ask_optional_color(prompt):
    """Ask for a color; empty = None; invalid → re-ask."""
    while True:
        val = input(prompt).strip()
        if not val:
            return None
        if _is_valid_color_string(val):
            return val
        print("  Invalid color. Use a named color (e.g. 'red') or hex like '#ff0000'. Try again.")

def _resolve_symbol(sym):
    """Return the exact symbol as in df_nodes if present, else None."""
    s = sym.strip()
    candidates = [s, s.capitalize(), s.upper(), s.lower()]
    for cand in candidates:
        if cand in symbols_raw:
            return cand
    return None

def _resolve_feature(feat):
    """Return the exact feature string as in df_edges if present, else None."""
    f = feat.strip()
    if f in features_raw:
        return f
    lower_map = {x.lower(): x for x in features_raw}
    return lower_map.get(f.lower(), None)

def interaction_plot(nodes_csv_dir: str, 
                     edges_csv_dir: str, 
                     style: dict|None = None,
                     file_name: str = 'interaction_plot'
                     # save_dir = FILE_DIR,
                    ) -> dict:
    """Visuzalize Protein-ligand interactions given complex results from `get_protein_ligand_interaction` tool"""
    # 1) Reads raw node/edge CSVs for a complex.
    # 2) Reindexes ligand/protein atoms into a compact, contiguous index space.
    # 5) Saves the plot as an HTML file.

    # -------------------------------------------------------------------------
    # 1) READ & REINDEX
    # -------------------------------------------------------------------------
    df_nodes = pd.read_csv(nodes_csv_dir)
    df_edges = pd.read_csv(edges_csv_dir)
    df_nodes, df_edges = _reindex(df_nodes, df_edges)

    # Get all node types and edge types
    node_types = {str(s).strip() for s in df_nodes["symbol"].dropna().unique()}
    edge_types = {str(f).strip() for f in df_edges["feature"].dropna().unique()}
    
    # -------------------------------------------------------------------------
    # 2) ASK USER FOR VISUAL OPTIONS (with validation)
    # -------------------------------------------------------------------------
    if style is None:
        style = {
            # Basic
            'show_legend' : True,
            'show_atom_labels' : True,
            'background_color' : "white", #"#f5f5f5",
            # Aromatic rings
            'aromatic_ring_color' : "#d3d3d3",
            'aromatic_ring_opacity' : 0.30,
            # Atoms
            'atom_color' : DEFAULT_ATOM_COLORS,
            'atom_global_color' : None,
            # Bonds
            'ligand_ligand_bond_color' : "#555555",
            'peptide_bond_color' : "#cccccc",
            'proten_ligand_bond_color' : "#0066ff",
            'edge_feature_color' : DEFAULT_EDGE_FEATURE_COLORS,
            'bond_global_color' : None,
            # Radii
            'atom_radius' : 0.25,
            'bond_radius' : 0.10,
            
            
        }

    
    file_name = os.path.join(FILE_DIR, file_name + '.html')
    
    # -------------------------------------------------------------------------
    # 3) MAIN PLOTTING LOGIC
    # -------------------------------------------------------------------------

    # Sanitize df
    # df_edges["index1"] = df_edges["index1"].apply(_parse_idx)
    # df_edges["index2"] = df_edges["index2"].apply(_parse_idx)
    df_nodes[["x", "y", "z"]] = df_nodes[["x", "y", "z"]].astype(float)

    # get aromatic cycles
    arom_edges = df_edges[df_edges["feature"] == "aromatic"]
    
    G_arom = nx.Graph()
    for _, e in arom_edges.iterrows():
        srcs = e["index1"] if isinstance(e["index1"], list) else [e["index1"]]
        dsts = e["index2"] if isinstance(e["index2"], list) else [e["index2"]]
        for s in srcs:
            for d in dsts:
                G_arom.add_edge(s, d)

    aromatic_cycles = nx.cycle_basis(G_arom)

    # colors for bond TYPES
    FEATURE_COLORS = style['edge_feature_color']
    TYPE_COLORS = {
        "ligand-ligand": style['ligand_ligand_bond_color'],
        "protein-protein": style['peptide_bond_color'],
        "ligand-protein": style['proten_ligand_bond_color'],
    }

    fig = go.Figure()

    # aromatic ring fill
    for cycle in aromatic_cycles:
        ring_points = df_nodes.loc[cycle, ["x", "y", "z"]].values
        res = _create_ring_fill_mesh(ring_points)
        if res is None:
            continue
        xs, ys, zs, i, j, k = res
        fig.add_trace(go.Mesh3d(
            x=xs, y=ys, z=zs,
            i=i, j=j, k=k,
            opacity=style['aromatic_ring_opacity'],
            color=style['aromatic_ring_color'],
            showscale=False,
            name="aromatic ring",
            lighting=dict(ambient=0.5, diffuse=0.8, specular=0.3, roughness=0.5),
            showlegend=False,
        ))

    # atoms
    for idx, node in df_nodes.iterrows():
        center = np.array([node["x"], node["y"], node["z"]], dtype=float)
        atom_label = str(node.get("atom_name", node["symbol"]))

        mol_type = str(node.get("molecule", "")).lower()
        if mol_type == "protein":
            color = style['atom_global_color'] if style['atom_global_color'] is not None else "#e0c6ff"
            radius_here = style['atom_radius'] * PROTEIN_NODE_RADIUS_FACTOR
        else:
            color = _get_atom_color(node["symbol"])
            radius_here = style['atom_radius']

        xs, ys, zs, i, j, k = _create_sphere_mesh(
            center, radius=radius_here, n_lat=12, n_lon=24
        )

        customdata = np.array([atom_label] * len(xs))

        fig.add_trace(go.Mesh3d(
            x=xs, y=ys, z=zs,
            i=i, j=j, k=k,
            color=color,
            opacity=1.0,
            flatshading=False,
            showscale=False,
            customdata=customdata,
            hovertemplate=(
                "%{customdata}<br>"
                "x: %{x:.3f}<br>"
                "y: %{y:.3f}<br>"
                "z: %{z:.3f}<extra></extra>"
            ),
            lighting=dict(
                ambient=0.4, diffuse=0.9,
                specular=0.5, roughness=0.4, fresnel=0.1
            ),
            name=atom_label,
            showlegend=False,
        ))

    if style['show_atom_labels']:
        fig.add_trace(go.Scatter3d(
            x=df_nodes["x"],
            y=df_nodes["y"],
            z=df_nodes["z"],
            mode="text",
            text=df_nodes["symbol"].astype(str),
            textposition="top center",
            showlegend=False,
        ))

    # bonds
    for _, e in df_edges.iterrows():
        feat_raw = str(e["feature"]).strip()
        feat = feat_raw
        bond_type_raw = str(e["type"]).lower()

        srcs = e["index1"] if isinstance(e["index1"], list) else [e["index1"]]
        dsts = e["index2"] if isinstance(e["index2"], list) else [e["index2"]]

        for s in srcs:
            for d in dsts:
                if s not in df_nodes.index or d not in df_nodes.index:
                    continue

                p0 = df_nodes.loc[s, ["x", "y", "z"]].values.astype(float)
                p1 = df_nodes.loc[d, ["x", "y", "z"]].values.astype(float)

                # default style
                radius_here = style['bond_radius']
                dashed = False
                double_bond = False

                # base color from bond TYPE or feature
                base_color = TYPE_COLORS.get(
                    bond_type_raw,
                    FEATURE_COLORS.get(feat, "#555555")
                )

                name = f"{bond_type_raw} ({feat})" if feat and feat.lower() != "nan" else bond_type_raw

                # refine style based on type / feature
                if bond_type_raw == "ligand-ligand":
                    # double bonds as two parallel tubes
                    if feat == "double":
                        double_bond = True

                elif bond_type_raw == "protein-protein":
                    radius_here = style['bond_radius'] * PEPTIDE_RADIUS_FACTOR

                elif bond_type_raw == "ligand-protein":
                    dashed = True
                    radius_here = style['bond_radius'] * PL_RADIUS_FACTOR
                    if feat == "ligand_hbonf_acceptor":
                        name = "H-bond acceptor (ligand-protein)"
                    elif feat == "ligand_hbonf_donnor":
                        name = "H-bond donor (ligand-protein)"
                    else:
                        name = "Ligand-protein interaction"

                # feature-specific color (subtype) overrides type color if available
                color = base_color
                if feat in FEATURE_COLORS:
                    color = FEATURE_COLORS[feat]

                # global bond color override (highest priority)
                if style['bond_global_color'] is not None:
                    color = style['bond_global_color']

                if double_bond:
                    offset_vec = _double_offset_vector(
                        p0, p1, scale=style['bond_radius'] * 1.5
                    )
                    for sign in (+1, -1):
                        _add_cylinder_trace(
                            fig,
                            p0 + sign * offset_vec,
                            p1 + sign * offset_vec,
                            color=color,
                            radius=radius_here * 0.9,
                            opacity=1.0,
                            name="Ligand-ligand double bond",
                        )
                elif dashed:
                    _add_dashed_cylinders(
                        fig,
                        p0,
                        p1,
                        color=color,
                        radius=radius_here,
                        opacity=1.0,
                        name=name,
                        num_dashes=PL_DASHES,
                        duty_cycle=PL_DUTY,
                    )
                else:
                    _add_cylinder_trace(
                        fig,
                        p0,
                        p1,
                        color=color,
                        radius=radius_here,
                        opacity=1.0,
                        name=name,
                    )

    # scene / layout
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False, showgrid=False, zeroline=False),
            yaxis=dict(visible=False, showgrid=False, zeroline=False),
            zaxis=dict(visible=False, showgrid=False, zeroline=False),
            aspectmode="data",
            bgcolor=style['background_color'],
        ),
        paper_bgcolor=style['background_color'],
        margin=dict(l=0, r=0, t=0, b=0),
        scene_camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)),
    )

    # full-screen-ish
    fig.update_layout(autosize=True, width=1600, height=900)
    config = {
        "scrollZoom": True,
        "displaylogo": False,
        "responsive": True,
    }

    # legend (uses type and key subtype colors you may have overridden)
    if style['show_legend']:
        legend_items = [
            ("Ligand–ligand (single/double/aromatic)",
             dict(color=TYPE_COLORS.get("ligand-ligand", "#555555"), width=8, dash="solid")),
            ("Peptide (protein–protein)",
             dict(color=TYPE_COLORS.get("protein-protein", "#cccccc"), width=5, dash="solid")),
            ("H-bond acceptor (prot–lig)",
             dict(color=FEATURE_COLORS.get("ligand_hbonf_acceptor", "#555555"), width=5, dash="dash")),
            ("H-bond donor (prot–lig)",
             dict(color=FEATURE_COLORS.get("ligand_hbonf_donnor", "#555555"), width=5, dash="dash")),
            ("hydrophobic (prot–lig)",
             dict(color=FEATURE_COLORS.get("hydrophobic", "#555555"), width=5, dash="dash")),
            ("ligand_halogen_donnor (prot–lig)",
             dict(color=FEATURE_COLORS.get("ligand_halogen_donnor", "#555555"), width=5, dash="dash")),
            ("ligand_ionic_negative (prot–lig)",
             dict(color=FEATURE_COLORS.get("ionic_ligand_negative", "#555555"), width=5, dash="dash")),
            ("waterbridge (prot–lig)",
             dict(color=FEATURE_COLORS.get("waterbridge", "#555555"), width=5, dash="dash")),
        ]

        for name, line_style in legend_items:
            fig.add_trace(go.Scatter3d(
                x=[None], y=[None], z=[None],
                mode="lines",
                line=line_style,
                name=name,
                showlegend=True,
            ))

        fig.update_layout(
            showlegend=True,
            legend=dict(
                title="Interactions",
                itemsizing="constant",
                orientation="v",
                x=0.01,
                y=0.99,
                font=dict(color="black"),
            ),
        )
    else:
        fig.update_layout(showlegend=False)

    # show full-window
    # fig.show(config=config)

    # save HTML if requested
    fig.write_html(file_name, config=config)

    return {
        'message' : f"interaction plot is saved at {file_name}.",
        'html_file' : file_name
    }