import os
foldX_rel_path = "../software/FoldX"
foldX_abs_path = os.path.abspath(foldX_rel_path)
os.environ["FOLDX_BINARY"] = foldX_abs_path

from pyfoldx.structure import Structure

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import gemmi
from collections import defaultdict
import warnings
from admet_ai import ADMETModel
from .dock_utilities import *

# Turn off rdkit warning
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# from dotenv import load_dotenv
# load_dotenv(f"envs/env_{st.session_state['user']}")
FILE_DIR = os.getenv("FILE_DIR")  
# FILE_DIR = 'files/'

## ========= helper fnctions =============
# def obabel(*args):
#     obabel_loc = 'obabel'   # /project/drug_discovery/drug_discovery_agent/env_chem_agent/bin/
#     command = obabel_loc + ' ' + args[0]
#     # command = 'timeout 1000s ' + command     # only for linux
#     os.system(command)

def _add_tag_to_path(filepath, tag):
    base, ext = os.path.splitext(filepath)
    return f"{base}_{tag}{ext}"


def remove_counterions(smiles):
    """Break the smiles into fragments and choose the one with most heavy atoms"""
    mol = Chem.MolFromSmiles(smiles)
    frags = Chem.GetMolFrags(mol, asMols=True)
    parent = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    return Chem.MolToSmiles(parent)

    
def extract_native_ligand_from_pdb(name: str, 
                                   input_file: str, 
                                   # save_dir: str = FILE_DIR
                                  ):
    save_dir = FILE_DIR
    Keyword = name
    true_id = Keyword[-3:] if len(Keyword) > 3 else Keyword
    
    ntv_pdb = Keyword.upper() + '.pdb'
    ntv_pdb_nFile = os.path.join(save_dir, ntv_pdb)
    extract_entity(input_file, ntv_pdb_nFile, [Keyword, 'HETATM'])
    extract_chains(ntv_pdb_nFile)
    ntv_nFiles = sorted(glob.glob(save_dir + '/' + Keyword + '_*.pdb'))
    try:
        ntv_mol = get_molblock(true_id, save_dir)
        ntv_smiles = Chem.MolToSmiles(ntv_mol)
        correct_bond_order(ntv_nFiles, ntv_mol)
    except:
        ntv_smiles = None
        print(f"The molecule {name} not availabe at RCSB")
    return {
        'Message' : 'Native ligands are extracted',
        'SMILES' : ntv_smiles,
        'save_dir' : ntv_nFiles,
    }


def optimize_protein_foldx(pdb_file: str):
    
    st = Structure('my_protein', path=pdb_file)

    # Calculate initial energy
    initial_energy = st.getTotalEnergy()
    # Get energy per residue
    # df_res_initial = st.getResiduesEnergy()
    
    # Repair the structure
    print("Optimization started. this may take some time ...")
    st_repaired = st.repair()
    print("Done!")

    # Calculate the repaired energy
    repaired_energy = st_repaired.getTotalEnergy()
    # Get energy per residue
    # df_res_repaired = st_repaired.getResiduesEnergy()
    
    # Save the repaired structure
    output_file = _add_tag_to_path(pdb_file, 'repaired')            
    st_repaired.toPdbFile(output_file)

    return {
        'message' : 'Done!',
        'initial_energy' : initial_energy['total'].item(),
        'repaired_energy' : repaired_energy['total'].item(),
        'saved_file' : output_file
    }


def protonate_protein(pdb_file: str,
                      ph: float = 7.2,
                      force_field: str = 'amber'
                     ):
    """Protonate and optimize protein based on propka and pdb2pqr"""
    # force_field = ['AMBER' ,'CHARMM' ,'PARSE', 'TYL06', 'PEOEPB', 'SWANSON']
    
    ph_str = str(ph).replace(".", "")             # 7.2 -> 72
    output_file = _add_tag_to_path(pdb_file, f'ph{ph_str}')    
    response = os.system(f"pdb2pqr --ff={force_field.upper()} --with-ph={ph} --pdb-output {output_file} {pdb_file} {FILE_DIR}/output_pqr.pqr")

    if response == 0:
        return{
            'message' : 'Successful!',
            'ph' : ph,
            'saved_file' : output_file
        }
    else:
        return{
            'message' : 'Not Successful!',
            'ph' : ph,
            'saved_file' : None
        }

## ========= Main functions ==============

# 1. Ligand analysis =====================
def smiles_to_3d(ligand_smiles: str, 
                 ligand_name: str = 'my_ligand', 
                 # save_dir: str = FILE_DIR,
                 method: str = 'obabel',
                 force_field: str = 'MMFF94',
                 convergence_criteria: str = '0.00001',
                 maximum_steps: int = 10000,
                 verbos: bool = True,
                ):
    """Convert molecular SMILES notation to sdf and pdb files by ligand energy minimization. It automatically removes counter ions if any"""

    save_dir = FILE_DIR
    
    ligand_sdf_dir = os.path.join(save_dir, ligand_name + '.sdf')
    ligand_pdb_dir = os.path.join(save_dir, ligand_name + '.pdb')
    # ligand_cif_dir = os.path.join(save_dir, ligand_name + '.cif')
    log_dir = os.path.join(save_dir, ligand_name + '_obabel.log')

    # Remove counterions if any
    ligand_smiles = remove_counterions(ligand_smiles)
    
    # generate sdf file
    obabel(f"""-:{'"'+ligand_smiles+'"'} -O {ligand_sdf_dir} --title {ligand_name} --gen3d \
    --best --minimize --ff {force_field} --steps {maximum_steps} --sd \
    --crit {convergence_criteria} > {log_dir} 2>&1""")

    # sdf to pdb
    obabel(f"-i sdf {ligand_sdf_dir} -o pdb -O {ligand_pdb_dir}")
    # obabel(f"-i sdf {ligand_sdf_dir} -o cif -O {ligand_cif_dir}")

    # check the log file
    with open(log_dir, 'r') as file:
        file_content = file.read()
        if file_content == '1 molecule converted\n':   
            flag = 'successful'
            if verbos:
                print("Conformers saved to:")
                print(ligand_sdf_dir)
                print(ligand_pdb_dir)
        else:
            flag = 'error'
            print(f'There is an error, for {ligand_name}, checkout the log file:')
            print(log_dir)


    # Second attempt with gen2d
    if flag == 'error':
        print("Trying gen2d ...")
        # generate sdf file
        obabel(f"""-:{'"'+ligand_smiles+'"'} -O {ligand_sdf_dir} --title {ligand_name} --gen2d \
        --best --minimize --ff {force_field} --steps {maximum_steps} --sd \
        --crit {convergence_criteria} > {log_dir} 2>&1""")
        # sdf to pdb
        obabel(f"-i sdf {ligand_sdf_dir} -o pdb -O {ligand_pdb_dir}")

        # check the log file
        with open(log_dir, 'r') as file:
            file_content = file.read()
            if file_content == '1 molecule converted\n':   
                flag = 'successful'
                if verbos:
                    print("Conformers saved to:")
                    print(ligand_sdf_dir)
                    print(ligand_pdb_dir)
            else:
                flag = 'error'
                print(f'There is an error, for {ligand_name}, checkout the log file:')
                print(log_dir)
            
    m = stereoisomers(ligand_smiles, check_only = True)
    warning = None
    if m['message'].split(':')[0] == '0':
        warning = 'The information of stereocenters is not specifid in this smiles. Therefore, obabel generates a random conformation. Use `stereoisomers` to generate all possible stereocenters'
        
    if flag == 'successful':
        return {
            'message' : 'SMILES is successfully converted into the 3D struture.',
            'next_steps' : ['Protonation and energy minimization', 'Prepare ligand', 'plot'],
            'ligand_sdf_file' : ligand_sdf_dir,
            'ligand_pdb_file' : ligand_pdb_dir,
            'warning' : warning
        }
    else:
        return {
            'message' : 'SMILES is not convertable into a valid 3D struture',
            'ligand_sdf_file' : None,
            'ligand_pdb_file' : None,
            'warning' : None
        }


def protonate_and_optimize_ligand(input_file: str, 
                                  ph: float|None = 7.2, 
                                  # save_dir: str = FILE_DIR,
                                  force_field: str = 'MMFF94',
                                  convergence_criteria: str = '0.00001',
                                  maximum_steps: int = 10000
                                 ):

    """Add hydrogens to molecules based on a given ph and perfrom energy minimization using obabel"""

    save_dir = FILE_DIR
    
    # PDB2PQR/PropKa for ptotein

    ph_str = str(ph).replace(".", "")             # 7.2 -> 72
    output_file_noH = _add_tag_to_path(input_file, 'noH')
    output_file_ph = _add_tag_to_path(input_file, f'ph{ph_str}')
    output_file_ph_min = _add_tag_to_path(input_file, f'ph{ph_str}_min') 
    log_dir = os.path.join(save_dir, 'obabel_optimize.log')

    if ph:
        # Remove hydrogens
        obabel(f"{input_file} -O {output_file_noH} -d")
        # Add hydrogens with the given ph
        obabel(f"{output_file_noH} -O {output_file_ph} -p {ph}")
        input_file = output_file_ph
    
    # Perform energy minimization
    obabel(f"""{input_file} -O {output_file_ph_min} \
    --best --minimize --ff {force_field} --steps {maximum_steps} --sd \
    --crit {convergence_criteria} > {log_dir} 2>&1""")

    # check the log file
    with open(log_dir, 'r') as file:
        file_content = file.read()
    if '1 molecule converted' in file_content:
        flag = 'successful'
        return {
            'message' : 'Successful!',
            'ph' : ph,
            'saved_file' : output_file_ph_min
        }
    else:
        return {
            'message' : 'Not Successful!',
            'ph' : ph,
            'saved_file' : None,
        }

        
def cxsmiles2smiles(cxsmiles: str, remove_counterions: bool = True):

    """Convert CXSMILES to canonical SMILES"""
    
    # Convert into cononical smiles
    mol = Chem.MolFromSmiles(cxsmiles)
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

    # Remove counterions if any
    if remove_counterions:
        canonical_smiles = remove_counterions(canonical_smiles)

    return canonical_smiles


def stereoisomers(smiles: str, check_only: bool = False):

    """Check and generate different stereoisomers for SMILES"""

    # Check if the information of stereocenters are there
    def _check_stereocenters(smiles):
        return True if smiles.find('@') > -1 else False

    # Extract stereocenters
    def _get_stereoisomers(smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True)
        isomers = list(EnumerateStereoisomers(mol, options=opts))
        smile_set = []
        for isomer in isomers:
            smile_set.append(Chem.MolToSmiles(isomer, isomericSmiles=True, canonical=True))
        return smile_set
        
    if check_only:
        if _check_stereocenters(smiles):
            message = '1:The information of stereocenters is specifid in this smiles. '
        else:
            message = '0:The information of stereocenters is not specifid in this smiles. '
        return {
            'message' : message
        }
    else:
        message = ''
        if _check_stereocenters(smiles):
            message = 'The information of stereocenters is already specifid in this smiles. '
            
        ## Generate all possible isoforms
        smiles_list = [smiles]
        stereoisomeric_smiles = _get_stereoisomers(smiles)
        smiles_list.extend(stereoisomeric_smiles)
        name_list = [f'stereoisomer_{i}' for i in range(len(smiles_list))]

        if len(smiles_list) == 1:
            message += 'No other stereoisomers were found.'
        else:
            message += f'{len(smiles_list) - 1} other stereoisomers was found.'
            
        return {
            'message' : message,
            'stereoisomeric_smiles' : smiles_list,
            'stereoisomeric_names' : name_list 
        }

    
def calculate_ligand_info(smiles: str = None, sdf_file: str = None) -> dict:
    """Extract molecular properties from a ligand provided as SMILES or SDF file."""
    
    if not smiles and not sdf_file:
        raise ValueError("Please provide either a SMILES string or an SDF file path.")
    if smiles and sdf_file:
        raise ValueError("Provide only one input: either SMILES or SDF file, not both.")

    mol = None
    validate_smiles = None
    error = None
    warning = None

    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=True)
            validate_smiles = mol is not None
            warning = "No pH is taken into account, implicit hydrogens based on valence rules and the atom’s formal charge are assigned. Use `protonate_and_optimize_ligand` to protonate at your desired pH"
            if not validate_smiles:
                error = "Invalid SMILES string."
        except Exception as e:
            error = str(e)
            validate_smiles = False

    elif sdf_file:
        try:
            suppl = Chem.SDMolSupplier(sdf_file, removeHs=False, sanitize=True)
            mols = [m for m in suppl if m is not None]
            if not mols:
                error = "No valid molecules found in SDF file."
            else:
                mol = mols[0]  # take the first molecule
                validate_smiles = None  # not applicable
        except Exception as e:
            error = str(e)

    if mol is None:
        return {"validate_smiles": validate_smiles, "error": error}

    # --- Compute properties ---
    props = {
        "message" : "",
        "validate_smiles": validate_smiles,
        "Molecular weight": Descriptors.MolWt(mol),
        "Molecular formula": rdMolDescriptors.CalcMolFormula(mol),
        "clogP": Crippen.MolLogP(mol),
        "Molar refractivity": Crippen.MolMR(mol),
        "Topological polar surface area (TPSA)": rdMolDescriptors.CalcTPSA(mol),
        "Number of atoms": mol.GetNumAtoms(),
        "Number of bonds": mol.GetNumBonds(),
        "Number of rings": mol.GetRingInfo().NumRings(),
        "Number of rotatable bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Formal charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()),
        "Number of heavy atoms": rdMolDescriptors.CalcNumHeavyAtoms(mol),
        "Num of H Acceptors": Lipinski.NumHAcceptors(mol),
        "Num of H Donors": Lipinski.NumHDonors(mol),
        "Number of hydrogens": sum(a.GetTotalNumHs() for a in mol.GetAtoms()),
        "warning" : warning
    }

    return props


def admet_predict(smiles: str) -> dict:
    """ADMET prediction based on ADMET-AI"""
    
    def ADMETModel_init_silent():
    
        import io
        import contextlib
        
        # Suppress torch FutureWarning
        warnings.filterwarnings(
            "ignore",
            category=FutureWarning,
            message="You are using `torch.load` with `weights_only=False`.*",
        )
        
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            model = ADMETModel()
            
        return model
    
    def ADMETModel_predict_silent(model, smiles):
    
        import io
        import contextlib
        
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            preds = model.predict(smiles=smiles)
            
        return preds
    
    
    model = ADMETModel_init_silent()
    preds = ADMETModel_predict_silent(model = model, smiles=smiles)
    
    return preds
    

def smiles_pattern_search(patterns: list, 
                          pattern_type: str = 'smiles',
                          search_name: str = 'my_search'):
    """Serach in PDBbind data for co-crystaizd ligands with a specific patterns given a list of SMILES or SMARTS"""
    
    pdbbind_dir = "../data/PDBbind.csv"
    df_smiles = pd.read_csv(pdbbind_dir)

    # Convert to patterns to mol
    if pattern_type.lower == 'smarts':
        pattern_mols = [Chem.MolFromSmarts(p) for p in patterns]
    else:
        pattern_mols = [Chem.MolFromSmiles(p) for p in patterns]
    
    df_match = pd.DataFrame(columns=['pdb_id', 'smiles', 'pattern'])
    for i, row in df_smiles.iterrows():
        smiles = row['smiles']
        pdb_id = row['pdb_id']
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        for p, pat_mol in zip(patterns, pattern_mols):
            if pat_mol is not None:
                if mol.HasSubstructMatch(pat_mol):
                    df_match.loc[len(df_match)] = [pdb_id, smiles, p]
    if len(df_match):
        file_dir = os.path.join(FILE_DIR, search_name + '.csv')
        message = f"The search result is saved at {file_dir}"
        df_match.to_csv(file_dir)
        # df_match.to_excel("file_dir.xlsx")
    else:
        message = "Nothing found"
        file_dir = None
    
    return {
        "message"  : message,
        "file_dir" : file_dir,
        "next_steps" : "Visualization",
    }



# 2. Protein analysis =====================

def extract_pdb_components(protein_pdb_dir: str):
    """Extract protein chains and native ligands in a pdb file and save in separate files"""

    # Extract protein chains
    protein_extracted = extract_proteins(protein_pdb_dir)
    for chain in protein_extracted.keys():
        tmp = protein_extracted[chain]
        protein_extracted[chain] = os.path.join(os.path.dirname(protein_pdb_dir), f"{tmp}.pdb")

    # Find native mols
    native_ligands = list_ligands(protein_pdb_dir, prnt=False)
    print(f"\n * Native molecules: {native_ligands} \n")
    
    # Extract native mols
    print("Extract native molecules \n")
    ligand_extracted = dict()
    for name in native_ligands:
        tmp = extract_native_ligand_from_pdb(name, protein_pdb_dir)
        ligand_extracted[name] = {
            'SMILES' : tmp['SMILES'],
            'path' : tmp['save_dir']
        }
    return {
        'message' : 'pdb components are extracted!',
        'next_steps' : ['Protonation and energy minimization', 'Prepare protein', 'plot'],
        'proteins' : protein_extracted,
        'ligands'  : ligand_extracted
    }


def protonate_and_optimize_protein(pdb_file: str,
                                   ph: float|None = 7.2,
                                   force_field: str = 'amber',
                                   foldx_repair: bool = False
                                  ):
    """Protonate and optimize protein based on propka, pdb2pqr, and foldX (FoldX may take a long time)"""
    if ph:
        out_ph = protonate_protein(pdb_file, ph, force_field)
        pdb_file = out_ph['saved_file']

    if foldx_repair:
        out_opt = optimize_protein_foldx(pdb_file)
        out_opt['ph'] = ph
        return out_opt
    else:
        return out_ph


# 3. Protein and ligand =====================

def obabel(arg: str = "--help"):
    """Run obabel given a desired argument"""
    os.system(f"obabel {arg}")


def extract_pdb_info(path):
    """Data extractor for key crystallographic/experimental metadata + chains and HETATM info. Works with PDB or mmCIF. Prefers mmCIF metadata."""
    def _ok(x):
        return x and str(x).strip() not in (".", "?", "")
    def _try(vals):
        for v in vals:
            if _ok(v): return str(v).strip()
        return None
    def _get(db, key):
        try:
            return db.find_value(key)
        except Exception:
            return None
    def _first_in_loop(db, cat, item):
        try:
            t = db.find(cat)
            if t and t.ok():
                col = t.find_column(item)
                if col != -1 and t.row_count() > 0:
                    return t.value(0, col)
        except Exception:
            pass
        return None

    ext = os.path.splitext(path)[1].lower()
    st = gemmi.read_structure(path)
    st.setup_entities()
    model = st[0] if len(st) else None

    # --- count atoms, chains, HETATMs --------------------------------------
    total_atoms = 0
    chain_names = []
    het_by_chain = {}
    if model:
        for ch in model:
            chain_names.append(ch.name)
            ligs = defaultdict(int)
            for res in ch:
                total_atoms += len(res)
                # if res.is_water():
                #     continue
                # if getattr(res, "het_flag", " ") != " ":    # get ATOM and HETATM
                if getattr(res, "het_flag", "\x00") == "H":   # get only HETATM
                    ligs[res.name] += 1
            het_by_chain[ch.name] = dict(sorted(ligs.items()))
    chain_names = sorted(set(chain_names))

    # --- fallback cell/space group from coordinates ------------------------
    sg_fallback = getattr(st.spacegroup, "hm", None) if getattr(st, "spacegroup", None) else None
    cell = st.cell
    cell_fallback = {
        "a": round(cell.a, 3), "b": round(cell.b, 3), "c": round(cell.c, 3),
        "alpha": round(cell.alpha, 2), "beta": round(cell.beta, 2), "gamma": round(cell.gamma, 2)
    }

    # --- mmCIF metadata -----------------------------------------------------
    cif = None
    if ext in (".cif", ".mmcif", ".mcif"):
        try:
            cif = gemmi.cif.read(path).sole_block()
        except Exception:
            cif = None

    exp_method = None
    resolution = None
    space_group = None
    cell_dims = None
    radiation_source = None
    wavelength = None
    refinement_software = []
    rfree = None
    rwork = None
    keywords = None
    total_entities = None
    polymer_mw_kda = None

    if cif:
        exp_method = _get(cif, "_exptl.method")
        resolution = _try([_get(cif, "_refine.ls_d_res_high"), _get(cif, "_refine_hist.d_res_high")])
        space_group = _try([_get(cif, "_symmetry.space_group_name_H-M"), _get(cif, "_space_group.name_H-M")])

        a = _get(cif, "_cell.length_a"); b = _get(cif, "_cell.length_b"); c = _get(cif, "_cell.length_c")
        al = _get(cif, "_cell.angle_alpha"); be = _get(cif, "_cell.angle_beta"); ga = _get(cif, "_cell.angle_gamma")
        if all(map(_ok, [a,b,c,al,be,ga])):
            cell_dims = {
                "a": round(float(a), 3), "b": round(float(b), 3), "c": round(float(c), 3),
                "alpha": round(float(al), 2), "beta": round(float(be), 2), "gamma": round(float(ga), 2)
            }

        source = _try([_first_in_loop(cif, "_diffrn_source", "source"),
                       _first_in_loop(cif, "_diffrn_source", "type")])
        beamline = _try([_first_in_loop(cif, "_diffrn_source", "pdbx_synchrotron_beamline"),
                         _first_in_loop(cif, "_synchrotron", "beamline")])
        radiation_source = f"{source} ({beamline})" if source and beamline else source

        wavelength = _try([_first_in_loop(cif, "_diffrn_radiation_wavelength", "wavelength"),
                           _get(cif, "_diffrn_radiation.wavelength")])

        tbl = cif.find("_software")
        if tbl and tbl.ok():
            cls_col = tbl.find_column("classification")
            name_col = tbl.find_column("name")
            for i in range(tbl.row_count()):
                cls_ = tbl.value(i, cls_col) if cls_col != -1 else ""
                name_ = tbl.value(i, name_col) if name_col != -1 else ""
                if name_ and cls_ and any(k in cls_.lower() for k in
                                          ["refinement", "scaling", "reduction", "phasing", "building"]):
                    refinement_software.append(name_.strip())

        rfree = _get(cif, "_refine.ls_R_factor_R_free")
        rwork = _get(cif, "_refine.ls_R_factor_R_work")
        keywords = _try([_get(cif, "_struct_keywords.text"), _get(cif, "_struct_keywords.pdbx_keywords")])
        total_entities = len(st.entities)
        mw = 0
        for e in st.entities:
            if e.entity_type == gemmi.EntityType.POLYMER:
                mw += 110 * e.get_polymer().length()
        polymer_mw_kda = round(mw/1000.0, 2) if mw else None

    # --- fallback parsing for PDB headers -----------------------------------
    elif ext == ".pdb":
        with open(path, "r") as f:
            for line in f:
                if line.startswith("EXPDTA"):
                    exp_method = line[10:].strip()
                elif line.startswith("REMARK   2 RESOLUTION."):
                    # Example: "REMARK   2 RESOLUTION.  3.18 ANGSTROMS."
                    try:
                        resolution = float(line.split()[3])
                    except Exception:
                        pass
                elif line.startswith("REMARK   2") and "ANGSTROMS" in line:
                    # generic catch
                    parts = [p for p in line.split() if p.replace('.', '', 1).isdigit()]
                    if parts:
                        try:
                            resolution = float(parts[0])
                        except Exception:
                            pass

    # --- assemble final dict ------------------------------------------------
    return {
        "Structure and Crystallization": {
            "Experimental Method": exp_method or "N/A",
            "Resolution (Å)": float(resolution) if _ok(resolution) else None,
            "Space Group": space_group or sg_fallback or "N/A",
            "Cell Dimensions": cell_dims or cell_fallback
        },
        "Experimental and Refinement Details": {
            "Radiation Source": radiation_source or "N/A",
            "Wavelength (Å)": float(wavelength) if _ok(wavelength) else None,
            "Refinement Software": refinement_software,
            "R-Free": float(rfree) if _ok(rfree) else None,
            "R-Work": float(rwork) if _ok(rwork) else None
        },
        "Key Biological and Molecular Features": {
            "Keywords": keywords,
            "Entity Information": {
                "Total entity count": total_entities,
                "Total deposited atoms": total_atoms,
                "Polymer molecular weight (kDa, approx.)": polymer_mw_kda
            }
        },
        "Chains": chain_names,
        "HETATM per chain": het_by_chain
    }
    