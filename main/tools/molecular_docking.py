import os
import time
import pandas as pd
import pymol2
from .dock_utilities import *
from .generate_graph import *
from .protein_ligand_basics import *

FILE_DIR = os.getenv("FILE_DIR")

## ============== Helper functions ======================
def binding_box(pdb_file, constraint_file, scale: float = 1.1):
    """Make a binding box around the given pdb file"""
    GB = GridBox(pdb_file)
    center, bxsize = GB.labox(scale=scale)
    with open(constraint_file, 'w') as cfg:
        cfg.write(f'center_x = {center[0]}\n')
        cfg.write(f'center_y = {center[1]}\n')
        cfg.write(f'center_z = {center[2]}\n')
        cfg.write(f'\n')
        cfg.write(f'size_x = {bxsize[0]}\n')
        cfg.write(f'size_y = {bxsize[1]}\n')
        cfg.write(f'size_z = {bxsize[2]}\n')


def sdf2pdbqt(sdf_file: str,
              ligand_name: str,   # Later this should by defult be the same name as sdf
             ):
    """Given a ligand sdf file, prepare the pdbqt file required for docking with VINA family"""

    pdbqt_file = os.path.join(FILE_DIR, ligand_name + '.pdbqt')
    os.system(f"mk_prepare_ligand.py -i {sdf_file} -o {pdbqt_file} --rigid_macrocycles --bad_charge_ok") # ./

    if os.path.exists(pdbqt_file):        
        # Convert UNL to UNL A
        pdbqt_add_chid(pdbqt_file)
        return {
            'message' : 'The ligand pdbqt is prepared',
            'ligand_pdbqt_file' : pdbqt_file,
        }
    else:
        return {
            'message' : 'The ligand pdbqt file cannot be generated!',
            'ligand_pdbqt_file' : None,
        }

def make_protein_ligand_complex(
    protein_pdb_dir : str, 
    ligand_sdf_dir : str,
    complex_pdb_dir: str = 'my_complex.pdb',
    add_all_h : bool = False,
    ligand_name: str = 'LIG',
    ligand_chain: str = 'L'
):

    """Generate protein ligand complex file"""
    
    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        cmd.load(protein_pdb_dir, "protein")
        cmd.load(ligand_sdf_dir, "ligand")
    
        if add_all_h:
            cmd.h_add("protein")
            cmd.h_add("ligand")
    
        # Give ligand clean IDs
        cmd.alter("ligand", f"chain='{ligand_chain}'")
        cmd.alter("ligand", f"resn='{ligand_name}'")
        cmd.alter("ligand", "resi='1'")
    
        # Rename ligand atoms deterministically by adding numbers
        lig = cmd.get_model("ligand")
        for i, atm in enumerate(lig.atom, 1):
            cmd.alter(f"ligand and index {atm.index}", f"name='{atm.symbol}{i}'")
        cmd.sort()
    
        # cmd.save(ligand_sdf_dir_H, "ligand")
    
        # Create complex and write rich CONECTs
        cmd.create("complex", "protein or ligand")
        cmd.set("pdb_conect_all", 1)
        cmd.save(complex_pdb_dir, "complex")

        # Remove TER before ligand
        with open(complex_pdb_dir) as f:
            lines = f.readlines()
        
        with open(complex_pdb_dir, "w") as f:
            for line in lines:
                if line.startswith("TER"):
                    continue
                f.write(line)


def extract_native_ligand_from_pdb(
    name: str,
    input_file: str,
    retain_atom_numbers : bool = False,
    # save_dir: str = FILE_DIR,
):

    """
    Extract native ligand from pdb file and save it as pdb, sdf. Also obtain its SMILES. 
    (Update: there is an option to keep the atom numbers from pdb. This is needed for complex generation.)
    """

    save_dir = FILE_DIR
    
    os.makedirs(save_dir, exist_ok=True)
    ntv_sdf = os.path.join(save_dir, f"{name}.sdf")
    ntv_pdb = os.path.join(save_dir, f"{name}.pdb")

    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        cmd.load(input_file, "prot")
        cmd.select("lig", f"resn {name} and not polymer")
        # ntv_smiles = cmd.get_smiles("lig")

        cmd.save(ntv_sdf, "lig")
        
        if retain_atom_numbers:
            cmd.set("pdb_retain_ids", 1)
        cmd.save(ntv_pdb, "lig")
        if retain_atom_numbers:
            cmd.set("pdb_retain_ids", 0)

    # Read SDF with RDKit and get SMILES
    try:
        suppl = Chem.SDMolSupplier(ntv_sdf, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        ntv_smiles = Chem.MolToSmiles(mol, isomericSmiles=True) if mol is not None else None
    except:
        ntv_smiles = None 

    ntv_smiles = ""
    return {
        "Message": "Done!",
        "SMILES": ntv_smiles,
        "save_dir_sdf": ntv_sdf,
        "save_dir_pdb": ntv_pdb,
    }

## ============== Main functions ======================

def prepare_protein(pdb_file: str, 
                    protein_name: str = 'my_protein',
                    partial_charge: str = 'gasteiger',
                    project_name: str = 'my_docking',
                    docking_method: str = 'smina',
                    constraint_pdb_file: str | None = None,
                    constraint_box_scale: float = 1.1,
                    ph: float|None = None,
                    force_field: str = 'AMBER',
                    foldx_repair: bool = False
                   ):
    """Given a protein pdb file, prepare the protein query table. It performs protonation and energy minimization if specified. If method is a VINA family, it makes the binding box, removes water, adds partial charges and generates the pdbqt file."""

    pdbqt_file = os.path.join(FILE_DIR, f"{project_name}_target.pdbqt")
    constraint_file = os.path.join(FILE_DIR, f"{project_name}_constraint_file")
    df_protein_dir = os.path.join(FILE_DIR, f"{project_name}_df_protein.csv")

    if ph or foldx_repair:
        protein_ph = protonate_and_optimize_protein(pdb_file, ph, force_field, foldx_repair)
        pdb_file = protein_ph['saved_file']
        
    if docking_method in ['vina', 'smina', 'gnina']:
        # Create constraint file for blind search            
        pdb_file_box = constraint_pdb_file if constraint_pdb_file else pdb_file
        binding_box(pdb_file_box, constraint_file, constraint_box_scale)
    
        # Parameterise protein with Gasteiger charges and generate pdbqt
        # partial_charge = 'gasteiger' # ['gasteiger', 'mmff94', 'qeq', 'qtpie']
        # -xr: do not generate torsion tree (protein is ridgid)
        # Hint: when creating the pdbqt file, obabel does not work with pdb file that have het molecules.
        os.system(f"obabel {pdb_file} -xr -O {pdbqt_file} --partialcharge {partial_charge} > /dev/null 2>&1")
        # Addition of Partial Charge and polar hydrogens, and removal of non-polar hydrogens
        response = os.system(f"mk_prepare_receptor.py --pdbqt {pdbqt_file} -o {pdbqt_file} --skip_gpf")
        if response != 0:
            return {
            'message' : 'The protein could not be prepared',
            'constraint_file' : None,
            'protein_pdbqt_file' : None,
            'protein_df_file': None
            }
            
    # Add to df_protein
    df_protein = pd.DataFrame(
        {
            'protein_name' : [protein_name],
            'protein_sequence' : [''],
            'protein_pdb_file' : [pdb_file],
            'protein_pdbqt_file' : [pdbqt_file],
            'constraint_file' : [constraint_file]  
        }
    ) 
    df_protein.to_csv(df_protein_dir, index=False)

    return {
        'message' : 'Protein is prepared',
        'next_steps' : ['prepare_ligand', 'make_query_table', 'plot with box'],
        'constraint_file' : constraint_file,
        'protein_pdbqt_file' : pdbqt_file,
        'protein_df_file': df_protein_dir
    }


def prepare_ligand(smiles: str = None,
                   sdf_file: str = None,
                   ligand_name: str = 'my_ligand',
                   project_name: str = 'my_docking',
                   docking_method: str = 'smina',
                   ph: float|None = None,
                   force_field: str = 'MMFF94',
                   convergence_criteria: str = '0.00001',
                   maximum_steps: int = 10000,
                  ):
    """Given a ligabd SMILES or sdf file, prepare the ligand query table. It performs protonation and energy minimization if specified. If method is a VINA family, adds partial charges and generates the pdbqt file."""

    # at least either sdf_file or smiles should be given
    
    df_ligand_dir = os.path.join(FILE_DIR, f"{project_name}_{ligand_name}_df_ligand.csv")
    warning = None
    
    # smiles to 3d
    if smiles:
        ligand_3d = smiles_to_3d(ligand_smiles = smiles, 
                                 ligand_name = ligand_name, 
                                 # save_dir = FILE_DIR,
                                 force_field = force_field,
                                 convergence_criteria = convergence_criteria,
                                 maximum_steps = maximum_steps)
        if ligand_3d['ligand_sdf_file']:
            ligand_sdf_dir = ligand_3d['ligand_sdf_file']
            ligand_pdb_dir = ligand_3d['ligand_pdb_file']
            warning = ligand_3d['warning']
        else:
            ligand_sdf_dir = sdf_file
            ligand_pdb_dir = None
    else:
        ligand_sdf_dir = sdf_file
        ligand_pdb_dir = None

    if not ligand_sdf_dir:
        return {
        'message' : 'Ligand sdf file is not available',
        'ligand_df_file' : None
    }

    # Protonation and optimization
    if ph:
        ligand_ph = protonate_and_optimize_ligand(input_file = ligand_sdf_dir, 
                                                  ph = ph, 
                                                  # save_dir = FILE_DIR,
                                                  force_field = force_field,
                                                  convergence_criteria = convergence_criteria,
                                                  maximum_steps = maximum_steps)
        ligand_sdf_dir = ligand_ph['saved_file']
        
    
    if docking_method in ['vina', 'smina', 'gnina']:
        tmp = sdf2pdbqt(ligand_sdf_dir, ligand_name)
        if tmp['ligand_pdbqt_file']:
            ligand_pdbqt_dir = tmp['ligand_pdbqt_file']
            df_ligand = pd.DataFrame(
                {
                    'ligand_name' : [ligand_name], 
                    'ligand_smiles' : [smiles], 
                    'ligsnd_sdf_file' : [ligand_sdf_dir],
                    'ligand_pdb_file' :  [ligand_pdb_dir],
                    'ligand_pdbqt_file' : [ligand_pdbqt_dir]
                }
            )
        else:
            df_ligand = None
            
    elif docking_method == 'diffdock':
        df_ligand = pd.DataFrame(
                {
                    'ligand_description' : ligand_sdf_dir
                }
            )     

    if df_ligand is not None:
        df_ligand_dir = os.path.join(FILE_DIR, f"{project_name}_{ligand_name}_df_ligand.csv")
        df_ligand.to_csv(df_ligand_dir, index=False)
        return {
            'message' : 'ligand dataframe saved',
            'next_steps' : ['prepare_protein', 'make_query_table'],
            'ligand_df_file' : df_ligand_dir,
            'warning' : warning
        }
    else:
        df_ligand_dir = None
        return {
        'message' : 'ligand dataframe not prepared',
        'ligand_df_file' : None
        }


def make_query_table(protein_df_file: str,
                     ligand_df_file: str,
                     project_name: str = 'my_docking',
                     docking_method: str = 'smina',
                    ):
    """Make query table for molecular docking given a single protein df file and single of multiple ligabds df files"""

    df_protein = pd.read_csv(protein_df_file)
    df_ligand = pd.read_csv(ligand_df_file)
    
    df_protein_repeated = pd.concat([df_protein] * len(df_ligand), ignore_index=True)
    df_combined = pd.concat([df_protein_repeated, df_ligand.reset_index(drop=True)], axis=1)

    if docking_method == 'diffdock':
        df_combined['complex_name'] = df_combined['protein_name'] + '_' + df_combined['ligand_name']
        df_combined = df_combined[['complex_name', 'protein_pdb_file', 'ligsnd_sdf_file', 'protein_sequence']]
        df_combined.columns = ['complex_name', 'protein_path', 'ligand_description', 'protein_sequence']
            
    path = os.path.join(FILE_DIR, f"{project_name}_docking_query_{docking_method}.csv")
    df_combined.to_csv(path, index=False)

    return {
        'message': f"The docking query table for the project {project_name} for method {docking_method} is generated",
        'save_dir' : path,
        'next_steps' : ['run_molecular_docking'],
    }

# Perform Docking and get ligand pose sdf files

def run_molecular_docking(query_table_dir: str,
                          docking_method: str = 'smina',
                          exhaustiveness: str = '16',
                          project_name: str = 'my_docking',
                         ):
    """Perfom molecular docking using VINA, SMINA, GNINA or DiffDock"""

    # Parameters for Vina based methods
    n_cpu = '20'
    # exhaustiveness = '16' #[4, 8, 16, 32, 64, 128, 256]

    result_dir = os.path.join(FILE_DIR, project_name)
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    
    # Set directories
    df = pd.read_csv(query_table_dir, dtype=str)

    start = time.time()
    for i in range(len(df)):

        if method == 'diffdock':
            inference_steps = 10
            complex_name = df['protein_name'][i]
            result_dir_diffdock = os.path.join(result_dir, f'dock_{method}_{protein_name}')
            os.system(f"python -m ../software/DiffDock/inference --protein_ligand_csv {df_docking_diffdock} --out_dir {result_dir_diffdock} --inference_steps {inference_steps} --samples_per_complex 40 --actual_steps 18 --no_final_step_noise")
    
        else:
    
            protein_name = df['protein_name'][i]
            ligand_name = df['ligand_name'][i]
            dock_pdbqt_dir = os.path.join(result_dir, f'dock_{docking_method}_{protein_name}_{ligand_name}.pdbqt')
    
            # if os.path.exists(dock_pdbqt_dir):
            #     continue
                
            query_dict = {
                'protein_pdbqt_dir' : df['protein_pdbqt_file'][i],
                'ligand_pdbqt_dir' : df['ligand_pdbqt_file'][i],
                'constraint_file_dir' : df['constraint_file'][i],
                'output_pdbqt_dir' : dock_pdbqt_dir,
            }
                                                  
            run_docking(query_dict, method = docking_method, n_cpu = n_cpu, exh = exhaustiveness)
            # !echo {i} >> output.log
    
            # Extract pose sdf files
            ligand_pose_pdb_dir = os.path.join(result_dir, f'{docking_method}_{protein_name}_{ligand_name}_pose_.sdf')
            os.system(f'obabel {dock_pdbqt_dir} -O {ligand_pose_pdb_dir} -m > /dev/null 2>&1')
    
    end = time.time()
    time_min = (end - start) / 60

    try:
        os.system(f"mv output.log {FILE_DIR}")
    except:
        pass

    # get report
    scores = extract_score(dock_pdbqt_dir, docking_method)
    if scores:
        message = 'The docking result is ready'
    else:
        message = 'Could not find any conformations completely within the search space.'
    
    return {
        'message' : message,
        'next_steps' : ['protein_ligand_interaction'],
        'docked ligand sdf files' : ligand_pose_pdb_dir,
        'docking score' : scores,
        'time lapsed' : time_min,
    }


def get_protein_ligand_interaction(pdb_file: str,
                                   sdf_file: str,
                                   complex_name: str = 'my_complex',
                                   pocket_threshold: float|int = 6      # in Ångström
                                  ):

    """Calculate protein-ligand interaction given pdb file of a protein and sdf file of a ligad. In addition, generate: complex pdb file, protein pocket, and protein-ligand graph"""
    
    # Define file names to be saved
    complex_pdb_file = os.path.join(FILE_DIR, complex_name + '_complex.pdb')
    protein_pocket_file = os.path.join(FILE_DIR, complex_name + '_protein_pocket.csv')
    interaction_table_file = os.path.join(FILE_DIR, complex_name + '_interaction_table.csv')
    node_tabe_file = os.path.join(FILE_DIR, complex_name + '_node_list.csv')
    edge_list_file = os.path.join(FILE_DIR, complex_name + '_edge_list.csv')

    # Make complex pdb file
    make_protein_ligand_complex(pdb_file, sdf_file, complex_pdb_file)
    
    # Extract ligand from pdb
    tmp = extract_native_ligand_from_pdb(name = 'LIG', 
                                         input_file = complex_pdb_file, 
                                         retain_atom_numbers=True)
    lig_pdb_file = tmp['save_dir_pdb']
    lig_sdf_file = tmp['save_dir_sdf']
    
    # Get ligand atoms
    df_lignd_atom = ligand_atoms(lig_pdb_file)
    lig_atm_map = {i : j for i, j in zip(df_lignd_atom['index_in_pdb'], df_lignd_atom['index'])}
    
    # Get ligand bonds
    df_ligand_edge = ligand_bonds(lig_sdf_file)
    
    # Get protein pocket
    df_pocket = protein_pocket(complex_pdb_file, ligand_id = 'LIG', threshold = pocket_threshold)
    df_pocket['index_in_pocket'] = df_pocket.index + 1
    prot_res_map = {i : j for i, j in zip(df_pocket['index'], df_pocket['index_in_pocket'])}
    
    # Get protein residue connections in the pocket
    df_protein_edge = protein_residue_interaction(df_pocket)
    df_protein_edge['idx1_in_pocket'] = df_protein_edge['idx1'].map(prot_res_map)
    df_protein_edge['idx2_in_pocket'] = df_protein_edge['idx2'].map(prot_res_map)
    
    # Get protein-ligand interactions
    df_lig_prot = protein_ligand_interaction(complex_pdb_file, ligand_id=['LIG'])
    
    # In df_lig_prot, map pdb indices to new indices starting from 1
    if len(df_lig_prot):
        df_lig_prot['prot_res_in_pocket'] = df_lig_prot['prot_res_num'].map(prot_res_map)
        lig_atm_index = []
        for i in df_lig_prot['ligand_idx']:
            if isinstance(i, list):
                tmp = [lig_atm_map[j] for j in i]
                tmp = tmp[0] if len(tmp) == 1 else tmp
            else:
                tmp = lig_atm_map[i]
            lig_atm_index.append(tmp)
        df_lig_prot['lig_atm_index'] = lig_atm_index

        # Cobmine node features
    df_nodes = pd.DataFrame(columns=['index', 'molecule', 'symbol', 'x', 'y', 'z'])
    
    for _, row in df_lignd_atom.iterrows():
        df_nodes.loc[len(df_nodes)] = [
            row['index'],
            'ligand',
            row['symbol'],
            row['x'],
            row['y'],
            row['z'],
        ]
    
    for _, row in df_pocket.iterrows():
        df_nodes.loc[len(df_nodes)] = [
            row['index_in_pocket'],
            'protein',
            row['residue'],
            row['x'],
            row['y'],
            row['z'],
        ]
    
    
    # Cobmine edges
    df_edges = pd.DataFrame(columns=['index1', 'index2', 'type', 'feature'])
    
    for _, row in df_ligand_edge.iterrows():
        df_edges.loc[len(df_edges)] = [
            row['atom1_idx'],
            row['atom2_idx'],
            'ligand-ligand',
            row['bond_type'].lower(),
        ]
    
    for _, row in df_protein_edge.iterrows():
        df_edges.loc[len(df_edges)] = [
            row['idx1_in_pocket'],
            row['idx2_in_pocket'],
            'protein-protein',
            'peptide',
        ]
    
    for _, row in df_lig_prot.iterrows():
        if row['lig_atm_index'] != []:
            df_edges.loc[len(df_edges)] = [
                row['lig_atm_index'],
                row['prot_res_in_pocket'],
                'ligand-protein',
                row['feature'].lower(),
            ]

    # Save tables
    df_pocket.to_csv(protein_pocket_file, index=False)
    df_lig_prot.to_csv(interaction_table_file, index=False)
    df_nodes.to_csv(node_tabe_file, index=False)
    df_edges.to_csv(edge_list_file, index=False)

    return {
        'Message' : 'Done!',
        'complex_pdb_file' : complex_pdb_file,
        'interaction_table_file' : interaction_table_file,
        'protein_pocket_file' : protein_pocket_file,
        'node_tabe_file' : node_tabe_file,
        'edge_list_file' : edge_list_file,
        'next_steps' : ['read interaction_table_file', 'read protein_pocket_file', 'plot complex file by adding a different style for the interacting residues']
    }