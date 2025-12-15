import numpy as np
import json
import string
import requests
import re
import os
import time
from tqdm import tqdm
from Bio.PDB import MMCIFParser, PDBIO
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

af3_path = "/project/alphafold3/"
output_path = os.path.join(af3_path, "gpu_final_output")
ccd_path = "files"
FILES = "files"


def to_base(number, base):
    """
    Convert a number of based 10 to any given base
    """
    if not number:
        return [0]

    digits = []
    while number:
        digits.append(number % base)
        number //= base
    return list(reversed(digits))


def generate_labels(x):
    """
    Generate a vector of ID lables with its length being equal to the number of molecules x.
    x: number of molecules
    return ['A', 'B', ..., 'Z', 'BA', 'BB', ...]
    """
    x = x - 1
    letters = string.ascii_uppercase
    result = to_base(x, 26)
    result = [letters[i] for i in result]  
    return "".join(result)


def get_protein_seq(uniprot_id, verbos=True):
    # Get uniprot sequence
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'sequence' in data:
            protein_sequence = data['sequence']['value']
            if verbos:
                print(f'Sequence of {uniprot_id} is:\n{protein_sequence}\n')
                print(f'Length = {len(protein_sequence)}\n')
            return protein_sequence
        else:
            print(f'Sequence of {uniprot_id} is not availabe in uniprot\n')
    else:
        print(f'Error in uniprot id {uniprot_id}\n')


def find_msa_json_path(name, dir):
    name = name.lower()
    json_files = os.listdir(dir)
    if json_files.count(name):
        path = os.path.join(dir, f"{name}/{name}_data.json")
    else:
        path = None
    return path


def cif2pdb(cif_dir, pdb_dir):
    parser = MMCIFParser()
    structure = parser.get_structure('protein_structure', cif_dir)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_dir)


def delete_project(project_name, af3_path = "/project/alphafold3/"):
    cpu_input_path = os.path.join(af3_path, 'cpu_source_for_msa')
    file_name = f"{project_name}.json"
    file_name = os.path.join(cpu_input_path, file_name)
    os.remove(file_name)


def sdf_pdb_to_enriched_cif(sdf_path, output_cif_path, comp_id="LIG"):

    """
    Example usage:
    sdf_pdb_to_enriched_cif("path/to/ligand.sdf", "path/to/output.cif")
    """
    
    # Load molecule from SDF
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError("Failed to load molecule from SDF.")

    # Generate coordinates if missing
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

    # Get conformer and properties
    conf = mol.GetConformer()
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    weight = Descriptors.MolWt(mol)

    with open(output_cif_path, "w") as f:
        f.write(f"data_{comp_id}\n")
        f.write(f"_chem_comp.id {comp_id}\n")
        f.write("_chem_comp.name ?\n")
        f.write("_chem_comp.type non-polymer\n")
        f.write(f"_chem_comp.formula {formula}\n")
        f.write("_chem_comp.mon_nstd_parent_comp_id ?\n")
        f.write("_chem_comp.pdbx_synonyms ?\n")
        f.write(f"_chem_comp.formula_weight {weight:.2f}\n\n")

        # Atom block
        f.write("loop_\n")
        f.write("_chem_comp_atom.comp_id\n")
        f.write("_chem_comp_atom.atom_id\n")
        f.write("_chem_comp_atom.type_symbol\n")
        f.write("_chem_comp_atom.charge\n")
        f.write("_chem_comp_atom.pdbx_model_Cartn_x_ideal\n")
        f.write("_chem_comp_atom.pdbx_model_Cartn_y_ideal\n")
        f.write("_chem_comp_atom.pdbx_model_Cartn_z_ideal\n")

        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom_id = f"A{i+1}"
            symbol = atom.GetSymbol()
            charge = atom.GetFormalCharge()
            f.write(f"{comp_id} {atom_id} {symbol} {charge} {pos.x:.3f} {pos.y:.3f} {pos.z:.3f}\n")

        # Bond block
        f.write("\nloop_\n")
        f.write("_chem_comp_bond.atom_id_1\n")
        f.write("_chem_comp_bond.atom_id_2\n")
        f.write("_chem_comp_bond.value_order\n")
        f.write("_chem_comp_bond.pdbx_aromatic_flag\n")

        for bond in mol.GetBonds():
            a1 = f"A{bond.GetBeginAtomIdx() + 1}"
            a2 = f"A{bond.GetEndAtomIdx() + 1}"
            order = bond.GetBondTypeAsDouble()
            aromatic = 'Y' if bond.GetIsAromatic() else 'N'
            f.write(f"{a1} {a2} {order} {aromatic}\n")

    print(f"Enriched CIF file saved to: {output_cif_path}")


def make_af3_json(project_name, 
                  protein_dict, 
                  ligand_dict = {},
                  af3_path = "/project/alphafold3/", 
                  version=2,
                  seed_max = 10,
                  verbos = True
                 ):

    # Variables
    msa_json_path = os.path.join(af3_path, 'cpu_target_msa')
    cpu_input_path = os.path.join(af3_path, 'cpu_source_for_msa')
    if verbos:
        print(f'msa json path : {msa_json_path}')
        print(f'cpu input path : {cpu_input_path}')

    protein_name = protein_dict['protein_name']
    uniprot_id = protein_dict['uniprot_id']
    sequence = protein_dict['sequence']
    n_chain = protein_dict['n_chain']
    msa_json = [find_msa_json_path(i, msa_json_path) for i in protein_dict['protein_name']]
    if verbos:
        print(f'MSA json file : {msa_json}')

    # Get ids
    n_molecules = np.sum(n_chain) + np.sum([ligand_dict[k]['count'] for k in ligand_dict.keys()])
    n_molecules = int(n_molecules)
    ids = [generate_labels(i + 1) for i in range(n_molecules)]
    if verbos:
        print(f'number of molecules : {n_molecules}')

    # Make sequences dict for protein
    sequences = []
    m = 0
    for i, (id, seq, msa) in enumerate(zip(uniprot_id, sequence, msa_json)):
        # Get the number of chains
        n = n_chain[i]
        # Get sequence
        if seq == None:
            seq = get_protein_seq(id, verbos)
        # Get MSA if exists
        if msa == None:
            unpairedMsa = None
            pairedMsa = None
        else:
            with open(msa, 'r') as file:
                msa = json.load(file)
            unpairedMsa = msa['sequences'][0]['protein']['unpairedMsa']
            pairedMsa = msa['sequences'][0]['protein']['pairedMsa']
            print(f'MSA already exists for protein {protein_name}')
        # Fill in the dictionary
        seq_i = {
            'protein': {
                'id' : ids[m : n+m],
                'sequence' : seq,
                'unpairedMsa' : unpairedMsa,
                'pairedMsa' : pairedMsa,
            }
        }
        sequences.append(seq_i)
        m = m+n

    # Make sequences dict for ligands
    ccd_all = ""
    if len(ligand_dict):
        for ligand in ligand_dict.keys():
            n = ligand_dict[ligand]['count']
            type_lig = ligand_dict[ligand]['type']
            
            seq_i = {
                'ligand': {
                    'id' : ids[m : n+m],
                    type_lig : [ligand]
                }
            }
            sequences.append(seq_i)
            m = m+n
    
        # add ccd
        if type_lig == 'ccdCodes':
            with open(os.path.join(ccd_path, f"{ligand}.cif"), 'r') as file:
                ccd = file.read()
            ccd_all = ccd_all + ccd

    # Make final query dict
    if len(ligand_dict): 
        query_dict = {
            'name' : project_name,
            'sequences' : sequences,
            'modelSeeds': list(range(1, seed_max+1)),
            'dialect': "alphafold3",
            'userCCD' : ccd_all,
            'version': 3
        }
    else:
        query_dict = {
            'name' : project_name,
            'sequences' : sequences,
            'modelSeeds': list(range(1, seed_max+1)),
            'dialect': "alphafold3",
            'version': 2
        }
        

    # Write out the json file
    file_name = f"{project_name}.json"
    target_path = os.path.join(cpu_input_path, file_name)
    with open(target_path ,'w') as file:
        json.dump(query_dict, file, indent=4, separators=(',', ': '))

def run_af3(
    protein_sequence: list,
    n_chain_per_sequence: list|None = None,
    uniprot_id: list|None = None,
    ligand_sdf_dir : str|None = None,
    n_ligand : int = 1,
    project_name: str = 'my_project'
):
    
    seed_max = 1
    seed = 1

    if n_chain_per_sequence in None:
        n_chain_per_sequence = [1] * len(protein_sequence)
    if uniprot_id in None:
        uniprot_id = ['None'] * len(protein_sequence)
                
    protein_dict = {
    'protein_name' : [project_name],
    'uniprot_id'   : uniprot_id,
    'sequence'     : protein_sequence,
    'n_chain'      : n_chain_per_sequence,
    }

    if ligand_sdf_dir:
        ligand_sdf_dir = os.path.join(ccd_path, ligand_sdf_dir)
        ligand_cif_dir = ligand_sdf_dir.rsplit(".", 1)[0] + ".cif"
        sdf_pdb_to_enriched_cif(ligand_sdf_dir, ligand_cif_dir)
        ligand_dict = {
        'LIG' : {
            'type' : 'ccdCodes',
            'count': n_ligand,
        }
        }
    else:
        ligand_sdf = {}
        
    make_af3_json(project_name, protein_dict, ligand_dict, af3_path, seed_max = seed_max, verbos=True)

    # Wait for AF3
    protein_cif_dir = os.path.join(output_path, f"{project_name}/seed-{seed}_sample-0/{project_name}_seed-{seed}_sample-0_model.cif")
    while(not os.path.exists(protein_cif_dir)):
        time.sleep(2)
    print('The result is ready at:')
    print(protein_cif_dir)
    
    # Convert CIF to PDB
    protein_pdb_dir = os.path.join(FILES, project_name + 'pdb')
    cif2pdb(protein_cif_dir, protein_pdb_dir)

    return {
        'message' : f"The result is ready",
        'output_file' : protein_pdb_dir
    }