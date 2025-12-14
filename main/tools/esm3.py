import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import pickle
from esm.sdk import client
from esm.sdk.api import ESMProtein, GenerationConfig, SamplingConfig

FILE_DIR = os.getenv("FILE_DIR")
ESM3_TOKEN = os.getenv("ESM3_TOKEN")

# === utilities ===
def plk2protein(file_path):
    """Load protein object"""
    if os.path.exists(file_path):
        with open(file_path, "rb") as f:
            return pickle.load(f)
    return None

def protein2pkl(protein, file_path):
    """Save protein object"""
    with open(file_path, "wb") as f:
        pickle.dump(protein, f)

def pdb2protein(file_path):
    """Load protein object from pdb file"""
    if os.path.exists(file_path):
        protein = ESMProtein.from_pdb(file_path)
        protein.sequence = None
        return protein
    return None

def protein2pdb(protein, file_path):
    """Load protein object from pdb file"""
    protein.to_pdb(file_path)

def seq2protein(sequence):
    protein = ESMProtein(sequence = sequence)
    return protein

#===============

def run_esm3(
    protein_input: dict,
    task = str,
    # token = ESM3_TOKEN,
    protein_name = 'my_protein',
    model_name = 'esm3-large-2024-03'
    ) -> dict:

    """This function performs all ESM3 analysis, including protein folding, inverse folding, sequence completion, aminoacis annotation, etc."""
    
    token = ESM3_TOKEN
    
    protein_input_keys = ['pkl_file', 'pdb_file', 'sequence', 'sasa', 'function', 'residue_annotations']
    for i in protein_input_keys:
        if i not in protein_input:
            protein_input[i] = ""
    
    # Read protein
    protein = plk2protein(protein_input['pkl_file'])
    
    if protein is None:
        protein = pdb2protein(protein_input['pdb_file'])
    
    if protein is None:
        if len(protein_input['sequence']):
            protein = seq2protein(protein_input['sequence'])
        else:
            raise Exception('Protein is not defined')
    elif len(protein_input['sequence']):
        protein.sequence = protein_input['sequence']
    
    if protein is None:
        raise Exception('Protein is not defined')
    
        
    # Add protein attributes
    if len(protein_input['sasa']):
        protein.sasa = protein_input['sasa']
    if len(protein_input['function']):
        protein.sasa = protein_input['function']
    if len(protein_input['residue_annotations']):
        protein.sasa = protein_input['residue_annotations']

    # Build the model
    model = client(model=model_name, url="https://forge.evolutionaryscale.ai", token=token)
    
    # Run task
    protein_new = model.generate(protein, GenerationConfig(track=task, num_steps=8)) # temperature=0.7
    
    pkl_dir = f'{protein_name}.pkl'
    pkl_dir = os.path.join(FILE_DIR, pkl_dir)
    protein2pkl(protein_new, pkl_dir)
    
    if task == 'structure':
        pdb_dir = f'{protein_name}.pdb'
        pdb_dir = os.path.join(FILE_DIR, pdb_dir)
        protein2pdb(protein_new, pdb_dir)
    else:
        pdb_dir = ''

    protein_output = {
        'pkl_file' : pkl_dir,
        'pdb_file' : pdb_dir,
        'sequence' : protein_new.sequence,
        'sasa'     : protein_new.sasa,
        'function' : protein_new.function_annotations,
        'secondary_structure' : protein_new.secondary_structure,
        'residue_annotations' : '',
        'plddt' : protein_new.plddt,
        'ptm'   : protein_new.ptm,
    }
    
    return protein_output
