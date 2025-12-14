import os
import json
import requests
import pandas as pd

FILE_DIR = os.getenv("FILE_DIR")

def uniprot2seq(uniprot_id: str) -> str:
    """Get the sequence of the protein given the uniprot id"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'sequence' in data:
            protein_sequence = data['sequence']['value']
            return protein_sequence
        else:
            return f'Sequence of {uniprot_id} is not availabe in uniprot'
    else:
        return f'Error in uniprot id {uniprot_id}'


def get_pdb(pdb_id: str) -> dict:
    """Download the pdb file the protein given the pdb id"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        path = os.path.join(FILE_DIR, f"{pdb_id}.pdb")
        with open(path, "wb") as f:
            f.write(response.content)
        message = f"PDB file of {pdb_id} protein downloaded successfully at {path}."
        saved_file = path
    else:
        message = f"Failed to download file. HTTP status code: {response.status_code}"
        saved_file = None
    return {
        'message' : message,
        'saved_file' : saved_file,
        'next_steps' : ['extract_pdb_components', 'extract_pdb_info (e.g. resolution)', 'protonate_and_optimize_protein', 'prepare_protein']
    }


def read_csv(csv_file: str) -> str:
    """Read a csv file"""
    df = pd.read_csv(csv_file)
    return df.to_csv(index=False)

def write_text_file(txt: str, file_name: str) -> None:
    """Write plain text file; file_name should include the extension.*"""
    file_name = os.path.join(FILE_DIR, file_name)
    with open(file_name, "w", encoding="utf-8") as f:
        f.write(txt)

def read_text_file(file_name: str) -> str:
    """Read plain text file; file_name should include the extension.*"""
    file_name = os.path.join(FILE_DIR, file_name)
    with open(file_name, "r", encoding="utf-8") as f:
        return f.read()