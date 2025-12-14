## Automatic Docking Utilities

# Internal Modules
import os
import glob
import time
import shutil
import re

# External modules
import py3Dmol

# Data-related
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm

# Docking-related
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS, AllChem, Draw
# from spyrmsd import io, rmsd

# Binding interaction-related
from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport

RDLogger.DisableLog('rdApp.warning')
##############################
# AA Consntant and Bond Colour Dictionary

# Kyte and Doolittle Hydropathy Scale (1982)
AA_HB = {'ALA':  1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS':  2.5,
         'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE':  4.5,
         'LEU':  3.8, 'LYS': -3.9, 'MET':  1.9, 'PHE':  2.8, 'PRO': -1.6,
         'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL':  4.2}

# University of Calgary PI Scale
AA_PI = {'ALA':  6.0, 'ARG': 10.76, 'ASN': 5.41, 'ASP': 2.77, 'CYS': 5.07,
         'GLN': 5.65, 'GLU':  3.22, 'GLY': 5.97, 'HIS': 7.59, 'ILE': 6.02,
         'LEU': 5.98, 'LYS':  9.74, 'MET': 5.74, 'PHE': 5.48, 'PRO':  6.3,
         'SEC': 5.68, 'SER':  5.68, 'THR':  5.6, 'TRP': 5.89, 'TYR': 5.66,
         'VAL': 5.96}

BOND_COL = {'HYDROPHOBIC': ['0x59e382', 'GREEN'],
            'HBOND': ['0x59bee3', 'LIGHT BLUE'],
            'WATERBRIDGE': ['0x4c4cff', 'BLUE'],
            'SALTBRIDGE': ['0xefd033', 'YELLOW'],
            'PISTACKING': ['0xb559e3', 'PURPLE'],
            'PICATION': ['0xe359d8', 'VIOLET'],
            'HALOGEN': ['0x59bee3', 'LIGHT BLUE'],
            'METAL':['0xe35959', 'ORANGE']}


##############################
# Functions to clean and extract subunits of proteins in the pdb file

def oupt_parse(inpt_file: str) -> tuple:
    name = os.path.basename(inpt_file).split('.')[0]
    path = os.path.dirname(inpt_file)
    dirn = os.path.basename(path)
    return name, dirn, path

def ter_line(asn: str, resn: str, chid: str, resi: str) -> str:
    return f'TER   {asn}      {resn} {chid}{resi:>4} {" "*54}\n'

def extract_entity(inpt_file: str, oupt_file: str, keywords: list) -> None:
    oupt_name, oupt_dirn, oupt_path = oupt_parse(oupt_file)
    extract = []
    headers = ['ATOM', 'HETATM']
    with open(inpt_file, 'r') as inpt:
        for line in inpt:
            record = line.split()[0]
            if all(keyword in line for keyword in keywords) and \
               any(header in record for header in headers):
                extract.append(line)

    assert extract, 'Expected \'keywords\' in PDB lines; found 0 lines'

    with open(oupt_file, 'w') as oupt:
        for line in extract:
            oupt.write(line)
    # print(f'+ Entity extracted: {oupt_name}.pdb > {oupt_dirn} folder')

def extract_chains(inpt_file: str) -> None:
    inpt_name, oupt_dirn, oupt_path = oupt_parse(inpt_file)
    with open(inpt_file, 'r') as inpt:
        data = inpt.readlines()

    chid = sorted(set(line[21:22] for line in data))
    chid_leng = len(chid)
    chid_list = ', '.join(chid)
    print(f'+ Chains detected: {chid_leng} ({chid_list})')

    file_names = {}
    for id_ in chid:
        oupt_name = inpt_name + '_' + id_
        oupt_file = os.path.join(oupt_path, oupt_name + '.pdb')
        with open(oupt_file, 'w') as oupt:
            for line in data:
                if line[21:22] in id_:
                    oupt.write(line)
                    asn = f'{int(line[6:11])+1:>5}'
                    resn = line[17:20]
                    resi = line[22:26]
            oupt.write(ter_line(asn, resn, id_, resi))
            oupt.write('END')
        print(f'+ Chain extracted: {oupt_name}.pdb')
        file_names[id_] = oupt_name
    print(f'+ Whole Molecule: {inpt_name}.pdb')
    file_names['whole'] = inpt_name
    return file_names

def extract_proteins (inpt_file):
    # Separet the protein from the native ligand
    extract_entity(inpt_file, inpt_file[:-4] + "_Protein.pdb", ['ATOM'])
    # Separet protein subunits
    file_names = extract_chains(inpt_file[:-4] + "_Protein.pdb")
    return file_names
    
def list_ligands(inpt_file, prnt=True):
    chids = []
    with open(inpt_file, 'r') as inpt:
        for line in inpt:
            if (line[0:6] == "HETATM"):
                chids.append(line[17:20])
    
    chids = np.unique(chids)
    chids = [re.split('^ ', c)[-1] for c in chids]
    if (prnt):
        print(f'Native ligands: {chids}')
    else:
        return(chids)

def get_molblock(keyword: str, NTV_FLD) -> str:
    # url_path = 'https://files.rcsb.org/ligands/download/' + keyword + '_model.sdf'
    url_path = 'https://files.rcsb.org/ligands/download/'+ keyword +'_ideal.sdf'
    sdf_file = os.path.join(NTV_FLD, keyword + '.sdf')
    os.system(f'wget {url_path} -O {sdf_file} -q')
    molblock = [mol for mol in  Chem.SDMolSupplier(sdf_file) if mol is not None][0]
    os.remove(sdf_file)
    return molblock

def correct_bond_order(inpt_list: list, temp: Chem.rdchem.Mol) -> None:
    for inpt_file in inpt_list:
        targ = AllChem.MolFromPDBFile(inpt_file)
        cmol = AllChem.AssignBondOrdersFromTemplate(temp, targ)
        pdbb = Chem.MolToPDBBlock(cmol, flavor=4)
        with open(inpt_file, 'w') as oupt:
            oupt.write(pdbb)


def extract_native_ligand(name, inpt_file, save_dir):
    Keyword = name
    true_id = Keyword[-3:] if len(Keyword) > 3 else Keyword
    
    ntv_pdb = Keyword.upper() + '.pdb'
    ntv_pdb_nFile = os.path.join(save_dir, ntv_pdb)
    extract_entity(inpt_file, ntv_pdb_nFile, [Keyword, 'HETATM'])
    extract_chains(ntv_pdb_nFile)
    ntv_nFiles = sorted(glob.glob(save_dir + '/' + Keyword + '_*.pdb'))
    ntv_smiles = get_molblock(true_id, save_dir)
    try:
        correct_bond_order(ntv_nFiles, ntv_smiles)
    except:
        print("Error in correct_bond_order")
    return(ntv_smiles)


##############################
# Grid Box Calculation Methods

class GridBox:

    ranges = tuple[list[float], list[float], list[float]]
    coords = tuple[float, float, float]
    center_bxsize = tuple[tuple[float, float, float], tuple[float, float, float]]

    def __init__(self, inpt_file: str) -> None:
        self.inpt = open(inpt_file, 'r')
        self.data = self.inpt.read()
        self.cmol = Chem.MolFromPDBBlock(self.data, sanitize=False)
        self.conf = self.cmol.GetConformer()
        self.ntom = self.cmol.GetNumAtoms()
        self.inpt.close()

    def update_gridbox(self, mol_block: str) -> None:
        self.cmol = Chem.MolFromPDBBlock(mol_block, sanitize=False)
        self.conf = self.cmol.GetConformer()
        self.ntom = self.cmol.GetNumAtoms()

    def compute_coords(self) -> ranges:
        x_coord = [self.conf.GetAtomPosition(c).x for c in range(self.ntom)]
        y_coord = [self.conf.GetAtomPosition(c).y for c in range(self.ntom)]
        z_coord = [self.conf.GetAtomPosition(c).z for c in range(self.ntom)]
        return x_coord, y_coord, z_coord

    def compute_ranges(self) -> ranges:
        x, y, z = self.compute_coords()
        x_range = [min(x), max(x)]
        y_range = [min(y), max(y)]
        z_range = [min(z), max(z)]
        return x_range, y_range, z_range

    def compute_center(self, use_range: bool = True) -> coords:
        x, y, z = self.compute_ranges() if use_range else self.compute_coords()
        x_center = round(np.mean(x), 3)
        y_center = round(np.mean(y), 3)
        z_center = round(np.mean(z), 3)
        return x_center, y_center, z_center

    def generate_res_molblock(self, residues_list: list[str]) -> str:
        res_lines = [line for line in self.data.split('\n')
                     if line[22:26].lstrip() in residues_list
                     and 'END' not in line]
        res_block = '\n'.join(res_lines)
        return res_block

    def labox(self, scale: float = 2.0) -> coords:
        xr, yr, zr = self.compute_ranges()
        center = self.compute_center()
        
        x_size = abs(xr[0] - xr[1])
        y_size = abs(yr[0] - yr[1])
        z_size = abs(yr[0] - yr[1])
                       
        x_size = 1 if x_size < 1 else x_size
        y_size = 1 if y_size < 1 else y_size
        z_size = 1 if z_size < 1 else z_size
        
        bxsize = (round(x_size * scale, 3),
                  round(y_size * scale, 3),
                  round(z_size * scale, 3))
        return center, bxsize

    def eboxsize(self, gy_box_ratio: float = 0.23, modified: bool = False) -> center_bxsize:
        xc, yc, zc = self.compute_coords()
        center = self.compute_center(modified)
        distsq = [(x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2
                  for x, y, z in zip(xc, yc, zc)]
        bxsize = (round(np.sqrt(sum(distsq) / len(xc)) / gy_box_ratio, 3),) * 3
        return center, bxsize

    def autodock_grid(self) -> center_bxsize:
        xr, yr, zr = self.compute_ranges()
        center = self.compute_center()
        bxsize = (22.5, 22.5, 22.5)
        return center, bxsize

    def defined_by_res(self, residue_number: str, scale: float = 1.25) -> center_bxsize:
        res_list = residue_number.replace(',', ' ').split()
        res_block = self.generate_res_molblock(res_list)
        self.update_gridbox(res_block)
        return self.labox(scale=scale)


##############################
# Add chain ID "A" to the pdbqt file for ligand

def pdbqt_add_chid(inpt_file: str) -> None:
    with open(inpt_file, 'r') as inpt:
        data = inpt.read()
    new_data = data.replace('  UNL  ', '  UNL A')
    with open(inpt_file, 'w') as oupt:
        oupt.write(new_data)

##############################
# Create Report
# From pdbqt file for a single pose:

## vina
def vina_report(inpt_file: str) -> dict:

    def vn_scores(inpt_file: str) -> float:
        with open(inpt_file, 'r') as inpt_data:
            for line in inpt_data:
                if 'REMARK VINA RESULT' in line:
                    Affinity = float(line.split()[3])
                    rmsd_lb = float(line.split()[4])
                    rmsd_ub = float(line.split()[5])
        return Affinity, rmsd_lb, rmsd_ub

    def vn_report(name: str, scores: tuple) -> dict:
        return {'NAME': [name], 'Affinity': [scores[0]],
                'RMSD_LB': [scores[1]], 'RMSD_UB': [scores[2]]}

    inpt_name = os.path.basename(inpt_file).split('.')[0]
    scores = vn_scores(inpt_file)
    report = vn_report(inpt_name, scores)
    return report


def vina_report_complete(inpt_file: str) -> dict:

    def vn_scores(inpt_file: str) -> float:
        with open(inpt_file, 'r') as inpt_data:
            Affinity = []
            rmsd_lb = []
            rmsd_ub = []
            for line in inpt_data:
                if 'REMARK VINA RESULT' in line:
                    Affinity.append(float(line.split()[3]))
                    rmsd_lb.append(float(line.split()[4]))
                    rmsd_ub.append(float(line.split()[5]))
        return Affinity, rmsd_lb, rmsd_ub

    def vn_report(scores: tuple) -> dict:
        return {'Affinity': scores[0],
                'RMSD_LB': scores[1], 'RMSD_UB': scores[2]}

    scores = vn_scores(inpt_file)
    report = vn_report(scores)
    return report


## smina
def smina_report(inpt_file: str) -> dict:

    def vn_scores(inpt_file: str) -> float:
        with open(inpt_file, 'r') as inpt_data:
            for line in inpt_data:
                if 'REMARK minimizedAffinity' in line:
                    Affinity = float(line.split()[2])
        return Affinity

    def vn_report(name: str, scores: tuple) -> dict:
        return {'NAME': [name], 'Affinity': [scores]}

    inpt_name = os.path.basename(inpt_file).split('.')[0]
    scores = vn_scores(inpt_file)
    report = vn_report(inpt_name, scores)
    return report


def smina_report_complete(inpt_file: str) -> dict:

    def vn_scores(inpt_file: str) -> float:
        with open(inpt_file, 'r') as inpt_data:
            Affinity = []
            rmsd_lb = []
            rmsd_ub = []
            for line in inpt_data:
                if 'REMARK minimizedAffinity' in line:
                    Affinity.append(float(line.split()[2]))
        return Affinity

    def vn_report(scores: tuple) -> dict:
        return {'Affinity': scores}

    scores = vn_scores(inpt_file)
    report = vn_report(scores)
    return report


## gnina
def gnina_report(inpt_file: str) -> dict:

    def vn_scores(inpt_file: str) -> float:
        with open(inpt_file, 'r') as inpt_data:
            for line in inpt_data:
                if 'REMARK minimizedAffinity' in line:
                    Affinity = float(line.split()[2][0:-6])
                    CNNscore = float(line.split()[4][0:-6])
                    CNNaffinity = float(line.split()[6][0:-6])
        return Affinity, CNNscore, CNNaffinity

    def vn_report(name: str, scores: tuple) -> dict:
        return {'NAME': [name], 'Affinity': [scores[0]],
                'CNNscore': [scores[1]], 'CNNaffinity': [scores[2]]}

    inpt_name = os.path.basename(inpt_file).split('.')[0]
    scores = vn_scores(inpt_file)
    report = vn_report(inpt_name, scores)
    return report


## gnina
def gnina_report_complete(inpt_file: str) -> dict:

    def vn_scores(inpt_file: str) -> float:
        with open(inpt_file, 'r') as inpt_data:
            Affinity = []
            CNNscore = []
            CNNaffinity = []
            for line in inpt_data:
                if 'REMARK minimizedAffinity' in line:
                    Affinity.append(float(line.split()[2][0:-6]))
                    CNNscore.append(float(line.split()[4][0:-6]))
                    CNNaffinity.append(float(line.split()[6][0:-6]))
        return Affinity, CNNscore, CNNaffinity

    def vn_report(scores: tuple) -> dict:
        return {'Affinity': scores[0],
                'CNNscore': scores[1], 'CNNaffinity': scores[2]}

    scores = vn_scores(inpt_file)
    report = vn_report(scores)
    return report

#########################################
### Conformer

def select_conformer_sdf(sdf_in, sdf_out, conformer_number=0):
    counter = 0
    with open(sdf_in, 'r') as sdf_all, \
         open(sdf_out, 'w') as sdf_query:
        for line in sdf_all:
            if counter == conformer_number:
                sdf_query.write(line)
            if "$$$$" in line:
                counter += 1



###### Docking
def run_docking(query_dict, method = 'smina', n_cpu = 1, exh = 16):

    #########  make aliases
    # %alias vina Softwares/vina
    def vina(*args):
        command = 'software/vina ' + args[0]
        os.system(command)
        
    # %alias smina Softwares/smina.static
    def smina(*args):
        command = '../software/smina.static ' + args[0]
        os.system(command)
            
    # %alias gnina Softwares/gnina
    def gnina(*args):
        command = 'software/gnina ' + args[0]
        os.system(command)
    
    # %alias qvina Softwares/qvina-w
    def qvina(*args):
        command = 'software/qvina-w ' + args[0]
        os.system(command)
        
    # %alias qvina2 Softwares/qvina21
    def qvina2(*args):
        command = 'software/qvina21 ' + args[0]
        os.system(command)

    
    ######## Run docking
    Target_pdbqt_dir = query_dict['protein_pdbqt_dir']
    Ligand_pdbqt_dir = query_dict['ligand_pdbqt_dir']
    dock_output_pdbqt_dir = query_dict['output_pdbqt_dir']
    constraint_file_dir = query_dict['constraint_file_dir']
    dock_output_log_dir = 'output.log'

    start = time.time()
    
    if method == 'vina':
        start = time.time()
        vina(f""" --receptor {Target_pdbqt_dir} --ligand {Ligand_pdbqt_dir} \
        --out {dock_output_pdbqt_dir} --config {constraint_file_dir} --cpu {n_cpu} \
        --exhaustiveness {exh} --verbosity 2 | tee {dock_output_log_dir}""")

    if method == 'smina':
        smina(f""" --receptor {Target_pdbqt_dir} --ligand {Ligand_pdbqt_dir} \
        --out {dock_output_pdbqt_dir} --config {constraint_file_dir} --cpu {n_cpu} \
        --exhaustiveness {exh} --verbosity 2 | tee {dock_output_log_dir}""")
    
    if method == 'gnina':
        gnina(f""" --receptor {Target_pdbqt_dir} --ligand {Ligand_pdbqt_dir} \
        --out {dock_output_pdbqt_dir} --config {constraint_file_dir} --cpu {n_cpu} \
        --exhaustiveness {exh} --verbosity 2 --no_gpu| tee {dock_output_log_dir}""")                    #--no_gpu
    
    if method == 'qvina':
        qvina(f""" --receptor {Target_pdbqt_dir} --ligand {Ligand_pdbqt_dir} \
        --out {dock_output_pdbqt_dir} --config {constraint_file_dir} --cpu {n_cpu} \
        --exhaustiveness {exh} | tee {dock_output_log_dir}""")
    
    if method == 'qvina2':
        qvina2(f""" --receptor {Target_pdbqt_dir} --ligand {Ligand_pdbqt_dir} \
        --out {dock_output_pdbqt_dir} --config {constraint_file_dir} --cpu {n_cpu} \
        --exhaustiveness {exh} | tee {dock_output_log_dir}""")

    end = time.time()
    # print("output file generated:")
    # print(dock_output_pdbqt_dir)



######## extract score
def extract_score(docking_pdbqt_dir, method):

    try:
        ## Process output file
        if (method in ['vina', 'qvina', 'qvina2']):
            report = vina_report_complete(docking_pdbqt_dir)
        elif method == 'smina':
            report = smina_report_complete(docking_pdbqt_dir)
        elif method == 'gnina':
            report = gnina_report_complete(docking_pdbqt_dir)
        report = pd.DataFrame(report)
    
        if (method in ['vina', 'smina', 'qvina', 'qvina2']):
            score = min(report['Affinity'])
            score = {'best affinity' : f'{score}', 'all affinity' : report['Affinity']}
        elif (method == 'gnina'):
            score1 = min(report['Affinity'])
            score2 = max(report['CNNaffinity'])
            score = {'best affinity' : f'{score1}', 'best CNNaffinity' : f'{score2}',
                     'all affinity' : report['Affinity'], 'all CNNaffinity' : report['CNNaffinity']}
    except:
        with open(docking_pdbqt_dir, "r") as file:
            score = file.read()
            
    return score




#########################
def plot_ligand_2D(smiles_list):
    
    fig, axs = plt.subplots(len(smiles_list), 1, figsize=(12, 8))
    axs = axs.flatten()
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)

        mol_h = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol_h)
        AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())

        img = Draw.MolToImage(mol_h, size=(200, 200)) 

        # print(smiles)
        axs[i].imshow(img)
        axs[i].axis('off')

    plt.tight_layout()
    plt.show()


def plot_ligand_3D(smiles_list):
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)

        mol_h = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol_h)
        AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())

        mb = Chem.MolToMolBlock(mol_h)
        print(smiles)
        if mb.find(' H ') > 0:
            view = py3Dmol.view(width=300, height=300)
            view.addModel(mb, "mol")
            view.setStyle({'stick': {}})
            view.zoomTo()
            view.show()


#####################################################################################

def run_protein_ligand_interaction(complex_file):

    interaction_bonds = ['hydrophobic',
                     'hbond',
                     'waterbridge',
                     'saltbridge',
                     'pistacking',
                     'pication',
                     'halogen',
                     'metal']
    

    ## Perform plip analysis
    mol = PDBComplex()
    mol.load_pdb(complex_file)
    mol.analyze()
    
    ## get interaction_sets as a dict
    interactions = mol.interaction_sets['UNL:A:1']
    interaction_report = BindingSiteReport(interactions)
    
    ## Add interactions
    interaction_table = []
    for bond in interaction_bonds:
        bond_features = list(getattr(interaction_report, bond + '_features'))
        bond_info = list(getattr(interaction_report, bond + '_info'))
    
        if bond_info != []:
            # print(bond)
            bond_features = ['BONDTYPE'] + bond_features
            bond_info = [[bond] + list(info) for info in bond_info]
        
            interaction_line = [{feat: info_i for feat, info_i in zip(bond_features, info)} for info in bond_info]
            interaction_table.extend(interaction_line)
    
    interaction_df = pd.DataFrame(interaction_table)
    
    # For certain cases (e.g. hbond), the DIST col is NaN. For those cases, we replace the DIST with other distances availabe 
    # For this, we set a priority list
    distance_cols_to_replace = ['DIST_D-A', 'CENTDIST']
    for col in distance_cols_to_replace:
        if col in interaction_df.columns:
            interaction_df['DIST'] =  interaction_df['DIST'].fillna( interaction_df[col])

    return interaction_df





def interaction_csv2dict(inpt_file) -> dict:

    def s2f_dict(item: dict) -> dict:
        return {key: tuple(float(val) for val in value[1:-1].split(','))
                for key, value in item.items()}

    def b2c_dict(item: dict, usage) -> dict:
        usg_map = {'lbsp': 0, 'view': 1}
        return {key: BOND_COL[val.upper()][usg_map[usage]] for key, val in item.items()}

    inter_df = pd.read_csv(inpt_file)
    int_dict = inter_df.to_dict()
    int_dict['LIGCOO'] = s2f_dict(int_dict['LIGCOO'])
    int_dict['PROTCOO'] = s2f_dict(int_dict['PROTCOO'])
    int_dict['COLOR'] = b2c_dict(int_dict['BONDTYPE'], 'view')
    int_dict['COLORX'] = b2c_dict(int_dict['BONDTYPE'], 'lbsp')

    return int_dict



def show_interaction(viewer, interaction_dict, chain=None):

    residue_style = {
        'stick':
         {'colorscheme': 'orangeCarbon', 'radius': 0.15}}
    residue_label = {
        'alignment': 'bottomLeft',
        'showBackground': False,
        'inFront': True,
        'fontSize': 14,
        'fontColor': '0x000000',
        'screenOffset': {'x': 25, 'y': 25}}
    atom_label = {
        'alignment': 'bottomLeft',
        'showBackground': False,
        'inFront': True,
        'fontSize': 14,
        'fontColor': '0x000000',
        'screenOffset': {'x': 10, 'y': 10}}

    def find_midpoint(coords: list) -> tuple[float, float, float]:
        return tuple(round(coord, 3) for coord in np.mean(coords, axis=0))
    
    dist = interaction_dict['DIST'].values()
    bond = interaction_dict['BONDTYPE'].values()
    resn = list(interaction_dict['RESNR'].values())
    ligcoo = interaction_dict['LIGCOO'].values()
    prtcoo = interaction_dict['PROTCOO'].values()
    color = interaction_dict['COLORX'].values()

    if chain == None:
        viewer.addStyle({'and': [{'model': 0}, {'resi': resn}]}, residue_style)
        viewer.addResLabels({'and': [{'model': 0}, {'resi': resn}]}, residue_label)
    else:
        viewer.addStyle({'and': [{'model': 0, 'chain' : chain}, {'resi': resn}]}, residue_style)
        viewer.addResLabels({'and': [{'model': 0, 'chain' : chain}, {'resi': resn}]}, residue_label)
    
    for dis, col, lig, prt in zip(dist, color, ligcoo, prtcoo):
        mid = find_midpoint([lig, prt])
        viewer.addCylinder(
            {'start': {'x': lig[0], 'y': lig[1], 'z': lig[2]},
             'end': {'x': prt[0], 'y': prt[1], 'z': prt[2]},
             'radius': 0.05,
             'fromCap': 1,
             'toCap': 1,
             'color': col,
             'dashed': True})
        viewer.addLabel(
            str(dis) + ' Å',
            {'position': {'x': mid[0], 'y': mid[1], 'z': mid[2]},
             'alignment': 'bottomLeft',
             'inFront': False,
             'backgroundColor': col,
             'fontSize': 10,
             'screenOffset': {'x': 10, 'y': 10}})

    return viewer




def atom_match(loc, ligand_pdb_content):
    x, y, z = loc
    
    def dist(a, b):
        x1, y1, z1 = a
        x2, y2, z2 = b
        return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    
    d = []
    atom = []
    n_atom = []
    for line in ligand_pdb_content.split('\n'):
        if line[0:4] == 'ATOM':
            x_atom, y_atom, z_atom = re.findall("\d+.\d+", line)[0:3]
            x_atom, y_atom, z_atom = float(x_atom), float(y_atom), float(z_atom)
            atom.append (line[13])
            n_atom.append (re.findall("\d+", line)[0])
            d.append(dist((x, y, z), (x_atom, y_atom, z_atom)))
    d = np.array(d)
    arg_min_d = d.argmin()
    atom = atom[arg_min_d]
    n_atom = n_atom[arg_min_d]
    
    return {'atom': atom, 'atom_number': n_atom}






def smiles_to_3d(ligand_smiles, 
                 ligand_name, 
                 save_dir,
                 method = 'obabel',
                 Force_field = 'MMFF94',
                 Convergence_criteria = '0.00001',
                 Maximum_steps = 10000,
                 verbos = True):
    """
    SMILES to 3D energetically minimised ligand and conformer generation
    The result of this code block is a set of SD and PDB file for each conformer.
    
    method : 'obabel'
    Force_field: 'MMFF94', 'GAFF', 'Ghemical', 'MMFF94', 'MMFF94s', 'UFF'
    Convergence_criteria : '0.00001' '0.1', '0.01','0.001', '0.0001', '0.00001', '0.000001', '0.0000001'
    Maximum_steps : {min:1000, max:100000}
    
    """

    def obabel(*args):
        obabel_loc = '/home/byousefi/envs/chem/bin/obabel'
        command = obabel_loc + ' ' + args[0]
        command = 'timeout 1000s ' + command
        os.system(command)
    
    ligand_sdf_dir = os.path.join(save_dir, ligand_name + '.sdf')
    ligand_pdb_dir = os.path.join(save_dir, ligand_name + '.pdb')
    # ligand_cif_dir = os.path.join(save_dir, ligand_name + '.cif')
    log_dir = os.path.join(save_dir, ligand_name + '_obabel.log')

    # generate sdf file
    obabel(f"""-:{'"'+ligand_smiles+'"'} -O {ligand_sdf_dir} --title {ligand_name} --gen3d \
    --best --minimize --ff {Force_field} --steps {Maximum_steps} --sd \
    --crit {Convergence_criteria} > {log_dir} 2>&1""")

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
        print("Try gen2d")
        # generate sdf file
        obabel(f"""-:{'"'+ligand_smiles+'"'} -O {ligand_sdf_dir} --title {ligand_name} --gen2d \
        --best --minimize --ff {Force_field} --steps {Maximum_steps} --sd \
        --crit {Convergence_criteria} > {log_dir} 2>&1""")
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
            
    return 0 if flag == 'successful' else 1


def prepare_ligand(ligand_smiles_list,
                   ligand_sdf_list,
                   save_dir,
                   method = 'obabel',
                   Force_field = 'MMFF94',
                   Convergence_criteria = '0.00001',
                   Maximum_steps = 10000,
                   verbos = True,
                  ):
    """
    This functions reads ligands from a dataframe file, then
    1. converts them to 3d structure -> sdf and pdf files
    2. prepares them for autodoc_vina -> pdbqt file

    method : 'obabel'
    Force_field: 'MMFF94', 'GAFF', 'Ghemical', 'MMFF94', 'MMFF94s', 'UFF'
    Convergence_criteria : '0.00001' '0.1', '0.01','0.001', '0.0001', '0.00001', '0.000001', '0.0000001'
    Maximum_steps : {min:1000, max:100000}

    """
    # Read data and make ligand_df
    df_smiles_i = df_smiles.loc[df_smiles['pdb_id'] == pdb_id].copy()
    df_smiles_i.columns = ['protein_name', 'ligand_name', 'ligand_smiles']

    ligand_df['sdf_file'] = "NA"
    ligand_df['pdb_file'] = "NA"
    ligand_df['pdbqt_file'] = "NA"
    
    for index, ligand in tqdm(ligand_df.iterrows()):

        ligand_name = ligand['ligand_name']
        ligand_smiles = ligand['ligand_smiles']

        print(ligand_name)

        ligand_sdf_dir = os.path.join(save_dir, ligand_name + '.sdf')
        ligand_pdb_dir = os.path.join(save_dir, ligand_name + '.pdb')
        ligand_pdbqt_dir = os.path.join(save_dir, ligand_name + '.pdbqt')
        
        if not os.path.exists(ligand_sdf_dir):
            # Convert to 3D sdf and pdb
            print('===========')
            flag = smiles_to_3d(ligand_smiles, ligand_name, save_dir, method, Force_field, 
                                Convergence_criteria, Maximum_steps, verbos)
        else:
            flag = 0
            print('3D already exist')

        if flag == 0:
            ligand_df.loc[index, 'sdf_file'] = ligand_sdf_dir
            ligand_df.loc[index, 'pdb_file'] = ligand_pdb_dir
            
        if not os.path.exists(ligand_pdbqt_dir):
            if flag == 0:
                # Convert to pdbqt file
                # mk_prepare_ligand -> full address should be here
                # The command below will not run then
                os.system(f"mk_prepare_ligand.py -i {ligand_sdf_dir} -o {ligand_pdbqt_dir} > /dev/null 2>&1")
                ligand_df.loc[index, 'pdbqt_file'] = ligand_pdbqt_dir
        else:
            ligand_df.loc[index, 'pdbqt_file'] = ligand_pdbqt_dir
            print('PDBQT already exist')
    
    return ligand_df


##############################
# Generate protein-ligand complex pdb file 

def generate_cmpx_pdb(inpt_prot: str, inpt_lig: str, oupt_cmpx: str) -> None:

    """
    Adds proteind pdb lines with ['ATOM', 'CONECT', 'TER'] keywords form inpt_prot file
     and ligand pdb lines with ['ATOM', 'CONECT', 'END'] keywords form inpt_prot file
     to an empty pdb file that is oupt_cmpx file.
    """
    
    def write_line(line: str, keywords: list, oupt_file: str) -> None:
        header = line.split()[0]
        if header in keywords:
            oupt_file.write(line)

    def cmpx_writer() -> None:
        with open(oupt_cmpx, 'w') as oupt_file, \
             open(inpt_prot, 'r') as prot_file, \
             open(inpt_lig, 'r') as lig_file:
            for prot_line in prot_file:
                write_line(prot_line, prot_headers, oupt_file)
            for lig_line in lig_file:
                write_line(lig_line, lig_headers, oupt_file)

    prot_headers = ['ATOM', 'CONECT', 'TER']
    lig_headers = ['ATOM', 'CONECT', 'END']
    cmpx_writer()