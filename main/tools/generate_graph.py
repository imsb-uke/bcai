import os
import json
import pandas as pd
from rdkit import Chem
from Bio.Data import IUPACData
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport

# Global variables
protein_residues = list(IUPACData.protein_letters_3to1.keys())
protein_residues = [i.upper() for i in protein_residues]

# Functions
def ligand_atoms(ligand_file):

    """
    Extract ligand atoms from a sdf of pdb file.
    The output is a pandas dataframe
    """

    ext = ligand_file.split('.')[-1].lower()
    if ext == 'sdf':
        mol = Chem.SDMolSupplier(ligand_file, 
                                 removeHs=False, 
                                 sanitize=False,     # skip the clean‑up
                                 strictParsing=False
                                )[0]
    else:
        mol = Chem.MolFromPDBFile(ligand_file,
                                  removeHs=False,
                                  sanitize=False,
                                  proximityBonding=True
                                 )
    # Get 3d structure
    conf = mol.GetConformer()
    
    df = pd.DataFrame(columns=['index', 'index_in_pdb', 'symbol', 'x', 'y', 'z'])
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        pos = conf.GetAtomPosition(idx)
        x, y, z = pos.x, pos.y, pos.z
        atom_index = idx + 1
        if ext == 'pdb' and atom.GetPDBResidueInfo() is not None:          # new: index comes from pdb
            atom_index_in_pdb = atom.GetPDBResidueInfo().GetSerialNumber()
        else:
            atom_index_in_pdb = None
        df.loc[len(df)] = [atom_index, atom_index_in_pdb, symbol, x, y, z]
    return df

    
def ligand_bonds(ligand_file):

    """
    Extract atomic bonds in a ligand presented by a sdf or pdb file
    The output is a pandas dataframe
    """
    if ligand_file.split('.')[-1] == 'sdf':
        mol = Chem.SDMolSupplier(ligand_file, 
                                 removeHs=False, 
                                 sanitize=False,     # skip the clean‑up
                                 strictParsing=False
                                )[0]
    else:
        mol = Chem.MolFromPDBFile(ligand_file,
                                  removeHs=False,
                                  sanitize=False,
                                  proximityBonding=True
                                 )
    
    # Get 3d structure
    conf = mol.GetConformer()
    
    # Extract bonds and positions
    bonds_data = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atom1_index = atom1.GetIdx()
        atom2_index = atom2.GetIdx()
        atom1_pos = conf.GetAtomPosition(atom1_index)
        atom2_pos = conf.GetAtomPosition(atom2_index)
        bond_info = {
            "atom1_idx": atom1_index + 1,
            "atom2_idx": atom2_index + 1,
            "atom1_symbol": atom1.GetSymbol(),
            "atom2_symbol": atom2.GetSymbol(),
            "atom1_position": (atom1_pos.x, atom1_pos.y, atom1_pos.z),
            "atom2_position": (atom2_pos.x, atom2_pos.y, atom2_pos.z),
            "bond_type": str(bond.GetBondType()),
        }
        bonds_data.append(bond_info)
    
    df = pd.DataFrame(bonds_data)
    return df


def protein_interacting_residues(pdb_file, ligand_id='UNK'):

    """
    This code listst prorein residues that intect with a ligand.
    But it does not specify interacting ligand atoms
    For a more compregensive function, use protein_ligand_interaction()
    """
    
    interaction_types = [
        'hbond',
        'hydrophobic',
        'waterbridge',
        'saltbridge',
        'pistacking',
        'pication',
        'halogen',
        'metal'
    ]
    
    mol = PDBComplex()
    mol.load_pdb(pdb_file)
    mol.analyze()
    
    # Get ligand ids and loop over them
    interaction_all = []
    
    lig_id_set = str(mol).split('\n')[1:]
    for lig in lig_id_set:
        lig_id = lig.split(':')[0]
        if lig_id == ligand_id:                  # select a specific ligand
            ## Get interaction report
            interactions = mol.interaction_sets[lig]
            interaction_report = BindingSiteReport(interactions)
            
            ## Hint: to get each interaction, we have: 
            ## interaction_report.hbond_features for feature names and
            ## interaction_report.hbond_info for feature values.
    
            for i in interaction_types:
                interaction_features = getattr(interaction_report, f'{i}_features')
                interaction_info = getattr(interaction_report, f'{i}_info')
                
                if interaction_info != []:
                    ### Add interaction type i to the features
                    interaction_features = ['BONDTYPE'] + list(interaction_features)
                    interaction_info = [[i] + list(info) for info in interaction_info]
    
                    for info in interaction_info:
                        interaction_all.append({feat: info_i for feat, info_i in zip(interaction_features, info)})
            
    df_lig_prot = pd.DataFrame(interaction_all)
    return df_lig_prot


def protein_ligand_interaction(pdb_file, ligand_id=['UNK']):

    """
    Lists all protein ligand interactions
    """
        
    interaction_types = [
        'hbond',
        'hydrophobic',
        'waterbridge',
        'saltbridge',
        'pistacking',
        'pication',
        'halogen',
        'metal'
    ]
    
    mol = PDBComplex()
    mol.load_pdb(pdb_file)
    mol.analyze()
    # print(mol)      # lists all ligands inside the mol, we only need 'UNK'
    
    # Get ligand ids and loop over them
    result = []
    
    lig_id_set = str(mol).split('\n')[1:]
    for lig in lig_id_set:
        lig_id = lig.split(':')[0]
        if lig_id in ligand_id:                  # select a specific ligand
            ## Get interaction report
            interactions = mol.interaction_sets[lig]
            interactions_all = interactions.all_itypes
            
            for idx in interactions_all:
                try:
                    result_dict = get_ligand_protein_atoms_from_plip(idx)
                except:
                    continue
                result_dict['prot_chain'] = idx.reschain
                result_dict['prot_res_num'] = idx.resnr
                result_dict['prot_res_type'] = idx.restype
                result_dict['lig_id'] = idx.restype_l
                result.append(result_dict)
    
    df_lig_prot = pd.DataFrame(result)
    return df_lig_prot


def get_ligand_protein_atoms_from_plip(x):
    """
    x is a PLIP ineraction
    Here we assume that the protein atom indices come after ligand atoms
    """
    # print("======")
    # print(x)
    
    
    def get_type(idx1, idx2):
        idx1 = idx1[0] if type(idx1) == list else idx1
        idx2 = idx2[0] if type(idx2) == list else idx2
        return ('ligand', 'protein') if idx1 > idx2 else ('protein', 'ligand')  # ligand comes after protein
    
    itype = x.__class__.__name__

    lig_idx = []
    prot_idx = []
    water_idx = []
    feature = []

    if itype == 'hbond':
        idx1 = x.a_orig_idx     # hbonf acceptor
        idx2 = x.d_orig_idx     # hbond donor
        types = get_type(idx1, idx2)
        lig_idx = idx1 if types[0] == 'ligand' else idx2
        prot_idx = idx2 if types[0] == 'ligand' else idx1
        feature = 'hbond_ligand_acceptor' if types[0] == 'ligand' else 'hbond_ligand_donnor'
    
    elif itype == 'hydroph_interaction':
        lig_idx = x.ligatom_orig_idx
        prot_idx = x.bsatom_orig_idx
        feature = 'hydrophobic'

    elif itype == 'halogenbond':
        idx1 = x.acc_orig_idx   # acceptor
        idx2 = x.don_orig_idx   # donor
        types = get_type(idx1, idx2)
        lig_idx = idx1 if types[0] == 'ligand' else idx2
        prot_idx = idx2 if types[0] == 'ligand' else idx1
        feature = 'halogen_ligand_acceptor' if types[0] == 'ligand' else 'halogen_ligand_donnor'

    elif itype == 'waterbridge':
        water_idx = x.water_orig_idx  # water mol and not protein (always donor)
        acc_idx = x.a_orig_idx        # acceptor
        acc_type = x.a.residue.name
        if acc_type in protein_residues:
            prot_idx = acc_idx
        else:
            lig_idx = acc_idx
        feature = 'waterbridge'

    elif itype == 'saltbridge':
        idx1 = x.negative.atoms_orig_idx
        idx2 = x.positive.atoms_orig_idx
        types = get_type(idx1, idx2)
        lig_idx = idx1 if types[0] == 'ligand' else idx2
        prot_idx = idx2 if types[0] == 'ligand' else idx1
        feature = 'ionic_ligand_negative' if types[0] == 'ligand' else 'ionic_ligand_positive'

    elif itype == 'pistack':
        lig_idx = x.ligandring.atoms_orig_idx
        prot_idx = x.proteinring.atoms_orig_idx
        feature = itype

    elif itype == 'pication':
        ring_idx = x.ring.atoms_orig_idx
        charge_idx = x.charge.atoms_orig_idx
        types = get_type(ring_idx, charge_idx)
        lig_idx = ring_idx if types[0] == 'ligand' else charge_idx
        prot_idx = charge_idx if types[0] == 'ligand' else ring_idx
        feature = 'pication_ligand_ring' if types[0] == 'ligand' else 'pication_ligand_charge'
        
    elif itype == 'metal_complex':
        lig_idx = idx.metal_orig_idx
        prot_idx = idx.target_orig_idx
        feature = f'metalic_{idx.metal_type}'

    result = {
        'itype'       : itype,
        'ligand_idx'  : lig_idx,
        'protein_idx' : prot_idx,
        'water_idx'   : water_idx,
        'feature'     : feature,
    
    }
    return result


def protein_pocket(pdb_file, ligand_id = 'UNK', threshold = 6):

    """
    Get protein pocket
    Strategy 1:
    Get protein residues with at least one atom closer than a threshold (6 Ångström) to ligands (+ those interact with ligand)
    Strategy 2:
    Get protein residues facing the ligand (+ those interact with ligand)
    Here, for now, we implement Strategy 1 as it's more simple. 
    """

    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', pdb_file)
    model = structure[0]
    
    # Collect ligand atoms and protein atoms
    ligand_atoms = []
    protein_atoms = []
    for chain in model:
        for residue in chain:
            het_flag, seq_id, ins_code = residue.id
            if residue.get_resname() == ligand_id and het_flag != ' ':                # can be ligand, HO2, or aa
                ligand_atoms.extend(residue.get_atoms())
            # elif het_flag == ' ':
            elif is_aa(residue, standard=True):
                protein_atoms.extend(residue.get_atoms())
    
    # Build neighbor search index on protein atoms
    ns = NeighborSearch(protein_atoms)
    
    # Find residues with at least one atom within threshold of any ligand atom
    close_residues = set()
    for latom in ligand_atoms:
        center = latom.get_coord()
        neighbors = ns.search(center, threshold)
        for patom in neighbors:
            parent = patom.get_parent()
            chain_id = parent.get_parent().id
            het_flag, resseq, icode = parent.id
            # if het_flag == ' ':  # ensure it's a standard residue
            if is_aa(parent, standard=True):
                if 'CA' in parent:
                    resname = parent.get_resname()
                    x, y, z = parent['CA'].get_coord()
                    close_residues.add((chain_id, resname, resseq, x, y, z))
    
    df_pocket = pd.DataFrame(columns=['chain', 'residue', 'index', 'x', 'y', 'z'])
    for row in sorted(close_residues):
        df_pocket.loc[len(df_pocket)] = row

    return df_pocket


def protein_residue_interaction(df_pocket):

    df_pocket = df_pocket.sort_values('index').reset_index(drop=True)
    
    df_res_interaction = pd.DataFrame(columns=['chain1', 'chain2', 'res1', 'res2', 'idx1', 'idx2'])
    
    for i in df_pocket.index[:-1]:
        idx1 = df_pocket.loc[i, 'index']
        idx2 = df_pocket.loc[i+1, 'index']
        res1 = df_pocket.loc[i, 'residue']
        res2 = df_pocket.loc[i+1, 'residue']   
        chain1 = df_pocket.loc[i, 'chain']
        chain2 = df_pocket.loc[i+1, 'chain']
        
        if idx1+1 == idx2 and chain1 == chain2:
            df_res_interaction.loc[len(df_res_interaction)] = [chain1, chain2, res1, res2, idx1, idx2]
    return df_res_interaction


def generate_graph(pdb_file, sdf_file, pocket_threshold = 6):

    # Get ligand atoms
    df_lignd_atom = ligand_atoms(sdf_file)
    # Get ligand bonds
    df_ligand_edge = ligand_bonds(sdf_file)

    # Get protein pocket
    df_pocket = protein_pocket(pdb_file, ligand_id = 'UNK', threshold = pocket_threshold)
    # Correct protein numbers
    # The indices of protein residues should start from the end of ligand atom numbers and they should be sequental with no gap
    n_lig_atom = df_lignd_atom.shape[0]
    df_pocket['index_new'] = df_pocket.index + n_lig_atom + 1
    # We also make a maping dict
    prot_res_map = {i : j for i, j in zip(df_pocket['index'], df_pocket['index_new'])}

    # Get protein residue connections in the pocket
    df_protein_edge = protein_residue_interaction(df_pocket)
    df_protein_edge['idx1_new'] = df_protein_edge['idx1'].map(prot_res_map)
    df_protein_edge['idx2_new'] = df_protein_edge['idx2'].map(prot_res_map)

    # Get protein-ligand interactions
    df_lig_prot = protein_ligand_interaction(pdb_file)
    if len(df_lig_prot):
        df_lig_prot['idx_new'] = df_lig_prot['prot_res_num'].map(prot_res_map)
    else:
        print("No protein-ligand interaction found")
        return [],[]

    ## Now we have all the pieces**
    ### Node features
    ##### ```df_lignd_atom```
    ##### ```df_pocket```
    ### Edge list
    ##### ```df_ligand_edge```
    ##### ```df_protein_edge```
    ##### ```df_lig_prot```
        
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
            row['index_new'],
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
            row['idx1_new'],
            row['idx2_new'],
            'protein-protein',
            'peptide',
        ]
    
    for _, row in df_lig_prot.iterrows():
        if row['ligand_idx'] != []:
            df_edges.loc[len(df_edges)] = [
                row['ligand_idx'],
                row['idx_new'],
                'ligand-protein',
                row['feature'].lower(),
            ]

    return df_nodes, df_edges


if __name__ == "__main__":
    
    sdf_file = ("ligand.sdf")
    pdb_file = ("complex.pdb")

    pocket_threshold = 6  # in Ångström
    df_nodes, df_edges = generate_graph(pdb_file, sdf_file, pocket_threshold)