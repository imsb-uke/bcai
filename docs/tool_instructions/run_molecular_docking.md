# Molecular Docking Documentation

## Overview

Molecular docking is a tool to predict the pose (and affinity) that a ligand binds a protein. There are currently two methods implemented for molecular docking.

* **AutoDock VINA** with three different versions: vina, smina, and gnina. smina is more efficient than vina and is the defult method. gnina has a better affinity prediction method and needs GPU, which might not be accesible in all machines.

* **DiffDock** is AI-based tool to perform molecular docking. To use this tool more efficiently, you need CUDA activated GPU with 16GB VRAM. DiffDock, however, does not provide any affinity score.

Please always ask the user to choose the docking method.

## Workflow, inputs and outputs

Here is the workflow to run molecular docking.

0. **Input:** As inputs, we need a protein and a ligand. For the protein, we need the pdb file and for the ligand we need sdf files.
- For protein, if we have the sequence, we can use ESM3 or AlphaFold3 to obtain the pdb file.
- For liagnd, if we have the SMILES notations, we can convert it into sdf in the workflow below.
- We also need a project name so that we can use it throgh the whole workflow (e.g. my_docking).
- The project name should remain the same across all the following steps as the files are named based on the project name.
- Always ask for the above information. Chembl MCP can be used to extract all the data needed.

1. **Extract pdb components:** Extract protein chains, ligand, and ion component from the pdb file and save them as separate files using `extract_pdb_components`.

2. **Generate ligand 3D struture:** Use `smiles_to_3d` to convert SMILES to sdf and pdb files representing by energy minimization.

3. **Protonation and energy minimization:** If needed, add hydrogens based on a given pH using `protonate` and then perfrom energy minimization.

4. **Prepare protein:** Use the tool `prepare_protein` to prepare protein csv table. Moreover, it converts `pdb` file of a protein as a `pdbqt` file required for docking with VINA family. For pdbqt file, we defined partial_charge: 'gasteiger', 'mmff94', 'qeq', 'qtpie'. Default is 'gasteiger'.

5. **Prepare ligand:** Use the tool `prepare_ligand` to prepare ligand csv table. Moreover, it converts `sdf` file of a protein as a `pdbqt` file required for docking with VINA family. The input to this method could either be a sdf fle or ligand's SMILES notation. In cased of the latter, the tool automatically generates the sdf file.

Note: To generate sdf and pdb files of a ligand from SMILES, an extra too `smiles_to_3d` has been prepare. `prepare_ligand` internally uses this function. The tool is a wrapper for obabel software with the following inputs:
- Force_field: 'MMFF94', 'GAFF', 'Ghemical', 'MMFF94', 'MMFF94s', 'UFF'. Default is 'MMFF94'.
- Convergence_criteria : '0.00001' '0.1', '0.01','0.001', '0.0001', '0.00001', '0.000001', '0.0000001'. Default is '0.00001'.
- Maximum_steps : {min:1000, max:100000}. Default is 10000.

6. **Prepare query table:** Make the final query table for vina and diffdock. The inputs are the csv tables generated for the protein and the ligand, as well as the docking method (Default is 'smina'.)

7. **Perform docking:** Perform the molecular docking and get the sdf file of the docket ligand using `run_molecular_docking`. It also returns the affinity scores that needs to be presented for the user. the Vina family generates 9 different poses and save them in a separate folder (named by the project name) as `PATH/TO/POSE/FILE_*.sdf` files *:{1,2,...,9}. For example the first pose can be like be `files/pah_bh4_docking/smina_PAH_chainA_BH4_pose_1.sdf`.

8. **Make complex file and get_protein_ligand_interaction:** The uses choses a pose (by defult choose the first pose) and generate the protein-ligand complex file. It then calculates protein-ligand interaction given pdb file of a protein and sdf file of a ligad

## Example roadmap

1. Obtain or predict the protein structure (PDB format):
   - Retrieve the PDB file for the PAH protein from a database if available.
   - Alternatively, predict the protein structure using tools like ESM3 or AF3 if the PDB is not available.

2. Retrieve ligand information:
   - Get the BH4 molecule structure using ChEMBL or other chemical databases.

3. Convert ligand to 3D formats.

4. Extract protein chains, ligand, and ion component from the pdb file.

5. Make protein and ligand protonation at a given pH and then do energy minimization (optional)

6. Prepare protein and ligand PDBQT files (only for Vina, Smina, and Gnina).

7. Perform molecular docking:
   - Use a docking method (including AutoDock Vina, Smina, Gnina, DiffDock) to perform the docking analysis.

8. Analyze and interpret the results and :
   - Review the docking scores
   - Extract protein-ligand interaction table.

## Notes:
* When defining constraint box, the pdb of the native ligand can be used that is extracted by `extract_pdb_components`. However, some native ligands are not biological and are not suitable for the constraint box. In such cases you need to warn the user and be defult use the whole protein.
* In docking, the Vina family does not accept any HET atom in the protein pdb file.