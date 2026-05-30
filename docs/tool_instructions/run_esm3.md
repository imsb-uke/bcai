# ESM3 Tool Documentation

To complete the protein sequence, we need to provide the `protein_input` field in the request. Let's proceed by setting up this field correctly. Please hold on a moment.

## Overview

ESM3 is a protain language model designed to perform protein sequence analysis, including sequence completion, analyzing protein structures, and various protein properties based on given input data. An MCP tool named as ```run_esm3``` is provided to run ESM3 perfroming different tasks:

1. **Sequence:** Predicting or completeing the protein sequence.
2. **Structure:** Predicting the protein's 3d structure.
3. **Secondary Structure:** Predicting the secondary structure elements, like alpha helices and beta sheets.
4. **SASA (Solvent Accessible Surface Area):** Predicting the solvent-accessible surface area per residue in the protein.
5. **Function:** Predicting the function of each residue in the protein structure.
6. **Residue Annotations:** Providing detailed annotations for each residue, possibly including their role or importance.

## Input Structure

Important note: `protein_input` and `task` are the most important inputs. See examples below.

1. **protein_input:** a dictionary with the following keys
* pkl_file: A pickle file of the ESM protein object. If available, this will be used first.
* pdb_file: A PDB file of the predicted protein structure, if available. This is used if `pkl_file` is not provided.
* sequence: The protein sequence. If `pkl_file` and `pdb_file` are not provided, sequence is mandatory.
* sasa: Solvent Accessible Surface Area per residue. Optional.
* function: Function per residue. Optional.
* residue_annotations: Annotations per residue. Optional.

Note: If sequence is availabe, it will alaways be added to the protein, and it should match the length of the protein sequence.

2. **task:** Specifies the task to be performed. Options are 'sequence', 'structure', 'secondary_structure', 'sasa', 'function', or 'residue_annotations'.

3. **token**: Required for authentication if interacting with external services.
4. **protein_name**: (Optional) The name used to save files like `pkl_file` and `pdb_file`. Default is 'my_protein'.
5. **model_name**: (Optional) Model name, default is 'esm3-large-2024-03'.

## Outputs

The tools returns a dictionary with the following keys:
* pkl_file : the directory that the protein pkl file is stored. This s file is always automatically saved after each run allowing to perform additional tasks on the protein if the user needs.
* pdb_file : the directory that the protein pdb file is stored (only if the `task` is `Structure`)
* sequence : protein sequence
* sasa : protein sasa
* function : function annotations
* secondary_structure : protein secondary structure
* residue_annotations' : protein residue structure
* plddt' : predicted plddt
* ptm'   : predicted ptm

## Examples 

The quiried task defines the format of the `protein_input` dicitionary. Here we provide some examples:

### 3D struture prediction
```
{
protein_input : {
    'pkl_file' : '',
    'pdb_file' : '',
    'sequence' : 'AAATTTCCC',
    'sasa'     : '',
    'function' : '',
    'residue_annotations' : ''
},
task : 'structure'
}
```

### Completing a Sequence
To complete a protein sequence with missing parts denoted by underscores:
```
{
protein_input : {
    'pkl_file' : '',
    'pdb_file' : '',
    'sequence' : 'AAATTTC__C',
    'sasa'     : '',
    'function' : '',
    'residue_annotations' : ''
},
task : 'sequence'
}
```

### Inverse folding
Obtain protein sequence given a pdb file
```
{
protein_input : {
    'pkl_file' : '',
    'pdb_file' : 'files/protein.pdb',
    'sequence' : '',
    'sasa'     : '',
    'function' : '',
    'residue_annotations' : ''
}
task : 'sequence'
}
```

### Get solvent accessible surface area (SASA)
Here we consider three cases:
* Get SASA given sequence
```
{
protein_input : {
    'pkl_file' : '',
    'pdb_file' : '',
    'sequence' : 'AAATTTCCC',
    'sasa'     : '',
    'function' : '',
    'residue_annotations' : ''
},
task : 'sasa'
}
```
* Get SASA given pdb
```
{
protein_input : {
    'pkl_file' : '',
    'pdb_file' : 'files/protein.pdb',
    'sequence' : '',
    'sasa'     : '',
    'function' : '',
    'residue_annotations' : ''
},
task : 'sasa'
}
```
* Get SASA given the protein that has esm3 has already been performed on it for another task 
```
{
protein_input : {
    'pkl_file' : 'files/protein.pkl',
    'pdb_file' : '',
    'sequence' : '',
    'sasa'     : '',
    'function' : '',
    'residue_annotations' : ''
},
task : 'sasa'
}
```


## Error Handling

### Common Errors
- **Missing Fields**: Ensure that at least one of the required fields of the `protein_input` dictionary is provided.
- **Invalid Files**: Check for correct paths and formats for `pkl_file` or `pdb_file`.

## Best Practices

- Always validate your input against the specified format to minimize errors.
- Use the sequence field when no structure files (pkl or pdb) are available.
- Utilize provided task examples for various use cases to understand the expected output.
