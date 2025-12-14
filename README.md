# Instruction to install and run BioChemAIgent via Docker

## Run using the docker container
Download external MCP servers
```
cd mcp_external
git clone https://github.com/Augmented-Nature/PDB-MCP-Server.git
git clone https://github.com/Augmented-Nature/ChEMBL-MCP-Server.git
```

Build the image and run the conainer
```
docker build -t bcai-image .
docker run --rm -it -p 8501:8501 -v "$PWD":/workspace bcai-image bash
```

