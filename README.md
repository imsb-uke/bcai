# Instruction to install and run BioChemAIgent via Docker

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

Whithin the conainer, install external MCP servers
```
./mcp_external/install.sh
```

Finally, run BioChemAIgent
```
cd main
tmux new -s mysession
python client.py
ctrl+b  %
streamlit run BioChemAIgent_chat_UI.py --server.address 0.0.0.0 --server.port 8501
```