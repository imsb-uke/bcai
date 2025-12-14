#!/usr/bin/env bash
set -e

if [ ! -d "/workspace/mcp_external/PDB-MCP-Server" ]; then
	echo "Cloning repo PDB-MCP-Server.git..."
	git clone https://github.com/Augmented-Nature/PDB-MCP-Server.git /workspace/mcp_external/PDB-MCP-Server
	echo "DONE"
else
	echo "Directory already exists, skipping clone."
fi


if [ ! -d "/workspace/mcp_external/ChEMBL-MCP-Server" ]; then
	echo "Cloning repo https://github.com/Augmented-Nature/ChEMBL-MCP-Server.git..."
	git clone https://github.com/Augmented-Nature/ChEMBL-MCP-Server.git /workspace/mcp_external/ChEMBL-MCP-Server
	echo "DONE"
else
	echo "Directory already exists, skipping clone."
fi


cd /workspace/mcp_external
echo "Installing the repositories"
./install.sh
echo "DONE"

cd /workspace/main

echo "Activating mamba"
eval "$(mamba shell hook --shell bash)"
eval "$(mamba activate bcai)"
echo "DONE"

echo "starting the scripts"

python client.py & streamlit run BioChemAIgent_chat_UI.py & wait

exec "$@"
