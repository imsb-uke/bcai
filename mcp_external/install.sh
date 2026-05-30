#!/usr/bin/env bash
set -e

echo "=== Installing PDB MCP Server ==="
cd /workspace/mcp_external/PDB-MCP-Server
npm install
npm run build

echo "=== Installing ChEMBL MCP Server ==="
cd /workspace/mcp_external/ChEMBL-MCP-Server
npm install
npm run build

echo "=== MCP Installation Complete ==="