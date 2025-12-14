import os
from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP

load_dotenv("../env")
from tools import *

# Create an MCP server
mcp = FastMCP(name="BioChemAIgent")

for fn in TOOLS:
    mcp.tool()(fn)

# Run the server
if __name__ == "__main__":
    mcp.run(transport="stdio")