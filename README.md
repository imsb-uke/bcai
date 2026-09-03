This project is licensed under the PolyForm Noncommercial License 1.0.0 — see `LICENSE`. Free for academic, research, and other noncommercial use. Commercial use requires a separate commercial license. Contact: yousefi.bme@gmail.com, sbonn@uke.de

# BioChemAIgent: An AI-driven Protein Modeling and Docking Framework for Structure-Based Drug Discovery

Computational and AI-based methods have advanced drug discovery, yet most remain task-specific and require substantial expert integration. We introduce BioChemAIgent, an agentic framework that integrates state-of-the-art AI models with established computational chemistry tools to enable end-to-end small molecule analysis, protein modeling, and molecular docking, with a chatbot interface accessible online. To foster community engagement and extensibility, we established an openly accessible user interface alongside a centralized registry that consolidates resources for developing and integrating next-generation tools in drug discovery and structural biology.

## How to run

1. Fist clone external MCP servers

```
cd external-mcp-servers
git clone https://github.com/Augmented-Nature/PDB-MCP-Server.git
git clone https://github.com/Augmented-Nature/ChEMBL-MCP-Server.git
git clone https://github.com/Augmented-Nature/Augmented-Nature-UniProt-MCP-Server.git
```

2. Make user_data folder to writable from the docker
```
chmod 777 user_data
chmod -R 777 external-mcp-servers/
mv Augmented-Nature-UniProt-MCP-Server/ UniProt-MCP-Server/
mv docker_compose_deploy docker_compose
```

3. Build the docker compose and then run it. The UI is accewssible via URL: http://localhost:3000

```
docker compose build --no-cache api web
docker compose up -d api web
```

To stop the platform
```
docker compose stop api web
```

## Contributions
This repository provides a community-oriented registry specialized for the development of agentic systems for drug discovery and structural biology.
* Domain experts can bring forward ideas to be implemented;
* Developers can propose or refine MCP servers to broaden and improve the agent’s functionality.

1. Create a new issue
2. Select the template
   - `Add feature` to suggest new functionalities to be added,
   - `Add MCP server` to suggest an existing MCP server.
3. Make a description of your suggestion and publish the issue.


## Project status
BioChemAIgent is an open-source research framework for structure-based drug discovery.
The method is described in a bioRxiv preprint and has **not yet undergone peer review**.
The API may evolve as the framework is refined.


## Citation
If you use BioChemAIgent, please cite:

Yousefi, B., Laubach, N. C., Heins, S., Testa, L., Gersting, S. W., & Bonn, S. (2025). BioChemAIgent: An AI-driven Protein Modeling and Docking Framework for Structure-Based Drug Discovery. bioRxiv, 2025-12.


## License

This software is licensed under the PolyForm Noncommercial License 1.0.0 — see `LICENSE`.

Free for academic, research, and other noncommercial use. Commercial use — including
use in proprietary pipelines, SaaS products, or any revenue-generating activity —
requires a separate commercial license. Contact: yousefi.bme@gmail.com, sbonn@uke.de

SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0