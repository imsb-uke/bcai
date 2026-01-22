This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0).

# BioChemAIgent: An AI-driven Protein Modeling and Docking Framework for Structure-Based Drug Discovery

Computational and AI-based methods have advanced drug discovery, yet most remain task-specific and require substantial expert integration. We introduce BioChemAIgent, an agentic framework that integrates state-of-the-art AI models with established computational chemistry tools to enable end-to-end small molecule analysis, protein modeling, and molecular docking, with a chatbot interface accessible online. To foster community engagement and extensibility, we established an openly accessible user interface alongside a centralized registry that consolidates resources for developing and integrating next-generation tools in drug discovery and structural biology.

## Use BioChemAIgent
### Online user interface
The online user interface is availabe at [https://bcai.ims.bio/](https://bcai.ims.bio/)

### Local Installation and Running BioChemAIgent via Docker
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

**Use of external softwares:**
Instructions to install and run external softwares is provided in the `software` directory [./software](here) 

### User access tokens
* To use cloud-based LLM models, three options are availabe: OpenAI [https://openai.com/api/], Ollama [https://ollama.com/], OpenRouter [https://openrouter.ai/]. These access tokens should be added to the `env` file.
* The use of ESM3 [https://github.com/evolutionaryscale/esm] and AlphaFold3 [https://github.com/google-deepmind/alphafold3] requires user access tokens. These access tokens should be added to the `env` file.

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
