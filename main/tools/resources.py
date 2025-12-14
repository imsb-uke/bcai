def drug_discovery_protocol() -> str:
    """Provide the drug discovery protocol and steps to LLM as a md file"""
    with open ("../docs/drug_discovery_protocol.md", 'r') as file:
        text = file.read()
    return text


def collaborative_people() -> str:
    """Provide the collaborative people to LLM as a md file"""
    with open ("../docs/collaborative_people.md", 'r') as file:
        text = file.read()
    return text


def get_tools_doc(tool: str) -> str:
    """Provide the usage document for a given tool including `run_esm3`, `run_af3`, `run_vina`, `run_diffdock`"""
    path = f"../docs/tool_instructions/{tool}.md"
    try:
        with open (path, 'r') as file:
            text = file.read()
    except:
        text = 'wrong tool name'
    return text