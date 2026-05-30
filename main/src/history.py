import os
import json
import pickle
from datetime import datetime

def save_history(messages, history_file):
    """Save full history (with Python objects preserved)."""
    with open(history_file, "wb") as f:
        pickle.dump(messages, f)


def load_history(history_file):
    """Load full history, or empty list if none."""
    if os.path.exists(history_file):
        with open(history_file, "rb") as f:
            return pickle.load(f)
    return []


def message_to_dict(msg):
    """
    Convert a single ChatCompletionMessage into a a dictionary.
    """

    msg_dict = {
        "role": getattr(msg, "role", None),
        "content": getattr(msg, "content", None),
    }

    # Include optional fields if present
    for attr in [
        "refusal", "annotations", "audio",
        "function_call", "tool_calls",
        "name", "tool_name", "tool_id"
    ]:
        value = getattr(msg, attr, None)
        if value is not None:
            msg_dict[attr] = value

    return msg_dict
