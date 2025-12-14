import os
import time
import json
from dotenv import load_dotenv
import streamlit as st

load_dotenv("../env")
import src.utils_client as utils
import src.user as user
from src.history import save_history, load_history, message_to_dict

# -----------------------------------------------------------------------------
# Login page
# -----------------------------------------------------------------------------
st.set_page_config(page_title="BioChemAIgent Chat", page_icon="💊", layout="wide")
st.title("BioChemAIgent Chat")

user_file = os.environ["USER_FILE"]
if not os.path.exists(user_file):
    os.system(f"echo 'username,password,free_questions,n_questions' > {user_file}")

if "user" not in st.session_state:
    user.login()
    st.stop()

# -----------------------------------------------------------------------------
# Main page
# -----------------------------------------------------------------------------

# History and file directory
HISTORY_PARENT_DIR = "history"
SESSION_PARENT_DIR = "session"
os.environ["HISTORY_DIR"] = os.path.join(HISTORY_PARENT_DIR, st.session_state['user'])
os.environ["SESSION_DIR"] = os.path.join(SESSION_PARENT_DIR, st.session_state['user'])

HISTORY_DIR = os.getenv("HISTORY_DIR")
SESSION_DIR = os.getenv("SESSION_DIR")
FILE_DIR = os.getenv("FILE_DIR")

HISTORY_FILE = os.path.join(HISTORY_DIR, "chat_history.pkl")
SESSION_FILE = os.path.join(SESSION_DIR, "chat_session.json")

# ----------------------------
# Polling state for external LLM
# ----------------------------
if "history_mtime" not in st.session_state:
    if os.path.exists(HISTORY_FILE):
        st.session_state["history_mtime"] = os.path.getmtime(HISTORY_FILE)
    else:
        st.session_state["history_mtime"] = 0.0

if "waiting_for_response" not in st.session_state:
    st.session_state["waiting_for_response"] = False

def history_file_changed(path: str) -> bool:
    """Check if history file has changed since last check."""
    if not os.path.exists(path):
        return False

    last = st.session_state.get("history_mtime", 0.0)
    current = os.path.getmtime(path)

    if current != last and os.path.getsize(path) > 0:
        st.session_state["history_mtime"] = current
        return True

    return False
    
# ----------------------------
# Configurations and initiations
# ----------------------------
# List LLM models
llm_model_set = []
with open("../llm_models", 'r') as file:
    for l in file:
        l = l.strip()
        if l:
            llm_model_set.append(l)

# Set LLM configurations
if "model_full" not in st.session_state:
    model_full = os.getenv("LLM_MODEL")
    st.session_state['model_full'] = model_full

# if st.session_state.initialize:
#     st.session_state.initialize = False
#     st.write(" ")
#     st.write(" ")
#     # st.info("Example questions:")
#     st.write(" ")
#     st.write(" ")
#     utils.message_box("Visualize the caffeine molecule at pH = 7",    bg="#ffdddd", border="#ff6b6b", height=5)
#     utils.message_box("What are the PDB IDs for hemoglobin",  bg="#ddffdd", border="#68c96b", height=5)
#     utils.message_box("Blue message",   bg="#dde8ff", border="#6b8cff", height=5)
#     utils.message_box("Yellow message", bg="#fff8cc", border="#e6d04d", height=5)
#     utils.message_box("Orange message", bg="#ffe4cc", border="#ff9d4d", height=5)
#     utils.message_box("Purple message", bg="#f0ddff", border="#c68cff", height=5)
    
# Make folders
os.makedirs(HISTORY_DIR, exist_ok=True)
os.makedirs(SESSION_DIR, exist_ok=True)
os.makedirs(FILE_DIR, exist_ok=True)

# Set user session info
user_name = st.session_state['user']

# if 'OpenAI_API_token' not in st.session_state[user_name]:
# st.session_state[user_name]['OpenAI_API_token'] = os.getenv('OPENAI_API_KEY')
# st.session_state[user_name]['Ollama_API_token'] = os.getenv("OLLAMA_API_KEY")
# st.session_state[user_name]['OpenRouter_API_token'] = os.getenv("OPENROUTER_API_KEY")

session_info = {
    'user_name' : user_name,
    'history_file' : HISTORY_FILE,
    'file_dir' : FILE_DIR,
    'model_full' : st.session_state['model_full'],
    'token_provided' : st.session_state[user_name]['token_provided'],
    'openai_api_key' : st.session_state[user_name]['OpenAI_API_token'],
    'ollama_api_key' : st.session_state[user_name]['Ollama_API_token'],
    'openrouter_api_key' : st.session_state[user_name]['OpenRouter_API_token']
}

# ----------------------------
# Sidebar
# ----------------------------
with st.sidebar:

    user_name = st.session_state['user']
    if user_name in user.USER_FREE_QUESTIONS:
        n_q_left = user.USER_FREE_QUESTIONS[user_name] - st.session_state[user_name]['n_questions']
        st.sidebar.write(f"Logged in as: {user_name}")
        st.sidebar.write(f"Free questions left: {n_q_left}")
    
    st.title("BioChemAIgent")
    # st.subheader("LLM model")
    model_full = st.selectbox("Choose a LLM model", llm_model_set)
    st.session_state['model_full'] = model_full
    st.caption(f"Model: `{model_full}`")
    
    st.markdown("---")
    
    col1, col2, col3 = st.columns([2, 2, 2])
    with col1:
        logout = st.button("Logout")
    with col2:
        clear_chat = st.button("Clear chat")

    chat_name = st.text_input("Chat name", "my_chat").strip()

    col1, col2, col3 = st.columns([2, 2, 2])
    with col1:
        save_chat = st.button("Save chat")
    with col2:
        load_chat = st.button("Load chat")
    with col3:
        delete_chat = st.button("Delete chat")
        
    if logout:
        if "user" in st.session_state:
            del st.session_state["user"]
            st.rerun()

    if clear_chat:
        st.session_state.messages = []
        if os.path.exists(HISTORY_FILE):
            os.remove(HISTORY_FILE)
        st.rerun()

    if save_chat:
        if chat_name:
            history_name = chat_name + ".pkl"
            file_list = os.listdir(HISTORY_DIR)
            if history_name in file_list:
                history_name = chat_name + "_new.pkl"
            os.system(f"cp {HISTORY_FILE} {os.path.join(HISTORY_DIR, history_name)}")
            st.info(f"Saved to {history_name}")

    if load_chat:
        history_name = chat_name + ".pkl"
        src = os.path.join(HISTORY_DIR, history_name)
        if os.path.exists(src):
            os.system(f"cp {src} {HISTORY_FILE}")
            st.info(f"The chat is loaded from {history_name}")
        else:
            st.info(f"The chat history {history_name} does not exist!")

    if delete_chat:
        history_name = chat_name + ".pkl"
        src = os.path.join(HISTORY_DIR, history_name)
        if os.path.exists(src):
            os.system(f"rm {src}")
            
    st.subheader("Chat history")
    with st.container(height=300):
        # Get the history file list sorted by date
        file_list = sorted(os.listdir(HISTORY_DIR), 
                           key=lambda f: os.path.getmtime(os.path.join(HISTORY_DIR, f)), reverse=True)
        file_list = [i.split('.')[0] for i in file_list if i.split('.')[0] not in ['', 'chat_history']]
        for f in file_list:
            st.caption(f)

    st.markdown("---")


    st.subheader("Chat settings")

    plot_scale = st.slider("Plot zoom", 0, 100, 70) / 100
    plot_height = st.slider("Plot height", 0, 1000, 400)        
    
    # auto_confirm_flag = st.toggle("Auto confirm")
    # auto_confirm = False
    
    if st.button("Refresh"):
        st.rerun()
    
    st.markdown("---")

    st.subheader("Files")
    with st.container(height=300):
        # Get the file list sorted by date
        file_list = sorted(os.listdir(FILE_DIR), 
                           key=lambda f: os.path.getmtime(os.path.join(FILE_DIR, f)), reverse=True)
        for f in file_list:
            st.caption(f)
                
    st.markdown("---")

    st.subheader("Access tokens")
    user_name = st.session_state['user']
    openai_token = st.text_input("OpenAI API token", '').strip()
    ollama_token = st.text_input("Ollama API token", '').strip()
    openrouter_name = st.text_input("OpenRouter API token", '').strip()
    if st.button("Submit tokens"):
        st.session_state[user_name]['OpenAI_API_token'] = openai_token
        st.session_state[user_name]['Ollama_API_token'] = ollama_token
        st.session_state[user_name]['OpenRouter_API_token'] = openrouter_name
        st.session_state[user_name]['token_provided'] = True
        st.info(f"Tokes are submited")
        
# ----------------------------
# Display chat history
# ----------------------------
st.session_state.messages = load_history(HISTORY_FILE)

for msg in st.session_state.messages:
    role = msg["role"]
    content = msg["content"]
    if content:
        if role == "user":
            with st.chat_message("user"):
                content = utils.escape_brackets_outside_code(content)
                st.markdown(content)
        elif role == "assistant":
            with st.chat_message("assistant"):
                content = utils.escape_brackets_outside_code(content)
                st.markdown(content)
            if 'html_file' in globals() and html_file:
                st.components.v1.html(
                    f"""
                    <div style='transform:scale({plot_scale});transform-origin:top left;'>
                    {open(html_file).read()}
                    </div>
                    """, 
                    height=plot_height, scrolling=True
                )
                html_file = ''
        elif role == "tool":
            with st.chat_message("assistant", avatar="🛠️"):
                tool_name = msg.get('name', '')
                if tool_name:
                    st.markdown(f"Tool call: `{tool_name}`")
                    # Get html file to render
                    try:
                        if tool_name in [
                            'BioChemAIgent_MCP_Server__bcai__render_structures',
                            'BioChemAIgent_MCP_Server__bcai__interaction_plot'
                        ]:
                            content = utils.str2dict(msg['content'])
                            html_file = content['html_file']
                        else:
                            html_file = ''
                    except:
                        pass

# ----------------------------
# If waiting for the response, poll for file change
# ----------------------------
if st.session_state.get("waiting_for_response", False):
    if history_file_changed(HISTORY_FILE):
        # Response arrived: stop waiting and rerun to show it
        st.session_state["waiting_for_response"] = False
        st.rerun()
    else:
        # Still waiting – show info and check again shortly
        st.info("Waiting for Agent's response...")
        time.sleep(0.1)
        st.rerun()

# ----------------------------
# Get user input and process
# ----------------------------
user_input = st.chat_input("Ask the agent...")

if user_input:
    user_name = st.session_state['user']
    if user.user_question_counter(user_name) or st.session_state[user_name]['token_provided']:
        # Add to history
        save_history(load_history(HISTORY_FILE) + [{"role": "user", "content": user_input}], HISTORY_FILE)

        # Show immediately in UI
        with st.chat_message("user"):
            safe_input = utils.escape_brackets_outside_code(user_input)
            st.markdown(safe_input)

        # Trigger MCP client
        # 1. Save session state
        with open(SESSION_FILE, 'w') as f:
            json.dump(session_info, f, indent=4)

        # 2. Copy history file for the client
        os.system(f"cp {HISTORY_FILE} {os.path.join('work_space', user_name + '.pkl')}")

        # 3. Set waiting flags and record current mtime
        if os.path.exists(HISTORY_FILE):
            st.session_state["history_mtime"] = os.path.getmtime(HISTORY_FILE)
        else:
            st.session_state["history_mtime"] = 0.0

        st.session_state["waiting_for_response"] = True
        st.rerun()

    else:
        st.warning("Please provide an LLM token to continue ...")
