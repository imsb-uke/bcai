import os
import pandas as pd
import streamlit as st

USER_FILE = os.environ["USER_FILE"]
VALID_USERS = {}
USER_FREE_QUESTIONS = {}
N_QUESTIONS = {}
DEFULT_USER_FREE_QUESTIONS = 5


def read_user_file():
    global VALID_USERS, USER_FREE_QUESTIONS, N_QUESTIONS
    df = pd.read_csv(USER_FILE)
    VALID_USERS = {user : str(password) for user, password in zip(df['username'], df['password'])}
    USER_FREE_QUESTIONS = {user : int(free_questions) for user, free_questions in zip(df['username'], df['free_questions'])}
    N_QUESTIONS = {user : int(free_questions) for user, free_questions in zip(df['username'], df['n_questions'])}

def login():
    read_user_file()
    
    # Initialize state flag
    if "show_register" not in st.session_state:
        st.session_state.show_register = False

    col1, col2, col3 = st.columns([2, 2, 2])
    with col2: 
        st.markdown("### Login")

        with st.form("login_form"):
            username = st.text_input("Username").strip()
            password = st.text_input("Password", type="password").strip()
            submit = st.form_submit_button("Login")

        # Button outside form to switch to register mode
        new_user = st.button("New user")

        if submit:
            if username in VALID_USERS and VALID_USERS[username] == password:
                st.session_state["user"] = username
                st.session_state[username] = {
                    'n_questions' : N_QUESTIONS[username] if username in N_QUESTIONS else 0,  
                    'OpenAI_API_token' : os.getenv('OPENAI_API_KEY'),
                    'Ollama_API_token' : os.getenv("OLLAMA_API_KEY"),
                    'OpenRouter_API_token' : os.getenv("OPENROUTER_API_KEY"),
                    'ESM3_API_token' : None,
                    'AF3_API_token' : None,
                    'token_provided' : False,
                }
                st.success(f"Logged in as {username}")
                st.session_state.initialize = True
                st.rerun()
            else:
                st.error("Invalid credentials")

        # When 'New user' is clicked, turn on register mode
        if new_user:
            st.session_state.show_register = True

        register = False
        if st.session_state.show_register:
            st.markdown("### Register")
            with st.form("register_form"):
                reg_username = st.text_input("Set Username:").strip()
                reg_password = st.text_input("Set Password", type="password").strip()
                register = st.form_submit_button("Register")

        if register:
            print(f"New user is added {reg_username}")
            os.system(f"echo '{reg_username},{reg_password},{DEFULT_USER_FREE_QUESTIONS},0' >> {USER_FILE}")
            st.success("New user registered")
            st.session_state.show_register = False
            st.rerun()


def user_question_counter(username):
    df = pd.read_csv(USER_FILE)  
    # add a try except here
    if st.session_state[username]['n_questions'] < USER_FREE_QUESTIONS[username]:
        st.session_state[username]['n_questions'] += 1
        df.loc[df['username'] == username, 'n_questions'] = st.session_state[username]['n_questions']
        df.to_csv(USER_FILE, index=False)
        return True
    else:
        return False 
