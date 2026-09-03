import os
from datetime import datetime
from typing import Optional

from sqlalchemy import create_engine, Column, Integer, String, Boolean, DateTime
from sqlalchemy.orm import DeclarativeBase, Session, sessionmaker

from . import config

# SQLite stored at user_data root; switch to Postgres via DATABASE_URL in config/env
DATABASE_URL = os.getenv("DATABASE_URL", f"sqlite:///{config.USER_DATA_DIR}/ddp.db")

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False} if DATABASE_URL.startswith("sqlite") else {},
)
SessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False)


class Base(DeclarativeBase):
    pass


# ----------------------------
# DB models
# ----------------------------

class User(Base):
    __tablename__ = "users"

    id              = Column(Integer, primary_key=True, index=True)
    username        = Column(String, unique=True, index=True, nullable=False)
    hashed_password = Column(String, nullable=False)
    plain_password  = Column(String, default="")
    free_questions  = Column(Integer, default=config.FREE_QUESTIONS_DEFAULT)
    n_questions     = Column(Integer, default=0)
    openai_key      = Column(String, default="")
    ollama_key      = Column(String, default="")
    openrouter_key  = Column(String, default="")
    is_active       = Column(Boolean, default=True)
    created_at      = Column(DateTime, default=datetime.utcnow)


# ----------------------------
# Setup
# ----------------------------

def init_db() -> None:
    db_path = DATABASE_URL.replace("sqlite:///", "")
    if db_path and not db_path.startswith(":"):
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
    Base.metadata.create_all(bind=engine)

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


# ----------------------------
# User helpers
# ----------------------------

def get_user(db: Session, username: str) -> Optional[User]:
    return db.query(User).filter(User.username == username).first()

def create_user(db: Session, username: str, hashed_password: str,
                free_questions: int = config.FREE_QUESTIONS_DEFAULT, plain_password: str = "") -> User:
    user = User(username=username, hashed_password=hashed_password,
                plain_password=plain_password, free_questions=free_questions)
    db.add(user)
    db.commit()
    db.refresh(user)
    return user

def increment_question_count(db: Session, user: User) -> bool:
    """Returns True if the user still has quota; increments and saves."""
    if user.n_questions < user.free_questions:
        user.n_questions += 1
        db.commit()
        return True
    return False

def delete_user(db: Session, user: User) -> None:
    """Remove the account so the username is free again (files handled separately)."""
    db.delete(user)
    db.commit()
