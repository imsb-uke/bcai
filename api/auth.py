import os
from datetime import datetime, timedelta, timezone
from typing import Optional

import bcrypt
from jose import JWTError, jwt
from fastapi import Depends, HTTPException, status
from fastapi.security import OAuth2PasswordBearer
from sqlalchemy.orm import Session

from .db import get_db, get_user, User

# JWT settings — set JWT_SECRET in config/env before deploying
# Generate one with: python -c "import secrets; print(secrets.token_hex(32))"
SECRET_KEY  = os.getenv("JWT_SECRET", "change-me-before-deploying")
ALGORITHM   = "HS256"
TOKEN_HOURS = int(os.getenv("TOKEN_HOURS", "8"))

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/auth/login")


# ----------------------------
# Password helpers
# ----------------------------
def hash_password(plain: str) -> str:
    return bcrypt.hashpw(plain.encode(), bcrypt.gensalt()).decode()

def verify_password(plain: str, hashed: str) -> bool:
    return bcrypt.checkpw(plain.encode(), hashed.encode())


# ----------------------------
# JWT helpers
# ----------------------------
def create_token(username: str) -> str:
    expire = datetime.now(timezone.utc) + timedelta(hours=TOKEN_HOURS)
    return jwt.encode({"sub": username, "exp": expire}, SECRET_KEY, algorithm=ALGORITHM)

def decode_token(token: str) -> Optional[str]:
    """Returns username, or None if the token is invalid or expired."""
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        return payload.get("sub")
    except JWTError:
        return None


# ----------------------------
# FastAPI dependency
# ----------------------------
def current_user(
    token: str = Depends(oauth2_scheme),
    db: Session = Depends(get_db),
) -> User:
    """Inject into any endpoint that requires a logged-in user."""
    username = decode_token(token)
    if not username:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or expired token",
            headers={"WWW-Authenticate": "Bearer"},
        )
    user = get_user(db, username)
    if user is None or not user.is_active:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                            detail="User not found or inactive")
    return user
