import {
  createContext, useContext, useState, useEffect,
  useCallback, type ReactNode,
} from 'react';
import type { User } from './types';
import { getMe } from './api';

// ----------------------------
// Auth context — token + user, persisted to localStorage
// ----------------------------

interface AuthState {
  token  : string | null;
  user   : User   | null;
  login  : (token: string) => Promise<void>;
  logout : () => void;
}

const AuthCtx = createContext<AuthState>({
  token: null, user: null,
  login: async () => {}, logout: () => {},
});

export function AuthProvider({ children }: { children: ReactNode }) {
  const [token, setToken] = useState<string | null>(
    () => localStorage.getItem('ddp_token'),
  );
  const [user, setUser] = useState<User | null>(null);

  const logout = useCallback(() => {
    if (user?.username) {
      localStorage.removeItem(`ddp_jobs_${user.username}`);
      localStorage.removeItem(`ddp_session_${user.username}`);
    }
    localStorage.removeItem('ddp_token');
    setToken(null);
    setUser(null);
  }, [user]);

  const login = useCallback(async (t: string) => {
    localStorage.setItem('ddp_token', t);
    setToken(t);
    try {
      setUser(await getMe(t));
    } catch {
      logout();
    }
  }, [logout]);

  // Re-hydrate user on page refresh if a token is already in localStorage
  useEffect(() => {
    if (token && !user) {
      getMe(token).then(setUser).catch(logout);
    }
  }, []); // intentionally runs once on mount only

  return (
    <AuthCtx.Provider value={{ token, user, login, logout }}>
      {children}
    </AuthCtx.Provider>
  );
}

export function useAuth() {
  return useContext(AuthCtx);
}

// ----------------------------
// Theme — persisted in localStorage, applied as data-theme on <html>
// ----------------------------

export type Theme = 'dark' | 'light';

export function getTheme(): Theme {
  return (localStorage.getItem('ddp_theme') as Theme) ?? 'dark';
}

export function applyTheme(t: Theme) {
  document.documentElement.setAttribute('data-theme', t);
  localStorage.setItem('ddp_theme', t);
}

export function toggleTheme() {
  applyTheme(getTheme() === 'dark' ? 'light' : 'dark');
}

// Session list is now server-side (GET /sessions) — no localStorage needed.
