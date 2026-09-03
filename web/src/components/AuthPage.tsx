import { useState, type FormEvent } from 'react';
import { login as apiLogin, register as apiRegister } from '../api';
import { useAuth } from '../store';
import logoUrl from '../assets/logo.svg';
import '../styles/auth.css';

export default function AuthPage() {
  const { login } = useAuth();

  const [mode,     setMode]     = useState<'login' | 'register'>('login');
  const [username, setUsername] = useState('');
  const [password, setPassword] = useState('');
  const [error,    setError]    = useState('');
  const [busy,     setBusy]     = useState(false);

  async function submit(e: FormEvent) {
    e.preventDefault();
    setError('');
    setBusy(true);
    try {
      if (mode === 'register') {
        await apiRegister(username, password);
        // fall through to auto-login after registration
      }
      const token = await apiLogin(username, password);
      await login(token);
    } catch (err: unknown) {
      setError(err instanceof Error ? err.message : 'Something went wrong');
    } finally {
      setBusy(false);
    }
  }

  function switchMode(m: 'login' | 'register') {
    setMode(m);
    setError('');
  }

  return (
    <div className="auth-shell">
      <div className="auth-card">

        <div className="auth-logo">
          <img src={logoUrl} alt="" className="auth-logo-icon" />
          <div className="auth-logo-text-wrap">
            <span className="auth-logo-text">BioChemAIgent</span>
            <span className="auth-logo-sub">Agentic drug discovery platform</span>
          </div>
        </div>

        <div className="auth-tabs">
          <button
            className={mode === 'login' ? 'active' : ''}
            onClick={() => switchMode('login')}
          >
            Sign in
          </button>
          <button
            className={mode === 'register' ? 'active' : ''}
            onClick={() => switchMode('register')}
          >
            Register
          </button>
        </div>

        <form onSubmit={submit} className="auth-form">
          <label>
            Username
            <input
              type="text"
              value={username}
              onChange={e => setUsername(e.target.value)}
              autoFocus
              required
            />
          </label>
          <label>
            Password
            <input
              type="password"
              value={password}
              onChange={e => setPassword(e.target.value)}
              required
            />
          </label>

          {error && <div className="auth-error">{error}</div>}

          <button type="submit" className="auth-submit" disabled={busy}>
            {busy
              ? 'Please wait…'
              : mode === 'login' ? 'Sign in' : 'Create account'}
          </button>
        </form>

      </div>
    </div>
  );
}
