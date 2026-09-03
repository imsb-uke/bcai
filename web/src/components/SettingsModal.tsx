import { useState, useEffect, type FormEvent } from 'react';
import { useAuth }                            from '../store';
import { getKeys, updateKeys, deleteAccount }  from '../api';
import '../styles/modal.css';

interface Props { onClose: () => void }

export default function SettingsModal({ onClose }: Props) {
  const { token, logout } = useAuth();

  const [openai,     setOpenai]     = useState('');
  const [openrouter, setOpenrouter] = useState('');
  const [ollama,     setOllama]     = useState('');
  const [saved,      setSaved]      = useState(false);
  const [error,      setError]      = useState('');
  const [busy,       setBusy]       = useState(false);

  // Load masked current values on open
  useEffect(() => {
    if (!token) return;
    getKeys(token).then(k => {
      setOpenai(k.openai_key);
      setOpenrouter(k.openrouter_key);
      setOllama(k.ollama_key);
    }).catch(() => {});
  }, [token]);

  async function submit(e: FormEvent) {
    e.preventDefault();
    if (!token) return;
    setError('');
    setBusy(true);
    try {
      await updateKeys(token, {
        openai_key     : openai,
        openrouter_key : openrouter,
        ollama_key     : ollama,
      });
      setSaved(true);
      setTimeout(() => setSaved(false), 2000);
    } catch (err: unknown) {
      setError(err instanceof Error ? err.message : 'Save failed');
    } finally {
      setBusy(false);
    }
  }

  async function handleDeleteAccount() {
    if (!token) return;
    if (!window.confirm(
      'Delete your account? Your login is removed and the username freed, ' +
      'but your files are kept (folder renamed to <username>_deleted).'
    )) return;
    try {
      await deleteAccount(token);
      logout();   // clears token → returns to the login screen
    } catch (err: unknown) {
      setError(err instanceof Error ? err.message : 'Delete failed');
    }
  }

  return (
    // Click backdrop to close
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal-card" onClick={e => e.stopPropagation()}>

        <div className="modal-header">
          <span className="modal-title">Settings</span>
          <button className="modal-close" onClick={onClose}>✕</button>
        </div>

        <form onSubmit={submit} className="modal-form">

          <div className="modal-section-label">API Keys</div>
          <p className="modal-hint">
            Keys are stored in the database. With a key set, your account bypasses the free-question quota.
            Leave blank to keep the existing value.
          </p>

          <label>
            OpenAI key
            <input
              type        = "password"
              value       = {openai}
              onChange    = {e => setOpenai(e.target.value)}
              placeholder = "sk-…"
              autoComplete= "off"
            />
          </label>

          <label>
            OpenRouter key
            <input
              type        = "password"
              value       = {openrouter}
              onChange    = {e => setOpenrouter(e.target.value)}
              placeholder = "sk-or-…"
              autoComplete= "off"
            />
          </label>

          <label>
            Ollama key
            <input
              type        = "password"
              value       = {ollama}
              onChange    = {e => setOllama(e.target.value)}
              placeholder = "optional"
              autoComplete= "off"
            />
          </label>

          {error  && <div className="modal-error">{error}</div>}
          {saved  && <div className="modal-success">Saved ✓</div>}

          <button type="submit" className="modal-submit" disabled={busy}>
            {busy ? 'Saving…' : 'Save keys'}
          </button>

        </form>

        {/* Danger zone — account deletion keeps files, frees the username */}
        <div className="modal-danger">
          <div className="modal-section-label">Danger zone</div>
          <p className="modal-hint">
            Deleting your account removes your login and frees the username. Your files
            are kept (folder renamed to <code>&lt;username&gt;_deleted</code>).
          </p>
          <button type="button" className="modal-delete" onClick={handleDeleteAccount}>
            Delete account
          </button>
        </div>
      </div>
    </div>
  );
}
