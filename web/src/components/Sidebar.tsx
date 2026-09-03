import { useState }       from 'react';
import type { Session }   from '../types';
import { useAuth, toggleTheme, getTheme } from '../store';
import SettingsModal       from './SettingsModal';
import FeedbackModal       from './FeedbackModal';
import logoUrl             from '../assets/favicon.svg';
import '../styles/sidebar.css';

interface Props {
  sessions      : Session[];
  currentId     : string | undefined;
  models        : string[];
  selectedModel : string;
  onModelChange : (m: string) => void;
  onSelect      : (s: Session) => void;
  onNewChat     : () => void;
  onRename      : (s: Session) => void;
  onDelete      : (s: Session) => void;
}

export default function Sidebar({
  sessions, currentId, models, selectedModel,
  onModelChange, onSelect, onNewChat, onRename, onDelete,
}: Props) {
  const { user, logout } = useAuth();
  const [theme,         setTheme]         = useState(getTheme());
  const [showSettings,  setShowSettings]  = useState(false);
  const [showFeedback,  setShowFeedback]  = useState(false);

  function handleThemeToggle() {
    toggleTheme();
    setTheme(t => t === 'dark' ? 'light' : 'dark');
  }

  return (
    <aside className="sidebar">

      {/* Header: brand + controls */}
      <div className="sidebar-header">
        <span className="sidebar-brand">
          <img src={logoUrl} alt="" className="brand-logo" />
          BioChemAIgent
        </span>
        <div className="sidebar-header-actions">
          {/* Theme toggle */}
          <button
            className = "btn-icon"
            onClick   = {handleThemeToggle}
            title     = {theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode'}
          >
            {theme === 'dark' ? '☀' : '☾'}
          </button>
          {/* New chat — clears messages and starts a fresh session */}
          <button className="btn-icon" onClick={onNewChat} title="New conversation">
            +
          </button>
        </div>
      </div>

      {showSettings && <SettingsModal onClose={() => setShowSettings(false)} />}
      {showFeedback && (
        <FeedbackModal sessionId={currentId} onClose={() => setShowFeedback(false)} />
      )}

      {/* Model dropdown populated from GET /models */}
      <div className="sidebar-section">
        <span className="sidebar-label">Model</span>
        <select
          className="model-select"
          value={selectedModel}
          onChange={e => onModelChange(e.target.value)}
        >
          {/* Fallback option shown while the model list is still loading */}
          {models.length === 0 && (
            <option value={selectedModel}>{selectedModel || 'Loading…'}</option>
          )}
          {models.map(m => (
            <option key={m} value={m}>{m}</option>
          ))}
        </select>
      </div>

      {/* Session history */}
      <nav className="session-list">
        <span className="sidebar-label">Recent</span>
        {sessions.length === 0 && (
          <div className="session-empty">No conversations yet</div>
        )}
        {sessions.map(s => (
          <div
            key       = {s.id}
            className = {`session-item${s.id === currentId ? ' active' : ''}`}
            onClick   = {() => onSelect(s)}
            title     = {s.id}
          >
            {s.running && (
              <span className="session-running-dot" title="Still generating…" />
            )}
            <span className="session-title">{s.title || 'Untitled'}</span>
            <div className="session-actions">
              <button className="session-action" title="Rename"
                      onClick={e => { e.stopPropagation(); onRename(s); }}>✎</button>
              <button className="session-action session-action-del" title="Delete"
                      onClick={e => { e.stopPropagation(); onDelete(s); }}>✕</button>
            </div>
          </div>
        ))}
      </nav>

      {/* User info + sign-out */}
      {user && (
        <div className="sidebar-footer">
          <div className="user-info">
            <span className="user-name">{user.username}</span>
            <span className="user-quota">
              {user.n_questions} / {user.free_questions} queries
            </span>
          </div>
          <div className="footer-actions">
            <button className="btn-icon" onClick={() => setShowSettings(true)} title="Settings" style={{ fontSize: '22px' }}>
              ⚙
            </button>
            <button
              className = "btn-icon"
              onClick   = {() => setShowFeedback(true)}
              title     = "Rate this session"
            >
              ★
            </button>
            <button className="btn-icon" onClick={logout} title="Sign out" style={{ fontSize: '18px', marginLeft: 'auto' }}>
              ⏻
            </button>
          </div>
        </div>
      )}

    </aside>
  );
}
