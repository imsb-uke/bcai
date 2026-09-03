import { useState } from 'react';
import { useAuth }  from '../store';
import { submitFeedback } from '../api';
import '../styles/modal.css';

const OPTIONS = [
  'I did not get what I wanted',
  'The agent job was good but needed a lot of instructions',
  'The agent job was good with reasonable instructions',
  'The agent worked better than I expected',
] as const;

interface Props {
  sessionId : string | undefined;
  onClose   : () => void;
}

export default function FeedbackModal({ sessionId, onClose }: Props) {
  const { token }              = useAuth();
  const [selected, setSelected] = useState('');
  const [status,   setStatus]   = useState<'idle' | 'saving' | 'done' | 'error'>('idle');

  async function handleSubmit() {
    if (!selected || !token || !sessionId) return;
    setStatus('saving');
    try {
      await submitFeedback(token, sessionId, selected);
      setStatus('done');
    } catch {
      setStatus('error');
    }
  }

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal-card" onClick={e => e.stopPropagation()}>

        <div className="modal-header">
          <span className="modal-title">Session Feedback</span>
          <button className="modal-close" onClick={onClose}>✕</button>
        </div>

        <p className="modal-hint">
          How has this conversation gone <strong>so far</strong>?
          Your rating is saved privately and never sent to the AI.
        </p>

        {!sessionId && (
          <div className="modal-hint" style={{ color: 'var(--error)' }}>
            No active session — start or open a conversation first.
          </div>
        )}

        <div className="feedback-options" style={{ opacity: sessionId ? 1 : 0.4, pointerEvents: sessionId ? 'auto' : 'none' }}>
          {OPTIONS.map(opt => (
            <label key={opt} className={`feedback-option${selected === opt ? ' selected' : ''}`}>
              <input
                type    = "radio"
                name    = "feedback"
                value   = {opt}
                checked = {selected === opt}
                onChange= {() => { setSelected(opt); setStatus('idle'); }}
              />
              {opt}
            </label>
          ))}
        </div>

        {status === 'done' && (
          <div className="modal-success">Thank you — feedback recorded.</div>
        )}
        {status === 'error' && (
          <div className="modal-error">Something went wrong. Please try again.</div>
        )}

        <button
          className = "modal-submit"
          onClick   = {handleSubmit}
          disabled  = {!selected || !sessionId || status === 'saving' || status === 'done'}
        >
          {status === 'saving' ? 'Saving…' : status === 'done' ? 'Submitted ✓' : 'Submit'}
        </button>

      </div>
    </div>
  );
}
