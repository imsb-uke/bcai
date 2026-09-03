import { useEffect, useRef, useState, type FormEvent, type KeyboardEvent } from 'react';
import type { ChatMessage } from '../types';
import Message from './Message';
import logoUrl from '../assets/logo.svg';
import '../styles/chat.css';

interface Props {
  messages  : ChatMessage[];
  onSend    : (text: string) => void;
  onStop    : () => void;
  streaming : boolean;
}

export default function Chat({ messages, onSend, onStop, streaming }: Props) {
  const [input,  setInput]  = useState('');
  const bottomRef           = useRef<HTMLDivElement>(null);
  const textareaRef         = useRef<HTMLTextAreaElement>(null);

  // Scroll to bottom whenever new content arrives
  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  // Auto-grow the textarea as the user types
  useEffect(() => {
    const ta = textareaRef.current;
    if (!ta) return;
    ta.style.height = 'auto';
    ta.style.height = `${Math.min(ta.scrollHeight, 160)}px`;
  }, [input]);

  function submit(e: FormEvent) {
    e.preventDefault();
    const text = input.trim();
    if (!text || streaming) return;
    setInput('');
    onSend(text);
  }

  // Enter sends; Shift+Enter inserts a newline
  function handleKey(e: KeyboardEvent<HTMLTextAreaElement>) {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      submit(e as unknown as FormEvent);
    }
  }

  return (
    <section className="chat-panel">

      <div className="chat-messages">
        {messages.length === 0 && (
          <div className="chat-empty">
            <img src={logoUrl} alt="" className="chat-empty-icon" />
            <div className="chat-empty-title">BioChemAIgent</div>
            <div className="chat-empty-sub">Agentic drug discovery platform</div>
          </div>
        )}
        {messages.map(m => <Message key={m.id} msg={m} />)}
        <div ref={bottomRef} />
      </div>

      <form className="chat-input-bar" onSubmit={submit}>
        <textarea
          ref         = {textareaRef}
          className   = "chat-input"
          value       = {input}
          onChange    = {e => setInput(e.target.value)}
          onKeyDown   = {handleKey}
          placeholder = "Message the agent… (Enter to send, Shift+Enter for new line)"
          rows        = {1}
          disabled    = {streaming}
        />
        {streaming ? (
          <button
            type      = "button"
            className = "chat-send chat-stop"
            onClick   = {onStop}
            title     = "Stop generating"
          >
            ◼
          </button>
        ) : (
          <button
            type      = "submit"
            className = "chat-send"
            disabled  = {!input.trim()}
            title     = "Send (Enter)"
          >
            ↑
          </button>
        )}
      </form>

    </section>
  );
}
