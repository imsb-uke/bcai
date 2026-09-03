import { useState, useCallback, useRef, useEffect } from 'react';
import { useAuth } from '../store';
import { streamChat, attachChat, cancelChat, fetchModels, fetchHistory, fetchSessions, fetchNotifications, renameSession, deleteSession } from '../api';
import type { ChatMessage, Session, Job, JobNotification, SSEEvent } from '../types';
import Sidebar  from './Sidebar';
import Chat     from './Chat';
import JobPanel from './JobPanel';
import logoUrl  from '../assets/logo.svg';
import '../styles/layout.css';

// ----------------------------
// Helpers
// ----------------------------

let _counter = 0;
const newId = () => `m${++_counter}-${Date.now()}`;

function blankAssistant(): ChatMessage {
  return { id: newId(), role: 'assistant', content: '', toolCalls: [], streaming: true, error: false };
}

// Map a raw history message to a ChatMessage for display.
// role="tool" becomes a standalone tool-call row; user/assistant become bubbles
// with no tool badge (the separate tool rows render the calls).
function msgFromRaw(
  raw     : Record<string, unknown>,
  argsMap : Record<string, Record<string, unknown>> = {},
): ChatMessage {
  const role    = raw.role as string;
  const content = (raw.content as string | null) ?? '';

  if (role === 'tool') {
    const tcId = raw.tool_call_id as string | undefined;
    const args = tcId ? argsMap[tcId] : undefined;
    return {
      id       : newId(),
      role     : 'tool',
      content,
      toolCalls: [{ name: (raw.name as string) ?? 'tool', args }],
      streaming: false,
      error    : false,
    };
  }

  // Assistant/user bubbles carry no tool badge — tool rows render the calls
  return {
    id       : newId(),
    role     : role === 'user' ? 'user' : 'assistant',
    content,
    toolCalls: [],
    streaming: false,
    error    : false,
  };
}


// ----------------------------
// ChatPage — holds all runtime state
// ----------------------------

export default function ChatPage() {
  const { token, user } = useAuth();

  const [sessions,      setSessions]      = useState<Session[]>([]);
  const [sessionId,     setSessionId]     = useState<string>(() => crypto.randomUUID());
  const [messages,      setMessages]      = useState<ChatMessage[]>([]);
  const [models,        setModels]        = useState<string[]>([]);
  const [selectedModel, setSelectedModel] = useState<string>('');
  const [jobs,          setJobs]          = useState<Job[]>([]);
  const [notifQueue,    setNotifQueue]    = useState<JobNotification[]>([]);
  const [filesRefresh,  setFilesRefresh]  = useState(0);
  const [mobileNav,     setMobileNav]     = useState<'none' | 'sidebar' | 'panel'>('none');
  const abortRef          = useRef<AbortController | null>(null);
  const savedSessionRef   = useRef<string | null>(null);
  const jobsLoaded        = useRef(false);
  const sessionRestored   = useRef(false);

  // Rehydrate jobs and capture saved session once the user is known
  useEffect(() => {
    if (!user || jobsLoaded.current) return;
    jobsLoaded.current = true;
    try {
      const raw = localStorage.getItem(`ddp_jobs_${user.username}`);
      if (raw) setJobs(JSON.parse(raw));
    } catch {}
    savedSessionRef.current = localStorage.getItem(`ddp_session_${user.username}`);
  }, [user]);

  // Persist jobs to localStorage (only after initial load to avoid overwriting with [])
  useEffect(() => {
    if (!user || !jobsLoaded.current) return;
    localStorage.setItem(`ddp_jobs_${user.username}`, JSON.stringify(jobs));
  }, [jobs, user]);

  // Persist the active session (only after restore to avoid clearing the saved key)
  useEffect(() => {
    if (!user || !sessionRestored.current) return;
    if (sessionId) localStorage.setItem(`ddp_session_${user.username}`, sessionId);
    else           localStorage.removeItem(`ddp_session_${user.username}`);
  }, [sessionId, user]);

  // Restore the last open session once both user and sessions list are ready
  useEffect(() => {
    if (!user || sessions.length === 0 || sessionRestored.current) return;
    sessionRestored.current = true;
    const saved = savedSessionRef.current;
    if (!saved) return;
    const match = sessions.find(s => s.id === saved);
    if (match) switchSession(match);
  }, [user, sessions]);

  // Load models and session list once on mount
  const didInit = useRef(false);
  useEffect(() => {
    if (didInit.current || !token) return;
    didInit.current = true;

    fetchModels(token).then(({ models: ms, default: def }) => {
      setModels(ms);
      setSelectedModel(s => s || def);
    });
    refreshSessions();
  }, [token]);

  // Poll the session list every few seconds — keeps each session's "running"
  // badge live (including ones started from a different tab/device) and
  // catches turns that finished in the background while unwatched.
  useEffect(() => {
    if (!token) return;
    const id = setInterval(refreshSessions, 4000);
    return () => clearInterval(id);
  }, [token]);

  // Poll for job-completion notifications every 5 seconds
  useEffect(() => {
    if (!token) return;
    const id = setInterval(async () => {
      const items = await fetchNotifications(token);
      if (items.length === 0) return;
      setNotifQueue(q => [...q, ...items]);
      // Mark matching jobs as done in the panel
      setJobs(prev => prev.map(j => {
        const match = items.find(n => n.job_id === j.id);
        return match ? { ...j, status: match.status as Job['status'], result: match.message } : j;
      }));
    }, 5000);
    return () => clearInterval(id);
  }, [token]);

  // ----------------------------
  // Shared SSE consumption — used both for sending a new message and for
  // reattaching to a turn that's already running (see attachRunning below).
  // Both paths render into the same `messages` state the same way.
  // ----------------------------

  async function consumeStream(
    gen        : AsyncGenerator<SSEEvent>,
    base       : ChatMessage[],
    turn       : ChatMessage[],
    controller : AbortController,
  ) {
    // Only paint if this stream is still the active one (guards against
    // late tokens overwriting a session the user has switched away to).
    const render = () => {
      if (abortRef.current !== controller) return;
      setMessages([...base, ...turn.map(m => ({ ...m }))]);
    };

    render();   // show the spinner immediately

    try {
      for await (const event of gen) {
        const last = turn[turn.length - 1];

        if (event.type === 'token') {
          // Append into the current assistant bubble (create one if needed)
          if (last.role === 'assistant') {
            last.content += event.text;
          } else {
            turn.push({ ...blankAssistant(), content: event.text });
          }
          render();

        } else if (event.type === 'tool') {
          // Finalise current assistant; drop it if it never received text
          if (last.role === 'assistant') {
            if (last.content === '') turn.pop();
            else                     last.streaming = false;
          }
          // Tool row, then a fresh assistant bubble for the next round
          turn.push({
            id: newId(), role: 'tool', content: event.content,
            toolCalls: [{ name: event.name, args: event.args }],
            streaming: false, error: false,
          });
          turn.push(blankAssistant());
          render();

          // Register an async job from the real tool result (not a text regex)
          registerJobFromResult(event.name, event.content);
          setFilesRefresh(n => n + 1);

        } else if (event.type === 'done') {
          // Tokens already filled the bubble — just stop streaming
          const a = turn[turn.length - 1];
          if (a.role === 'assistant') {
            a.streaming = false;
            if (a.content === '') a.content = event.content;
          }
          render();
          setFilesRefresh(n => n + 1);

        } else if (event.type === 'error') {
          if (last.role === 'assistant') {
            last.streaming = false;
            last.error     = true;
            last.content   = event.message;
          }
          render();
        }
      }
    } catch (err) {
      const a = turn[turn.length - 1];
      if (a?.role === 'assistant') {
        a.streaming = false;
        if ((err as Error)?.name === 'AbortError') {
          // We only stopped watching — the run itself is unaffected unless
          // cancelStreaming() explicitly told the server to stop too. Show
          // whatever streamed so far without claiming it was cancelled.
          if (a.content === '') turn.pop();
        } else {
          a.error = true; a.content = String(err);
        }
      }
      render();
      abortRef.current = null;
      return;
    }
    abortRef.current = null;
    refreshSessions();
  }

  // ----------------------------
  // Send a message + stream the response
  // ----------------------------

  const sendMessage = useCallback(async (text: string) => {
    if (!token) return;
    if (abortRef.current) return;   // a response is already streaming — ignore

    const userMsg: ChatMessage = {
      id: newId(), role: 'user', content: text,
      toolCalls: [], streaming: false, error: false,
    };
    const history = [...messages, userMsg].map(m => ({ role: m.role, content: m.content }));

    const base = [...messages, userMsg];
    const turn : ChatMessage[] = [blankAssistant()];

    const controller = new AbortController();
    abortRef.current = controller;

    await consumeStream(streamChat(token, history, selectedModel, sessionId, controller.signal), base, turn, controller);
  }, [token, messages, selectedModel, sessionId]);

  // Reattach to a turn the sidebar shows as still running for this session —
  // no new message is sent, we just resume watching (replaying everything
  // emitted so far, then live events).
  async function attachRunning(id: string, base: ChatMessage[]) {
    if (!token) return;
    const turn: ChatMessage[] = [blankAssistant()];
    const controller = new AbortController();
    abortRef.current = controller;
    await consumeStream(attachChat(token, id, controller.signal), base, turn, controller);
  }

  // Stop watching only — the turn keeps running server-side. Used whenever
  // the view moves away from the session currently streaming.
  function stopWatching() {
    abortRef.current?.abort();
    abortRef.current = null;
  }

  // Actually stop the turn server-side — the explicit "Stop generating" action.
  function cancelStreaming() {
    if (token) void cancelChat(token, sessionId);
    stopWatching();
  }

  // Register a job by reading the real job_id from a tool result JSON
  function registerJobFromResult(toolName: string, content: string) {
    let data: Record<string, unknown>;
    try { data = JSON.parse(content); } catch { return; }
    const jobId = (data.job_id ?? data.id) as string | undefined;
    if (!jobId) return;                     // not an async job
    const status = (data.status as string) ?? 'running';
    setJobs(prev => {
      if (prev.some(j => j.id === jobId)) return prev;   // already tracked
      return [...prev, {
        id: jobId, tool: toolName.split('__').pop() ?? toolName,
        status: status as Job['status'], result: null,
        sessionId, notified: false,
      }];
    });
  }


  // ----------------------------
  // Session switching — load history from the backend, then reattach if
  // a turn is still running for it
  // ----------------------------

  async function switchSession(s: Session) {
    stopWatching();   // stop watching any in-flight stream — never stops it server-side
    setSessionId(s.id);
    setMessages([]);
    if (s.model) setSelectedModel(s.model);
    if (!token) return;

    const raw = await fetchHistory(token, s.id);
    // Build tool_call_id → parsed args map from assistant messages so tool rows
    // can display the actual arguments that were sent.
    const argsMap: Record<string, Record<string, unknown>> = {};
    for (const m of raw) {
      if (m.role === 'assistant' && Array.isArray(m.tool_calls)) {
        for (const tc of m.tool_calls as Array<Record<string, unknown>>) {
          const id   = tc.id as string | undefined;
          const fn   = tc.function as Record<string, unknown> | undefined;
          const argStr = fn?.arguments as string | undefined;
          if (id && argStr) {
            try { argsMap[id] = JSON.parse(argStr); } catch { /* skip */ }
          }
        }
      }
    }
    const base = raw
      .filter(m => m.role !== 'system')
      // Skip assistant messages that are pure tool-routing steps (no visible content)
      .filter(m => !(m.role === 'assistant' && !m.content && Array.isArray(m.tool_calls)))
      .map(m => msgFromRaw(m, argsMap));
    setMessages(base);

    if (s.running) {
      await attachRunning(s.id, base);
    }
  }

  function newChat() {
    stopWatching();   // stop watching any in-flight stream — never stops it server-side
    setSessionId(crypto.randomUUID());
    setMessages([]);
  }

  function refreshSessions() {
    if (!token) return;
    fetchSessions(token).then(list => {
      setSessions(list.map(s => ({
        id: s.id, title: s.title, model: s.model, created: s.modified, running: s.running,
      })));
    });
  }

  async function handleRename(s: Session) {
    if (!token) return;
    const title = window.prompt('Rename conversation:', s.title);
    if (!title || title === s.title) return;
    await renameSession(token, s.id, title);
    refreshSessions();
  }

  async function handleDelete(s: Session) {
    if (!token) return;
    if (!window.confirm(`Delete "${s.title}"? This cannot be undone.`)) return;
    await deleteSession(token, s.id);
    if (s.id === sessionId) newChat();
    refreshSessions();
  }


  // ----------------------------
  // "Tell agent" — send job result back into chat
  // ----------------------------

  const tellAgent = useCallback((job: Job) => {
    void sendMessage(`Job ${job.id} is complete. Result: ${job.result ?? 'done'}`);
    setJobs(prev => prev.map(j =>
      j.id === job.id ? { ...j, notified: true } : j,
    ));
  }, [sendMessage]);

  function dismissJob(job: Job) {
    setJobs(prev => prev.filter(j => j.id !== job.id));
  }

  function acknowledgeJob(job: Job) {
    setJobs(prev => prev.map(j =>
      j.id === job.id ? { ...j, notified: true } : j,
    ));
  }

  function dismissNotif() {
    setNotifQueue(q => q.slice(1));
  }

  function notifTellAgent(n: JobNotification) {
    void sendMessage(`Job ${n.job_id} (${n.type.split('__').pop()}) finished with status: ${n.status}.`);
    dismissNotif();
  }

  const isStreaming = messages.some(m => m.streaming);
  const topNotif    = notifQueue[0] ?? null;

  return (
    <div className={`app-shell nav-${mobileNav}`}>

      {/* Mobile-only top bar — toggles the side panels */}
      <div className="mobile-bar">
        <button className="mobile-btn" onClick={() => setMobileNav(n => n === 'sidebar' ? 'none' : 'sidebar')}>☰</button>
        <span className="mobile-title">
          <img src={logoUrl} alt="" className="brand-logo" />
          BioChemAIgent
        </span>
        <button className="mobile-btn" onClick={() => setMobileNav(n => n === 'panel' ? 'none' : 'panel')}>🗂</button>
      </div>

      {/* Backdrop closes any open panel on mobile */}
      {mobileNav !== 'none' && (
        <div className="mobile-backdrop" onClick={() => setMobileNav('none')} />
      )}

      {topNotif && (
        <div className="notif-banner">
          <span className="notif-msg">✅ {topNotif.message}</span>
          <div className="notif-actions">
            <button className="notif-btn-tell" onClick={() => notifTellAgent(topNotif)}>
              Tell agent ↑
            </button>
            <button className="notif-btn-dismiss" onClick={dismissNotif}>✕</button>
          </div>
          {notifQueue.length > 1 && (
            <span className="notif-more">(+{notifQueue.length - 1} more)</span>
          )}
        </div>
      )}
      <Sidebar
        sessions      = {sessions}
        currentId     = {sessionId}
        models        = {models}
        selectedModel = {selectedModel}
        onModelChange = {setSelectedModel}
        onSelect      = {s => { void switchSession(s); setMobileNav('none'); }}
        onNewChat     = {() => { newChat(); setMobileNav('none'); }}
        onRename      = {handleRename}
        onDelete      = {handleDelete}
      />
      <Chat
        messages  = {messages}
        onSend    = {sendMessage}
        onStop    = {cancelStreaming}
        streaming = {isStreaming}
      />
      <JobPanel
        jobs           = {jobs}
        onTellAgent    = {tellAgent}
        onDismissJob   = {dismissJob}
        onAcknowledge  = {acknowledgeJob}
        filesRefreshKey= {filesRefresh}
      />
    </div>
  );
}
