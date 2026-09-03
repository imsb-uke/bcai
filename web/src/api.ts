import type { User, SSEEvent, FileEntry, JobNotification } from './types';

// nginx proxies /api/* → backend at :8000, stripping the /api prefix
const BASE = '/api';

// ----------------------------
// Auth
// ----------------------------

export async function login(username: string, password: string): Promise<string> {
  // Backend uses OAuth2PasswordRequestForm — must be form-encoded, not JSON
  const res = await fetch(`${BASE}/auth/login`, {
    method : 'POST',
    headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
    body   : new URLSearchParams({ username, password }).toString(),
  });
  if (!res.ok) {
    const err = await res.json().catch(() => ({})) as Record<string, unknown>;
    throw new Error((err.detail as string) ?? 'Login failed');
  }
  const data = await res.json() as Record<string, unknown>;
  return data.access_token as string;
}

export async function register(username: string, password: string): Promise<void> {
  const res = await fetch(`${BASE}/auth/register`, {
    method : 'POST',
    headers: { 'Content-Type': 'application/json' },
    body   : JSON.stringify({ username, password }),
  });
  if (!res.ok) {
    const err = await res.json().catch(() => ({})) as Record<string, unknown>;
    throw new Error((err.detail as string) ?? 'Registration failed');
  }
}

export async function getMe(token: string): Promise<User> {
  const res = await fetch(`${BASE}/auth/me`, {
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) throw new Error('Session expired');
  return res.json() as Promise<User>;
}

// ----------------------------
// Models
// ----------------------------

export async function fetchModels(
  token: string,
): Promise<{ models: string[]; default: string }> {
  const res = await fetch(`${BASE}/models`, {
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) return { models: [], default: '' };
  return res.json() as Promise<{ models: string[]; default: string }>;
}

// ----------------------------
// Sessions + History
// ----------------------------

interface SessionInfo { id: string; title: string; model: string; modified: number; running: boolean }

export async function fetchSessions(token: string): Promise<SessionInfo[]> {
  const res = await fetch(`${BASE}/sessions`, {
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) return [];
  const data = await res.json() as { sessions: SessionInfo[] };
  return data.sessions ?? [];
}

export async function fetchHistory(
  token     : string,
  sessionId : string,
): Promise<Array<Record<string, unknown>>> {
  const res = await fetch(`${BASE}/history/${sessionId}`, {
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) return [];
  const data = await res.json() as { messages: Array<Record<string, unknown>> };
  return data.messages ?? [];
}

export async function renameSession(token: string, sessionId: string, title: string): Promise<void> {
  const res = await fetch(`${BASE}/history/${sessionId}`, {
    method : 'PATCH',
    headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${token}` },
    body   : JSON.stringify({ title }),
  });
  if (!res.ok) throw new Error('Rename failed');
}

export async function deleteSession(token: string, sessionId: string): Promise<void> {
  const res = await fetch(`${BASE}/history/${sessionId}`, {
    method : 'DELETE',
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) throw new Error('Delete failed');
}

// ----------------------------
// API keys
// ----------------------------

export async function getKeys(
  token: string,
): Promise<{ openai_key: string; openrouter_key: string; ollama_key: string }> {
  const res = await fetch(`${BASE}/auth/keys`, {
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) throw new Error('Failed to load keys');
  return res.json() as Promise<{ openai_key: string; openrouter_key: string; ollama_key: string }>;
}

export async function updateKeys(
  token : string,
  keys  : { openai_key?: string; openrouter_key?: string; ollama_key?: string },
): Promise<void> {
  const res = await fetch(`${BASE}/auth/keys`, {
    method : 'PATCH',
    headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${token}` },
    body   : JSON.stringify(keys),
  });
  if (!res.ok) throw new Error('Failed to update keys');
}

export async function deleteAccount(token: string): Promise<void> {
  const res = await fetch(`${BASE}/auth/account`, {
    method : 'DELETE',
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) throw new Error('Failed to delete account');
}

// ----------------------------
// Notifications
// ----------------------------

export async function fetchNotifications(token: string): Promise<JobNotification[]> {
  const res = await fetch(`${BASE}/notifications`, {
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) return [];
  const data = await res.json() as { notifications: JobNotification[] };
  return data.notifications ?? [];
}

// ----------------------------
// File manager
// ----------------------------

export async function listFiles(token: string): Promise<FileEntry[]> {
  const res = await fetch(`${BASE}/files/list`, {
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) return [];
  const data = await res.json() as { files: FileEntry[] };
  return data.files ?? [];
}

export async function uploadFile(token: string, file: File): Promise<void> {
  const form = new FormData();
  form.append('file', file);
  const res = await fetch(`${BASE}/files/upload`, {
    method : 'POST',
    headers: { Authorization: `Bearer ${token}` },
    body   : form,
  });
  if (!res.ok) throw new Error('Upload failed');
}

export async function deleteFile(token: string, path: string): Promise<void> {
  const res = await fetch(`${BASE}/files/${path}`, {
    method : 'DELETE',
    headers: { Authorization: `Bearer ${token}` },
  });
  if (!res.ok) throw new Error('Delete failed');
}

export function fileDownloadUrl(path: string): string {
  return `${BASE}/files/${path}`;
}

// ----------------------------
// Feedback
// ----------------------------

export async function submitFeedback(
  token     : string,
  sessionId : string,
  rating    : string,
): Promise<void> {
  const res = await fetch(`${BASE}/feedback/${sessionId}`, {
    method : 'POST',
    headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${token}` },
    body   : JSON.stringify({ rating }),
  });
  if (!res.ok) throw new Error('Failed to submit feedback');
}

// ----------------------------
// Chat — SSE via fetch + ReadableStream
// EventSource cannot send custom headers, so we parse the stream manually.
// Turns run decoupled from the HTTP request server-side (see turns.py /
// orchestrator.py): streamChat starts one (409s if already running for this
// session), attachChat only reattaches to watch one already in flight — both
// share the same frame-parsing loop since the wire format is identical.
// ----------------------------

async function* parseSSE(res: Response): AsyncGenerator<SSEEvent> {
  if (!res.ok) {
    const err = await res.json().catch(() => ({})) as Record<string, unknown>;
    yield { type: 'error', message: (err.detail as string) ?? `HTTP ${res.status}` };
    return;
  }

  const reader  = res.body!.getReader();
  const decoder = new TextDecoder();
  let   buffer  = '';

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;

    buffer += decoder.decode(value, { stream: true });

    // SSE frames are separated by double newlines: "data: {...}\n\n"
    const frames = buffer.split('\n\n');
    buffer = frames.pop() ?? '';   // keep the incomplete trailing frame

    for (const frame of frames) {
      const line = frame.trim();
      if (!line.startsWith('data:')) continue;
      try {
        yield JSON.parse(line.slice(5).trim()) as SSEEvent;
      } catch {
        // malformed frame — skip silently
      }
    }
  }
}

export function streamChat(
  token     : string,
  messages  : Array<{ role: string; content: string }>,
  model     : string,
  sessionId : string,
  signal?   : AbortSignal,
): AsyncGenerator<SSEEvent> {
  return (async function* () {
    const res = await fetch(`${BASE}/chat`, {
      method : 'POST',
      headers: {
        'Content-Type': 'application/json',
        Authorization : `Bearer ${token}`,
      },
      body  : JSON.stringify({ messages, model, session_id: sessionId }),
      signal,
    });
    yield* parseSSE(res);
  })();
}

// Reattach to a turn already running server-side for this session — no new
// message, just resume watching (replays everything emitted so far, then
// live events). Used after a session switch, logout/login, dropped
// connection, or another tab.
export function attachChat(
  token     : string,
  sessionId : string,
  signal?   : AbortSignal,
): AsyncGenerator<SSEEvent> {
  return (async function* () {
    const res = await fetch(`${BASE}/chat/stream/${sessionId}`, {
      headers: { Authorization: `Bearer ${token}` },
      signal,
    });
    yield* parseSSE(res);
  })();
}

// Actually stop the turn server-side (not just stop watching it).
export async function cancelChat(token: string, sessionId: string): Promise<void> {
  await fetch(`${BASE}/chat/${sessionId}/cancel`, {
    method : 'POST',
    headers: { Authorization: `Bearer ${token}` },
  });
}
