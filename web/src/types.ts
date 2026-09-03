// ----------------------------
// Message / chat
// ----------------------------

export type MessageRole = 'user' | 'assistant' | 'tool';

export interface ToolCallEvent {
  name : string;
  args?: Record<string, unknown>;
}

export interface ChatMessage {
  id        : string;
  role      : MessageRole;
  content   : string;
  toolCalls : ToolCallEvent[];
  streaming : boolean;   // true while tokens are still arriving
  error     : boolean;
}

// ----------------------------
// Sessions
// ----------------------------

export interface Session {
  id      : string;
  title   : string;   // custom title, or first ~60 chars of the first user message
  model   : string;   // LLM used (restored on switch); '' if never recorded
  created : number;   // history file mtime (epoch seconds)
  running : boolean;  // a turn is currently in flight server-side for this session
}

// ----------------------------
// User
// ----------------------------

export interface User {
  username       : string;
  free_questions : number;
  n_questions    : number;
  is_active      : boolean;
}

// ----------------------------
// SSE event shapes from POST /chat
// ----------------------------

export type SSEEvent =
  | { type: 'token';   text: string    }
  | { type: 'tool';    name: string; content: string; args?: Record<string, unknown> }
  | { type: 'done';    content: string }
  | { type: 'error';   message: string };

// ----------------------------
// Files
// ----------------------------

export interface FileEntry {
  path     : string;   // relative to user's files/ dir
  name     : string;
  size     : number;
  modified : number;
}

// ----------------------------
// Job notifications
// ----------------------------

export interface JobNotification {
  job_id  : string;
  type    : string;
  status  : string;
  message : string;
}

// ----------------------------
// Async jobs
// ----------------------------

export type JobStatus = 'running' | 'done' | 'failed' | 'cancelled';

export interface Job {
  id        : string;
  tool      : string;
  status    : JobStatus;
  result    : string | null;
  sessionId : string;    // which chat session started this job
  notified  : boolean;   // true once "Tell agent" has been clicked
}
