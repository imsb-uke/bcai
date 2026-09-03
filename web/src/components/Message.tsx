import ReactMarkdown        from 'react-markdown';
import remarkGfm            from 'remark-gfm';
import type { ChatMessage } from '../types';
import { useAuth }          from '../store';
import ToolCall             from './ToolCall';
import HtmlViewer, { isHtmlTool } from './HtmlViewer';
import MoleculeSpinner      from './MoleculeSpinner';

// Escape [ and ] outside code spans so SMILES strings (e.g. [C@@H], [Na+])
// are not misread as markdown link labels by remark-gfm.
function escapeBrackets(text: string): string {
  return text.split(/(`[^`]*`)/).map(part =>
    (part.startsWith('`') && part.endsWith('`'))
      ? part
      : part.replace(/\[/g, '\\[').replace(/\]/g, '\\]')
  ).join('');
}

interface Props { msg: ChatMessage }

export default function Message({ msg }: Props) {
  const { token } = useAuth();

  // Tool-role messages: show the call indicator + optional inline HTML viewer
  if (msg.role === 'tool') {
    const tc = msg.toolCalls[0];
    return (
      <div className="msg-row msg-assistant">
        <div className="msg-tools">
          <ToolCall name={tc?.name ?? 'tool'} args={tc?.args} result={msg.content} />
        </div>
        {tc && isHtmlTool(tc.name) && token && (
          <HtmlViewer toolName={tc.name} content={msg.content} token={token} />
        )}
      </div>
    );
  }

  const isUser = msg.role === 'user';

  return (
    <div className={`msg-row ${isUser ? 'msg-user' : 'msg-assistant'}`}>

      <div className={`msg-bubble${msg.error ? ' msg-error' : ''}`}>
        {msg.streaming && msg.content === '' ? (
          <MoleculeSpinner />
        ) : isUser ? (
          // User text stays verbatim — no markdown parsing
          <>
            {msg.content.split('\n').map((line, i, arr) => (
              <span key={i}>{line}{i < arr.length - 1 && <br />}</span>
            ))}
          </>
        ) : (
          <>
            <div className="md-body">
              <ReactMarkdown remarkPlugins={[remarkGfm]}>{escapeBrackets(msg.content)}</ReactMarkdown>
            </div>
            {msg.streaming && <span className="streaming-cursor" aria-hidden="true" />}
          </>
        )}
      </div>

    </div>
  );
}
