interface Props {
  name    : string;
  args?   : Record<string, unknown>;
  result? : string;
}

// Tool rows are only created after a tool finishes, so this always shows the
// completed (✓) state. Clicking expands to show the full tool name, actual
// args, and the raw result returned by the tool.
export default function ToolCall({ name, args, result }: Props) {
  // Strip the server-key prefix the engine adds (e.g. "BioChemAIgent__run_screening" → "run_screening")
  const display = name.includes('__') ? name.split('__').pop()! : name;

  return (
    <details className="tool-call tool-done">
      <summary className="tool-summary">
        <span className="tool-icon">✓</span>
        <span className="tool-name">{display}</span>
      </summary>
      <div className="tool-detail">
        <code className="tool-full-name">{name}</code>
        {args && (
          <pre className="tool-args">{JSON.stringify(args, null, 2)}</pre>
        )}
        {result && (
          <pre className="tool-result">{result}</pre>
        )}
      </div>
    </details>
  );
}
