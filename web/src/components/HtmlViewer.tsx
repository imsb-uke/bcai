import { useState, useEffect } from 'react';

// Tools whose result should be rendered as an inline HTML viewer
const HTML_TOOLS = ['render_structures', 'interaction_plot'];

interface Props {
  toolName : string;
  content  : string;   // raw JSON string from the tool result
  token    : string;
}

export function isHtmlTool(toolName: string): boolean {
  const bare = toolName.includes('__') ? toolName.split('__').pop()! : toolName;
  return HTML_TOOLS.some(t => bare.includes(t));
}

export default function HtmlViewer({ toolName, content, token }: Props) {
  const [html,   setHtml]   = useState<string | null>(null);
  const [zoom,   setZoom]   = useState(70);
  const [height, setHeight] = useState(400);
  const [error,  setError]  = useState('');

  useEffect(() => {
    // Parse html_file path from tool result JSON
    let htmlFile = '';
    try {
      const parsed = JSON.parse(content);
      htmlFile = parsed.html_file ?? parsed.output ?? '';
    } catch {
      htmlFile = '';
    }
    if (!htmlFile) { setError('No HTML file in tool result'); return; }

    // Strip everything up to and including /files/ to get the relative path
    const marker  = '/files/';
    const idx     = htmlFile.indexOf(marker);
    const relPath = idx >= 0 ? htmlFile.slice(idx + marker.length) : htmlFile;

    fetch(`/api/files/${relPath}`, {
      headers: { Authorization: `Bearer ${token}` },
    })
      .then(r => { if (!r.ok) throw new Error(`HTTP ${r.status}`); return r.text(); })
      .then(setHtml)
      .catch(e => setError(String(e)));
  }, [content, token]);

  if (error) return <div className="html-viewer-error">{error}</div>;
  if (!html)  return <div className="html-viewer-loading">Loading viewer…</div>;

  function openInNewTab() {
    const blob = new Blob([html!], { type: 'text/html' });
    const url  = URL.createObjectURL(blob);
    window.open(url, '_blank');
    // blob URLs are kept alive as long as the tab is open; no need to revoke
  }

  return (
    <div className="html-viewer">
      <div className="html-viewer-controls">
        <label>
          Zoom&nbsp;
          <input type="range" min={20} max={100} value={zoom}
                 onChange={e => setZoom(+e.target.value)} />
          &nbsp;{zoom}%
        </label>
        <label>
          Height&nbsp;
          <input type="range" min={200} max={1000} step={50} value={height}
                 onChange={e => setHeight(+e.target.value)} />
          &nbsp;{height}px
        </label>
        <button className="html-viewer-newtab" onClick={openInNewTab} title="Open in new tab">
          ↗
        </button>
      </div>
      <iframe
        className    = "html-viewer-frame"
        style        = {{ height }}
        srcDoc       = {`<div style="transform:scale(${zoom/100});transform-origin:top left;width:${Math.round(10000/zoom)}%">${html}</div>`}
        sandbox      = "allow-scripts"
        title        = {toolName}
      />
    </div>
  );
}
