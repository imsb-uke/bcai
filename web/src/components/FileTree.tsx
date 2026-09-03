import { useState } from 'react';
import type { FileEntry } from '../types';

// ----------------------------
// Tree model — build a nested structure from flat relative paths
// ----------------------------

interface TreeNode {
  name     : string;
  path     : string;            // full relative path (files only)
  isDir    : boolean;
  size     : number;
  modified : number;
  children : Record<string, TreeNode>;
}

function buildTree(files: FileEntry[]): TreeNode {
  const root: TreeNode = { name: '', path: '', isDir: true, size: 0, modified: 0, children: {} };
  for (const f of files) {
    const parts = f.path.split('/');
    let node = root;
    parts.forEach((part, i) => {
      const isLeaf = i === parts.length - 1;
      if (!node.children[part]) {
        node.children[part] = {
          name    : part,
          path    : parts.slice(0, i + 1).join('/'),
          isDir   : !isLeaf,
          size    : isLeaf ? f.size : 0,
          modified: isLeaf ? f.modified : 0,
          children: {},
        };
      }
      // bubble up the most recent modified time to parent dirs
      if (f.modified > node.children[part].modified) {
        node.children[part].modified = f.modified;
      }
      node = node.children[part];
    });
  }
  return root;
}

function formatSize(bytes: number): string {
  if (bytes < 1024)        return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

// File extensions that can be displayed as plain text
const TEXT_EXTS = new Set([
  'txt', 'log', 'csv', 'tsv', 'json', 'yaml', 'yml',
  'py', 'sh', 'md', 'fasta', 'fa', 'pdb', 'sdf', 'mol2',
  'cif', 'mmcif', 'xml', 'toml', 'ini', 'cfg',
]);

function fileType(name: string): 'html' | 'text' | 'binary' {
  const ext = name.split('.').pop()?.toLowerCase() ?? '';
  if (ext === 'html' || ext === 'htm') return 'html';
  if (TEXT_EXTS.has(ext))              return 'text';
  return 'binary';
}

// ----------------------------
// File viewer modal
// ----------------------------

interface ViewerProps {
  path    : string;
  name    : string;
  content : string | null;   // null = loading, '' = error/unsupported
  isHtml  : boolean;
  onClose : () => void;
}

function FileViewer({ path, name, content, isHtml, onClose }: ViewerProps) {
  function openNewTab() {
    if (!content) return;
    const mime = isHtml ? 'text/html' : 'text/plain';
    const blob = new Blob([content], { type: mime });
    const url  = URL.createObjectURL(blob);
    window.open(url, '_blank');
  }

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal-card file-viewer-card" onClick={e => e.stopPropagation()}>
        <div className="modal-header">
          <span className="modal-title" title={path}>{name}</span>
          <div style={{ display: 'flex', gap: 6 }}>
            {content && (
              <button className="html-viewer-newtab" onClick={openNewTab} title="Open in new tab">↗</button>
            )}
            <button className="modal-close" onClick={onClose}>✕</button>
          </div>
        </div>

        {content === null && <pre className="file-viewer-content">Loading…</pre>}
        {content === ''   && <pre className="file-viewer-content">Cannot display this file.</pre>}
        {content && isHtml && (
          <iframe
            className = "file-viewer-frame"
            srcDoc    = {content}
            sandbox   = "allow-scripts"
            title     = {name}
          />
        )}
        {content && !isHtml && (
          <pre className="file-viewer-content">{content}</pre>
        )}
      </div>
    </div>
  );
}

// ----------------------------
// Components
// ----------------------------

export type SortBy = 'date' | 'name';

function sortNodes(nodes: TreeNode[], sortBy: SortBy): TreeNode[] {
  return nodes.sort((a, b) =>
    Number(b.isDir) - Number(a.isDir) ||
    (sortBy === 'name'
      ? a.name.localeCompare(b.name)
      : b.modified - a.modified),
  );
}

interface Props {
  files      : FileEntry[];
  sortBy     : SortBy;
  token      : string;
  onDownload : (f: FileEntry) => void;
  onDelete   : (path: string) => void;
}

export default function FileTree({ files, sortBy, token, onDownload, onDelete }: Props) {
  const tree = buildTree(files);
  const top  = sortNodes(Object.values(tree.children), sortBy);

  const [viewer, setViewer] = useState<{ path: string; name: string; content: string | null; isHtml: boolean } | null>(null);

  async function openFile(path: string, name: string) {
    const type = fileType(name);
    if (type === 'binary') {
      setViewer({ path, name, content: '', isHtml: false });
      return;
    }
    const isHtml = type === 'html';
    setViewer({ path, name, content: null, isHtml });
    try {
      const res  = await fetch(`/api/files/${path}`, {
        headers: { Authorization: `Bearer ${token}` },
      });
      const text = res.ok ? await res.text() : '';
      setViewer({ path, name, content: text, isHtml });
    } catch {
      setViewer({ path, name, content: '', isHtml });
    }
  }

  return (
    <>
      <div className="file-tree">
        {top.map(node => (
          <TreeRow key={node.path} node={node} depth={0} sortBy={sortBy}
                   onDownload={onDownload} onDelete={onDelete} onOpen={openFile} />
        ))}
      </div>

      {viewer && (
        <FileViewer
          path={viewer.path} name={viewer.name} content={viewer.content}
          isHtml={viewer.isHtml} onClose={() => setViewer(null)}
        />
      )}
    </>
  );
}

interface RowProps {
  node       : TreeNode;
  depth      : number;
  sortBy     : SortBy;
  onDownload : (f: FileEntry) => void;
  onDelete   : (path: string) => void;
  onOpen     : (path: string, name: string) => void;
}

function TreeRow({ node, depth, sortBy, onDownload, onDelete, onOpen }: RowProps) {
  const [open, setOpen] = useState(true);

  if (node.isDir) {
    const kids = sortNodes(Object.values(node.children), sortBy);
    return (
      <>
        <div className="file-row file-dir" onClick={() => setOpen(o => !o)}>
          <span className="file-caret">{open ? '▾' : '▸'}</span>
          <span className="file-name">{node.name}</span>
        </div>
        {open && (
          <div className="file-children">
            {kids.map(k => (
              <TreeRow key={k.path} node={k} depth={depth + 1} sortBy={sortBy}
                       onDownload={onDownload} onDelete={onDelete} onOpen={onOpen} />
            ))}
          </div>
        )}
      </>
    );
  }

  return (
    <div className="file-row" title={node.path}>
      <span className="file-name file-name-link" onClick={() => onOpen(node.path, node.name)}>
        {node.name}
      </span>
      <span className="file-size">{formatSize(node.size)}</span>
      <div className="file-actions">
        <button className="file-btn"
                onClick={() => onDownload({ path: node.path, name: node.name, size: node.size, modified: 0 })}
                title="Download">↓</button>
        <button className="file-btn file-btn-del"
                onClick={() => onDelete(node.path)} title="Delete">✕</button>
      </div>
    </div>
  );
}
