import { useState, useEffect, useRef } from 'react';
import { useAuth }                      from '../store';
import { listFiles, uploadFile, deleteFile, fileDownloadUrl } from '../api';
import type { Job, FileEntry }          from '../types';
import type { SortBy }                  from './FileTree';
import JobCard                          from './JobCard';
import FileTree                         from './FileTree';
import '../styles/jobpanel.css';

interface Props {
  jobs            : Job[];
  onTellAgent     : (job: Job) => void;
  onDismissJob    : (job: Job) => void;
  onAcknowledge   : (job: Job) => void;
  filesRefreshKey : number;
}

export default function JobPanel({ jobs, onTellAgent, onDismissJob, onAcknowledge, filesRefreshKey }: Props) {
  const { token }                       = useAuth();
  const [files,     setFiles]           = useState<FileEntry[]>([]);
  const [uploading, setUploading]       = useState(false);
  const [sortBy,    setSortBy]          = useState<SortBy>('date');
  const fileInputRef                    = useRef<HTMLInputElement>(null);
  const running = jobs.filter(j => j.status === 'running').length;

  function refreshFiles() {
    if (!token) return;
    listFiles(token).then(setFiles);
  }

  useEffect(() => { refreshFiles(); }, [token, filesRefreshKey]);

  async function handleUpload(e: React.ChangeEvent<HTMLInputElement>) {
    const file = e.target.files?.[0];
    if (!file || !token) return;
    setUploading(true);
    try {
      await uploadFile(token, file);
      refreshFiles();
    } finally {
      setUploading(false);
      if (fileInputRef.current) fileInputRef.current.value = '';
    }
  }

  async function handleDelete(path: string) {
    if (!token) return;
    await deleteFile(token, path);
    refreshFiles();
  }

  async function handleDownload(f: FileEntry) {
    if (!token) return;
    const res  = await fetch(fileDownloadUrl(f.path), {
      headers: { Authorization: `Bearer ${token}` },
    });
    const blob = await res.blob();
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement('a');
    a.href     = url;
    a.download = f.name;
    a.click();
    URL.revokeObjectURL(url);
  }

  return (
    <aside className="job-panel">
      {/* Hidden file input — outside any visible layout element to avoid rendering artifacts */}
      <input ref={fileInputRef} type="file" hidden onChange={handleUpload} />

      {/* Files section */}
      <div className="job-panel-header">
        <span className="job-panel-title">Files</span>
      </div>
      <div className="job-panel-toolbar">
        <div className="sort-group">
          <span className="sort-label">Sort by</span>
          <div style={{ display: 'flex', gap: 4 }}>
            <button className={`btn-sort ${sortBy === 'date' ? 'btn-sort-active' : ''}`}
                    onClick={() => setSortBy('date')} title="Sort by date">Date</button>
            <button className={`btn-sort ${sortBy === 'name' ? 'btn-sort-active' : ''}`}
                    onClick={() => setSortBy('name')} title="Sort by name">Name</button>
          </div>
        </div>
        <div style={{ display: 'flex', gap: 4 }}>
          <button className="btn-icon" onClick={refreshFiles} title="Refresh">⟳</button>
          <button className="btn-icon" onClick={() => fileInputRef.current?.click()}
                  title="Upload file" disabled={uploading}>
            {uploading ? '…' : '↑'}
          </button>
        </div>
      </div>

      <div className="file-list">
        {files.length === 0 ? (
          <div className="file-empty">No files yet. Generated structures and results will appear here.</div>
        ) : (
          <FileTree files={files} sortBy={sortBy} token={token ?? ''} onDownload={handleDownload} onDelete={handleDelete} />
        )}
      </div>

      {/* Jobs section */}
      <div className="job-panel-header job-panel-header-jobs">
        <span className="job-panel-title">Jobs</span>
        {running > 0 && (
          <span className="job-count">{running} running</span>
        )}
      </div>

      <div className="job-list">
        {jobs.length === 0 ? (
          <div className="job-empty">
            Async jobs (screening, AF3, ESM3) appear here when the agent starts them.
          </div>
        ) : (
          jobs.map(j => (
            <JobCard key={j.id} job={j} onTellAgent={onTellAgent} onDismiss={onDismissJob} onAcknowledge={onAcknowledge} />
          ))
        )}
      </div>

    </aside>
  );
}
