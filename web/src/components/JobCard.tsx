import type { Job } from '../types';

interface Props {
  job           : Job;
  onTellAgent   : (job: Job) => void;
  onDismiss     : (job: Job) => void;
  onAcknowledge : (job: Job) => void;
}

const STATUS_LABEL: Record<string, string> = {
  running  : 'Running',
  done     : 'Done',
  failed   : 'Failed',
  cancelled: 'Cancelled',
};

export default function JobCard({ job, onTellAgent, onDismiss, onAcknowledge }: Props) {
  return (
    <div className={`job-card job-${job.status}`}>

      <div className="job-header">
        <span className="job-tool">{job.tool}</span>
        <span className={`job-badge badge-${job.status}`}>
          {STATUS_LABEL[job.status] ?? job.status}
        </span>
        {job.status !== 'running' && (
          <button className="job-dismiss" onClick={() => onDismiss(job)} title="Dismiss">✕</button>
        )}
      </div>

      {/* Truncated ID — full ID visible in title on hover */}
      <div className="job-id" title={job.id}>
        {job.id.length > 20 ? `${job.id.slice(0, 20)}…` : job.id}
      </div>

      {job.result && (
        <div className="job-result">{job.result}</div>
      )}

      {/* "Tell agent" posts the result back into the active chat session */}
      {job.status === 'done' && !job.notified && (
        <div className="job-actions">
          <button className="btn-tell-agent" onClick={() => onTellAgent(job)}>
            Tell agent ↑
          </button>
          <button className="btn-acknowledge" onClick={() => onAcknowledge(job)} title="Agent already knows">
            Got it
          </button>
        </div>
      )}

    </div>
  );
}
