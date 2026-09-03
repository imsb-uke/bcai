import '../styles/chat.css';

export default function MoleculeSpinner() {
  return (
    <div className="mol-spinner-wrap" aria-label="Agent is processing…">
      <svg
        className = "mol-spinner"
        viewBox   = "0 0 220 150"
        width     = "88"
        height    = "60"
        fill      = "none"
        xmlns     = "http://www.w3.org/2000/svg"
      >
        {/* Hexagon 1 (left ring) */}
        <polygon
          points       = "78,39 110,57 110,93 78,111 46,93 46,57"
          stroke       = "#378ADD"
          strokeWidth  = "2.5"
          strokeLinejoin="round"
        />
        {/* Hexagon 2 (right ring) */}
        <polygon
          points       = "142,39 174,57 174,93 142,111 110,93 110,57"
          stroke       = "#378ADD"
          strokeWidth  = "2.5"
          strokeLinejoin="round"
        />

        {/* Aromatic circles */}
        <circle className="ring"   cx="78"  cy="75" r="20" stroke="#85B7EB" strokeWidth="1.5" />
        <circle className="ring ring-r" cx="142" cy="75" r="20" stroke="#85B7EB" strokeWidth="1.5" />

        {/* Substituents — left side */}
        <line x1="46" y1="57"  x2="26" y2="46"  stroke="#378ADD" strokeWidth="2" />
        <circle className="dot d1" cx="21" cy="43"  r="5" fill="#378ADD" />
        <line x1="46" y1="93"  x2="26" y2="104" stroke="#378ADD" strokeWidth="2" />
        <circle className="dot d2" cx="21" cy="107" r="5" fill="#378ADD" />

        {/* Substituents — right side */}
        <line x1="174" y1="57"  x2="194" y2="46"  stroke="#378ADD" strokeWidth="2" />
        <circle className="dot d3" cx="199" cy="43"  r="5" fill="#378ADD" />
        <line x1="174" y1="93"  x2="194" y2="104" stroke="#378ADD" strokeWidth="2" />
        <circle className="dot d4" cx="199" cy="107" r="5" fill="#378ADD" />

        {/* Substituents — top and bottom */}
        <line x1="78"  y1="39"  x2="78"  y2="19"  stroke="#378ADD" strokeWidth="2" />
        <circle className="dot d5" cx="78"  cy="14"  r="5" fill="#85B7EB" />
        <line x1="142" y1="111" x2="142" y2="131" stroke="#378ADD" strokeWidth="2" />
        <circle className="dot d6" cx="142" cy="136" r="5" fill="#85B7EB" />
      </svg>
      <span className="mol-spinner-label">Processing…</span>
    </div>
  );
}
