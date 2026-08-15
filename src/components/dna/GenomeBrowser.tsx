"use client";

import { useMemo } from "react";
import { Oligo } from "@/lib/dna/types";

const BASE_COLORS: Record<string, string> = {
  A: "#22c55e", // green
  C: "#3b82f6", // blue
  G: "#eab308", // yellow
  T: "#ef4444", // red
};

const BASE_TEXT_COLORS: Record<string, string> = {
  A: "#15803d",
  C: "#1e40af",
  G: "#a16207",
  T: "#b91c1c",
};

interface GenomeBrowserProps {
  oligos: Oligo[];
  forwardPrimer: string;
  reversePrimer: string;
  primerLength: number;
  selectedOligo: number | null;
  onSelectOligo: (idx: number) => void;
  /** Mutated reads to overlay (optional). If provided, shows mutations as highlights. */
  mutatedReads?: { index: number; sequence: string; original: string }[];
  maxOligosToShow?: number;
}

/**
 * Visual genome browser: displays oligos as horizontal tracks with color-coded
 * bases. Click an oligo to see its detailed sequence.
 */
export function GenomeBrowser({
  oligos,
  forwardPrimer,
  reversePrimer,
  primerLength,
  selectedOligo,
  onSelectOligo,
  mutatedReads,
  maxOligosToShow = 50,
}: GenomeBrowserProps) {
  const visibleOligos = oligos.slice(0, maxOligosToShow);

  // Per-position GC content (across all visible oligos)
  const gcStats = useMemo(() => {
    if (visibleOligos.length === 0) return [];
    const oligoLen = visibleOligos[0].sequence.length;
    const stats: { pos: number; gc: number; a: number; c: number; g: number; t: number }[] = [];
    for (let pos = 0; pos < oligoLen; pos++) {
      let gc = 0;
      let a = 0, c = 0, g = 0, t = 0;
      let total = 0;
      for (const oligo of visibleOligos) {
        const base = oligo.sequence[pos];
        if (base === "G" || base === "C") gc++;
        if (base === "A") a++;
        else if (base === "C") c++;
        else if (base === "G") g++;
        else if (base === "T") t++;
        total++;
      }
      stats.push({ pos, gc: gc / total, a, c, g, t });
    }
    return stats;
  }, [visibleOligos]);

  if (oligos.length === 0) {
    return (
      <div className="flex h-64 items-center justify-center text-muted-foreground border rounded-lg">
        No oligos to display. Encode a file to see the genome browser.
      </div>
    );
  }

  const oligoLen = oligos[0].sequence.length;

  return (
    <div className="space-y-4">
      {/* GC content track */}
      <div>
        <div className="text-xs font-medium mb-1 text-muted-foreground">
          GC Content per Position (across {visibleOligos.length} oligos)
        </div>
        <GcChart stats={gcStats} primerLength={primerLength} />
      </div>

      {/* Oligo tracks */}
      <div>
        <div className="text-xs font-medium mb-2 text-muted-foreground">
          Oligo Tracks ({oligos.length} total{oligos.length > maxOligosToShow ? `, showing first ${maxOligosToShow}` : ""})
        </div>
        <div className="space-y-1 max-h-[500px] overflow-y-auto pr-2 border rounded-lg p-2 bg-muted/30">
          {visibleOligos.map((oligo) => (
            <OligoTrack
              key={oligo.index}
              oligo={oligo}
              forwardPrimer={forwardPrimer}
              reversePrimer={reversePrimer}
              primerLength={primerLength}
              selected={selectedOligo === oligo.index}
              onClick={() => onSelectOligo(oligo.index)}
            />
          ))}
        </div>
      </div>

      {/* Selected oligo detail */}
      {selectedOligo !== null && oligos[selectedOligo] && (
        <OligoDetail
          oligo={oligos[selectedOligo]}
          forwardPrimer={forwardPrimer}
          reversePrimer={reversePrimer}
          primerLength={primerLength}
        />
      )}
    </div>
  );
}

function OligoTrack({
  oligo,
  forwardPrimer,
  reversePrimer,
  primerLength,
  selected,
  onClick,
}: {
  oligo: Oligo;
  forwardPrimer: string;
  reversePrimer: string;
  primerLength: number;
  selected: boolean;
  onClick: () => void;
}) {
  const seq = oligo.sequence;
  return (
    <button
      onClick={onClick}
      className={`w-full text-left rounded p-1 transition-colors ${
        selected ? "bg-primary/15 ring-1 ring-primary" : "hover:bg-muted"
      }`}
    >
      <div className="flex items-center gap-2">
        <span className="text-[10px] font-mono text-muted-foreground w-12 shrink-0">
          #{oligo.index.toString().padStart(4, "0")}
        </span>
        <div className="flex-1 overflow-hidden">
          <div className="flex h-4 gap-[1px]">
            {seq.split("").map((base, i) => {
              const isPrimer = i < primerLength || i >= seq.length - primerLength;
              return (
                <div
                  key={i}
                  className="flex-1 min-w-[2px] rounded-sm"
                  style={{
                    backgroundColor: BASE_COLORS[base] ?? "#888",
                    opacity: isPrimer ? 0.4 : 1,
                  }}
                  title={`pos ${i}: ${base}${isPrimer ? " (primer)" : ""}`}
                />
              );
            })}
          </div>
        </div>
        <span className="text-[10px] font-mono text-muted-foreground w-16 shrink-0 text-right">
          GC {(oligo.gc * 100).toFixed(0)}%
        </span>
      </div>
    </button>
  );
}

function OligoDetail({
  oligo,
  forwardPrimer,
  reversePrimer,
  primerLength,
}: {
  oligo: Oligo;
  forwardPrimer: string;
  reversePrimer: string;
  primerLength: number;
}) {
  const seq = oligo.sequence;
  const payload = seq.slice(primerLength, seq.length - primerLength);

  // Group into 10-base chunks for readability
  const chunks: string[] = [];
  for (let i = 0; i < payload.length; i += 10) {
    chunks.push(payload.slice(i, i + 10));
  }

  return (
    <div className="border rounded-lg p-4 space-y-3 bg-card">
      <div className="flex items-center justify-between">
        <h4 className="font-semibold">Oligo #{oligo.index} Detail</h4>
        <div className="flex gap-3 text-xs text-muted-foreground">
          <span>Length: {oligo.length} nt</span>
          <span>GC: {(oligo.gc * 100).toFixed(1)}%</span>
          <span>Max homopolymer: {oligo.maxHomopolymer}</span>
          <span>Seed: {oligo.seed}</span>
        </div>
      </div>

      <div>
        <div className="text-xs text-muted-foreground mb-1">Forward primer (5&apos; end, {primerLength} nt):</div>
        <div className="font-mono text-xs break-all bg-muted/50 p-2 rounded">
          {colorizeBases(forwardPrimer)}
        </div>
      </div>

      <div>
        <div className="text-xs text-muted-foreground mb-1">Payload ({payload.length} nt):</div>
        <div className="font-mono text-xs break-all bg-muted/50 p-2 rounded leading-relaxed">
          {chunks.map((chunk, i) => (
            <span key={i}>
              {colorizeBases(chunk)}
              {(i + 1) % 4 === 0 && i < chunks.length - 1 ? "\n" : " "}
            </span>
          ))}
        </div>
      </div>

      <div>
        <div className="text-xs text-muted-foreground mb-1">Reverse primer (3&apos; end, {primerLength} nt):</div>
        <div className="font-mono text-xs break-all bg-muted/50 p-2 rounded">
          {colorizeBases(reversePrimer)}
        </div>
      </div>

      <div className="flex gap-2 text-xs">
        <BaseLegend />
      </div>
    </div>
  );
}

function colorizeBases(seq: string): React.ReactNode {
  return seq.split("").map((base, i) => (
    <span key={i} style={{ color: BASE_TEXT_COLORS[base] ?? "#666" }}>
      {base}
    </span>
  ));
}

function BaseLegend() {
  return (
    <div className="flex gap-3">
      {Object.entries(BASE_COLORS).map(([base, color]) => (
        <div key={base} className="flex items-center gap-1">
          <div className="w-3 h-3 rounded" style={{ backgroundColor: color }} />
          <span className="font-mono">{base}</span>
        </div>
      ))}
    </div>
  );
}

function GcChart({
  stats,
  primerLength,
}: {
  stats: { pos: number; gc: number }[];
  primerLength: number;
}) {
  if (stats.length === 0) return null;
  const maxLen = stats.length;
  const width = 800;
  const height = 60;
  const barWidth = width / maxLen;

  return (
    <div className="border rounded-lg p-2 bg-card overflow-x-auto">
      <svg width={width} height={height} className="w-full">
        {/* 50% line */}
        <line
          x1={0}
          y1={height / 2}
          x2={width}
          y2={height / 2}
          stroke="#94a3b8"
          strokeDasharray="4 2"
          strokeWidth={0.5}
        />
        {/* 40% and 60% bounds */}
        <line x1={0} y1={height * 0.4} x2={width} y2={height * 0.4} stroke="#cbd5e1" strokeWidth={0.5} />
        <line x1={0} y1={height * 0.6} x2={width} y2={height * 0.6} stroke="#cbd5e1" strokeWidth={0.5} />

        {/* GC bars */}
        {stats.map((s, i) => {
          const barH = s.gc * height;
          const isPrimer = i < primerLength || i >= maxLen - primerLength;
          return (
            <rect
              key={i}
              x={i * barWidth}
              y={height - barH}
              width={Math.max(barWidth - 0.5, 1)}
              height={barH}
              fill={isPrimer ? "#94a3b8" : "#0ea5e9"}
              opacity={isPrimer ? 0.4 : 0.9}
            />
          );
        })}

        {/* Labels */}
        <text x={4} y={10} fontSize={9} fill="#64748b">100%</text>
        <text x={4} y={height / 2 + 3} fontSize={9} fill="#64748b">50%</text>
        <text x={4} y={height - 2} fontSize={9} fill="#64748b">0%</text>
      </svg>
    </div>
  );
}
