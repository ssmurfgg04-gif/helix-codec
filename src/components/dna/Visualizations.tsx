"use client";

import { useMemo } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { SequencingRead } from "@/lib/dna/simulate";

/**
 * Q-Score Visualization Component
 *
 * Displays per-base Phred quality scores as a color-coded bar chart.
 * Green = Q30+ (high confidence), Yellow = Q10-30, Red = Q<10 (erasure candidate).
 */
export function QScoreViewer({ read }: { read: SequencingRead | null }) {
  if (!read || !read.quality) {
    return (
      <div className="text-xs text-muted-foreground p-4 border rounded-lg">
        No Q-scores available. Run a simulation to generate per-base quality scores.
      </div>
    );
  }

  const q = read.quality;
  const maxLen = Math.min(q.length, 200); // show first 200 bases

  const colorForQ = (quality: number): string => {
    if (quality >= 30) return "#22c55e"; // green
    if (quality >= 20) return "#84cc16"; // lime
    if (quality >= 10) return "#eab308"; // yellow
    if (quality >= 5) return "#f97316"; // orange
    return "#ef4444"; // red
  };

  const lowQ = Array.from(q).filter((v) => v < 10).length;
  const meanQ = Array.from(q).reduce((s, v) => s + v, 0) / q.length;

  return (
    <div className="space-y-3">
      <div className="grid grid-cols-4 gap-3 text-xs">
        <Stat label="Length" value={`${q.length} bases`} />
        <Stat label="Mean Q" value={meanQ.toFixed(1)} />
        <Stat label="Low Q (<10)" value={`${lowQ} (${((lowQ / q.length) * 100).toFixed(1)}%)`} />
        <Stat label="Accuracy" value={`${(100 * (1 - Math.pow(10, -meanQ / 10))).toFixed(2)}%`} />
      </div>

      {/* Q-score bar chart */}
      <div className="border rounded-lg p-3 bg-card overflow-x-auto">
        <div className="flex items-end h-32 gap-[1px] min-w-[600px]">
          {Array.from(q.slice(0, maxLen)).map((quality, i) => (
            <div
              key={i}
              className="flex-1 min-w-[2px] rounded-t-sm transition-all hover:opacity-80"
              style={{
                height: `${(quality / 40) * 100}%`,
                backgroundColor: colorForQ(quality),
              }}
              title={`pos ${i}: Q${quality} (${read.sequence[i]})`}
            />
          ))}
        </div>
        <div className="flex justify-between text-[10px] text-muted-foreground mt-1">
          <span>pos 0</span>
          <span>pos {maxLen - 1}</span>
        </div>
      </div>

      {/* Sequence with color */}
      <div className="border rounded-lg p-3 bg-card overflow-x-auto">
        <div className="font-mono text-[10px] break-all leading-relaxed">
          {read.sequence.slice(0, maxLen).split("").map((base, i) => (
            <span
              key={i}
              style={{ color: colorForQ(q[i]) }}
              title={`pos ${i}: ${base} Q${q[i]}`}
            >
              {base}
            </span>
          ))}
        </div>
      </div>

      {/* Legend */}
      <div className="flex gap-3 text-xs">
        <LegendItem color="#22c55e" label="Q30+ (99.9%)" />
        <LegendItem color="#84cc16" label="Q20-30 (99%)" />
        <LegendItem color="#eab308" label="Q10-20 (90%)" />
        <LegendItem color="#f97316" label="Q5-10 (68%)" />
        <LegendItem color="#ef4444" label="Q<5 (32%)" />
      </div>
    </div>
  );
}

/**
 * Merkle Tree Graph Visualization
 *
 * Renders the SHA-256 Merkle tree as an SVG graph with nodes and edges.
 */
export function MerkleTreeGraph({
  leaves,
  nodes,
  root,
  depth,
}: {
  leaves: string[];
  nodes: string[][];
  root: string;
  depth: number;
}) {
  if (!root) {
    return (
      <div className="text-xs text-muted-foreground p-4">
        No Merkle tree available.
      </div>
    );
  }

  const width = 800;
  const nodeRadius = 12;
  const levelHeight = 50;
  const totalHeight = (depth + 1) * levelHeight + 20;

  return (
    <div className="border rounded-lg p-3 bg-card overflow-x-auto">
      <svg width={width} height={totalHeight} className="w-full">
        {/* Render tree levels from bottom (leaves) to top (root) */}
        {nodes.map((level, levelIdx) => {
          const y = totalHeight - 20 - levelIdx * levelHeight;
          const isLeafLevel = levelIdx === 0;
          const isRootLevel = levelIdx === nodes.length - 1;

          return level.slice(0, 16).map((hash, i) => {
            const x = ((i + 0.5) / Math.min(level.length, 16)) * width;
            const shortHash = hash.slice(0, 6);

            // Draw edge to parent
            if (!isRootLevel && nodes[levelIdx + 1]) {
              const parentIdx = Math.floor(i / 2);
              const parentLevel = nodes[levelIdx + 1];
              if (parentIdx < parentLevel.length) {
                const parentX = ((parentIdx + 0.5) / Math.min(parentLevel.length, 16)) * width;
                const parentY = totalHeight - 20 - (levelIdx + 1) * levelHeight;
                return (
                  <g key={`${levelIdx}-${i}`}>
                    <line x1={x} y1={y} x2={parentX} y2={parentY} stroke="#cbd5e1" strokeWidth={0.5} />
                    <circle cx={x} cy={y} r={nodeRadius} fill={isRootLevel ? "#0ea5e9" : isLeafLevel ? "#22c55e" : "#94a3b8"} />
                    <text x={x} y={y + 3} fontSize={7} fill="white" textAnchor="middle" fontWeight="bold">
                      {shortHash.slice(0, 4)}
                    </text>
                  </g>
                );
              }
            }

            return (
              <g key={`${levelIdx}-${i}`}>
                <circle cx={x} cy={y} r={nodeRadius} fill={isRootLevel ? "#0ea5e9" : isLeafLevel ? "#22c55e" : "#94a3b8"} />
                <text x={x} y={y + 3} fontSize={7} fill="white" textAnchor="middle" fontWeight="bold">
                  {shortHash.slice(0, 4)}
                </text>
              </g>
            );
          });
        })}

        {/* Root hash display */}
        <text x={width / 2} y={12} fontSize={10} fill="#0ea5e9" textAnchor="middle" fontWeight="bold">
          Root: {root.slice(0, 24)}...
        </text>
      </svg>

      <div className="flex gap-3 text-xs mt-2">
        <LegendItem color="#22c55e" label="Leaves (chunks)" />
        <LegendItem color="#94a3b8" label="Internal nodes" />
        <LegendItem color="#0ea5e9" label="Merkle root" />
      </div>
    </div>
  );
}

/**
 * Mutation Hotspot Heatmap
 *
 * Shows per-position error rate across all oligos as a heatmap.
 */
export function MutationHeatmap({
  perOligoStats,
  oligoCount,
}: {
  perOligoStats: {
    index: number;
    avgSubstitutions: number;
    avgInsertions: number;
    avgDeletions: number;
    readCount: number;
  }[];
  oligoCount: number;
}) {
  const data = useMemo(() => {
    return Array.from({ length: oligoCount }, (_, i) => {
      const stats = perOligoStats.find((s) => s.index === i);
      if (!stats || stats.readCount === 0) return { index: i, errorRate: 0, dropped: true };
      const totalErrors = stats.avgSubstitutions + stats.avgInsertions + stats.avgDeletions;
      const errorRate = totalErrors / 200; // per-base rate (200nt oligo)
      return { index: i, errorRate, dropped: false };
    });
  }, [perOligoStats, oligoCount]);

  const maxRate = Math.max(...data.map((d) => d.errorRate), 0.001);

  const colorForRate = (rate: number, dropped: boolean): string => {
    if (dropped) return "#1e293b"; // dark = dropped
    const normalized = rate / maxRate;
    if (normalized < 0.2) return "#22c55e"; // green
    if (normalized < 0.4) return "#84cc16";
    if (normalized < 0.6) return "#eab308";
    if (normalized < 0.8) return "#f97316";
    return "#ef4444"; // red
  };

  return (
    <div className="space-y-3">
      <div className="border rounded-lg p-3 bg-card">
        <div className="text-xs text-muted-foreground mb-2">Per-oligo error rate (avg errors / 200nt)</div>
        <div className="flex flex-wrap gap-[2px]">
          {data.map((d) => (
            <div
              key={d.index}
              className="w-6 h-6 rounded-sm transition-all hover:scale-110 cursor-default"
              style={{ backgroundColor: colorForRate(d.errorRate, d.dropped) }}
              title={`Oligo ${d.index}: ${d.dropped ? "DROPPED" : `${(d.errorRate * 100).toFixed(1)}% error rate`}`}
            />
          ))}
        </div>
      </div>

      <div className="flex gap-3 text-xs">
        <LegendItem color="#22c55e" label="<2%" />
        <LegendItem color="#eab308" label="2-4%" />
        <LegendItem color="#f97316" label="4-6%" />
        <LegendItem color="#ef4444" label=">6%" />
        <LegendItem color="#1e293b" label="Dropped" />
      </div>
    </div>
  );
}

function Stat({ label, value }: { label: string; value: string }) {
  return (
    <div className="space-y-0.5">
      <div className="text-muted-foreground">{label}</div>
      <div className="font-mono font-bold">{value}</div>
    </div>
  );
}

function LegendItem({ color, label }: { color: string; label: string }) {
  return (
    <div className="flex items-center gap-1">
      <div className="w-3 h-3 rounded" style={{ backgroundColor: color }} />
      <span>{label}</span>
    </div>
  );
}
