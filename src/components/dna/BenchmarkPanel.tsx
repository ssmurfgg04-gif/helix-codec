"use client";

import { useState, useCallback } from "react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";
import { Progress } from "@/components/ui/progress";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { useToast } from "@/hooks/use-toast";
import { Loader2, Gauge, TrendingDown, TrendingUp } from "lucide-react";
import { CodecConfig, DEFAULT_CONFIG } from "@/lib/dna/types";

interface BenchmarkPoint {
  label: string;
  errorRate: number;
  coverage: number;
  oligoCount: number;
  encoded: boolean;
  decoded: boolean;
  hashMatch: boolean;
  encodeMs: number;
  decodeMs: number;
  netDensityBitsPerNt: number;
  oligosRecovered: number;
  oligosErased: number;
  readsTotal: number;
}

interface BenchmarkSummary {
  maxErrorRateRecovered: number;
  minCoverageRecovered: number;
  avgDensityBitsPerNt: number;
  avgEncodeMs: number;
  avgDecodeMs: number;
}

interface BenchmarkResultData {
  points: BenchmarkPoint[];
  summary: BenchmarkSummary;
}

export function BenchmarkPanel() {
  const [config] = useState<CodecConfig>(DEFAULT_CONFIG);
  const [payloadSize, setPayloadSize] = useState(2048);
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState(0);
  const [errorResult, setErrorResult] = useState<BenchmarkResultData | null>(null);
  const [coverageResult, setCoverageResult] = useState<BenchmarkResultData | null>(null);
  const { toast } = useToast();

  const runBenchmark = useCallback(
    async (type: "error" | "coverage") => {
      setRunning(true);
      setProgress(10);
      try {
        const res = await fetch("/api/dna", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            op: type === "error" ? "benchmark-error" : "benchmark-coverage",
            payloadSize,
            config,
          }),
        });
        setProgress(60);
        if (!res.ok) {
          const err = await res.json();
          throw new Error(err.error ?? "Benchmark failed");
        }
        const json = await res.json();
        setProgress(100);
        if (type === "error") {
          setErrorResult(json.benchmark as BenchmarkResultData);
        } else {
          setCoverageResult(json.benchmark as BenchmarkResultData);
        }
        toast({
          title: "Benchmark complete",
          description: `Max error rate recovered: ${((type === "error" ? json.benchmark.summary.maxErrorRateRecovered : 0.003) * 100).toFixed(2)}%`,
        });
      } catch (e) {
        toast({ title: "Benchmark failed", description: (e as Error).message, variant: "destructive" });
      } finally {
        setRunning(false);
        setTimeout(() => setProgress(0), 1000);
      }
    },
    [payloadSize, config, toast],
  );

  return (
    <div className="space-y-4">
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Gauge className="h-4 w-4" /> Benchmark Suite
          </CardTitle>
          <CardDescription>
            Measure recovery success rate vs. error rate and coverage depth. Uses a synthetic payload
            with mixed entropy.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2">
            <div className="flex justify-between text-sm">
              <Label>Payload size: {payloadSize.toLocaleString()} bytes</Label>
              <span className="text-muted-foreground">synthetic test data</span>
            </div>
            <Slider
              min={512}
              max={10240}
              step={512}
              value={[payloadSize]}
              onValueChange={([v]) => setPayloadSize(v)}
            />
          </div>

          <Tabs defaultValue="error">
            <TabsList className="grid w-full grid-cols-2">
              <TabsTrigger value="error">
                <TrendingUp className="h-3 w-3 mr-1" /> Error Rate Sweep
              </TabsTrigger>
              <TabsTrigger value="coverage">
                <TrendingDown className="h-3 w-3 mr-1" /> Coverage Sweep
              </TabsTrigger>
            </TabsList>
            <TabsContent value="error" className="space-y-3 pt-3">
              <p className="text-xs text-muted-foreground">
                Encodes the payload once, then decodes at increasing error rates (0% to 5%) with
                20x coverage. Reports the maximum error rate at which recovery still succeeds.
              </p>
              <Button onClick={() => runBenchmark("error")} disabled={running} className="w-full">
                {running ? (
                  <>
                    <Loader2 className="h-4 w-4 mr-2 animate-spin" /> Running...
                  </>
                ) : (
                  <>
                    <Gauge className="h-4 w-4 mr-2" /> Run Error Sweep
                  </>
                )}
              </Button>
              {errorResult && <BenchmarkChart result={errorResult} xKey="errorRate" xLabel="Error Rate" />}
            </TabsContent>
            <TabsContent value="coverage" className="space-y-3 pt-3">
              <p className="text-xs text-muted-foreground">
                Encodes the payload once, then decodes at increasing coverage depths (5x to 30x) with
                Illumina error profile. Reports the minimum coverage needed for recovery.
              </p>
              <Button onClick={() => runBenchmark("coverage")} disabled={running} className="w-full">
                {running ? (
                  <>
                    <Loader2 className="h-4 w-4 mr-2 animate-spin" /> Running...
                  </>
                ) : (
                  <>
                    <Gauge className="h-4 w-4 mr-2" /> Run Coverage Sweep
                  </>
                )}
              </Button>
              {coverageResult && <BenchmarkChart result={coverageResult} xKey="coverage" xLabel="Coverage" />}
            </TabsContent>
          </Tabs>

          {progress > 0 && <Progress value={progress} className="h-1" />}

          {(errorResult || coverageResult) && (
            <div className="grid gap-3 md:grid-cols-3">
              {errorResult && (
                <>
                  <SummaryCard
                    label="Max error recovered"
                    value={`${(errorResult.summary.maxErrorRateRecovered * 100).toFixed(2)}%`}
                  />
                  <SummaryCard
                    label="Avg density"
                    value={`${errorResult.summary.avgDensityBitsPerNt.toFixed(3)} bits/nt`}
                  />
                  <SummaryCard
                    label="Avg decode time"
                    value={`${errorResult.summary.avgDecodeMs.toFixed(0)} ms`}
                  />
                </>
              )}
            </div>
          )}
        </CardContent>
      </Card>
    </div>
  );
}

function SummaryCard({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded-lg border p-3 bg-muted/30">
      <div className="text-xs text-muted-foreground">{label}</div>
      <div className="font-mono text-lg font-bold">{value}</div>
    </div>
  );
}

function BenchmarkChart({
  result,
  xKey,
  xLabel,
}: {
  result: BenchmarkResultData;
  xKey: "errorRate" | "coverage";
  xLabel: string;
}) {
  const points = result.points;
  if (points.length === 0) return null;

  const width = 700;
  const height = 200;
  const padding = { top: 20, right: 20, bottom: 30, left: 40 };
  const chartW = width - padding.left - padding.right;
  const chartH = height - padding.top - padding.bottom;

  const xValues = points.map((p) => p[xKey]);
  const xMin = Math.min(...xValues);
  const xMax = Math.max(...xValues);
  const xRange = xMax - xMin || 1;

  const xScale = (v: number) => padding.left + ((v - xMin) / xRange) * chartW;
  const yScale = (success: boolean) => padding.top + (success ? 0.1 : 0.9) * chartH;

  // Build line for oligos recovered / total
  const recoveredLine = points.map((p, i) => {
    const x = xScale(p[xKey]);
    const y = padding.top + chartH - (p.oligosRecovered / p.oligoCount) * chartH;
    return `${i === 0 ? "M" : "L"} ${x} ${y}`;
  }).join(" ");

  return (
    <div className="border rounded-lg p-2 bg-card overflow-x-auto">
      <svg width={width} height={height} className="w-full">
        {/* Axes */}
        <line x1={padding.left} y1={padding.top} x2={padding.left} y2={padding.top + chartH} stroke="#94a3b8" />
        <line x1={padding.left} y1={padding.top + chartH} x2={padding.left + chartW} y2={padding.top + chartH} stroke="#94a3b8" />

        {/* Grid lines */}
        {[0, 0.5, 1].map((v) => {
          const y = padding.top + chartH - v * chartH;
          return (
            <g key={v}>
              <line x1={padding.left} y1={y} x2={padding.left + chartW} y2={y} stroke="#e2e8f0" strokeDasharray="2 2" />
              <text x={padding.left - 4} y={y + 3} fontSize={9} fill="#64748b" textAnchor="end">
                {(v * 100).toFixed(0)}%
              </text>
            </g>
          );
        })}

        {/* Recovered line */}
        <path d={recoveredLine} fill="none" stroke="#0ea5e9" strokeWidth={2} />

        {/* Success/fail points */}
        {points.map((p, i) => {
          const x = xScale(p[xKey]);
          const y = yScale(p.hashMatch);
          return (
            <g key={i}>
              <circle
                cx={x}
                cy={y}
                r={5}
                fill={p.hashMatch ? "#22c55e" : "#ef4444"}
                stroke="#fff"
                strokeWidth={1}
              />
              <text x={x} y={padding.top + chartH + 15} fontSize={9} fill="#64748b" textAnchor="middle">
                {xKey === "errorRate" ? `${(p.errorRate * 100).toFixed(1)}%` : `${p.coverage}x`}
              </text>
            </g>
          );
        })}

        {/* X-axis label */}
        <text x={padding.left + chartW / 2} y={height - 4} fontSize={10} fill="#64748b" textAnchor="middle">
          {xLabel}
        </text>

        {/* Legend */}
        <g transform={`translate(${padding.left + 10}, ${padding.top + 5})`}>
          <circle cx={0} cy={0} r={4} fill="#22c55e" />
          <text x={8} y={3} fontSize={9} fill="#64748b">Recovered</text>
          <circle cx={70} cy={0} r={4} fill="#ef4444" />
          <text x={78} y={3} fontSize={9} fill="#64748b">Failed</text>
        </g>
      </svg>

      {/* Table */}
      <div className="mt-3 max-h-48 overflow-y-auto">
        <table className="w-full text-xs">
          <thead className="sticky top-0 bg-card">
            <tr className="text-left">
              <th className="p-1">{xLabel}</th>
              <th className="p-1">Reads</th>
              <th className="p-1">Oligos recovered</th>
              <th className="p-1">Erased</th>
              <th className="p-1">Decode ms</th>
              <th className="p-1">Result</th>
            </tr>
          </thead>
          <tbody>
            {points.map((p, i) => (
              <tr key={i} className="border-t">
                <td className="p-1 font-mono">
                  {xKey === "errorRate" ? `${(p.errorRate * 100).toFixed(2)}%` : `${p.coverage}x`}
                </td>
                <td className="p-1 font-mono">{p.readsTotal}</td>
                <td className="p-1 font-mono">
                  {p.oligosRecovered}/{p.oligoCount}
                </td>
                <td className="p-1 font-mono">{p.oligosErased}</td>
                <td className="p-1 font-mono">{p.decodeMs}</td>
                <td className="p-1">
                  <span className={p.hashMatch ? "text-green-600" : "text-red-600"}>
                    {p.hashMatch ? "OK" : "FAIL"}
                  </span>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
