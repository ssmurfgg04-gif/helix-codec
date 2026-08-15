"use client";

import { useState, useCallback, useRef, useEffect } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { Droplets, Play, RotateCcw, CheckCircle2, XCircle, Loader2 } from "lucide-react";

interface DropletInfo {
  seed: number;
  degree: number;
  sources: number[];
  payload: string; // hex for display
}

interface PeelingStep {
  iteration: number;
  chunksRemaining: number;
  degree1Droplets: number;
  recoveredChunk: number | null;
  totalChunks: number;
}

interface FountainResult {
  droplets: DropletInfo[];
  peelingSteps: PeelingStep[];
  decoded: boolean;
  K: number;
  overhead: number;
  encodeTimeMs: number;
  decodeTimeMs: number;
}

export function FountainPanel() {
  const [payload, setPayload] = useState("Hello, DNA Fountain! This is a test of the LT code peeling decoder.");
  const [chunkSize, setChunkSize] = useState(4);
  const [overhead, setOverhead] = useState(0.2);
  const [result, setResult] = useState<FountainResult | null>(null);
  const [running, setRunning] = useState(false);
  const [currentStep, setCurrentStep] = useState(0);
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const runFountain = useCallback(async () => {
    setRunning(true);
    setResult(null);
    setCurrentStep(0);

    // Small delay to let UI update
    await new Promise((r) => setTimeout(r, 50));

    try {
      const data = new TextEncoder().encode(payload);
      const K = Math.max(1, Math.ceil(data.length / chunkSize));
      const numDroplets = Math.ceil(K * (1 + overhead));

      // Simple PRNG (xorshift32)
      let prngState = 42 >>> 0;
      const rng = () => {
        prngState ^= prngState << 13;
        prngState ^= prngState >>> 17;
        prngState ^= prngState << 5;
        prngState = prngState >>> 0;
        return prngState / 0x100000000;
      };
      const rngInt = (max: number) => Math.floor(rng() * max);

      // Build chunks
      const chunks: Uint8Array[] = [];
      for (let i = 0; i < K; i++) {
        const chunk = new Uint8Array(chunkSize);
        for (let b = 0; b < chunkSize; b++) {
          const idx = i * chunkSize + b;
          chunk[b] = idx < data.length ? data[idx] : 0;
        }
        chunks.push(chunk);
      }

      // RSD degree distribution (simplified)
      const sampleDegree = (): number => {
        const rho = (d: number) => d === 1 ? 1 / K : 1 / (d * (d - 1));
        const c = 0.1;
      const delta = 0.5;
      const S = c * Math.log(K / delta) * Math.sqrt(K);
      const KOverS = Math.floor(K / S);
      const tau = (d: number) => {
        if (d <= KOverS - 1) return S / (K * d);
        if (d === KOverS) return (S * Math.log(S / delta)) / K;
        return 0;
      };
      let Z = 0;
      const mu: number[] = [];
      for (let d = 1; d <= K; d++) {
        mu.push(rho(d) + tau(d));
        Z += mu[d - 1];
      }
      const r = rng();
      let cum = 0;
      for (let d = 1; d <= K; d++) {
        cum += mu[d - 1] / Z;
        if (r <= cum) return d;
      }
      return K;
      };

      // Generate droplets
      const encodeStart = performance.now();
      const droplets: DropletInfo[] = [];
      for (let i = 0; i < numDroplets; i++) {
        const seed = Math.floor(rng() * 0xffffffff);
        // Use seed to derive degree and sources
        let s = seed >>> 0;
        const seedRng = () => {
          s ^= s << 13; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
          return s / 0x100000000;
        };
        const degree = sampleDegree();
        const sources = new Set<number>();
        while (sources.size < degree) {
          sources.add(Math.floor(seedRng() * K));
        }
        const srcArr = Array.from(sources).sort((a, b) => a - b);
        // XOR chunks
        const payloadBytes = new Uint8Array(chunkSize);
        for (const idx of srcArr) {
          for (let b = 0; b < chunkSize; b++) {
            payloadBytes[b] ^= chunks[idx][b];
          }
        }
        droplets.push({
          seed,
          degree,
          sources: srcArr,
          payload: Array.from(payloadBytes).map((b) => b.toString(16).padStart(2, "0")).join(""),
        });
      }
      const encodeTimeMs = performance.now() - encodeStart;

      // Peeling decoder with step tracking
      const decodeStart = performance.now();
      const recovered: (Uint8Array | null)[] = new Array(K).fill(null);
      const remaining = new Set<number>();
      for (let i = 0; i < K; i++) remaining.add(i);

      const workDroplets = droplets.map((d) => ({
        ...d,
        workPayload: droplets.find((dd) => dd.seed === d.seed)!.payload
          .match(/.{2}/g)!.map((h) => parseInt(h, 16)),
        remainingSources: new Set(d.sources),
      }));

      const peelingSteps: PeelingStep[] = [];
      let iter = 0;
      let progress = true;

      while (progress && remaining.size > 0) {
        progress = false;
        const deg1: typeof workDroplets = [];
        for (const wd of workDroplets) {
          if (wd.remainingSources.size === 1) deg1.push(wd);
        }

        peelingSteps.push({
          iteration: iter,
          chunksRemaining: remaining.size,
          degree1Droplets: deg1.length,
          recoveredChunk: null,
          totalChunks: K,
        });

        for (const wd of deg1) {
          if (wd.remainingSources.size !== 1) continue;
          const chunkIdx = Array.from(wd.remainingSources)[0];
          if (recovered[chunkIdx] !== null) {
            wd.remainingSources.delete(chunkIdx);
            continue;
          }
          recovered[chunkIdx] = new Uint8Array(workDroplets.find((d) => d.seed === wd.seed)!.workPayload);
          remaining.delete(chunkIdx);
          progress = true;

          for (const other of workDroplets) {
            if (other === wd) continue;
            if (other.remainingSources.has(chunkIdx)) {
              for (let b = 0; b < chunkSize; b++) {
                other.workPayload[b] ^= recovered[chunkIdx]![b];
              }
              other.remainingSources.delete(chunkIdx);
            }
          }
          wd.remainingSources.clear();

          peelingSteps.push({
            iteration: iter,
            chunksRemaining: remaining.size,
            degree1Droplets: deg1.length - 1,
            recoveredChunk: chunkIdx,
            totalChunks: K,
          });
        }
        iter++;
      }

      const decodeTimeMs = performance.now() - decodeStart;
      const decoded = remaining.size === 0;

      setResult({
        droplets: droplets.slice(0, 50), // Show first 50
        peelingSteps,
        decoded,
        K,
        overhead: numDroplets / K - 1,
        encodeTimeMs,
        decodeTimeMs,
      });
    } catch (e) {
      console.error("Fountain error:", e);
    } finally {
      setRunning(false);
    }
  }, [payload, chunkSize, overhead]);

  // Draw the Tanner graph visualization
  const drawGraph = useCallback(
    (step: number) => {
      const canvas = canvasRef.current;
      if (!canvas || !result) return;
      const ctx = canvas.getContext("2d");
      if (!ctx) return;

      const w = canvas.width;
      const h = canvas.height;
      ctx.clearRect(0, 0, w, h);

      const K = result.K;
      const currentStepData = result.peelingSteps[step];
      if (!currentStepData) return;

      // Determine which chunks are recovered at this step
      const recoveredSet = new Set<number>();
      for (let i = 0; i <= step; i++) {
        const s = result.peelingSteps[i];
        if (s.recoveredChunk !== null) recoveredSet.add(s.recoveredChunk);
      }

      // Draw chunks (left side)
      const chunkRadius = Math.min(20, 300 / K);
      const chunkSpacing = (h - 40) / K;
      const chunkX = 60;

      for (let i = 0; i < K; i++) {
        const y = 20 + i * chunkSpacing;
        const isRecovered = recoveredSet.has(i);
        ctx.beginPath();
        ctx.arc(chunkX, y, chunkRadius, 0, 2 * Math.PI);
        ctx.fillStyle = isRecovered ? "#10b981" : "#e2e8f0";
        ctx.fill();
        ctx.strokeStyle = isRecovered ? "#059669" : "#94a3b8";
        ctx.lineWidth = 2;
        ctx.stroke();
        ctx.fillStyle = isRecovered ? "#fff" : "#475569";
        ctx.font = `${Math.max(8, chunkRadius)}px monospace`;
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(String(i), chunkX, y);
      }

      // Draw droplets (right side) - show first 20
      const numShowDroplets = Math.min(20, result.droplets.length);
      const dropletSpacing = (h - 40) / numShowDroplets;
      const dropletX = w - 60;

      for (let i = 0; i < numShowDroplets; i++) {
        const d = result.droplets[i];
        const y = 20 + i * dropletSpacing;
        const isDeg1 = d.degree === 1;
        ctx.beginPath();
        ctx.arc(dropletX, y, chunkRadius * 0.8, 0, 2 * Math.PI);
        ctx.fillStyle = isDeg1 ? "#f59e0b" : "#3b82f6";
        ctx.fill();
        ctx.strokeStyle = isDeg1 ? "#d97706" : "#2563eb";
        ctx.lineWidth = 1.5;
        ctx.stroke();
        ctx.fillStyle = "#fff";
        ctx.font = `${Math.max(7, chunkRadius * 0.7)}px monospace`;
        ctx.fillText(`d${d.degree}`, dropletX, y);
      }

      // Draw edges
      ctx.strokeStyle = "#cbd5e1";
      ctx.lineWidth = 0.5;
      for (let i = 0; i < numShowDroplets; i++) {
        const d = result.droplets[i];
        const dy = 20 + i * dropletSpacing;
        for (const src of d.sources) {
          if (src >= K) continue;
          const sy = 20 + src * chunkSpacing;
          ctx.beginPath();
          ctx.moveTo(chunkX + chunkRadius, sy);
          ctx.lineTo(dropletX - chunkRadius * 0.8, dy);
          ctx.stroke();
        }
      }
    },
    [result],
  );

  // Redraw canvas when step or result changes
  useEffect(() => {
    if (result) {
      drawGraph(currentStep);
    }
  }, [currentStep, result, drawGraph]);

  return (
    <div className="space-y-4">
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Droplets className="h-5 w-5" />
            DNA Fountain Encoder / Decoder
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2">
            <Label htmlFor="payload">Payload (text)</Label>
            <textarea
              id="payload"
              className="w-full min-h-[80px] rounded-md border border-input bg-background px-3 py-2 text-sm font-mono"
              value={payload}
              onChange={(e) => setPayload(e.target.value)}
              disabled={running}
            />
          </div>
          <div className="grid grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label htmlFor="chunkSize">Chunk Size (bytes): {chunkSize}</Label>
              <Input
                id="chunkSize"
                type="range"
                min={2}
                max={32}
                value={chunkSize}
                onChange={(e) => setChunkSize(Number(e.target.value))}
                disabled={running}
              />
            </div>
            <div className="space-y-2">
              <Label htmlFor="overhead">Overhead: {(overhead * 100).toFixed(0)}%</Label>
              <Input
                id="overhead"
                type="range"
                min={0}
                max={1}
                step={0.05}
                value={overhead}
                onChange={(e) => setOverhead(Number(e.target.value))}
                disabled={running}
              />
            </div>
          </div>
          <div className="flex gap-2">
            <Button onClick={runFountain} disabled={running}>
              {running ? (
                <>
                  <Loader2 className="h-4 w-4 mr-2 animate-spin" />
                  Running...
                </>
              ) : (
                <>
                  <Play className="h-4 w-4 mr-2" />
                  Run Fountain Code
                </>
              )}
            </Button>
            {result && (
              <Button
                variant="outline"
                onClick={() => {
                  setResult(null);
                  setCurrentStep(0);
                }}
              >
                <RotateCcw className="h-4 w-4 mr-2" />
                Reset
              </Button>
            )}
          </div>
        </CardContent>
      </Card>

      {result && (
        <>
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                Results
                {result.decoded ? (
                  <CheckCircle2 className="h-5 w-5 text-green-500" />
                ) : (
                  <XCircle className="h-5 w-5 text-red-500" />
                )}
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                <div>
                  <p className="text-xs text-muted-foreground">Data Chunks (K)</p>
                  <p className="text-2xl font-bold">{result.K}</p>
                </div>
                <div>
                  <p className="text-xs text-muted-foreground">Droplets</p>
                  <p className="text-2xl font-bold">{result.droplets.length}{result.droplets.length === 50 ? "+" : ""}</p>
                </div>
                <div>
                  <p className="text-xs text-muted-foreground">Overhead</p>
                  <p className="text-2xl font-bold">{(result.overhead * 100).toFixed(0)}%</p>
                </div>
                <div>
                  <p className="text-xs text-muted-foreground">Decode Time</p>
                  <p className="text-2xl font-bold">{result.decodeTimeMs.toFixed(1)}ms</p>
                </div>
              </div>
              <div className="mt-4 flex gap-2">
                <Badge variant={result.decoded ? "default" : "destructive"}>
                  {result.decoded ? "✅ Decoded Successfully" : "❌ Decoding Failed"}
                </Badge>
                <Badge variant="secondary">Peeling Steps: {result.peelingSteps.filter((s) => s.recoveredChunk !== null).length}</Badge>
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle>Tanner Graph Visualization</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="flex justify-between items-center mb-2 text-xs text-muted-foreground">
                <span>← Data Chunks (green = recovered)</span>
                <span>Droplets (orange = degree 1) →</span>
              </div>
              <canvas
                ref={canvasRef}
                width={600}
                height={400}
                className="w-full border rounded-md bg-background"
              />
              <div className="mt-4 space-y-2">
                <div className="flex items-center justify-between">
                  <Label>Peeling Step: {currentStep} / {result.peelingSteps.length - 1}</Label>
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={() => {
                      const next = Math.min(currentStep + 1, result.peelingSteps.length - 1);
                      setCurrentStep(next);
                    }}
                    disabled={currentStep >= result.peelingSteps.length - 1}
                  >
                    Next Step →
                  </Button>
                </div>
                <Progress value={(currentStep / Math.max(1, result.peelingSteps.length - 1)) * 100} />
                {result.peelingSteps[currentStep] && (
                  <div className="text-xs space-y-1 mt-2 p-2 bg-muted/30 rounded">
                    <p>Iteration: {result.peelingSteps[currentStep].iteration}</p>
                    <p>Chunks remaining: {result.peelingSteps[currentStep].chunksRemaining} / {result.K}</p>
                    <p>Degree-1 droplets: {result.peelingSteps[currentStep].degree1Droplets}</p>
                    {result.peelingSteps[currentStep].recoveredChunk !== null && (
                      <p className="text-green-600 font-medium">
                        → Recovered chunk #{result.peelingSteps[currentStep].recoveredChunk}
                      </p>
                    )}
                  </div>
                )}
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle>Droplets (first 50)</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="space-y-1 max-h-64 overflow-y-auto text-xs font-mono">
                {result.droplets.map((d, i) => (
                  <div key={i} className="flex items-center gap-2 border-b border-border/50 py-1">
                    <span className="text-muted-foreground w-8">#{i}</span>
                    <Badge variant={d.degree === 1 ? "default" : "secondary"} className="w-12">
                      d={d.degree}
                    </Badge>
                    <span className="text-muted-foreground">
                      src=[{d.sources.slice(0, 8).join(",")}{d.sources.length > 8 ? "..." : ""}]
                    </span>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>
        </>
      )}
    </div>
  );
}
