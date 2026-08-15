"use client";

import { useState, useCallback } from "react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { useToast } from "@/hooks/use-toast";
import { Loader2, FlaskConical, Activity, AlertTriangle } from "lucide-react";
import { EncodedFile } from "@/lib/dna/types";
import { SimulationResult } from "@/lib/dna/simulate";

export interface SimulatePanelProps {
  encoded: EncodedFile;
  onSimulated: (result: SimulationResult, config: SimConfig) => void;
}

export interface SimConfig {
  substitutionRate: number;
  insertionRate: number;
  deletionRate: number;
  coverage: number;
  dropoutRate: number;
  seed: number;
}

const PRESETS: Record<string, { label: string; config: SimConfig; description: string }> = {
  real2024: {
    label: "Real 2024",
    config: { substitutionRate: 0.025, insertionRate: 0.016, deletionRate: 0.082, coverage: 25, dropoutRate: 0.02, seed: 42 },
    description: "Preuss et al. 2026 measured rates — del 8.2%, sub 2.5%, ins 1.6%",
  },
  illumina: {
    label: "Illumina",
    config: { substitutionRate: 0.001, insertionRate: 0.0005, deletionRate: 0.001, coverage: 20, dropoutRate: 0, seed: 42 },
    description: "Chandak 2018 — sub-dominant (~0.25% total)",
  },
  nanopore: {
    label: "Nanopore (ONT)",
    config: { substitutionRate: 0.02, insertionRate: 0.03, deletionRate: 0.04, coverage: 15, dropoutRate: 0.05, seed: 42 },
    description: "R10.4.1 — indel-heavy (~9% total)",
  },
  pacbio: {
    label: "PacBio",
    config: { substitutionRate: 0.005, insertionRate: 0.05, deletionRate: 0.03, coverage: 10, dropoutRate: 0.02, seed: 42 },
    description: "HiFi (CCS) — insertion-dominant (~8.5% total)",
  },
  high: {
    label: "High Stress",
    config: { substitutionRate: 0.02, insertionRate: 0.015, deletionRate: 0.015, coverage: 30, dropoutRate: 0.05, seed: 42 },
    description: "5% total error, 30x coverage",
  },
  clean: {
    label: "Perfect",
    config: { substitutionRate: 0, insertionRate: 0, deletionRate: 0, coverage: 1, dropoutRate: 0, seed: 42 },
    description: "No errors (sanity check)",
  },
};

export function SimulatePanel({ encoded, onSimulated }: SimulatePanelProps) {
  const [preset, setPreset] = useState<keyof typeof PRESETS>("real2024");
  const [config, setConfig] = useState<SimConfig>(PRESETS.real2024.config);
  const [simulating, setSimulating] = useState(false);
  const [progress, setProgress] = useState(0);
  const { toast } = useToast();

  const applyPreset = (key: keyof typeof PRESETS) => {
    setPreset(key);
    setConfig(PRESETS[key].config);
  };

  const handleSimulate = useCallback(async () => {
    setSimulating(true);
    setProgress(20);
    try {
      const res = await fetch("/api/dna", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          op: "simulate",
          encoded,
          mutationConfig: config,
        }),
      });
      setProgress(70);
      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.error ?? "Simulate failed");
      }
      const json = await res.json();
      setProgress(100);
      onSimulated(json.simulation as SimulationResult, config);
      toast({
        title: "Simulation complete",
        description: `${json.simulation.totalReads} reads, ${json.simulation.totalErrors} errors, ${json.simulation.droppedOligos.length} dropped`,
      });
    } catch (e) {
      toast({ title: "Simulation failed", description: (e as Error).message, variant: "destructive" });
    } finally {
      setSimulating(false);
      setTimeout(() => setProgress(0), 1000);
    }
  }, [encoded, config, onSimulated, toast]);

  const totalErrorRate = config.substitutionRate + config.insertionRate + config.deletionRate;
  const errorPct = (totalErrorRate * 100).toFixed(2);

  return (
    <div className="space-y-4">
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <FlaskConical className="h-4 w-4" /> Sequencing Error Simulator
          </CardTitle>
          <CardDescription>
            Model real-world DNA synthesis and sequencing errors. Each oligo is &quot;read&quot; multiple
            times (coverage) with independent errors.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {/* Presets */}
          <div className="space-y-2">
            <Label>Error Profile Presets</Label>
            <div className="grid grid-cols-2 gap-2 md:grid-cols-3">
              {(Object.keys(PRESETS) as (keyof typeof PRESETS)[]).map((key) => (
                <Button
                  key={key}
                  variant={preset === key ? "default" : "outline"}
                  size="sm"
                  onClick={() => applyPreset(key)}
                  className="justify-start"
                >
                  {PRESETS[key].label}
                </Button>
              ))}
            </div>
            <p className="text-xs text-muted-foreground">{PRESETS[preset].description}</p>
          </div>

          <Tabs defaultValue="rates">
            <TabsList className="grid w-full grid-cols-2">
              <TabsTrigger value="rates">Error Rates</TabsTrigger>
              <TabsTrigger value="coverage">Coverage & Dropout</TabsTrigger>
            </TabsList>
            <TabsContent value="rates" className="space-y-4 pt-4">
              <SliderRow
                label="Substitution rate"
                value={config.substitutionRate}
                min={0}
                max={0.05}
                step={0.0005}
                format={(v) => `${(v * 100).toFixed(3)}%`}
                onChange={(v) => setConfig({ ...config, substitutionRate: v })}
              />
              <SliderRow
                label="Insertion rate"
                value={config.insertionRate}
                min={0}
                max={0.05}
                step={0.0005}
                format={(v) => `${(v * 100).toFixed(3)}%`}
                onChange={(v) => setConfig({ ...config, insertionRate: v })}
              />
              <SliderRow
                label="Deletion rate"
                value={config.deletionRate}
                min={0}
                max={0.05}
                step={0.0005}
                format={(v) => `${(v * 100).toFixed(3)}%`}
                onChange={(v) => setConfig({ ...config, deletionRate: v })}
              />
              <div className="flex items-center gap-2 rounded-lg border p-3 bg-muted/30">
                <Activity className="h-4 w-4 text-primary" />
                <span className="text-sm">Total error rate:</span>
                <Badge variant={totalErrorRate > 0.02 ? "destructive" : "secondary"}>{errorPct}%</Badge>
                {totalErrorRate > 0.02 && (
                  <span className="text-xs text-destructive flex items-center gap-1">
                    <AlertTriangle className="h-3 w-3" /> May exceed RS capacity
                  </span>
                )}
              </div>
            </TabsContent>
            <TabsContent value="coverage" className="space-y-4 pt-4">
              <SliderRow
                label="Coverage depth"
                value={config.coverage}
                min={1}
                max={50}
                step={1}
                format={(v) => `${v}x`}
                onChange={(v) => setConfig({ ...config, coverage: v })}
              />
              <SliderRow
                label="Dropout rate (lost oligos)"
                value={config.dropoutRate}
                min={0}
                max={0.2}
                step={0.01}
                format={(v) => `${(v * 100).toFixed(1)}%`}
                onChange={(v) => setConfig({ ...config, dropoutRate: v })}
              />
              <SliderRow
                label="Random seed"
                value={config.seed}
                min={0}
                max={1000}
                step={1}
                format={(v) => v.toString()}
                onChange={(v) => setConfig({ ...config, seed: v })}
              />
            </TabsContent>
          </Tabs>

          <Button onClick={handleSimulate} disabled={simulating} className="w-full">
            {simulating ? (
              <>
                <Loader2 className="h-4 w-4 mr-2 animate-spin" /> Simulating...
              </>
            ) : (
              <>
                <FlaskConical className="h-4 w-4 mr-2" /> Run Simulation
              </>
            )}
          </Button>
          {progress > 0 && <Progress value={progress} className="h-1" />}
        </CardContent>
      </Card>
    </div>
  );
}

function SliderRow({
  label,
  value,
  min,
  max,
  step,
  format,
  onChange,
}: {
  label: string;
  value: number;
  min: number;
  max: number;
  step: number;
  format: (v: number) => string;
  onChange: (v: number) => void;
}) {
  return (
    <div className="space-y-2">
      <div className="flex justify-between text-sm">
        <Label>{label}</Label>
        <span className="font-mono text-muted-foreground">{format(value)}</span>
      </div>
      <Slider min={min} max={max} step={step} value={[value]} onValueChange={([v]) => onChange(v)} />
    </div>
  );
}

export function SimulationStatsCard({ result }: { result: SimulationResult }) {
  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Activity className="h-4 w-4" /> Simulation Results
        </CardTitle>
      </CardHeader>
      <CardContent>
        <div className="grid grid-cols-2 gap-3 text-sm md:grid-cols-4">
          <Stat label="Total reads" value={result.totalReads.toLocaleString()} />
          <Stat label="Avg read length" value={`${result.avgReadLength.toFixed(0)} nt`} />
          <Stat label="Total errors" value={result.totalErrors.toLocaleString()} />
          <Stat label="Dropped oligos" value={result.droppedOligos.length.toString()} />
          <Stat label="Simulation time" value={`${result.simulationTimeMs} ms`} />
        </div>
      </CardContent>
    </Card>
  );
}

function Stat({ label, value }: { label: string; value: string }) {
  return (
    <div className="space-y-0.5">
      <div className="text-xs text-muted-foreground">{label}</div>
      <div className="font-mono text-sm">{value}</div>
    </div>
  );
}
