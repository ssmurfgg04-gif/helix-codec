"use client";

import { useState, useCallback } from "react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { Switch } from "@/components/ui/switch";
import { Label } from "@/components/ui/label";
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from "@/components/ui/accordion";
import { useToast } from "@/hooks/use-toast";
import { Loader2, FolderSearch, CheckCircle2, XCircle, Download, Zap } from "lucide-react";
import { EncodedFile, CodecConfig } from "@/lib/dna/types";
import { SimulationResult } from "@/lib/dna/simulate";

export interface DecodePanelProps {
  encoded: EncodedFile;
  config: CodecConfig;
  simulation: SimulationResult | null;
}

export interface DecodeResultData {
  data: string | null;
  hash: string;
  hashMatches: boolean;
  stats: {
    totalReads: number;
    readsUsed: number;
    clustersFormed: number;
    oligosRecovered: number;
    oligosErased: number;
    oligosFailedInnerRS: number;
    oligosFailedOuterRS: number;
    consensusSuccessRate: number;
    decodeTimeMs: number;
  };
  perOligo: {
    index: number;
    readCount: number;
    consensusLength: number;
    crcPassed: boolean;
    innerRS: { corrected: number; success: boolean };
    seed: number;
    payloadBytes: string;
    isParity: boolean;
  }[];
}

export function DecodePanel({ encoded, config, simulation }: DecodePanelProps) {
  const [decoding, setDecoding] = useState(false);
  const [progress, setProgress] = useState(0);
  const [result, setResult] = useState<DecodeResultData | null>(null);
  const [useSoftInfo, setUseSoftInfo] = useState(true);
  const { toast } = useToast();

  const handleDecode = useCallback(async () => {
    if (!simulation) {
      toast({ title: "No simulation", description: "Run a simulation first.", variant: "destructive" });
      return;
    }
    setDecoding(true);
    setProgress(20);
    try {
      const res = await fetch("/api/dna", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          op: "decode",
          reads: simulation.reads,
          metadata: encoded.metadata,
          config,
          forwardPrimer: encoded.forwardPrimer,
          reversePrimer: encoded.reversePrimer,
          useSoftInfo,
        }),
      });
      setProgress(70);
      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.error ?? "Decode failed");
      }
      const json = await res.json();
      setProgress(100);
      setResult(json.decoded as DecodeResultData);
      toast({
        title: json.decoded.hashMatches ? "Recovery successful" : "Recovery failed",
        description: json.decoded.hashMatches
          ? "SHA-256 hash matches original."
          : "Hash mismatch — data corrupted.",
        variant: json.decoded.hashMatches ? "default" : "destructive",
      });
    } catch (e) {
      toast({ title: "Decode failed", description: (e as Error).message, variant: "destructive" });
    } finally {
      setDecoding(false);
      setTimeout(() => setProgress(0), 1000);
    }
  }, [simulation, encoded, config, toast]);

  const handleDownload = useCallback(() => {
    if (!result?.data) return;
    const bytes = Uint8Array.from(atob(result.data), (c) => c.charCodeAt(0));
    const blob = new Blob([bytes], { type: "application/octet-stream" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `recovered_${encoded.metadata.fileName}`;
    a.click();
    URL.revokeObjectURL(url);
  }, [result, encoded.metadata.fileName]);

  return (
    <div className="space-y-4">
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <FolderSearch className="h-4 w-4" /> Recovery Engine
          </CardTitle>
          <CardDescription>
            Cluster reads, build consensus, apply Reed-Solomon, and recover the original file.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="flex items-center justify-between rounded-lg border p-3 bg-muted/30">
            <div>
              <Label htmlFor="softinfo" className="flex items-center gap-2">
                <Zap className="h-3 w-3" /> Soft-information passthrough
              </Label>
              <p className="text-xs text-muted-foreground mt-1">
                Use Q-score erasure hints (2x RS capacity). Banal 2026 style.
              </p>
            </div>
            <Switch id="softinfo" checked={useSoftInfo} onCheckedChange={setUseSoftInfo} />
          </div>

          <Button onClick={handleDecode} disabled={decoding || !simulation} className="w-full">
            {decoding ? (
              <>
                <Loader2 className="h-4 w-4 mr-2 animate-spin" /> Decoding...
              </>
            ) : (
              <>
                <FolderSearch className="h-4 w-4 mr-2" /> Decode Reads
              </>
            )}
          </Button>
          {progress > 0 && <Progress value={progress} className="h-1" />}

          {result && (
            <>
              <div className="rounded-lg border p-4 space-y-3">
                <div className="flex items-center justify-between">
                  <span className="text-sm font-medium">Recovery Status</span>
                  {result.hashMatches ? (
                    <Badge className="bg-green-600 hover:bg-green-600">
                      <CheckCircle2 className="h-3 w-3 mr-1" /> Verified
                    </Badge>
                  ) : (
                    <Badge variant="destructive">
                      <XCircle className="h-3 w-3 mr-1" /> Failed
                    </Badge>
                  )}
                </div>
                <div className="grid grid-cols-2 gap-3 text-sm md:grid-cols-4">
                  <Stat label="Reads used" value={`${result.stats.readsUsed}/${result.stats.totalReads}`} />
                  <Stat label="Clusters" value={result.stats.clustersFormed.toString()} />
                  <Stat label="Oligos recovered" value={result.stats.oligosRecovered.toString()} />
                  <Stat label="Oligos erased" value={result.stats.oligosErased.toString()} />
                  <Stat label="Inner RS failures" value={result.stats.oligosFailedInnerRS.toString()} />
                  <Stat label="Outer RS failures" value={result.stats.oligosFailedOuterRS.toString()} />
                  <Stat
                    label="Consensus rate"
                    value={`${(result.stats.consensusSuccessRate * 100).toFixed(1)}%`}
                  />
                  <Stat label="Decode time" value={`${result.stats.decodeTimeMs} ms`} />
                </div>
                <div className="text-xs text-muted-foreground">
                  SHA-256: <code className="font-mono">{result.hash.slice(0, 32)}...</code>
                </div>
                {result.data && (
                  <Button onClick={handleDownload} variant="outline" size="sm">
                    <Download className="h-3 w-3 mr-1" /> Download recovered file
                  </Button>
                )}
              </div>

              <Accordion type="single" collapsible>
                <AccordionItem value="per-oligo">
                  <AccordionTrigger className="text-sm">
                    Per-oligo recovery details ({result.perOligo.length} oligos)
                  </AccordionTrigger>
                  <AccordionContent>
                    <div className="max-h-96 overflow-y-auto">
                      <table className="w-full text-xs">
                        <thead className="sticky top-0 bg-background">
                          <tr className="text-left">
                            <th className="p-1">#</th>
                            <th className="p-1">Type</th>
                            <th className="p-1">Reads</th>
                            <th className="p-1">Cons. len</th>
                            <th className="p-1">CRC</th>
                            <th className="p-1">Inner RS</th>
                            <th className="p-1">Seed</th>
                          </tr>
                        </thead>
                        <tbody>
                          {result.perOligo.map((p) => (
                            <tr key={p.index} className="border-t">
                              <td className="p-1 font-mono">{p.index}</td>
                              <td className="p-1">
                                <Badge variant="outline" className="text-[10px]">
                                  {p.isParity ? "PARITY" : "DATA"}
                                </Badge>
                              </td>
                              <td className="p-1 font-mono">{p.readCount}</td>
                              <td className="p-1 font-mono">{p.consensusLength}</td>
                              <td className="p-1">
                                {p.crcPassed ? (
                                  <CheckCircle2 className="h-3 w-3 text-green-600" />
                                ) : (
                                  <XCircle className="h-3 w-3 text-red-600" />
                                )}
                              </td>
                              <td className="p-1 font-mono">
                                {p.innerRS.success ? `+${p.innerRS.corrected}` : "FAIL"}
                              </td>
                              <td className="p-1 font-mono">{p.seed}</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </AccordionContent>
                </AccordionItem>
              </Accordion>
            </>
          )}
        </CardContent>
      </Card>
    </div>
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
