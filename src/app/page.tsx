"use client";

import { useState, useCallback } from "react";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Dna, Microscope, FlaskConical, Gauge, FileText, Package, Hexagon, Droplets } from "lucide-react";
import { EncodedFile, CodecConfig, DEFAULT_CONFIG } from "@/lib/dna/types";
import { SimulationResult } from "@/lib/dna/simulate";
import { BioArchive } from "@/lib/dna/bioarchive";
import { EncodePanel, EncodeStatsCard, EncodeStats } from "@/components/dna/EncodePanel";
import { SimulatePanel, SimulationStatsCard, SimConfig } from "@/components/dna/SimulatePanel";
import { DecodePanel } from "@/components/dna/DecodePanel";
import { BenchmarkPanel } from "@/components/dna/BenchmarkPanel";
import { GenomeBrowser } from "@/components/dna/GenomeBrowser";
import { SpecPanel } from "@/components/dna/SpecPanel";
import { ArchivePanel } from "@/components/dna/ArchivePanel";
import { HolographicPanel } from "@/components/dna/HolographicPanel";
import { FountainPanel } from "@/components/dna/FountainPanel";

export default function Home() {
  const [encoded, setEncoded] = useState<EncodedFile | null>(null);
  const [encodeStats, setEncodeStats] = useState<EncodeStats | null>(null);
  const [simulation, setSimulation] = useState<SimulationResult | null>(null);
  const [simConfig, setSimConfig] = useState<SimConfig | null>(null);
  const [selectedOligo, setSelectedOligo] = useState<number | null>(null);
  const [archive, setArchive] = useState<BioArchive | null>(null);

  const handleEncoded = useCallback((enc: EncodedFile, stats: EncodeStats) => {
    setEncoded(enc);
    setEncodeStats(stats);
    setSimulation(null);
    setSimConfig(null);
    setArchive(null);
    setSelectedOligo(0);
  }, []);

  const handleSimulated = useCallback((sim: SimulationResult, cfg: SimConfig) => {
    setSimulation(sim);
    setSimConfig(cfg);
  }, []);

  const handleArchiveCreated = useCallback((arc: BioArchive) => {
    setArchive(arc);
  }, []);

  return (
    <div className="min-h-screen bg-gradient-to-br from-background via-background to-emerald-950/20">
      <header className="border-b bg-background/80 backdrop-blur sticky top-0 z-40">
        <div className="container mx-auto px-4 py-3 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="flex h-9 w-9 items-center justify-center rounded-lg bg-primary text-primary-foreground">
              <Dna className="h-5 w-5" />
            </div>
            <div>
              <h1 className="text-lg font-bold leading-none">Helix Codec</h1>
              <p className="text-xs text-muted-foreground">Biological Archival File Format</p>
            </div>
          </div>
          <div className="flex items-center gap-2">
            <Badge variant="outline" className="text-xs">v2.0</Badge>
            <Badge variant="outline" className="text-xs hidden md:flex">
              <Dna className="h-3 w-3 mr-1" /> BioArchive · Merkle · XChaCha20 · Holographic Sharding
            </Badge>
          </div>
        </div>
      </header>

      <main className="container mx-auto px-4 py-6 space-y-6">
        {!encoded && (
          <Card className="border-primary/30 bg-gradient-to-br from-primary/5 to-emerald-950/10">
            <CardContent className="pt-6">
              <div className="space-y-3 max-w-3xl">
                <h2 className="text-2xl font-bold tracking-tight">
                  A mutation-aware biological archival file format.
                </h2>
                <p className="text-muted-foreground leading-relaxed">
                  Modern archival systems protect static bits against hardware failure.
                  <strong> BioArchive</strong> protects information against mutation, decay, and
                  generational drift — designed for a medium that mutates, degrades, replicates,
                  and has biochemical constraints. Includes Reed-Solomon ECC, Merkle integrity trees,
                  XChaCha20-Poly1305 encryption, lifecycle policies, audit/scrubbing, generational
                  lineage, and a novel Holographic DNA Sharding codec for 1.5x redundancy.
                </p>
                <div className="flex flex-wrap gap-2 pt-2">
                  <Badge variant="secondary">BioArchive v1 container</Badge>
                  <Badge variant="secondary">SHA-256 Merkle tree</Badge>
                  <Badge variant="secondary">XChaCha20-Poly1305</Badge>
                  <Badge variant="secondary">RS(38,30) inner + outer erasure</Badge>
                  <Badge variant="secondary">Holographic sharding (1.5x)</Badge>
                  <Badge variant="secondary">Generational lineage</Badge>
                  <Badge variant="secondary">Audit &amp; scrubbing</Badge>
                  <Badge variant="secondary">Lifecycle policies</Badge>
                </div>
              </div>
            </CardContent>
          </Card>
        )}

        <Tabs defaultValue="studio" className="space-y-4">
          <TabsList className="grid w-full grid-cols-2 md:grid-cols-8 h-auto">
            <TabsTrigger value="studio" className="flex flex-col items-center gap-1 py-2">
              <Dna className="h-4 w-4" />
              <span className="text-xs">Studio</span>
            </TabsTrigger>
            <TabsTrigger value="browser" disabled={!encoded} className="flex flex-col items-center gap-1 py-2">
              <Microscope className="h-4 w-4" />
              <span className="text-xs">Browser</span>
            </TabsTrigger>
            <TabsTrigger value="mutate" disabled={!encoded} className="flex flex-col items-center gap-1 py-2">
              <FlaskConical className="h-4 w-4" />
              <span className="text-xs">Recover</span>
            </TabsTrigger>
            <TabsTrigger value="archive" disabled={!encoded} className="flex flex-col items-center gap-1 py-2">
              <Package className="h-4 w-4" />
              <span className="text-xs">Archive</span>
            </TabsTrigger>
            <TabsTrigger value="holographic" className="flex flex-col items-center gap-1 py-2">
              <Hexagon className="h-4 w-4" />
              <span className="text-xs">Holographic</span>
            </TabsTrigger>
            <TabsTrigger value="fountain" className="flex flex-col items-center gap-1 py-2">
              <Droplets className="h-4 w-4" />
              <span className="text-xs">Fountain</span>
            </TabsTrigger>
            <TabsTrigger value="benchmark" className="flex flex-col items-center gap-1 py-2">
              <Gauge className="h-4 w-4" />
              <span className="text-xs">Benchmark</span>
            </TabsTrigger>
            <TabsTrigger value="spec" className="flex flex-col items-center gap-1 py-2">
              <FileText className="h-4 w-4" />
              <span className="text-xs">Spec</span>
            </TabsTrigger>
          </TabsList>

          <TabsContent value="studio" className="space-y-4">
            <EncodePanel onEncoded={handleEncoded} />
            {encodeStats && <EncodeStatsCard stats={encodeStats} />}
            {encoded && (
              <Card>
                <CardHeader>
                  <CardTitle className="flex items-center gap-2">
                    <Dna className="h-4 w-4" /> Encoded Oligos
                    <Badge variant="secondary">{encoded.oligos.length}</Badge>
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="space-y-2 max-h-96 overflow-y-auto">
                    {encoded.oligos.map((oligo) => (
                      <div key={oligo.index} className="text-xs border rounded p-2 bg-muted/30">
                        <div className="flex justify-between mb-1">
                          <span className="font-mono">Oligo #{oligo.index}</span>
                          <span className="text-muted-foreground">
                            GC {(oligo.gc * 100).toFixed(1)}% · maxHp {oligo.maxHomopolymer} · seed {oligo.seed}
                          </span>
                        </div>
                        <div className="font-mono text-[10px] break-all leading-relaxed">
                          <span className="text-blue-600 dark:text-blue-400">{oligo.sequence.slice(0, 20)}</span>
                          <span>{oligo.sequence.slice(20, oligo.sequence.length - 20)}</span>
                          <span className="text-blue-600 dark:text-blue-400">{oligo.sequence.slice(oligo.sequence.length - 20)}</span>
                        </div>
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>
            )}
          </TabsContent>

          <TabsContent value="browser" className="space-y-4">
            {encoded && (
              <GenomeBrowser
                oligos={encoded.oligos}
                forwardPrimer={encoded.forwardPrimer}
                reversePrimer={encoded.reversePrimer}
                primerLength={DEFAULT_CONFIG.primerLength}
                selectedOligo={selectedOligo}
                onSelectOligo={setSelectedOligo}
              />
            )}
          </TabsContent>

          <TabsContent value="mutate" className="space-y-4">
            {encoded && (
              <>
                <SimulatePanel encoded={encoded} onSimulated={handleSimulated} />
                {simulation && <SimulationStatsCard result={simulation} />}
                <DecodePanel encoded={encoded} config={DEFAULT_CONFIG} simulation={simulation} />
              </>
            )}
          </TabsContent>

          <TabsContent value="archive" className="space-y-4">
            {encoded && <ArchivePanel encoded={encoded} onArchiveCreated={handleArchiveCreated} />}
          </TabsContent>

          <TabsContent value="holographic" className="space-y-4">
            <HolographicPanel />
          </TabsContent>

          <TabsContent value="fountain" className="space-y-4">
            <FountainPanel />
          </TabsContent>

          <TabsContent value="benchmark" className="space-y-4">
            <BenchmarkPanel />
          </TabsContent>

          <TabsContent value="spec" className="space-y-4">
            <SpecPanel />
          </TabsContent>
        </Tabs>
      </main>

      <footer className="border-t mt-12 bg-background/50">
        <div className="container mx-auto px-4 py-4 text-xs text-muted-foreground flex flex-wrap items-center justify-between gap-2">
          <div>
            Helix Codec v2.0 · BioArchive · Mutation-aware biological archival file format
          </div>
          <div className="flex items-center gap-3">
            <span>Inspired by Goldman 2013, Erlich &amp; Zielinski 2017, HEDGES 2020, Shamir 1979, ZFS</span>
          </div>
        </div>
      </footer>
    </div>
  );
}
