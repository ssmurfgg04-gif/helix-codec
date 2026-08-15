"use client";

import { useState, useCallback } from "react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import { Slider } from "@/components/ui/slider";
import { Switch } from "@/components/ui/switch";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { useToast } from "@/hooks/use-toast";
import { Loader2, Upload, FileText, Dna, Download } from "lucide-react";
import { CodecConfig, DEFAULT_CONFIG, EncodedFile, computeLayout } from "@/lib/dna/types";

export interface EncodePanelProps {
  onEncoded: (encoded: EncodedFile, stats: EncodeStats) => void;
}

export interface EncodeStats {
  rawSize: number;
  compressedSize: number;
  oligoCount: number;
  payloadBytesPerOligo: number;
  netDensityBitsPerNt: number;
  overheadPercent: number;
  screeningRetries: number;
  encodeTimeMs: number;
}

export function EncodePanel({ onEncoded }: EncodePanelProps) {
  const [config, setConfig] = useState<CodecConfig>(DEFAULT_CONFIG);
  const [fileName, setFileName] = useState("hello.txt");
  const [textContent, setTextContent] = useState(
    "Hello, DNA! This is a test of synthetic DNA data storage. The quick brown fox jumps over the lazy dog. ".repeat(
      8,
    ),
  );
  const [fileBytes, setFileBytes] = useState<Uint8Array | null>(null);
  const [encoding, setEncoding] = useState(false);
  const [progress, setProgress] = useState(0);
  const { toast } = useToast();

  const layout = computeLayout(config);

  const handleFileUpload = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const file = e.target.files?.[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        const bytes = new Uint8Array(reader.result as ArrayBuffer);
        setFileBytes(bytes);
        setFileName(file.name);
        setTextContent("");
        toast({
          title: "File loaded",
          description: `${file.name} (${bytes.length.toLocaleString()} bytes)`,
        });
      };
      reader.readAsArrayBuffer(file);
    },
    [toast],
  );

  const handleEncode = async () => {
    let data: Uint8Array;
    let contentType: string;

    if (fileBytes) {
      data = fileBytes;
      contentType = "application/octet-stream";
    } else if (textContent) {
      data = new TextEncoder().encode(textContent);
      contentType = "text/plain";
    } else {
      toast({ title: "No input", description: "Please provide text or upload a file.", variant: "destructive" });
      return;
    }

    setEncoding(true);
    setProgress(10);
    try {
      const base64 = btoa(String.fromCharCode(...data));
      setProgress(30);
      const res = await fetch("/api/dna", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          op: "encode",
          data: base64,
          fileName,
          contentType,
          config,
        }),
      });
      setProgress(70);
      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.error ?? "Encode failed");
      }
      const json = await res.json();
      setProgress(100);
      onEncoded(json.encoded as EncodedFile, json.stats as EncodeStats);
      toast({
        title: "Encoded successfully",
        description: `${json.stats.oligoCount} oligos, ${(json.stats.netDensityBitsPerNt).toFixed(2)} bits/nt`,
      });
    } catch (e) {
      toast({ title: "Encode failed", description: (e as Error).message, variant: "destructive" });
    } finally {
      setEncoding(false);
      setTimeout(() => setProgress(0), 1000);
    }
  };

  return (
    <div className="grid gap-4 md:grid-cols-2">
      {/* Input */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Upload className="h-4 w-4" /> Input Data
          </CardTitle>
          <CardDescription>Provide text or upload a file to encode into DNA.</CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2">
            <Label htmlFor="file-name">File name</Label>
            <Input
              id="file-name"
              value={fileName}
              onChange={(e) => setFileName(e.target.value)}
              placeholder="hello.txt"
            />
          </div>

          <div className="space-y-2">
            <Label htmlFor="file-upload">Upload file (binary)</Label>
            <Input id="file-upload" type="file" onChange={handleFileUpload} />
            {fileBytes && (
              <p className="text-xs text-muted-foreground">
                Loaded {fileBytes.length.toLocaleString()} bytes
              </p>
            )}
          </div>

          <div className="space-y-2">
            <Label htmlFor="text-input">Or enter text directly</Label>
            <Textarea
              id="text-input"
              value={textContent}
              onChange={(e) => {
                setTextContent(e.target.value);
                if (fileBytes) setFileBytes(null);
              }}
              placeholder="Type text to encode..."
              className="min-h-[120px] font-mono text-xs"
            />
            <p className="text-xs text-muted-foreground">
              {new TextEncoder().encode(textContent).length.toLocaleString()} bytes
            </p>
          </div>
        </CardContent>
      </Card>

      {/* Configuration */}
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Dna className="h-4 w-4" /> Codec Configuration
          </CardTitle>
          <CardDescription>Tune oligo length, error correction, and constraints.</CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2">
            <div className="flex justify-between text-sm">
              <Label>Oligo length: {config.oligoLength} nt</Label>
              <span className="text-muted-foreground">total per strand</span>
            </div>
            <Slider
              min={100}
              max={300}
              step={20}
              value={[config.oligoLength]}
              onValueChange={([v]) => setConfig({ ...config, oligoLength: v })}
            />
          </div>

          <div className="space-y-2">
            <div className="flex justify-between text-sm">
              <Label>Primer length: {config.primerLength} nt</Label>
              <span className="text-muted-foreground">each end</span>
            </div>
            <Slider
              min={10}
              max={30}
              step={2}
              value={[config.primerLength]}
              onValueChange={([v]) => setConfig({ ...config, primerLength: v })}
            />
          </div>

          <div className="space-y-2">
            <div className="flex justify-between text-sm">
              <Label>Inner RS parity: {config.innerParityBytes} bytes</Label>
              <span className="text-muted-foreground">per oligo</span>
            </div>
            <Slider
              min={4}
              max={20}
              step={2}
              value={[config.innerParityBytes]}
              onValueChange={([v]) => setConfig({ ...config, innerParityBytes: v })}
            />
          </div>

          <div className="space-y-2">
            <div className="flex justify-between text-sm">
              <Label>Outer RS parity: {Math.round(config.outerParityRatio * 100)}%</Label>
              <span className="text-muted-foreground">of data oligos</span>
            </div>
            <Slider
              min={0}
              max={0.5}
              step={0.05}
              value={[config.outerParityRatio]}
              onValueChange={([v]) => setConfig({ ...config, outerParityRatio: v })}
            />
          </div>

          <div className="flex items-center justify-between">
            <Label htmlFor="compress">DEFLATE compression</Label>
            <Switch
              id="compress"
              checked={config.compress}
              onCheckedChange={(v) => setConfig({ ...config, compress: v })}
            />
          </div>

          <div className="rounded-lg border p-3 space-y-1 text-xs bg-muted/30">
            <div className="font-medium mb-1">Computed Layout</div>
            <div className="flex justify-between">
              <span className="text-muted-foreground">Payload per oligo</span>
              <span className="font-mono">{layout.payloadBytes} bytes ({layout.payloadBytes * 4} nt)</span>
            </div>
            <div className="flex justify-between">
              <span className="text-muted-foreground">Address</span>
              <span className="font-mono">{layout.addressBytes} bytes ({layout.addressBytes * 4} nt)</span>
            </div>
            <div className="flex justify-between">
              <span className="text-muted-foreground">Inner RS parity</span>
              <span className="font-mono">{layout.innerParityBytes} bytes ({layout.innerParityBytes * 4} nt)</span>
            </div>
            <div className="flex justify-between">
              <span className="text-muted-foreground">CRC-16</span>
              <span className="font-mono">{layout.crcBytes} bytes ({layout.crcBytes * 4} nt)</span>
            </div>
            <div className="flex justify-between border-t pt-1 mt-1">
              <span className="text-muted-foreground">Total inner</span>
              <span className="font-mono">{layout.totalInnerBytes} bytes ({layout.totalInnerBytes * 4} nt)</span>
            </div>
          </div>

          <Button onClick={handleEncode} disabled={encoding} className="w-full">
            {encoding ? (
              <>
                <Loader2 className="h-4 w-4 mr-2 animate-spin" /> Encoding...
              </>
            ) : (
              <>
                <Dna className="h-4 w-4 mr-2" /> Encode to DNA
              </>
            )}
          </Button>
          {progress > 0 && <Progress value={progress} className="h-1" />}
        </CardContent>
      </Card>
    </div>
  );
}

export function EncodeStatsCard({ stats }: { stats: EncodeStats }) {
  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <FileText className="h-4 w-4" /> Encoding Results
        </CardTitle>
      </CardHeader>
      <CardContent>
        <div className="grid grid-cols-2 gap-3 text-sm md:grid-cols-4">
          <Stat label="Original size" value={`${stats.rawSize.toLocaleString()} B`} />
          <Stat label="Compressed" value={`${stats.compressedSize.toLocaleString()} B`} />
          <Stat label="Oligos" value={stats.oligoCount.toString()} />
          <Stat label="Payload/oligo" value={`${stats.payloadBytesPerOligo} B`} />
          <Stat
            label="Net density"
            value={`${stats.netDensityBitsPerNt.toFixed(3)} bits/nt`}
            highlight
          />
          <Stat label="Overhead" value={`${stats.overheadPercent.toFixed(1)}%`} />
          <Stat label="Screening retries" value={stats.screeningRetries.toString()} />
          <Stat label="Encode time" value={`${stats.encodeTimeMs} ms`} />
        </div>
      </CardContent>
    </Card>
  );
}

function Stat({ label, value, highlight }: { label: string; value: string; highlight?: boolean }) {
  return (
    <div className="space-y-0.5">
      <div className="text-xs text-muted-foreground">{label}</div>
      <div className={`font-mono text-sm ${highlight ? "text-primary font-bold" : ""}`}>{value}</div>
    </div>
  );
}
