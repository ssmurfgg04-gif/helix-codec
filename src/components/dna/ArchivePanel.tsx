"use client";

import { useState, useCallback } from "react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { Switch } from "@/components/ui/switch";
import { Progress } from "@/components/ui/progress";
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from "@/components/ui/accordion";
import { useToast } from "@/hooks/use-toast";
import { Loader2, Package, TreePine, Shield, Clock, Download } from "lucide-react";
import { EncodedFile } from "@/lib/dna/types";
import { BioArchive } from "@/lib/dna/bioarchive";

export interface ArchivePanelProps {
  encoded: EncodedFile;
  onArchiveCreated: (archive: BioArchive) => void;
}

export function ArchivePanel({ encoded, onArchiveCreated }: ArchivePanelProps) {
  const [useEncryption, setUseEncryption] = useState(false);
  const [password, setPassword] = useState("");
  const [useLifecycle, setUseLifecycle] = useState(true);
  const [retention, setRetention] = useState("100y");
  const [storageClass, setStorageClass] = useState("deep_bio_archive");
  const [creating, setCreating] = useState(false);
  const [progress, setProgress] = useState(0);
  const [archive, setArchive] = useState<BioArchive | null>(null);
  const { toast } = useToast();

  const handleCreate = useCallback(async () => {
    setCreating(true);
    setProgress(20);
    try {
      const encryption = useEncryption && password
        ? {
            cipher: "xchacha20-poly1305" as const,
            kdf: "argon2id" as const,
            keyId: "pending",
            salt: "",
            nonce: "",
          }
        : { cipher: "none" as const, kdf: "none" as const };

      const lifecycle = useLifecycle
        ? {
            retention,
            storageClass: storageClass as "hot" | "warm" | "cold" | "deep_bio_archive",
            replicationTarget: 3,
            migrationInterval: "10y",
            decayPolicy: "repair_if_mutation_gt_5%",
          }
        : undefined;

      const res = await fetch("/api/dna", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          op: "archive",
          encoded,
          encryption,
          lifecycle,
        }),
      });
      setProgress(70);
      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.error ?? "Archive creation failed");
      }
      const json = await res.json();
      setProgress(100);
      setArchive(json.archive as BioArchive);
      onArchiveCreated(json.archive as BioArchive);
      toast({
        title: "BioArchive created",
        description: `Archive ID: ${json.archive.manifest.archiveId}`,
      });
    } catch (e) {
      toast({ title: "Archive creation failed", description: (e as Error).message, variant: "destructive" });
    } finally {
      setCreating(false);
      setTimeout(() => setProgress(0), 1000);
    }
  }, [encoded, useEncryption, password, useLifecycle, retention, storageClass, onArchiveCreated, toast]);

  const handleDownloadManifest = useCallback(() => {
    if (!archive) return;
    const blob = new Blob([JSON.stringify(archive.manifest, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `manifest_${archive.manifest.archiveId}.json`;
    a.click();
    URL.revokeObjectURL(url);
  }, [archive]);

  const handleDownloadArchive = useCallback(() => {
    if (!archive) return;
    const blob = new Blob([JSON.stringify(archive, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `archive_${archive.manifest.archiveId}.bioarc.json`;
    a.click();
    URL.revokeObjectURL(url);
  }, [archive]);

  return (
    <div className="space-y-4">
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Package className="h-4 w-4" /> Create BioArchive Container
          </CardTitle>
          <CardDescription>
            Wrap the encoded oligos into a self-describing archive with manifest, Merkle tree,
            per-chunk checksums, and optional encryption + lifecycle policy.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="flex items-center justify-between">
            <div>
              <Label htmlFor="enc">Encrypt with XChaCha20-Poly1305</Label>
              <p className="text-xs text-muted-foreground">Authenticated encryption (HKDF-SHA256 key derivation)</p>
            </div>
            <Switch id="enc" checked={useEncryption} onCheckedChange={setUseEncryption} />
          </div>
          {useEncryption && (
            <div className="space-y-2">
              <Label htmlFor="pwd">Password</Label>
              <Input
                id="pwd"
                type="password"
                value={password}
                onChange={(e) => setPassword(e.target.value)}
                placeholder="Enter encryption password..."
              />
            </div>
          )}

          <div className="flex items-center justify-between">
            <div>
              <Label htmlFor="life">Lifecycle policy</Label>
              <p className="text-xs text-muted-foreground">Retention, storage class, migration interval</p>
            </div>
            <Switch id="life" checked={useLifecycle} onCheckedChange={setUseLifecycle} />
          </div>
          {useLifecycle && (
            <div className="grid grid-cols-2 gap-3">
              <div className="space-y-2">
                <Label htmlFor="ret">Retention</Label>
                <Input id="ret" value={retention} onChange={(e) => setRetention(e.target.value)} placeholder="100y" />
              </div>
              <div className="space-y-2">
                <Label htmlFor="sc">Storage class</Label>
                <select
                  id="sc"
                  value={storageClass}
                  onChange={(e) => setStorageClass(e.target.value)}
                  className="w-full rounded-md border px-3 py-2 text-sm bg-background"
                >
                  <option value="hot">Hot</option>
                  <option value="warm">Warm</option>
                  <option value="cold">Cold</option>
                  <option value="deep_bio_archive">Deep Bio Archive</option>
                </select>
              </div>
            </div>
          )}

          <Button onClick={handleCreate} disabled={creating} className="w-full">
            {creating ? (
              <>
                <Loader2 className="h-4 w-4 mr-2 animate-spin" /> Creating archive...
              </>
            ) : (
              <>
                <Package className="h-4 w-4 mr-2" /> Create BioArchive
              </>
            )}
          </Button>
          {progress > 0 && <Progress value={progress} className="h-1" />}
        </CardContent>
      </Card>

      {archive && (
        <>
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Package className="h-4 w-4" /> Archive Manifest
                <Badge variant="secondary">{archive.manifest.format}</Badge>
              </CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="grid grid-cols-2 gap-3 text-sm md:grid-cols-3">
                <Stat label="Archive ID" value={archive.manifest.archiveId} mono />
                <Stat label="Created" value={new Date(archive.manifest.created).toLocaleString()} />
                <Stat label="Original file" value={archive.manifest.payload.originalName} />
                <Stat label="Size" value={`${archive.manifest.payload.sizeBytes.toLocaleString()} B`} />
                <Stat label="SHA-256" value={`${archive.manifest.payload.sha256.slice(0, 16)}...`} mono />
                <Stat label="MIME" value={archive.manifest.payload.mime} />
              </div>

              <div className="grid grid-cols-2 gap-3 text-sm md:grid-cols-4">
                <Stat label="Data chunks" value={archive.manifest.chunking.dataChunks.toString()} />
                <Stat label="Parity chunks" value={archive.manifest.chunking.parityChunks.toString()} />
                <Stat label="Total chunks" value={archive.manifest.chunking.totalChunks.toString()} />
                <Stat label="Overhead ratio" value={`${archive.manifest.ecc.overheadRatio.toFixed(2)}x`} />
              </div>

              <div className="grid grid-cols-2 gap-3 text-sm md:grid-cols-3">
                <Stat
                  label="Encryption"
                  value={archive.manifest.encryption.cipher === "none" ? "None" : "XChaCha20-Poly1305"}
                  badge={archive.manifest.encryption.cipher === "none" ? "secondary" : "default"}
                />
                <Stat label="Inner ECC" value={archive.manifest.ecc.inner} />
                <Stat label="Outer ECC" value={archive.manifest.ecc.outer} />
              </div>

              {archive.manifest.lifecycle && (
                <div className="rounded-lg border p-3 space-y-2 bg-muted/30">
                  <div className="flex items-center gap-2 text-sm font-medium">
                    <Clock className="h-4 w-4" /> Lifecycle Policy
                  </div>
                  <div className="grid grid-cols-2 gap-2 text-xs md:grid-cols-4">
                    <Stat label="Retention" value={archive.manifest.lifecycle.retention} />
                    <Stat label="Storage class" value={archive.manifest.lifecycle.storageClass} />
                    <Stat label="Replication" value={`${archive.manifest.lifecycle.replicationTarget}x`} />
                    <Stat label="Migration" value={archive.manifest.lifecycle.migrationInterval} />
                  </div>
                  <div className="text-xs text-muted-foreground">
                    Decay policy: {archive.manifest.lifecycle.decayPolicy}
                  </div>
                </div>
              )}

              <div className="flex gap-2">
                <Button onClick={handleDownloadManifest} variant="outline" size="sm">
                  <Download className="h-3 w-3 mr-1" /> Download Manifest
                </Button>
                <Button onClick={handleDownloadArchive} variant="outline" size="sm">
                  <Download className="h-3 w-3 mr-1" /> Download .bioarc
                </Button>
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <TreePine className="h-4 w-4" /> Merkle Integrity Tree
              </CardTitle>
              <CardDescription>SHA-256 Merkle tree over all chunk payloads.</CardDescription>
            </CardHeader>
            <CardContent>
              <MerkleTreeView archive={archive} />
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Shield className="h-4 w-4" /> Chunk Index ({archive.chunks.length} chunks)
              </CardTitle>
              <CardDescription>FASTA-like chunk records with barcodes, addresses, and checksums.</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="max-h-96 overflow-y-auto">
                <table className="w-full text-xs">
                  <thead className="sticky top-0 bg-card">
                    <tr className="text-left">
                      <th className="p-1">ID</th>
                      <th className="p-1">Type</th>
                      <th className="p-1">Barcode</th>
                      <th className="p-1">Checksum</th>
                      <th className="p-1">GC</th>
                      <th className="p-1">MaxHp</th>
                      <th className="p-1">Seed</th>
                    </tr>
                  </thead>
                  <tbody>
                    {archive.chunks.map((c) => (
                      <tr key={c.chunkId} className="border-t hover:bg-muted/50">
                        <td className="p-1 font-mono">{c.chunkId}</td>
                        <td className="p-1">
                          <Badge variant="outline" className="text-[10px]">
                            {c.type.toUpperCase()}
                          </Badge>
                        </td>
                        <td className="p-1 font-mono text-[10px]">{c.barcode}</td>
                        <td className="p-1 font-mono text-[10px]">{c.checksum.slice(0, 12)}...</td>
                        <td className="p-1 font-mono">{(c.gc * 100).toFixed(0)}%</td>
                        <td className="p-1 font-mono">{c.maxHomopolymer}</td>
                        <td className="p-1 font-mono">{c.seed}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
        </>
      )}
    </div>
  );
}

function MerkleTreeView({ archive }: { archive: BioArchive }) {
  const tree = archive.merkleTree;
  return (
    <div className="space-y-3">
      <div className="grid grid-cols-2 gap-3 text-sm md:grid-cols-4">
        <Stat label="Algorithm" value="sha256" />
        <Stat label="Leaves" value={tree.leaves.length.toString()} />
        <Stat label="Depth" value={tree.depth.toString()} />
        <Stat label="Root" value={`${tree.root.slice(0, 16)}...`} mono />
      </div>
      <Accordion type="single" collapsible>
        <AccordionItem value="tree">
          <AccordionTrigger className="text-sm">
            View tree levels ({tree.nodes.length} levels)
          </AccordionTrigger>
          <AccordionContent>
            <div className="space-y-2 max-h-96 overflow-y-auto">
              {tree.nodes.map((level, i) => (
                <div key={i} className="space-y-1">
                  <div className="text-xs font-medium text-muted-foreground">
                    Level {i} ({level.length} nodes{i === 0 ? " — leaves" : i === tree.nodes.length - 1 ? " — root" : ""})
                  </div>
                  <div className="flex flex-wrap gap-1">
                    {level.slice(0, 32).map((hash, j) => (
                      <span
                        key={j}
                        className="font-mono text-[10px] px-1.5 py-0.5 rounded bg-muted"
                        title={hash}
                      >
                        {hash.slice(0, 8)}
                      </span>
                    ))}
                    {level.length > 32 && (
                      <span className="text-[10px] text-muted-foreground">+{level.length - 32} more</span>
                    )}
                  </div>
                </div>
              ))}
            </div>
          </AccordionContent>
        </AccordionItem>
      </Accordion>
    </div>
  );
}

function Stat({
  label,
  value,
  mono,
  badge,
}: {
  label: string;
  value: string;
  mono?: boolean;
  badge?: "default" | "secondary";
}) {
  return (
    <div className="space-y-0.5">
      <div className="text-xs text-muted-foreground">{label}</div>
      {badge ? (
        <Badge variant={badge}>{value}</Badge>
      ) : (
        <div className={`text-sm ${mono ? "font-mono" : ""} truncate`}>{value}</div>
      )}
    </div>
  );
}
