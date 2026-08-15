"use client";

import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Dna, Layers, Shield, GitBranch, FileCode, Database, Package, Hexagon, TreePine, Clock } from "lucide-react";
import { DEFAULT_CONFIG, computeLayout } from "@/lib/dna/types";

const layout = computeLayout(DEFAULT_CONFIG);

export function SpecPanel() {
  return (
    <div className="space-y-4">
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <FileCode className="h-4 w-4" /> Helix Codec v1 — Format Specification
          </CardTitle>
          <CardDescription>
            A biological archival file format for synthetic DNA data storage.
          </CardDescription>
        </CardHeader>
        <CardContent className="prose prose-sm dark:prose-invert max-w-none space-y-4">
          <div>
            <h3 className="text-base font-semibold flex items-center gap-2">
              <Dna className="h-4 w-4 text-primary" /> Oligo Structure
            </h3>
            <p className="text-sm text-muted-foreground">
              Each oligo is a single synthetic DNA strand of <code className="font-mono">{DEFAULT_CONFIG.oligoLength} nt</code>,
              consisting of a forward primer, an inner payload block, and a reverse primer.
            </p>
            <pre className="bg-muted/50 rounded p-3 text-xs overflow-x-auto font-mono">
{`5' [FWD PRIMER ${DEFAULT_CONFIG.primerLength}nt] [ADDRESS 16nt] [PAYLOAD ${layout.payloadBytes * 4}nt] [INNER RS PARITY ${layout.innerParityBytes * 4}nt] [CRC-16 8nt] [REV PRIMER ${DEFAULT_CONFIG.primerLength}nt] 3'
                              └─────────── inner block: ${layout.totalInnerBytes * 4}nt = ${layout.totalInnerBytes} bytes ───────────┘`}
            </pre>
          </div>

          <div className="grid gap-4 md:grid-cols-2">
            <SpecSection
              icon={<Layers className="h-4 w-4" />}
              title="Encoding Pipeline"
              steps={[
                "Read input file as bytes.",
                "DEFLATE compress (level 9, pako).",
                "Compute SHA-256 hash of original data.",
                `Split into ${layout.payloadBytes}-byte payload chunks (one per data oligo).`,
                `Apply outer RS(${DEFAULT_CONFIG.oligoLength}, k) across oligos — ${Math.round(DEFAULT_CONFIG.outerParityRatio * 100)}% parity.`,
                "Per oligo: build inner block = address(4B) + payload + inner RS parity + CRC-16.",
                "Apply fixed whitening XOR to address (breaks up zero patterns).",
                "Map bytes → DNA via direct 2-bit mapping (00=A, 01=C, 10=G, 11=T).",
                "Screen for GC content (40–60%) and homopolymer (max 3).",
                "If screen fails: XOR payload with seed-derived keystream, re-encode. Try up to 64 seeds.",
                "Prepend forward primer, append reverse primer.",
              ]}
            />
            <SpecSection
              icon={<Shield className="h-4 w-4" />}
              title="Decoding Pipeline"
              steps={[
                "Trim primers from each read (Hamming distance ≤ 2).",
                "Convert DNA → bytes (2-bit inverse mapping).",
                "Unwhiten address; extract 3-byte index for clustering.",
                "Cluster reads by oligo index.",
                "Per-cluster: column-wise plurality consensus on DNA.",
                "Decode consensus DNA → bytes.",
                "Verify CRC-16; flag failures.",
                `Apply inner RS(${layout.addressBytes + layout.payloadBytes + layout.innerParityBytes}, ${layout.addressBytes + layout.payloadBytes}) to correct residual errors.`,
                "Unwhiten address; extract seed; reverse XOR-with-seed on payload.",
                "Outer RS with erasure decoding: recover missing oligos.",
                "Concatenate payloads in index order; trim to fileSize.",
                "DEFLATE decompress.",
                "Verify SHA-256 hash matches metadata.",
              ]}
            />
          </div>

          <div className="grid gap-4 md:grid-cols-2">
            <SpecSection
              icon={<GitBranch className="h-4 w-4" />}
              title="Error Correction"
              steps={[
                `Inner RS: RS(${layout.addressBytes + layout.payloadBytes + layout.innerParityBytes}, ${layout.addressBytes + layout.payloadBytes}) = ${layout.innerParityBytes} parity bytes per oligo.`,
                `  → Corrects up to ${Math.floor(layout.innerParityBytes / 2)} unknown errors per oligo.`,
                `  → Or up to ${layout.innerParityBytes} erasures per oligo (if positions known).`,
                `Outer RS: RS(n, k) across oligos, with ${Math.round(DEFAULT_CONFIG.outerParityRatio * 100)}% parity.`,
                "  → Recovers missing oligos (synthesis dropouts) via erasure decoding.",
                "  → Uses pure-erasure Forney algorithm (positions known from address index).",
                "CRC-16/CCITT-FALSE per oligo: detects residual errors after RS.",
                "SHA-256 file hash: end-to-end integrity verification.",
              ]}
            />
            <SpecSection
              icon={<Database className="h-4 w-4" />}
              title="Biological Constraints"
              steps={[
                `GC content: 40–60% (target 45–55%).`,
                `Max homopolymer run: 3 nt.`,
                `Oligo length: ${DEFAULT_CONFIG.oligoLength} nt (within synthesis sweet spot 150–230 nt).`,
                `Primer length: ${DEFAULT_CONFIG.primerLength} nt each end (standard PCR primer length).`,
                "Address whitening: fixed XOR pattern (0x1B, 0x4B, 0x24, 0x6D) breaks up zero patterns for small oligo indices.",
                "Seed-based re-encoding: xorshift32 keystream derived from 1-byte seed, up to 64 retries per oligo.",
                "Direct 2-bit mapping: 2 bits/base theoretical density (vs. 1.58 bits/base for Goldman trit coding).",
              ]}
            />
          </div>

          <div>
            <h3 className="text-base font-semibold flex items-center gap-2">
              <Package className="h-4 w-4 text-primary" /> BioArchive Container (v1)
            </h3>
            <p className="text-sm text-muted-foreground mb-2">
              A self-describing archive container wrapping the encoded oligos. Inspired by ZFS
              (checksums everywhere), object storage (manifests), PAR2 (recovery blocks), and
              Git/IPFS (content-addressed Merkle trees).
            </p>
            <pre className="bg-muted/50 rounded p-3 text-xs overflow-x-auto font-mono">
{`archive.bioarc/
  manifest.json       # format, payload, chunking, ECC, bio encoding, recovery, lifecycle
  merkle.json         # SHA-256 Merkle tree (leaves = per-chunk payload hashes)
  chunks/
    chunk_000001.dna  # FASTA-like: >chunk_000001 barcode=... address=... checksum=...
    chunk_000002.dna
    ...
  metadata.enc        # (optional) XChaCha20-Poly1305 encrypted metadata
  lineage.json        # generational history (mutation rates, repair events)
  recovery_report.json # audit / scrubbing results`}
            </pre>
            <div className="grid gap-4 md:grid-cols-2 mt-3">
              <SpecSection
                icon={<TreePine className="h-4 w-4" />}
                title="Merkle Integrity Tree"
                steps={[
                  "Each chunk's payload is hashed with SHA-256 → leaf.",
                  "Leaves are paired and hashed up to a root (binary tree).",
                  "Root hash becomes the archive ID (first 16 hex chars).",
                  "Verification: recompute leaves, compare to stored tree.",
                  "Merkle proofs: O(log n) proof for any single chunk.",
                  "Tamper detection: any single byte change propagates to root.",
                ]}
              />
              <SpecSection
                icon={<Shield className="h-4 w-4" />}
                title="Encryption Layer"
                steps={[
                  "Cipher: XChaCha20-Poly1305 (authenticated, 24-byte nonce).",
                  "KDF: HKDF-SHA256 (Argon2id recommended for password hardening).",
                  "Pipeline: compress → encrypt → chunk → ECC → DNA encode.",
                  "Metadata can be encrypted separately (metadata.enc).",
                  "DNA sequence reveals nothing about filename, structure, or plaintext.",
                  "Key ID (SHA-256 of derived key) for key management without exposing key.",
                ]}
              />
              <SpecSection
                icon={<Clock className="h-4 w-4" />}
                title="Lifecycle & Audit"
                steps={[
                  "Retention policy (e.g. 100y), storage class (hot/warm/cold/deep).",
                  "Replication target (e.g. 3x physical copies).",
                  "Migration interval (e.g. 10y — regenerate from healthy copy).",
                  "Decay policy (e.g. repair if mutation > 5%).",
                  "Audit/scrubbing: verify checksums, detect mutation hotspots.",
                  "Recovery probability: (healthy + repaired) / total chunks.",
                ]}
              />
              <SpecSection
                icon={<GitBranch className="h-4 w-4" />}
                title="Generational Lineage"
                steps={[
                  "Each replication event = one generation.",
                  "Track: observed mutation rate, repair events, chunk health.",
                  "Parent hash links to previous generation's Merkle root.",
                  "Repair log: per-chunk corrections (consensus, inner ECC, outer ECC).",
                  "Recovery probability tracked across generations.",
                  "Enables 'evolution' analysis of archive health over time.",
                ]}
              />
            </div>
          </div>

          <div>
            <h3 className="text-base font-semibold flex items-center gap-2">
              <Hexagon className="h-4 w-4 text-primary" /> Holographic DNA Sharding Codec
            </h3>
            <p className="text-sm text-muted-foreground mb-2">
              A novel erasure code where each shard carries a polynomial projection of the entire
              dataset. Achieves 100% recovery at 1.5x redundancy (vs. 2x for traditional RS).
            </p>
            <div className="grid gap-4 md:grid-cols-2">
              <SpecSection
                icon={<Hexagon className="h-4 w-4" />}
                title="Holographic Encoding"
                steps={[
                  "Split data into blocks of K bytes.",
                  "Per block: build polynomial P(x) = sum data[i] * x^i over GF(256).",
                  "Evaluate P at N distinct points (x = 1..N) → N shard bytes.",
                  "Bijective Feistel shuffle: spread adjacent bytes across shards.",
                  "Each shard contains evaluations from ALL blocks (holographic).",
                  "Losing any shard degrades all blocks equally, not specific data.",
                ]}
              />
              <SpecSection
                icon={<Hexagon className="h-4 w-4" />}
                title="Holographic Decoding"
                steps={[
                  "Take any K of N available shards.",
                  "Per block: solve Vandermonde system via GF(256) Gaussian elimination.",
                  "Lagrange interpolation recovers polynomial coefficients = original bytes.",
                  "Reverse Feistel shuffle to restore byte order.",
                  "Threshold: recovery succeeds iff shards available >= K.",
                  "Max loss tolerance: (N - K) / N = 33% at 1.5x overhead.",
                ]}
              />
            </div>
          </div>

          <div>
            <h3 className="text-base font-semibold">Parameters</h3>
            <div className="grid gap-2 md:grid-cols-2 mt-2">
              <ParamTable
                title="Codec Parameters"
                rows={[
                  ["Version", "1"],
                  ["Oligo length", `${DEFAULT_CONFIG.oligoLength} nt`],
                  ["Primer length", `${DEFAULT_CONFIG.primerLength} nt × 2`],
                  ["Address", `${layout.addressBytes} B (3B index + 1B seed)`],
                  ["Payload per oligo", `${layout.payloadBytes} B (${layout.payloadBytes * 4} nt)`],
                  ["Inner RS parity", `${layout.innerParityBytes} B (RS ${layout.addressBytes + layout.payloadBytes + layout.innerParityBytes},${layout.addressBytes + layout.payloadBytes})`],
                  ["Outer RS parity", `${Math.round(DEFAULT_CONFIG.outerParityRatio * 100)}% of data oligos`],
                  ["CRC", "CRC-16/CCITT-FALSE (2 B)"],
                  ["Compression", "DEFLATE (pako, level 9)"],
                  ["Hash", "SHA-256"],
                ]}
              />
              <ParamTable
                title="Mutation Model Defaults"
                rows={[
                  ["Illumina substitution", "1e-3 per position"],
                  ["Illumina insertion", "5e-4 per position"],
                  ["Illumina deletion", "1e-3 per position"],
                  ["Illumina coverage", "20×"],
                  ["Nanopore total", "~9% (indel-heavy)"],
                  ["PacBio total", "~8.5% (insertion-heavy)"],
                  ["Dropout rate", "0–5% (synthesis failure)"],
                  ["Coverage range", "1× to 50×"],
                  ["PRNG", "xorshift32 (deterministic with seed)"],
                ]}
              />
            </div>
          </div>

          <div>
            <h3 className="text-base font-semibold">Theoretical Density</h3>
            <div className="space-y-2 text-sm">
              <div className="flex items-center justify-between border-b pb-1">
                <span>Shannon ceiling (theoretical)</span>
                <Badge variant="outline">2.000 bits/nt</Badge>
              </div>
              <div className="flex items-center justify-between border-b pb-1">
                <span>DNA Fountain (Erlich &amp; Zielinski 2017)</span>
                <Badge variant="outline">1.570 bits/nt</Badge>
              </div>
              <div className="flex items-center justify-between border-b pb-1">
                <span>Goldman 2013 (rotational trit)</span>
                <Badge variant="outline">0.830 bits/nt</Badge>
              </div>
              <div className="flex items-center justify-between border-b pb-1">
                <span>Helix Codec v1 (direct 2-bit, with overhead)</span>
                <Badge variant="secondary">~0.6–1.2 bits/nt</Badge>
              </div>
              <div className="flex items-center justify-between">
                <span>Helix Codec v1 (raw, before overhead)</span>
                <Badge variant="secondary">2.000 bits/nt</Badge>
              </div>
            </div>
            <p className="text-xs text-muted-foreground mt-2">
              Density is reduced by: primers ({Math.round((DEFAULT_CONFIG.primerLength * 2 / DEFAULT_CONFIG.oligoLength) * 100)}%),
              address ({Math.round((layout.addressBytes / layout.totalInnerBytes) * 100)}%),
              inner RS ({Math.round((layout.innerParityBytes / layout.totalInnerBytes) * 100)}%),
              outer RS ({Math.round(DEFAULT_CONFIG.outerParityRatio * 100)}%),
              CRC ({Math.round((layout.crcBytes / layout.totalInnerBytes) * 100)}%).
            </p>
          </div>

          <div>
            <h3 className="text-base font-semibold">References</h3>
            <ul className="text-xs space-y-1 text-muted-foreground list-disc pl-4">
              <li>Goldman et al. (2013). &quot;Towards practical, high-capacity, low-maintenance information storage in synthesized DNA.&quot; Nature 494:77-80.</li>
              <li>Erlich &amp; Zielinski (2017). &quot;DNA Fountain enables a robust and efficient storage architecture.&quot; Science 355:6328.</li>
              <li>Organick et al. (2018). &quot;Random access in large-scale DNA data storage.&quot; Nature Biotechnology 36:242-248.</li>
              <li>Press, Jones et al. (2020). &quot;HEDGES error-correcting code for DNA storage corrects insertions and deletions.&quot; PNAS 117:31.</li>
              <li>Chandak et al. (2018). &quot;Improved read/write cost tradeoff in a DNA-based storage system using rotating fountain codes.&quot; ISMB.</li>
              <li>Wikiversity: &quot;Reed-Solomon codes for coders&quot; — implementation reference.</li>
              <li>Rizzo (1997). &quot;Effective Erasure Codes for Reliable Computer Communication Protocols.&quot; ACM SIGCOMM CCR.</li>
            </ul>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}

function SpecSection({
  icon,
  title,
  steps,
}: {
  icon: React.ReactNode;
  title: string;
  steps: string[];
}) {
  return (
    <div className="rounded-lg border p-3 space-y-2">
      <div className="flex items-center gap-2 font-medium text-sm">
        {icon} {title}
      </div>
      <ol className="text-xs space-y-1 list-decimal pl-4 text-muted-foreground">
        {steps.map((s, i) => (
          <li key={i} className="leading-relaxed">{s}</li>
        ))}
      </ol>
    </div>
  );
}

function ParamTable({ title, rows }: { title: string; rows: [string, string][] }) {
  return (
    <div className="rounded-lg border p-3 space-y-2">
      <div className="font-medium text-sm">{title}</div>
      <div className="space-y-1">
        {rows.map(([k, v]) => (
          <div key={k} className="flex justify-between text-xs">
            <span className="text-muted-foreground">{k}</span>
            <span className="font-mono">{v}</span>
          </div>
        ))}
      </div>
    </div>
  );
}
