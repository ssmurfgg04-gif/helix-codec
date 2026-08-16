/**
 * Round-trip benchmark for helix-codec: encode → decode on real and synthetic datasets.
 *
 * Tests multiple Illumina codec configs across diverse data types
 * (genomic, text, random). Uses coverage=10 (Illumina-typical) with zero
 * channel errors for round-trip verification.
 *
 * Run: npx tsx scripts/roundtrip-bench.ts
 *
 * KNOWN ISSUES:
 *   1. LDPC "auto" decoder has ~5-8% per-oligo failure rate even at zero
 *      channel errors (hard→BP fallback doesn't always converge). Outer RS
 *      must have sufficient parity to recover these erasures.
 *      → Fix: use outerParityRatio ≥ 0.10 for reliable round-trip at scale.
 *   2. Nanopore + useConvolutionalInner: encode fails with "offset is out
 *      of bounds" — CRC marker insertion makes conv-encoded data larger
 *      than layout expects (codec.ts line 431). Bug in layout calculation.
 *   3. Nanopore without conv inner: decode fails with "require is not
 *      defined" — WASM module loading issue in ESM context.
 *   4. coverage=1 round-trips fail at scale: single-read decode can't
 *      resolve LDPC ambiguities. Use coverage ≥ 5 for reliable decode.
 */

import { readFile } from "fs/promises";
import { join } from "path";
import { gunzipSync } from "zlib";
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { computeDensity } from "../src/lib/dna/presets";
import type { CodecConfig } from "../src/lib/dna/types";
import { createHash } from "crypto";

// ─── Config presets to benchmark ───────────────────────────────────

const CONFIGS: Record<string, CodecConfig> = {
  // Standard Illumina — 300nt oligos, constrained mapping, 8B LDPC parity, 10% outer RS
  "constrained-300nt": {
    oligoLength: 300,
    primerLength: 12,
    innerCode: "ldpc",
    ldpcDecoder: "auto",
    mappingMode: "constrained",
    innerParityBytes: 8,
    outerParityRatio: 0.1,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 1,
    interleaveDepth: 0,
    channel: "illumina",
    lowCoverageTrigger: 5,
  },
  // High-density Illumina — 700nt oligos, 8B parity, 10% outer RS
  "hd-700nt": {
    oligoLength: 700,
    primerLength: 12,
    innerCode: "ldpc",
    ldpcDecoder: "auto",
    mappingMode: "constrained",
    innerParityBytes: 8,
    outerParityRatio: 0.1,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 1,
    interleaveDepth: 0,
    channel: "illumina",
    lowCoverageTrigger: 5,
  },
  // Ultra-high-density Illumina — 1100nt oligos, 8B parity, 10% outer RS
  "uhd-1100nt": {
    oligoLength: 1100,
    primerLength: 12,
    innerCode: "ldpc",
    ldpcDecoder: "auto",
    mappingMode: "constrained",
    innerParityBytes: 8,
    outerParityRatio: 0.1,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 1,
    interleaveDepth: 0,
    channel: "illumina",
    lowCoverageTrigger: 5,
  },
  // Yin-Yang coding — 2.0 bits/nt, homopolymer-free by construction
  "yinyang-300nt": {
    oligoLength: 300,
    primerLength: 12,
    innerCode: "ldpc",
    ldpcDecoder: "auto",
    mappingMode: "yinyang",
    innerParityBytes: 8,
    outerParityRatio: 0.1,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 1,
    interleaveDepth: 0,
    channel: "illumina",
    lowCoverageTrigger: 5,
  },
};

// ─── Dataset definitions ───────────────────────────────────────────

interface Dataset {
  name: string;
  data: () => Promise<Uint8Array>;
  description: string;
}

const TEST_DATA_DIR = join(import.meta.dirname ?? ".", "..", "test-data");

function sha256hex(data: Uint8Array): string {
  return createHash("sha256").update(data).digest("hex");
}

const DATASETS: Dataset[] = [
  {
    name: "random-1kb",
    description: "Random binary (1 KB)",
    data: async () => {
      const buf = await readFile(join(TEST_DATA_DIR, "random_1kb.bin"));
      return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
    },
  },
  {
    name: "random-10kb",
    description: "Random binary (10 KB)",
    data: async () => {
      const buf = await readFile(join(TEST_DATA_DIR, "random_10kb.bin"));
      return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
    },
  },
  {
    name: "pride-and-prejudice",
    description: "Gutenberg text (~723KB)",
    data: async () => {
      const buf = await readFile(join(TEST_DATA_DIR, "pride_and_prejudice.txt"));
      return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
    },
  },
  {
    name: "random-100kb",
    description: "Random binary (100 KB)",
    data: async () => {
      const buf = await readFile(join(TEST_DATA_DIR, "random_100kb.bin"));
      return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
    },
  },
  {
    name: "random-1mb",
    description: "Random binary (1 MB)",
    data: async () => {
      const buf = await readFile(join(TEST_DATA_DIR, "random_1mb.bin"));
      return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
    },
  },
  {
    name: "ecoli-k12",
    description: "E. coli K-12 genome (~4.6MB FASTA)",
    data: async () => {
      const gz = await readFile(join(TEST_DATA_DIR, "ecoli.fa.gz"));
      const text = gunzipSync(gz);
      return new Uint8Array(text.buffer, text.byteOffset, text.byteLength);
    },
  },
];

// ─── Result type ───────────────────────────────────────────────────

interface BenchResult {
  dataset: string;
  config: string;
  inputBytes: number;
  encodeMs: number;
  oligoCount: number;
  decodeMs: number;
  success: boolean;
  netDensityBitsPerNt: number;
  theoreticalDensityBitsPerNt: number;
  compressedBytes: number;
  screeningRetries: number;
  oligosRecovered: number;
  oligosErased: number;
  failInner: number;
  failOuter: number;
}

// ─── Helpers ───────────────────────────────────────────────────────

function formatBytes(bytes: number): string {
  if (bytes < 1024) return `${bytes}B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)}KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)}MB`;
}

const COVERAGE = 10;

// ─── Main ──────────────────────────────────────────────────────────

async function main() {
  console.log("╔══════════════════════════════════════════════════════════════════╗");
  console.log("║   helix-codec round-trip benchmark  encode → decode (clean)     ║");
  console.log(`║   Coverage = ${COVERAGE}× (Illumina-typical, zero channel errors)          ║`);
  console.log("╚══════════════════════════════════════════════════════════════════╝\n");

  const results: BenchResult[] = [];

  for (const ds of DATASETS) {
    console.log(`\n── Dataset: ${ds.name} (${ds.description}) ──`);

    let data: Uint8Array;
    try {
      data = await ds.data();
      console.log(`  Loaded: ${formatBytes(data.length)} (SHA256: ${sha256hex(data).slice(0, 16)}…)`);
    } catch (e: any) {
      console.warn(`  ⚠ Failed to load: ${e.message}. Skipping.`);
      continue;
    }

    for (const [configName, cfg] of Object.entries(CONFIGS)) {
      process.stdout.write(`  → ${configName} … `);

      try {
        // Encode
        const encodeResult = await encodeFile(data, cfg, {
          fileName: `${ds.name}.bin`,
          contentType: "application/octet-stream",
        });

        const encodeMs = encodeResult.stats.encodeTimeMs;
        const oligoCount = encodeResult.stats.oligoCount;
        const netDensity = encodeResult.stats.netDensityBitsPerNt;
        const compressedBytes = encodeResult.stats.compressedSize;
        const screeningRetries = encodeResult.stats.screeningRetries;

        // Simulate clean channel with coverage
        const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_CLEAN, coverage: COVERAGE, simulator: "basic" });

        // Decode
        const t0 = Date.now();
        const decodeResult = await decodeReads(
          sim.reads,
          encodeResult.encoded.metadata,
          cfg,
          encodeResult.encoded.forwardPrimer,
          encodeResult.encoded.reversePrimer,
        );
        const decodeMs = Date.now() - t0;

        const success = decodeResult.hashMatches;
        const theoreticalDensity = computeDensity(cfg, "total");

        if (success) {
          console.log(`✓ enc ${encodeMs}ms | dec ${decodeMs}ms | ${oligoCount} oligos | ${netDensity.toFixed(3)} b/nt`);
        } else {
          console.log(`✗ HASH MISMATCH (rec ${decodeResult.stats.oligosRecovered}/${oligoCount} erased ${decodeResult.stats.oligosErased} failI ${decodeResult.stats.oligosFailedInnerRS} failO ${decodeResult.stats.oligosFailedOuterRS})`);
        }

        results.push({
          dataset: ds.name,
          config: configName,
          inputBytes: data.length,
          encodeMs,
          oligoCount,
          decodeMs,
          success,
          netDensityBitsPerNt: netDensity,
          theoreticalDensityBitsPerNt: theoreticalDensity,
          compressedBytes,
          screeningRetries,
          oligosRecovered: decodeResult.stats.oligosRecovered,
          oligosErased: decodeResult.stats.oligosErased,
          failInner: decodeResult.stats.oligosFailedInnerRS,
          failOuter: decodeResult.stats.oligosFailedOuterRS,
        });
      } catch (e: any) {
        console.log(`✗ ERROR: ${e.message?.slice(0, 100)}`);
        results.push({
          dataset: ds.name,
          config: configName,
          inputBytes: data.length,
          encodeMs: -1, oligoCount: -1, decodeMs: -1, success: false,
          netDensityBitsPerNt: -1,
          theoreticalDensityBitsPerNt: computeDensity(cfg, "total"),
          compressedBytes: -1, screeningRetries: -1,
          oligosRecovered: -1, oligosErased: -1, failInner: -1, failOuter: -1,
        });
      }
    }
  }

  // ─── Summary table ────────────────────────────────────────────────

  console.log("\n\n╔══════════════════════════════════════════════════════════════════╗");
  console.log("║                       BENCHMARK RESULTS                           ║");
  console.log("╚══════════════════════════════════════════════════════════════════╝\n");

  console.table(
    results.map((r) => ({
      dataset: r.dataset,
      config: r.config,
      size: formatBytes(r.inputBytes),
      "enc(ms)": r.encodeMs,
      oligos: r.oligoCount,
      "dec(ms)": r.decodeMs,
      pass: r.success ? "✓" : "✗",
      "net(b/nt)": r.netDensityBitsPerNt > 0 ? r.netDensityBitsPerNt.toFixed(3) : "—",
      "thy(b/nt)": r.theoreticalDensityBitsPerNt.toFixed(3),
      erased: r.oligosErased,
      retries: r.screeningRetries,
    })),
  );

  // ─── Pass/fail summary ────────────────────────────────────────────
  const passed = results.filter((r) => r.success).length;
  const failed = results.filter((r) => !r.success && r.encodeMs >= 0).length;
  const errored = results.filter((r) => r.encodeMs === -1).length;
  console.log(`\nTotal: ${results.length} | ✓ Passed: ${passed} | ✗ Failed: ${failed} | ✗ Errored: ${errored}`);

  // ─── Density summary by config ────────────────────────────────────
  console.log("\n── Net density by config (avg over passing datasets) ──");
  for (const configName of Object.keys(CONFIGS)) {
    const cr = results.filter((r) => r.config === configName && r.success);
    const all = results.filter((r) => r.config === configName);
    if (cr.length === 0) {
      console.log(`  ${configName}: no successful results`);
      continue;
    }
    const avgDen = cr.reduce((s, r) => s + r.netDensityBitsPerNt, 0) / cr.length;
    const avgEnc = cr.reduce((s, r) => s + r.encodeMs, 0) / cr.length;
    const avgDec = cr.reduce((s, r) => s + r.decodeMs, 0) / cr.length;
    console.log(`  ${configName}: avg ${avgDen.toFixed(3)} b/nt | enc ${avgEnc.toFixed(0)}ms | dec ${avgDec.toFixed(0)}ms | ${cr.length}/${all.length} pass`);
  }

  // ─── Throughput summary ───────────────────────────────────────────
  console.log("\n── Throughput by config (avg MB/s, passing datasets) ──");
  for (const configName of Object.keys(CONFIGS)) {
    const cr = results.filter((r) => r.config === configName && r.success);
    if (cr.length === 0) continue;
    const avgEncTput = cr.reduce((s, r) => s + r.inputBytes / (r.encodeMs / 1000), 0) / cr.length;
    const avgDecTput = cr.reduce((s, r) => s + r.inputBytes / (r.decodeMs / 1000), 0) / cr.length;
    console.log(`  ${configName}: enc ${(avgEncTput / 1e6).toFixed(2)} MB/s | dec ${(avgDecTput / 1e6).toFixed(2)} MB/s`);
  }

  // ─── Erasure rate analysis ────────────────────────────────────────
  console.log("\n── Per-oligo erasure rate (erased / total oligos) ──");
  for (const configName of Object.keys(CONFIGS)) {
    const cr = results.filter((r) => r.config === configName && r.oligoCount > 0);
    if (cr.length === 0) continue;
    const rates = cr.map(r => r.oligosErased / r.oligoCount);
    const avgRate = rates.reduce((s, r) => s + r, 0) / rates.length;
    console.log(`  ${configName}: avg ${(avgRate * 100).toFixed(1)}% (range ${(Math.min(...rates) * 100).toFixed(1)}%–${(Math.max(...rates) * 100).toFixed(1)}%)`);
  }

  if (failed + errored > 0) {
    console.log(`\n⚠ ${failed + errored} result(s) had issues — see details above.`);
  }

  // ─── Known issues ─────────────────────────────────────────────────
  console.log("\n── Known Issues ──");
  console.log("  1. LDPC auto-decoder ~5-8% per-oligo failure rate at zero errors.");
  console.log("     Requires outerParityRatio ≥ 0.10 for reliable recovery at scale.");
  console.log("  2. Nanopore + conv inner: encode bug (CRC markers exceed layout).");
  console.log("  3. Nanopore without conv inner: decode bug (require() in ESM).");
  console.log("  4. coverage=1 fails at scale — use coverage ≥ 5.");
}

main().catch((e) => {
  console.error("FATAL:", e);
  process.exit(1);
});
