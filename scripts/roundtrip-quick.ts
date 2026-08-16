/**
 * Quick smoke test for roundtrip (encode→decode) — small datasets only.
 * Run: npx tsx scripts/roundtrip-quick.ts
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { computeDensity } from "../src/lib/dna/presets";
import type { CodecConfig } from "../src/lib/dna/types";
import { readFile } from "fs/promises";
import { join } from "path";
import { createHash } from "crypto";

const TEST_DATA_DIR = join(import.meta.dirname ?? ".", "..", "test-data");

function sha256(data: Uint8Array): string {
  return createHash("sha256").update(data).digest("hex");
}

function formatBytes(bytes: number): string {
  if (bytes < 1024) return `${bytes}B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)}KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)}MB`;
}

const CONFIGS: Record<string, CodecConfig> = {
  "default-constrained": {
    oligoLength: 300,
    primerLength: 12,
    innerCode: "ldpc",
    ldpcDecoder: "auto",
    mappingMode: "constrained",
    innerParityBytes: 8,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 5,
    interleaveDepth: 0,
    channel: "illumina",
    lowCoverageTrigger: 5,
  },
  "v55-density": {
    oligoLength: 700,
    primerLength: 12,
    innerCode: "ldpc",
    ldpcDecoder: "auto",
    mappingMode: "constrained",
    innerParityBytes: 8,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 5,
    interleaveDepth: 0,
    channel: "illumina",
    lowCoverageTrigger: 5,
  },
  "nanopore-v61": {
    oligoLength: 150,
    primerLength: 12,
    innerCode: "ldpc",
    ldpcDecoder: "auto",
    mappingMode: "constrained",
    innerParityBytes: 10,
    outerParityRatio: 0.5,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 5,
    interleaveDepth: 0,
    channel: "nanopore",
    lowCoverageTrigger: 3,
    useConvolutionalInner: false,
  },
};

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
  screeningRetries: number;
}

async function runOne(
  data: Uint8Array,
  dsName: string,
  configName: string,
  cfg: CodecConfig,
): Promise<BenchResult> {
  const theoreticalDensity = computeDensity(cfg, "total");
  try {
    const encodeResult = await encodeFile(data, cfg, {
      fileName: `${dsName}.bin`,
      contentType: "application/octet-stream",
    });
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
    const t0 = Date.now();
    const decodeResult = await decodeReads(
      sim.reads,
      encodeResult.encoded.metadata,
      cfg,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );
    const decodeMs = Date.now() - t0;
    return {
      dataset: dsName,
      config: configName,
      inputBytes: data.length,
      encodeMs: encodeResult.stats.encodeTimeMs,
      oligoCount: encodeResult.stats.oligoCount,
      decodeMs,
      success: decodeResult.hashMatches,
      netDensityBitsPerNt: encodeResult.stats.netDensityBitsPerNt,
      theoreticalDensityBitsPerNt: theoreticalDensity,
      screeningRetries: encodeResult.stats.screeningRetries,
    };
  } catch (e: any) {
    console.error(`    ✗ ERROR: ${e.message?.slice(0, 200)}`);
    return {
      dataset: dsName,
      config: configName,
      inputBytes: data.length,
      encodeMs: -1,
      oligoCount: -1,
      decodeMs: -1,
      success: false,
      netDensityBitsPerNt: -1,
      theoreticalDensityBitsPerNt: theoreticalDensity,
      screeningRetries: -1,
    };
  }
}

async function main() {
  console.log("╔══════════════════════════════════════════════════════════════╗");
  console.log("║   helix-codec round-trip benchmark: encode → decode (clean)  ║");
  console.log("╚══════════════════════════════════════════════════════════════╝\n");

  const results: BenchResult[] = [];

  // Datasets (load incrementally)
  const datasets: Array<{ name: string; desc: string; load: () => Promise<Uint8Array> }> = [
    {
      name: "random-1kb",
      desc: "Random binary 1KB",
      load: async () => {
        const buf = await readFile(join(TEST_DATA_DIR, "random_1kb.bin"));
        return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
      },
    },
    {
      name: "random-10kb",
      desc: "Random binary 10KB",
      load: async () => {
        const buf = await readFile(join(TEST_DATA_DIR, "random_10kb.bin"));
        return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
      },
    },
    {
      name: "random-100kb",
      desc: "Random binary 100KB",
      load: async () => {
        const buf = await readFile(join(TEST_DATA_DIR, "random_100kb.bin"));
        return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
      },
    },
    {
      name: "pride-and-prejudice",
      desc: "Gutenberg text ~723KB",
      load: async () => {
        const buf = await readFile(join(TEST_DATA_DIR, "pride_and_prejudice.txt"));
        return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
      },
    },
    {
      name: "random-1mb",
      desc: "Random binary 1MB",
      load: async () => {
        const buf = await readFile(join(TEST_DATA_DIR, "random_1mb.bin"));
        return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
      },
    },
  ];

  for (const ds of datasets) {
    console.log(`\n── Dataset: ${ds.name} (${ds.desc}) ──`);
    let data: Uint8Array;
    try {
      data = await ds.load();
      console.log(`  Loaded: ${formatBytes(data.length)} (SHA256: ${sha256(data).slice(0, 16)}…)`);
    } catch (e: any) {
      console.warn(`  ⚠ Failed: ${e.message}. Skipping.`);
      continue;
    }

    for (const [configName, cfg] of Object.entries(CONFIGS)) {
      process.stdout.write(`  → ${configName} … `);
      const r = await runOne(data, ds.name, configName, cfg);
      results.push(r);
      if (r.success) {
        console.log(`✓ enc ${r.encodeMs}ms | dec ${r.decodeMs}ms | ${r.oligoCount} oligos | ${r.netDensityBitsPerNt.toFixed(3)} b/nt`);
      } else if (r.encodeMs === -1) {
        console.log(`✗ ERROR`);
      } else {
        console.log(`✗ HASH MISMATCH`);
      }
    }
  }

  // Summary table
  console.log("\n\n╔══════════════════════════════════════════════════════════════╗");
  console.log("║                     BENCHMARK RESULTS                        ║");
  console.log("╚══════════════════════════════════════════════════════════════╝\n");

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
      retries: r.screeningRetries,
    })),
  );

  const passed = results.filter((r) => r.success).length;
  const failed = results.filter((r) => !r.success && r.encodeMs >= 0).length;
  const errored = results.filter((r) => r.encodeMs === -1).length;
  console.log(`\nTotal: ${results.length} | ✓ Passed: ${passed} | ✗ Failed: ${failed} | ✗ Errored: ${errored}`);

  // Density by config
  console.log("\n── Avg net density by config ──");
  for (const configName of Object.keys(CONFIGS)) {
    const cr = results.filter((r) => r.config === configName && r.success);
    if (cr.length === 0) { console.log(`  ${configName}: no successful results`); continue; }
    const avgDen = cr.reduce((s, r) => s + r.netDensityBitsPerNt, 0) / cr.length;
    const avgEnc = cr.reduce((s, r) => s + r.encodeMs, 0) / cr.length;
    const avgDec = cr.reduce((s, r) => s + r.decodeMs, 0) / cr.length;
    console.log(`  ${configName}: avg ${avgDen.toFixed(3)} b/nt | enc ${avgEnc.toFixed(0)}ms | dec ${avgDec.toFixed(0)}ms`);
  }

  if (failed + errored > 0) {
    console.log(`\n⚠ ${failed + errored} result(s) had issues.`);
  }
}

main().catch((e) => {
  console.error("FATAL:", e);
  process.exit(1);
});
