/**
 * Real-Dataset Validation Suite for helix-codec.
 *
 * Tests encoding → clean decode roundtrip against REAL biological datasets
 * (no synthetic/fake data). Loads FASTA files, encodes to DNA oligos,
 * simulates reads, decodes, and verifies lossless recovery.
 *
 * Tiers:
 *   Small  — SARS-CoV-2 (~30KB), UniProt protein (~1KB)     [seconds]
 *   Medium — E. coli K-12 (~4.6MB), Yeast (~12MB)           [minutes]
 *   Large  — Human chr21 (~47MB)                             [tens of minutes]
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import {
  V51_DEFAULT_CONFIG,
  ULTIMATE_V55_DENSITY_CONFIG,
  ULTIMATE_V63_HD_CONFIG,
  ULTIMATE_V59_FAST_ENCODE_CONFIG,
} from "../src/lib/dna/presets";
import type { CodecConfig } from "../src/lib/dna/types";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// ─── Helpers ────────────────────────────────────────────────────────────────

/** Extract sequence data from a FASTA file (strip headers, join lines). */
function parseFasta(content: string): string {
  const lines = content.split("\n");
  const seqLines: string[] = [];
  for (const line of lines) {
    if (line.startsWith(">")) continue;
    if (line.trim() === "") continue;
    seqLines.push(line.trim());
  }
  return seqLines.join("");
}

/** Convert a string to UTF-8 bytes. */
function seqToBytes(seq: string): Uint8Array {
  return new Uint8Array(Buffer.from(seq, "utf-8"));
}

/** Read a FASTA file and return the sequence as bytes. */
function loadFasta(filePath: string): { name: string; bytes: Uint8Array; seqLen: number } {
  const content = fs.readFileSync(filePath, "utf-8");
  const seq = parseFasta(content);
  const bytes = seqToBytes(seq);
  const name = path.basename(filePath);
  return { name, bytes, seqLen: seq.length };
}

// ─── Test Presets ───────────────────────────────────────────────────────────

const ALL_PRESETS: { name: string; config: CodecConfig }[] = [
  { name: "v51-default (300nt)", config: V51_DEFAULT_CONFIG },
  { name: "v59-fast (300nt SRT)", config: ULTIMATE_V59_FAST_ENCODE_CONFIG },
  { name: "v55-density (700nt)", config: ULTIMATE_V55_DENSITY_CONFIG },
  { name: "v63-hd (1100nt)", config: ULTIMATE_V63_HD_CONFIG },
];

// ─── Test Result ────────────────────────────────────────────────────────────

interface TestResult {
  dataset: string;
  preset: string;
  dataSize: number;
  oligoCount: number;
  netDensity: number;
  encodeMs: number;
  decodeMs: number;
  screeningRetries: number;
  allConstraintsSatisfied: boolean;
  roundtripOk: boolean;
  hashOk: boolean;
  gcRange: { min: number; max: number };
  maxHomopolymer: number;
  error?: string;
}

// ─── Single Dataset Test ────────────────────────────────────────────────────

async function testDatasetWithPreset(
  datasetName: string,
  data: Uint8Array,
  preset: { name: string; config: CodecConfig },
): Promise<TestResult> {
  const cfg = { ...preset.config };

  // Encode
  const encodeResult = await encodeFile(data, cfg, {
    fileName: datasetName,
    contentType: "application/octet-stream",
  });

  const encoded = encodeResult.encoded;
  const stats = encodeResult.stats;

  // Validate all oligos satisfy constraints
  let allOk = true;
  let gcMin = 1, gcMax = 0;
  let maxHp = 0;
  for (const oligo of encoded.oligos) {
    const gc = oligo.gc;
    if (gc < gcMin) gcMin = gc;
    if (gc > gcMax) gcMax = gc;
    if (oligo.maxHomopolymer > maxHp) maxHp = oligo.maxHomopolymer;
    if (gc < cfg.constraints.gcMin || gc > cfg.constraints.gcMax) {
      allOk = false;
    }
    if (oligo.maxHomopolymer > cfg.constraints.maxHomopolymer) {
      allOk = false;
    }
  }

  // Simulate clean reads (30× coverage, no errors)
  const simResult = simulate(encoded.oligos, {
    ...PRESET_CLEAN,
    coverage: 30,
    simulator: "basic", // use basic for speed
  });

  // Decode
  const decodeStart = Date.now();
  const decodeResult = await decodeReads(
    simResult.reads,
    encoded.metadata,
    cfg,
    encoded.forwardPrimer,
    encoded.reversePrimer,
  );
  const decodeMs = Date.now() - decodeStart;

  // Verify roundtrip
  const decoded = decodeResult.data;
  let roundtripOk = false;
  if (decoded && decoded.length === data.length) {
    roundtripOk = true;
    for (let i = 0; i < data.length; i++) {
      if (decoded[i] !== data[i]) {
        roundtripOk = false;
        break;
      }
    }
  }

  // Verify hash
  const hashOk = decodeResult.hashMatches;

  return {
    dataset: datasetName,
    preset: preset.name,
    dataSize: data.length,
    oligoCount: stats.oligoCount,
    netDensity: stats.netDensityBitsPerNt,
    encodeMs: stats.encodeTimeMs,
    decodeMs,
    screeningRetries: stats.screeningRetries,
    allConstraintsSatisfied: allOk,
    roundtripOk,
    hashOk,
    gcRange: { min: gcMin, max: gcMax },
    maxHomopolymer: maxHp,
  };
}

// ─── Run Tier ───────────────────────────────────────────────────────────────

async function runTier(
  tierName: string,
  datasets: { name: string; bytes: Uint8Array }[],
  presets: { name: string; config: CodecConfig }[],
): Promise<{ results: TestResult[]; allPassed: boolean }> {
  console.log(`\n${"=".repeat(70)}`);
  console.log(`  ${tierName}`);
  console.log(`${"=".repeat(70)}`);

  const results: TestResult[] = [];
  let allPassed = true;

  for (const ds of datasets) {
    console.log(`\n  DATASET: ${ds.name} (${(ds.bytes.length / 1024).toFixed(1)} KB, ${ds.bytes.length} bytes)`);

    for (const preset of presets) {
      const label = `${preset.name}`;
      process.stdout.write(`  [RUN] ${label}: encoding...`);
      try {
        const result = await testDatasetWithPreset(ds.name, ds.bytes, preset);
        results.push(result);

        const ok = result.allConstraintsSatisfied && result.roundtripOk && result.hashOk;
        if (!ok) allPassed = false;

        const status = ok ? "PASS" : "FAIL";
        console.log(`\r  [${status}] ${label}: ` +
          `${result.oligoCount} oligos, ` +
          `density=${result.netDensity.toFixed(3)} b/nt, ` +
          `enc=${result.encodeMs}ms, dec=${result.decodeMs}ms, ` +
          `retries=${result.screeningRetries}, ` +
          `GC=[${result.gcRange.min.toFixed(2)},${result.gcRange.max.toFixed(2)}], ` +
          `maxHp=${result.maxHomopolymer}, ` +
          `roundtrip=${result.roundtripOk ? "OK" : "FAIL"}, ` +
          `hash=${result.hashOk ? "OK" : "FAIL"}`);

        if (!result.allConstraintsSatisfied) {
          console.log(`    WARNING: CONSTRAINT VIOLATION — GC or homopolymer out of range`);
        }
        if (!result.roundtripOk) {
          console.log(`    WARNING: ROUNDTRIP MISMATCH — decoded data differs from original`);
        }
        if (!result.hashOk) {
          console.log(`    WARNING: HASH MISMATCH — SHA-256 mismatch`);
        }
      } catch (e: any) {
        console.log(`\r  [FAIL] ${label}: ERROR — ${e.message?.slice(0, 200)}`);
        results.push({
          dataset: ds.name,
          preset: preset.name,
          dataSize: ds.bytes.length,
          oligoCount: 0,
          netDensity: 0,
          encodeMs: 0,
          decodeMs: 0,
          screeningRetries: 0,
          allConstraintsSatisfied: false,
          roundtripOk: false,
          hashOk: false,
          gcRange: { min: 0, max: 0 },
          maxHomopolymer: 0,
          error: e.message?.slice(0, 200),
        });
        allPassed = false;
      }
    }
  }

  return { results, allPassed };
}

// ─── Main ───────────────────────────────────────────────────────────────────

async function main() {
  console.log("=".repeat(70));
  console.log("  helix-codec Real-Dataset Validation Suite");
  console.log("  Testing against REAL biological data — no synthetic/fake data");
  console.log("=".repeat(70));

  const datasetsDir = path.resolve(__dirname, "..", "datasets");
  const allResults: TestResult[] = [];
  let overallPass = true;

  // ─── Small Tier ────────────────────────────────────────────────────────
  const smallDir = path.join(datasetsDir, "small");
  const smallDatasets: { name: string; bytes: Uint8Array }[] = [];

  const sarsFile = path.join(smallDir, "sars-cov-2.fa");
  if (fs.existsSync(sarsFile)) {
    const { name, bytes } = loadFasta(sarsFile);
    smallDatasets.push({ name: `SARS-CoV-2 (${name})`, bytes });
    console.log(`  Loaded SARS-CoV-2: ${bytes.length} bytes`);
  } else {
    console.log(`  WARNING: SARS-CoV-2 FASTA not found at ${sarsFile}`);
  }

  const uniprotFile = path.join(smallDir, "uniprot-p00533.fa");
  if (fs.existsSync(uniprotFile)) {
    const { name, bytes } = loadFasta(uniprotFile);
    smallDatasets.push({ name: `UniProt EGFR (${name})`, bytes });
    console.log(`  Loaded UniProt EGFR: ${bytes.length} bytes`);
  }

  if (smallDatasets.length > 0) {
    const { results, allPassed } = await runTier(
      "SMALL — Sanity/Correctness Tests",
      smallDatasets,
      ALL_PRESETS,
    );
    allResults.push(...results);
    if (!allPassed) overallPass = false;
  }

  // ─── Medium Tier ───────────────────────────────────────────────────────
  const mediumDir = path.join(datasetsDir, "medium");
  const mediumDatasets: { name: string; bytes: Uint8Array }[] = [];

  const ecoliFile = path.join(mediumDir, "ecoli-k12.fa");
  if (fs.existsSync(ecoliFile)) {
    const { name, bytes } = loadFasta(ecoliFile);
    mediumDatasets.push({ name: `E.coli K-12 (${name})`, bytes });
    console.log(`  Loaded E.coli K-12: ${bytes.length} bytes`);
  }

  const yeastFile = path.join(mediumDir, "yeast.fa");
  if (fs.existsSync(yeastFile)) {
    const { name, bytes } = loadFasta(yeastFile);
    mediumDatasets.push({ name: `Yeast S288C (${name})`, bytes });
    console.log(`  Loaded Yeast S288C: ${bytes.length} bytes`);
  }

  if (mediumDatasets.length > 0) {
    const mediumPresets = ALL_PRESETS.filter(p =>
      p.name.includes("v51") || p.name.includes("v59")
    );
    const { results, allPassed } = await runTier(
      "MEDIUM — Realistic Single-Genome Tests",
      mediumDatasets,
      mediumPresets,
    );
    allResults.push(...results);
    if (!allPassed) overallPass = false;
  }

  // ─── Large Tier ────────────────────────────────────────────────────────
  const largeDir = path.join(datasetsDir, "large");
  const largeDatasets: { name: string; bytes: Uint8Array }[] = [];

  const chr21File = path.join(largeDir, "chr21.fa");
  if (fs.existsSync(chr21File)) {
    const { name, bytes } = loadFasta(chr21File);
    largeDatasets.push({ name: `Human chr21 (${name})`, bytes });
    console.log(`  Loaded Human chr21: ${bytes.length} bytes`);
  }

  if (largeDatasets.length > 0) {
    const largePresets = ALL_PRESETS.filter(p => p.name.includes("v51"));
    const { results, allPassed } = await runTier(
      "LARGE — Stress Test Scaling",
      largeDatasets,
      largePresets,
    );
    allResults.push(...results);
    if (!allPassed) overallPass = false;
  }

  // ─── Summary ───────────────────────────────────────────────────────────
  console.log(`\n${"=".repeat(70)}`);
  console.log("  SUMMARY");
  console.log(`${"=".repeat(70)}`);

  const byDataset = new Map<string, TestResult[]>();
  for (const r of allResults) {
    if (!byDataset.has(r.dataset)) byDataset.set(r.dataset, []);
    byDataset.get(r.dataset)!.push(r);
  }

  let passCount = 0;
  let failCount = 0;
  for (const [ds, results] of byDataset) {
    console.log(`\n  ${ds}`);
    for (const r of results) {
      const ok = r.allConstraintsSatisfied && r.roundtripOk && r.hashOk;
      if (ok) passCount++; else failCount++;
      const status = ok ? "PASS" : "FAIL";
      console.log(`    [${status}] ${r.preset}: density=${r.netDensity.toFixed(3)} b/nt, ` +
        `enc=${r.encodeMs}ms, dec=${r.decodeMs}ms, ` +
        `roundtrip=${r.roundtripOk}, hash=${r.hashOk}` +
        (r.error ? `, error=${r.error}` : ""));
    }
  }

  console.log(`\n  Total: ${passCount} passed, ${failCount} failed, ${passCount + failCount} tests`);

  // Save results as JSON
  const resultsPath = path.join(datasetsDir, "real-dataset-results.json");
  fs.writeFileSync(resultsPath, JSON.stringify(allResults, null, 2));
  console.log(`  Results saved to ${resultsPath}`);

  const overall = failCount === 0;
  console.log(`\n  ${overall ? "ALL TESTS PASSED" : "SOME TESTS FAILED"}`);
  if (!overall) process.exit(1);
}

main().catch((e) => {
  console.error("Fatal error:", e);
  process.exit(1);
});
