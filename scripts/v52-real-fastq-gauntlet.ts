/**
 * v52 REAL FASTQ VALIDATION GAUNTLET
 *
 * Runs the v52 codec against REAL sequencing reads downloaded from ENA:
 *
 *   1. ERR1797975 (Erlich 2017, Illumina HiSeq, 152nt reads)
 *   2. SRR6831225 (Organick 2018 subset of SRP135605, Illumina, 230nt reads)
 *   3. DRR421226  (Takahashi 2024 Nanopore, 35nt short reads)
 *
 * For each dataset, we:
 *   a. Parse the FASTQ to extract real reads (sequences + Q-scores)
 *   b. Encode a known test payload with Helix using a primer that matches
 *      a subsequence found in the real reads
 *   c. Attempt to decode using the real reads (no synthetic noise)
 *   d. Report recovery rate
 *
 * NOTE: The real reads are from the original Erlich/Organick/Takahashi
 * experiments, which used DIFFERENT primers and DIFFERENT payload than what
 * Helix would encode. So we can't actually decode Helix-encoded data from
 * these reads. Instead, this script:
 *   - Parses the FASTQ to extract REAL error profiles (sub/ins/del rates)
 *   - Re-simulates the same number of reads at the same error profile
 *     against Helix-encoded data
 *   - Reports whether Helix can recover the data
 *
 * This gives us a "realistic error profile" validation rather than a
 * synthetic-noise validation.
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { ULTIMATE_NANOPORE_V52_CONFIG, V51_DEFAULT_CONFIG, computeDensity } from "../src/lib/dna/presets";
import * as crypto from "crypto";
import * as fs from "fs";
import * as zlib from "zlib";

const REPORT = (msg: string) => console.log(`[v52-real-fastq] ${msg}`);

// --- FASTQ parsing ---

interface FastqRead {
  readId: string;
  sequence: string;
  quality: string;  // ASCII-encoded Q-scores (Phred+33)
}

/**
 * Parse a FASTQ file (possibly gzipped) into an array of reads.
 * Stops after `maxReads` reads to bound memory.
 *
 * Handles truncated gzip files (e.g., partial downloads) by using
 * Z_SYNC_FLUSH instead of Z_FINISH.
 */
function parseFastq(path: string, maxReads: number = 5000): FastqRead[] {
  const buf = fs.readFileSync(path);
  let text: string;
  if (path.endsWith(".gz")) {
    // Use Z_SYNC_FLUSH to handle truncated gzip files (partial downloads)
    try {
      text = zlib.gunzipSync(buf, {
        finishFlush: zlib.constants.Z_SYNC_FLUSH,
      }).toString("utf8");
    } catch (e) {
      // If even sync-flush fails, try inflateRawSync with no header
      try {
        // Strip the 10-byte gzip header and try raw inflate
        text = zlib.inflateRawSync(buf.subarray(10), {
          finishFlush: zlib.constants.Z_SYNC_FLUSH,
        }).toString("utf8");
      } catch (e2) {
        // Last resort: empty
        text = "";
      }
    }
  } else {
    text = buf.toString("utf8");
  }
  const lines = text.split("\n");
  const reads: FastqRead[] = [];
  for (let i = 0; i + 3 < lines.length && reads.length < maxReads; i += 4) {
    const header = lines[i];
    if (!header.startsWith("@")) continue;
    const sequence = lines[i + 1];
    const plus = lines[i + 2];
    const quality = lines[i + 3];
    if (!plus.startsWith("+") || sequence.length !== quality.length) continue;
    reads.push({
      readId: header.slice(1).trim(),
      sequence: sequence.trim(),
      quality: quality.trim(),
    });
  }
  return reads;
}

/**
 * Estimate the substitution rate from Q-scores.
 * Average Q-score → error probability: p_err = 10^(-Q/10).
 */
function estimateSubstitutionRate(reads: FastqRead[]): {
  avgQ: number;
  subRate: number;
  avgLen: number;
  readCount: number;
} {
  if (reads.length === 0) {
    return { avgQ: 0, subRate: 0, avgLen: 0, readCount: 0 };
  }
  let totalQ = 0;
  let totalLen = 0;
  for (const r of reads) {
    for (let i = 0; i < r.quality.length; i++) {
      totalQ += r.quality.charCodeAt(i) - 33;
    }
    totalLen += r.quality.length;
  }
  const avgQ = totalQ / totalLen;
  const subRate = Math.pow(10, -avgQ / 10);
  const avgLen = totalLen / reads.length;
  return { avgQ, subRate, avgLen, readCount: reads.length };
}

/**
 * Estimate indel rate by checking the read length distribution.
 * If reads have varying lengths (vs. the modal length), the variance
 * indicates indel rate.
 */
function estimateIndelRate(reads: FastqRead[]): { insRate: number; delRate: number; modeLen: number } {
  if (reads.length === 0) return { insRate: 0, delRate: 0, modeLen: 0 };
  // Find mode length
  const lenCounts = new Map<number, number>();
  for (const r of reads) {
    lenCounts.set(r.sequence.length, (lenCounts.get(r.sequence.length) ?? 0) + 1);
  }
  let modeLen = 0;
  let modeCount = 0;
  for (const [len, cnt] of lenCounts) {
    if (cnt > modeCount) {
      modeCount = cnt;
      modeLen = len;
    }
  }
  // Count reads that deviate from mode length
  let shorter = 0;
  let longer = 0;
  let total = 0;
  for (const r of reads) {
    if (r.sequence.length < modeLen) shorter++;
    else if (r.sequence.length > modeLen) longer++;
    total++;
  }
  // Rough estimate: shorter = deletions, longer = insertions
  // (This is a rough heuristic — true indel estimation requires alignment)
  const delRate = shorter / total / modeLen;
  const insRate = longer / total / modeLen;
  return { insRate, delRate, modeLen };
}

// --- Validation ---

async function validateDataset(
  name: string,
  fastqPath: string,
  payloadSize: number,
  config: any,
  coverageLevels: number[],
  useUltra: boolean = false,
): Promise<void> {
  REPORT(`\n=== ${name} ===`);
  REPORT(`FASTQ: ${fastqPath}`);

  if (!fs.existsSync(fastqPath)) {
    REPORT(`❌ FASTQ file not found, skipping`);
    return;
  }

  // Parse FASTQ
  const reads = parseFastq(fastqPath, 5000);
  REPORT(`Parsed ${reads.length} reads`);
  if (reads.length === 0) {
    REPORT(`❌ No reads parsed, skipping`);
    return;
  }

  // Estimate error profile
  const subInfo = estimateSubstitutionRate(reads);
  const indelInfo = estimateIndelRate(reads);
  REPORT(`Real error profile:`);
  REPORT(`  avg Q-score: ${subInfo.avgQ.toFixed(2)}`);
  REPORT(`  avg read length: ${subInfo.avgLen.toFixed(1)} nt (mode: ${indelInfo.modeLen})`);
  REPORT(`  estimated sub rate: ${(subInfo.subRate * 100).toFixed(3)}%`);
  REPORT(`  estimated ins rate: ${(indelInfo.insRate * 100).toFixed(3)}%`);
  REPORT(`  estimated del rate: ${(indelInfo.delRate * 100).toFixed(3)}%`);

  // Encode a known payload
  const payload = crypto.randomBytes(payloadSize);
  const tEnc = Date.now();
  const enc = await encodeFile(
    Buffer.from(payload),
    config,
    { fileName: `${name}.bin`, contentType: "application/octet-stream" },
  );
  const encMs = Date.now() - tEnc;
  const density = computeDensity(config, "total");
  REPORT(`Encoded ${payloadSize}B → ${enc.encoded.oligos.length} oligos in ${encMs}ms`);
  REPORT(`Density: ${density.toFixed(3)} bits/nt (total-oligo)`);

  // Simulate reads at the REAL error profile
  const realNoise = {
    substitutionRate: Math.max(subInfo.subRate, 0.0005),  // floor at 0.05%
    insertionRate: Math.max(indelInfo.insRate, 0.0001),
    deletionRate: Math.max(indelInfo.delRate, 0.0001),
    dropoutRate: 0.02,
    coverage: coverageLevels[0],
    seed: 42,
  };

  REPORT(`\nCov | Reads  | Recovery | Time(ms) | OligosOK | Status`);
  REPORT(`----|--------|----------|----------|----------|--------`);

  for (const cov of coverageLevels) {
    const sim = simulate(enc.encoded.oligos, { ...realNoise, coverage: cov, seed: 42 });
    const t0 = Date.now();
    let dec;
    try {
      if (useUltra) {
        dec = await decodeReadsUltra(
          sim.reads,
          enc.encoded.metadata,
          config,
          enc.encoded.forwardPrimer,
          enc.encoded.reversePrimer,
        );
      } else {
        const jsDec = await decodeReads(
          sim.reads,
          enc.encoded.metadata,
          config,
          enc.encoded.forwardPrimer,
          enc.encoded.reversePrimer,
          true,
        );
        dec = {
          data: jsDec.data,
          hash: jsDec.hash,
          hashMatches: jsDec.hashMatches,
          stats: jsDec.stats,
          perOligo: jsDec.perOligo,
        };
      }
    } catch (e: any) {
      REPORT(`  ${cov}x | ${sim.totalReads.toString().padStart(6)} | ---      | ---      | ---      | [decode failed: ${e.message?.slice(0, 60)}]`);
      continue;
    }
    const ms = Date.now() - t0;
    const pass = dec.hashMatches;
    const pct = ((dec.stats.oligosRecovered / enc.encoded.oligos.length) * 100).toFixed(2);
    REPORT(`  ${cov}x | ${sim.totalReads.toString().padStart(6)} | ${pct.padStart(6)}% | ${ms.toString().padStart(8)} | ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length} | ${pass ? "✅ PASS" : "❌ FAIL"}`);
  }
}

async function main() {
  REPORT("╔══════════════════════════════════════════════════════════════╗");
  REPORT("║  v52 REAL FASTQ VALIDATION GAUNTLET                         ║");
  REPORT("║  Erlich / Organick / Takahashi — REAL error profiles        ║");
  REPORT("╚══════════════════════════════════════════════════════════════╝");

  const FASTQ_DIR = __dirname + "/../fastq-data";

  // === Dataset 1: Erlich 2017 (Illumina, ERR1797975) ===
  // Real Illumina HiSeq reads, ~152nt, low error rate
  await validateDataset(
    "Erlich 2017 (ERR1797975, Illumina)",
    `${FASTQ_DIR}/ERR1797975_1.partial.fastq.gz`,
    64 * 1024,  // 64KB payload
    {
      ...V51_DEFAULT_CONFIG,
      oligoLength: 300,
      mappingMode: "direct" as const,
      innerCode: "ldpc" as const,
      innerParityBytes: 8,
      outerParityRatio: 0.15,
      lowCoverageTrigger: 5,
      channel: "illumina" as const,
    },
    [3, 5, 10],
    true,
  );

  // === Dataset 2: Organick 2018 (Illumina, SRR6831225 from SRP135605) ===
  // Real Illumina reads, ~230nt
  await validateDataset(
    "Organick 2018 (SRR6831225 ⊂ SRP135605, Illumina)",
    `${FASTQ_DIR}/SRR6831225_1.partial.fastq.gz`,
    64 * 1024,  // 64KB payload (subset)
    {
      ...V51_DEFAULT_CONFIG,
      oligoLength: 300,
      mappingMode: "direct" as const,
      innerCode: "ldpc" as const,
      innerParityBytes: 8,
      outerParityRatio: 0.20,
      lowCoverageTrigger: 5,
      channel: "illumina" as const,
    },
    [5, 10],
    true,
  );

  // === Dataset 3: Takahashi 2024 (Nanopore, DRR421226) ===
  // Real Nanopore reads, short (35nt) — high indel rate
  await validateDataset(
    "Takahashi 2024 (DRR421226, Nanopore)",
    `${FASTQ_DIR}/DRR421226.partial.fastq.gz`,
    8 * 1024,  // 8KB payload (v52 conv inner has lower density)
    {
      ...ULTIMATE_NANOPORE_V52_CONFIG,
      // Override for v52 full pipeline
      useConvolutionalInner: true,
      channel: "nanopore" as const,
      lowCoverageTrigger: 5,
    },
    [15, 25],
    false,  // JS path (v52 conv inner requires JS)
  );

  REPORT("\n╔══════════════════════════════════════════════════════════════╗");
  REPORT("║  REAL FASTQ VALIDATION COMPLETE                             ║");
  REPORT("╚══════════════════════════════════════════════════════════════╝");
}

main().catch((e) => {
  console.error("FATAL:", e);
  process.exit(1);
});
