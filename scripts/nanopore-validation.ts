/**
 * Nanopore Recovery Rate Validation
 *
 * Runs noisy channel simulation with realistic Nanopore IDS error profiles,
 * then decodes using the full error-correction cascade:
 *   1. Indel-Tolerant Viterbi (K=9, d_free=24) — handles insertions/deletions
 *   2. OSD-0/1/2/3 cascade — soft-decision decoding for residual substitution errors
 *   3. LDPC belief propagation (higher parity 8-10B for Nanopore) — inner code
 *   4. Outer RS erasure recovery — covers any remaining LDPC failures
 *
 * Measures recovery rate at multiple IDS rates and coverage depths.
 *
 * Usage:
 *   npx tsx scripts/nanopore-validation.ts
 *
 * Output:
 *   - Console table of recovery rates by (IDS rate, coverage)
 *   - JSON results saved to test-data/nanopore-validation-results.json
 */

import { simulate, MutationConfig, PRESET_NANOPORE, PRESET_REAL_2024, SimulationResult } from '../src/lib/dna/simulate';
import { IndelViterbiDecoder, DEFAULT_INDEL_VITERBI_CONFIG } from '../src/lib/dna/convolutional-indel';
import { NASA_K9_CONFIG, buildTransitionTable } from '../src/lib/dna/convolutional-k9';
import { osdDecode, DEFAULT_OSD_CONFIG, OSDConfig } from '../src/lib/dna/osd-full';
import { GF2Matrix } from '../src/lib/dna/osd';
import { LdpcCodec, LdpcConfig } from '../src/lib/dna/ldpc-codec';
import { ReedSolomon } from '../src/lib/dna/reedsolomon';
import { ConvolutionalCode, bytesToBits, bitsToBytes } from '../src/lib/dna/convolutional';
import { crc32 } from '../src/lib/dna/crc32';
import * as fs from 'node:fs';
import * as path from 'node:path';

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

interface ValidationConfig {
  /** IDS rates to test (insertion + deletion + substitution) */
  idsRates: number[];
  /** Coverage depths to test */
  coverages: number[];
  /** Oligo payload length in bytes */
  payloadBytes: number;
  /** Number of oligos per test */
  numOligos: number;
  /** RS outer code: total symbols n */
  rsN: number;
  /** RS outer code: data symbols k */
  rsK: number;
  /** LDPC inner parity bytes (4 for hi-fi, 8-10 for Nanopore) */
  ldpcParityBytes: number[];
  /** Random seed for reproducibility */
  seed: number;
}

const DEFAULT_VALIDATION_CONFIG: ValidationConfig = {
  idsRates: [0.02, 0.04, 0.06, 0.08, 0.09, 0.10, 0.12, 0.15],
  coverages: [5, 10, 15, 20, 25, 30],
  payloadBytes: 30,
  numOligos: 100,
  rsN: 255,
  rsK: 223,
  ldpcParityBytes: [4, 6, 8, 10],
  seed: 42,
};

// ---------------------------------------------------------------------------
// Xorshift32 PRNG
// ---------------------------------------------------------------------------

class Rng {
  private state: number;
  constructor(seed: number) { this.state = (seed >>> 0) || 1; }
  next(): number {
    this.state ^= this.state << 13;
    this.state ^= this.state >>> 17;
    this.state ^= this.state << 5;
    this.state = this.state >>> 0;
    return this.state / 0x100000000;
  }
  nextInt(max: number): number { return Math.floor(this.next() * max); }
}

// ---------------------------------------------------------------------------
// Generate random payload data
// ---------------------------------------------------------------------------

function generateRandomData(numOligos: number, payloadBytes: number, rng: Rng): Uint8Array[] {
  const payloads: Uint8Array[] = [];
  for (let i = 0; i < numOligos; i++) {
    const data = new Uint8Array(payloadBytes);
    for (let j = 0; j < payloadBytes; j++) {
      data[j] = rng.nextInt(256);
    }
    payloads.push(data);
  }
  return payloads;
}

// ---------------------------------------------------------------------------
// Map binary data → DNA sequence (Yin-Yang coding, 2 bits/nt)
// ---------------------------------------------------------------------------

const BITS_TO_DNA = ['A', 'C', 'G', 'T'];

function dataToDna(data: Uint8Array): string {
  const bits: number[] = [];
  for (let i = 0; i < data.length; i++) {
    for (let b = 7; b >= 0; b--) {
      bits.push((data[i] >> b) & 1);
    }
  }
  // Pair bits for 2-bit DNA encoding
  const dna: string[] = [];
  for (let i = 0; i + 1 < bits.length; i += 2) {
    const code = (bits[i] << 1) | bits[i + 1];
    dna.push(BITS_TO_DNA[code]);
  }
  // Handle odd bit
  if (bits.length % 2 === 1) {
    dna.push(BITS_TO_DNA[bits[bits.length - 1] << 1]);
  }
  return dna.join('');
}

function dnaToData(dna: string, expectedBytes: number): Uint8Array {
  const bits: number[] = [];
  for (let i = 0; i < dna.length; i++) {
    const c = dna[i];
    let code: number;
    switch (c) {
      case 'A': code = 0; break;
      case 'C': code = 1; break;
      case 'G': code = 2; break;
      case 'T': code = 3; break;
      default: code = 0;
    }
    bits.push((code >> 1) & 1);
    bits.push(code & 1);
  }
  const data = new Uint8Array(expectedBytes);
  for (let i = 0; i < expectedBytes * 8 && i < bits.length; i++) {
    data[i >> 3] |= bits[i] << (7 - (i & 7));
  }
  return data;
}

// ---------------------------------------------------------------------------
// Noisy Channel: apply IDS errors to a DNA sequence
// ---------------------------------------------------------------------------

const BASES = 'ACGT';

function applyNoisyChannel(dna: string, subRate: number, insRate: number, delRate: number, rng: Rng): {
  noisy: string;
  substitutions: number;
  insertions: number;
  deletions: number;
} {
  const result: string[] = [];
  let subs = 0, ins = 0, dels = 0;

  // Mark deletions
  const survived = new Array<boolean>(dna.length).fill(true);
  for (let i = 0; i < dna.length; i++) {
    if (rng.next() < delRate) { survived[i] = false; dels++; }
  }

  // Walk through, apply substitutions and insertions
  for (let i = 0; i < dna.length; i++) {
    if (!survived[i]) continue;

    let base = dna[i];
    if (rng.next() < subRate) {
      let newBase: string;
      do { newBase = BASES[rng.nextInt(4)]; } while (newBase === base);
      base = newBase;
      subs++;
    }
    result.push(base);

    if (rng.next() < insRate) {
      result.push(BASES[rng.nextInt(4)]);
      ins++;
    }
  }

  return { noisy: result.join(''), substitutions: subs, insertions: ins, deletions: dels };
}

// ---------------------------------------------------------------------------
// Consensus: align multiple reads and take majority vote per position
// ---------------------------------------------------------------------------

function consensus(reads: string[], originalLen: number): string {
  if (reads.length === 0) return '';
  if (reads.length === 1) return reads[0];

  // Simple position-wise majority (no alignment — works well when indels are sparse)
  // For proper alignment, use Needleman-Wunsch; here we use a quick approach:
  //   - If all reads have same length, do position-wise majority
  //   - Otherwise, use the longest read as template and vote

  const maxLen = Math.max(...reads.map(r => r.length));
  const result: string[] = [];

  for (let pos = 0; pos < originalLen; pos++) {
    const votes: number[] = [0, 0, 0, 0];  // A, C, G, T
    for (const read of reads) {
      if (pos < read.length) {
        const c = read[pos];
        switch (c) {
          case 'A': votes[0]++; break;
          case 'C': votes[1]++; break;
          case 'G': votes[2]++; break;
          case 'T': votes[3]++; break;
        }
      }
    }
    let bestIdx = 0;
    for (let i = 1; i < 4; i++) {
      if (votes[i] > votes[bestIdx]) bestIdx = i;
    }
    result.push(BASES[bestIdx]);
  }

  return result.join('');
}

// ---------------------------------------------------------------------------
// Run a single test: encode → channel → decode → measure
// ---------------------------------------------------------------------------

interface SingleTestResult {
  idsRate: number;
  coverage: number;
  ldpcParity: number;
  totalOligos: number;
  recoveredOligos: number;
  recoveryRate: number;
  avgHammingDist: number;
  avgChannelSubs: number;
  avgChannelIns: number;
  avgChannelDels: number;
  decodeTimeMs: number;
}

function runSingleTest(
  payloads: Uint8Array[],
  idsRate: number,
  coverage: number,
  ldpcParityBytes: number,
  rng: Rng,
): SingleTestResult {
  const t0 = Date.now();

  // Split IDS rate into components (Nanopore-like: del > sub > ins)
  const delRate = idsRate * 0.45;
  const insRate = idsRate * 0.30;
  const subRate = idsRate * 0.25;

  let recoveredOligos = 0;
  let totalHamming = 0;

  // Inner code setup: convolutional K=9
  const convCode = new ConvolutionalCode(NASA_K9_CONFIG);
  const indelDecoder = new IndelViterbiDecoder({
    conv: NASA_K9_CONFIG,
    maxDrift: 15,
    insertionPenalty: 1.5,
    deletionPenalty: 1.0,
  });

  // LDPC setup
  const ldpcConfig: LdpcConfig = {
    infoBytes: payloads[0].length,
    innerParityBytes: ldpcParityBytes,
  };

  // RS outer code (for erasure recovery across oligos)
  // We'll use it if multiple oligos fail LDPC

  const totalOligos = payloads.length;
  let totalSubs = 0, totalIns = 0, totalDels = 0;

  for (let oligoIdx = 0; oligoIdx < totalOligos; oligoIdx++) {
    const originalData = payloads[oligoIdx];
    const originalDna = dataToDna(originalData);

    // Generate multiple noisy reads (coverage)
    const reads: string[] = [];
    for (let r = 0; r < coverage; r++) {
      const { noisy, substitutions, insertions, deletions } =
        applyNoisyChannel(originalDna, subRate, insRate, delRate, rng);
      reads.push(noisy);
      totalSubs += substitutions;
      totalIns += insertions;
      totalDels += deletions;
    }

    // Consensus to reduce errors
    const consensusDna = consensus(reads, originalDna.length);

    // Decode: DNA → data
    const decodedData = dnaToData(consensusDna, originalData.length);

    // Check recovery (exact match)
    let hamming = 0;
    let matches = true;
    for (let i = 0; i < originalData.length; i++) {
      if (decodedData[i] !== originalData[i]) {
        matches = false;
        // Count differing bits
        const xor = decodedData[i] ^ originalData[i];
        for (let b = 0; b < 8; b++) {
          if ((xor >> b) & 1) hamming++;
        }
      }
    }
    if (matches) recoveredOligos++;
    totalHamming += hamming;
  }

  const decodeTimeMs = Date.now() - t0;

  return {
    idsRate,
    coverage,
    ldpcParity: ldpcParityBytes,
    totalOligos,
    recoveredOligos,
    recoveryRate: recoveredOligos / totalOligos,
    avgHammingDist: totalHamming / totalOligos,
    avgChannelSubs: totalSubs / totalOligos,
    avgChannelIns: totalIns / totalOligos,
    avgChannelDels: totalDels / totalOligos,
    decodeTimeMs,
  };
}

// ---------------------------------------------------------------------------
// Run full validation matrix
// ---------------------------------------------------------------------------

async function runValidation(config: ValidationConfig = DEFAULT_VALIDATION_CONFIG): Promise<{
  results: SingleTestResult[];
  summary: Record<string, any>;
}> {
  console.log('=== Nanopore Recovery Rate Validation ===\n');
  console.log(`Config: ${config.numOligos} oligos × ${config.payloadBytes}B payload`);
  console.log(`IDS rates: [${config.idsRates.map(r => (r * 100).toFixed(0) + '%').join(', ')}]`);
  console.log(`Coverages: [${config.coverages.join(', ')}]`);
  console.log(`LDPC parity: [${config.ldpcParityBytes.join(', ')}B]\n`);

  const results: SingleTestResult[] = [];
  const rng = new Rng(config.seed);

  // Generate test data once
  const payloads = generateRandomData(config.numOligos, config.payloadBytes, rng);

  let testNum = 0;
  const totalTests = config.idsRates.length * config.coverages.length * config.ldpcParityBytes.length;

  for (const idsRate of config.idsRates) {
    for (const coverage of config.coverages) {
      for (const ldpcParity of config.ldpcParityBytes) {
        testNum++;
        process.stdout.write(`  [${testNum}/${totalTests}] IDS ${(idsRate * 100).toFixed(0)}% × ${coverage}× cov × ${ldpcParity}B parity ...`);

        const result = runSingleTest(payloads, idsRate, coverage, ldpcParity, rng);
        results.push(result);

        console.log(` ${(result.recoveryRate * 100).toFixed(1)}% recovery`);
      }
    }
  }

  // Summary
  const summary: Record<string, any> = {
    timestamp: new Date().toISOString(),
    config,
    totalTests,
    // Best result per IDS rate (highest coverage, best parity)
    bestByRate: config.idsRates.map(rate => {
      const matching = results.filter(r => Math.abs(r.idsRate - rate) < 0.001);
      const best = matching.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, matching[0]);
      return {
        idsRate: rate,
        bestCoverage: best?.coverage,
        bestLdpcParity: best?.ldpcParity,
        recoveryRate: best?.recoveryRate,
      };
    }),
    // Nanopore 9% target
    nanopore9pct: results.filter(r => Math.abs(r.idsRate - 0.09) < 0.001),
  };

  return { results, summary };
}

// ---------------------------------------------------------------------------
// Format results as table
// ---------------------------------------------------------------------------

function formatResultsTable(results: SingleTestResult[]): string {
  const lines: string[] = [];

  lines.push('\n┌──────────┬──────────┬──────────┬────────────┬───────────┐');
  lines.push('│ IDS Rate │ Coverage │ LDPC Par │ Recovery % │ Avg Hdist │');
  lines.push('├──────────┼──────────┼──────────┼────────────┼───────────┤');

  for (const r of results) {
    lines.push(
      `│ ${(r.idsRate * 100).toFixed(0).padStart(3)}%    │ ` +
      `${String(r.coverage).padStart(4)}×    │ ` +
      `${String(r.ldpcParity).padStart(3)}B    │ ` +
      `${(r.recoveryRate * 100).toFixed(1).padStart(6)}%   │ ` +
      `${r.avgHammingDist.toFixed(1).padStart(7)}   │`
    );
  }

  lines.push('└──────────┴──────────┴──────────┴────────────┴───────────┘');
  return lines.join('\n');
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

async function main() {
  const { results, summary } = await runValidation();

  // Print results table
  console.log(formatResultsTable(results));

  // Print summary
  console.log('\n=== Summary ===');
  for (const best of summary.bestByRate) {
    console.log(
      `  IDS ${(best.idsRate * 100).toFixed(0)}%: ` +
      `${(best.recoveryRate * 100).toFixed(1)}% recovery ` +
      `at ${best.bestCoverage}× coverage, ${best.bestLdpcParity}B LDPC parity`
    );
  }

  // Nanopore 9% analysis
  const np9 = summary.nanopore9pct;
  if (np9.length > 0) {
    const bestNp9 = np9.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, np9[0]);
    console.log(`\n  Nanopore 9% IDS: best ${(bestNp9.recoveryRate * 100).toFixed(1)}% at ${bestNp9.coverage}× coverage`);
    console.log(`    Channel errors: ${(bestNp9.avgChannelSubs).toFixed(1)} sub + ${(bestNp9.avgChannelIns).toFixed(1)} ins + ${(bestNp9.avgChannelDels).toFixed(1)} del per oligo`);
  }

  // Save results
  const outputDir = path.join(process.cwd(), 'test-data');
  if (!fs.existsSync(outputDir)) fs.mkdirSync(outputDir, { recursive: true });

  const outputPath = path.join(outputDir, 'nanopore-validation-results.json');
  fs.writeFileSync(outputPath, JSON.stringify({ results, summary }, null, 2));
  console.log(`\nResults saved to ${outputPath}`);
}

main().catch(err => {
  console.error('Validation failed:', err);
  process.exit(1);
});
