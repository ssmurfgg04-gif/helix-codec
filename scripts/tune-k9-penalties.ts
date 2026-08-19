/**
 * K=9 Viterbi Penalty Tuning + LDPC Parity Optimization
 *
 * Tests different (insertionPenalty, deletionPenalty) and LDPC parity byte
 * configurations against simulated Nanopore reads to find optimal settings.
 *
 * Uses the napi-rs native Viterbi addon for speed.
 *
 * Usage:
 *   npx tsx scripts/tune-k9-penalties.ts
 */
import * as path from "node:path";
import * as fs from "node:fs";

// Native Viterbi addon
const addonPath = path.resolve(process.cwd(), "rust/helix-dna-napi/target/release/libhelix_dna_napi.so");
let nativeAddon: any = null;
try {
  const mod = { exports: {} };
  (process as any).dlopen(mod, addonPath);
  nativeAddon = mod.exports;
  console.log(`[tune] Native Viterbi: ${nativeAddon.napiVersion()}`);
} catch {
  console.log(`[tune] Native Viterbi NOT available, using JS fallback`);
}

import { NASA_K9_CONFIG, buildTransitionTable } from "../src/lib/dna/convolutional-k9";
import { ConvolutionalCode, bytesToBits, bitsToBytes } from "../src/lib/dna/convolutional";
import { LDPCInnerCode, getCachedLDPCInner } from "../src/lib/dna/ldpc-codec";
import { ReedSolomon } from "../src/lib/dna/reedsolomon";
import { crc16Bytes } from "../src/lib/dna/crc16";

// PRNG
class Rng {
  private s: number;
  constructor(seed: number) { this.s = (seed >>> 0) || 1; }
  next(): number {
    this.s ^= this.s << 13; this.s ^= this.s >>> 17;
    this.s ^= this.s << 5; this.s = this.s >>> 0;
    return this.s / 0x100000000;
  }
  nextInt(max: number): number { return Math.floor(this.next() * max); }
}

const BASES = "ACGT";
const DNA_TO_BITS: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

// DNA <-> Bytes
function dnaToBytes(dna: string): Uint8Array {
  const bits: number[] = [];
  for (const c of dna) {
    const code = DNA_TO_BITS[c] ?? 0;
    bits.push((code >> 1) & 1);
    bits.push(code & 1);
  }
  const bytes = new Uint8Array(Math.floor(bits.length / 8));
  for (let b = 0; b < bytes.length * 8 && b < bits.length; b++) {
    bytes[b >> 3] |= bits[b] << (7 - (b & 7));
  }
  return bytes;
}

function bytesToDna(data: Uint8Array): string {
  const dna: string[] = [];
  for (const byte of data) {
    for (let bit = 7; bit >= 1; bit -= 2) {
      const code = ((byte >> bit) & 1) << 1 | ((byte >> (bit - 1)) & 1);
      dna.push(BASES[code]);
    }
  }
  return dna.join("");
}

// Nanopore channel
function applyNanoporeChannel(
  dna: string, subRate: number, insRate: number, delRate: number, rng: Rng,
): { noisy: string; subs: number; ins: number; dels: number } {
  const result: string[] = [];
  let subs = 0, ins = 0, dels = 0;
  for (let i = 0; i < dna.length; i++) {
    if (rng.next() < delRate) { dels++; continue; }
    let base = dna[i];
    if (rng.next() < subRate) {
      let nb: string;
      do { nb = BASES[rng.nextInt(4)]; } while (nb === base);
      base = nb; subs++;
    }
    result.push(base);
    if (rng.next() < insRate) {
      result.push(BASES[rng.nextInt(4)]);
      ins++;
    }
  }
  return { noisy: result.join(""), subs, ins, dels };
}

// Simple plurality consensus
function simpleConsensus(reads: string[], targetLen: number): string {
  if (reads.length === 0) return "";
  if (reads.length === 1) return reads[0].slice(0, targetLen);
  const result: string[] = [];
  for (let pos = 0; pos < targetLen; pos++) {
    const votes = [0, 0, 0, 0];
    for (const r of reads) {
      if (pos < r.length) {
        const idx = DNA_TO_BITS[r[pos]] ?? -1;
        if (idx >= 0) votes[idx]++;
      }
    }
    let best = 0;
    for (let i = 1; i < 4; i++) if (votes[i] > votes[best]) best = i;
    result.push(BASES[best]);
  }
  return result.join("");
}

interface TuneResult {
  insPen: number;
  delPen: number;
  maxDrift: number;
  ldpcParity: number;
  idsRate: number;
  coverage: number;
  recovered: number;
  totalOligos: number;
  recoveryRate: number;
  avgViterbiMs: number;
  avgDecodeMs: number;
}

function runTuneCombo(
  numOligos: number,
  payloadBytes: number,
  idsRate: number,
  coverage: number,
  ldpcParityBytes: number,
  insPen: number,
  delPen: number,
  maxDrift: number,
  rng: Rng,
): TuneResult {
  const delRate = idsRate * 0.45;
  const insRate = idsRate * 0.30;
  const subRate = idsRate * 0.25;
  const innerDataBytes = payloadBytes + ldpcParityBytes + 2; // +CRC-16

  // LDPC
  let ldpcCode: LDPCInnerCode | null = null;
  try { ldpcCode = getCachedLDPCInner(payloadBytes + ldpcParityBytes, payloadBytes); } catch {}

  // Convolutional code for encoding
  const convCode = new ConvolutionalCode(NASA_K9_CONFIG);
  const transTable = buildTransitionTable(NASA_K9_CONFIG);

  // Encode + decode each oligo
  let recovered = 0;
  let totalViterbiMs = 0;
  let totalDecodeMs = 0;

  for (let i = 0; i < numOligos; i++) {
    // Generate random payload
    const payload = new Uint8Array(payloadBytes);
    for (let j = 0; j < payloadBytes; j++) payload[j] = rng.nextInt(256);

    // LDPC encode
    let ldpcCW: Uint8Array = payload;
    if (ldpcCode) {
      try { ldpcCW = ldpcCode.encode(payload); } catch {}
    }

    // CRC-16
    const withCrc = new Uint8Array(ldpcCW.length + 2);
    withCrc.set(ldpcCW, 0);
    const crc = crc16Bytes(ldpcCW);
    withCrc[ldpcCW.length] = crc[0];
    withCrc[ldpcCW.length + 1] = crc[1];

    // Convolutional encode
    let convOut: Uint8Array;
    if (nativeAddon) {
      try { convOut = new Uint8Array(nativeAddon.convK9Encode(withCrc)); } catch { convOut = withCrc; }
    } else {
      // JS fallback encode
      const infoBits = bytesToBits(withCrc);
      const encodedBits = convCode.encode(infoBits);
      convOut = new Uint8Array(Math.ceil(encodedBits.length / 8));
      for (let b = 0; b < encodedBits.length; b++) {
        convOut[b >> 3] |= encodedBits[b] << (7 - (b & 7));
      }
    }

    // Map to DNA
    const dna = bytesToDna(convOut);

    // Generate noisy reads + consensus
    const reads: string[] = [];
    for (let r = 0; r < coverage; r++) {
      const { noisy } = applyNanoporeChannel(dna, subRate, insRate, delRate, rng);
      reads.push(noisy);
    }
    const consensus = simpleConsensus(reads, dna.length);
    const consBytes = dnaToBytes(consensus);

    // Viterbi decode
    const vitStart = Date.now();
    let afterConv: Uint8Array;
    if (nativeAddon) {
      try {
        const decoded = nativeAddon.viterbiK9Decode(
          Buffer.from(consBytes),
          { maxDrift, insertionPenalty: insPen, deletionPenalty: delPen, numInfoBits: innerDataBytes * 8 }
        );
        afterConv = new Uint8Array(decoded);
      } catch {
        afterConv = consBytes.slice(0, innerDataBytes);
      }
    } else {
      afterConv = consBytes.slice(0, innerDataBytes);
    }
    totalViterbiMs += Date.now() - vitStart;

    // LDPC decode
    const decStart = Date.now();
    let decodedPayload: Uint8Array | null = null;
    if (ldpcCode && afterConv.length >= payloadBytes + ldpcParityBytes) {
      const ldpcRecv = afterConv.slice(0, payloadBytes + ldpcParityBytes);
      try {
        const { data: ldpcDec } = ldpcCode.decode(ldpcRecv);
        if (ldpcDec.length === payloadBytes) decodedPayload = ldpcDec;
      } catch {}
    }
    totalDecodeMs += Date.now() - decStart;

    // CRC fallback
    if (!decodedPayload && afterConv.length >= innerDataBytes) {
      const dataPart = afterConv.slice(0, afterConv.length - 2);
      const computedCrc = crc16Bytes(dataPart);
      if (afterConv[afterConv.length - 2] === computedCrc[0] && afterConv[afterConv.length - 1] === computedCrc[1]) {
        decodedPayload = afterConv.slice(0, payloadBytes);
      }
    }

    // Check recovery
    if (decodedPayload && decodedPayload.length === payloadBytes) {
      let match = true;
      for (let b = 0; b < payloadBytes; b++) {
        if (decodedPayload[b] !== payload[b]) { match = false; break; }
      }
      if (match) recovered++;
    }
  }

  return {
    insPen, delPen, maxDrift, ldpcParity: ldpcParityBytes, idsRate, coverage,
    recovered, totalOligos: numOligos, recoveryRate: recovered / numOligos,
    avgViterbiMs: totalViterbiMs / numOligos,
    avgDecodeMs: totalDecodeMs / numOligos,
  };
}

async function main() {
  console.log("╔════════════════════════════════════════════════════════════════════╗");
  console.log("║   K=9 Viterbi Penalty Tuning + LDPC Parity Optimization         ║");
  console.log("║   Native Viterbi: " + (nativeAddon ? "✓ ENABLED" : "✗ JS fallback") + "                                        ║");
  console.log("╚════════════════════════════════════════════════════════════════════╝\n");

  const rng = new Rng(42);
  const numOligos = 20; // Reduced for faster sweep
  const payloadBytes = 30;

  // Penalty grid to test
  const penalties: { ins: number; del: number }[] = [
    { ins: 1.0, del: 1.0 },  // symmetric, low
    { ins: 1.2, del: 1.2 },  // symmetric, moderate-low
    { ins: 1.5, del: 1.5 },  // v4.1 default: balanced
    { ins: 1.8, del: 1.8 },  // symmetric, high
    { ins: 2.0, del: 2.0 },  // symmetric, very high
    { ins: 1.5, del: 1.0 },  // asymmetric (old buggy default)
    { ins: 1.2, del: 1.5 },  // asymmetric: lower ins, higher del
    { ins: 1.8, del: 1.5 },  // asymmetric: higher ins
  ];

  const maxDrifts = [10, 15]; // Reduced for speed
  const idsRates = [0.06, 0.09, 0.12];
  const coverages = [10, 20]; // Reduced
  const ldpcParities = [4, 8, 10]; // Reduced

  // Focus on Nanopore 9% IDS first — the critical operating point
  const results: TuneResult[] = [];
  let testNum = 0;

  console.log("Phase 1: Penalty sweep at 9% IDS, 20× coverage, 8B LDPC parity\n");

  for (const pen of penalties) {
    for (const md of maxDrifts) {
      testNum++;
      const label = `ins=${pen.ins} del=${pen.del} md=${md}`;
      process.stdout.write(`  [${testNum}] ${label}: `);

      const r = runTuneCombo(numOligos, payloadBytes, 0.09, 20, 8, pen.ins, pen.del, md, rng);
      results.push(r);

      console.log(`${(r.recoveryRate * 100).toFixed(1).padStart(6)}% (${r.recovered}/${r.totalOligos}) vit=${r.avgViterbiMs.toFixed(1)}ms`);
    }
  }

  // Find best penalty config
  const best = results.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, results[0]);
  console.log(`\n  Best at 9% IDS: ins=${best.insPen} del=${best.delPen} md=${best.maxDrift} → ${(best.recoveryRate * 100).toFixed(1)}%`);

  // Phase 2: LDPC parity sweep with best penalties
  console.log("\nPhase 2: LDPC parity sweep with best penalties\n");
  const parityResults: TuneResult[] = [];

  for (const ids of idsRates) {
    for (const cov of coverages) {
      for (const par of ldpcParities) {
        testNum++;
        const label = `IDS ${(ids*100).toFixed(0)}% × ${cov}× × ${par}B`;
        process.stdout.write(`  [${testNum}] ${label}: `);

        const r = runTuneCombo(numOligos, payloadBytes, ids, cov, par, best.insPen, best.delPen, best.maxDrift, rng);
        parityResults.push(r);

        console.log(`${(r.recoveryRate * 100).toFixed(1).padStart(6)}% vit=${r.avgViterbiMs.toFixed(1)}ms`);
      }
    }
  }

  // Summary
  console.log("\n" + "=".repeat(70));
  console.log("  PENALTY TUNING RESULTS");
  console.log("=".repeat(70));

  console.log("\n  Best penalty configuration:");
  console.log(`    insertionPenalty: ${best.insPen}`);
  console.log(`    deletionPenalty:  ${best.delPen}`);
  console.log(`    maxDrift:         ${best.maxDrift}`);
  console.log(`    Recovery at 9%:   ${(best.recoveryRate * 100).toFixed(1)}%`);
  console.log(`    Viterbi latency:  ${best.avgViterbiMs.toFixed(1)}ms/oligo`);

  // Per-IDS summary
  console.log("\n  Per-IDS recovery rates (best penalties, 20× coverage, 8B LDPC):");
  for (const ids of idsRates) {
    const matching = parityResults.filter(r => Math.abs(r.idsRate - ids) < 0.001 && r.coverage === 20 && r.ldpcParity === 8);
    if (matching.length > 0) {
      const bestForIds = matching.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, matching[0]);
      const status = bestForIds.recoveryRate >= 0.95 ? "✓" : bestForIds.recoveryRate >= 0.80 ? "~" : "✗";
      console.log(`    IDS ${(ids * 100).toFixed(0)}%: ${status} ${(bestForIds.recoveryRate * 100).toFixed(1)}%`);
    }
  }

  // LDPC parity recommendation
  console.log("\n  LDPC parity recommendation (at 9% IDS, 20× coverage):");
  const np9Results = parityResults.filter(r => Math.abs(r.idsRate - 0.09) < 0.001 && r.coverage === 20);
  for (const r of np9Results) {
    const status = r.recoveryRate >= 0.95 ? "✓" : r.recoveryRate >= 0.80 ? "~" : "✗";
    console.log(`    ${r.ldpcParity}B: ${status} ${(r.recoveryRate * 100).toFixed(1)}%`);
  }

  // Save results
  const outputDir = path.join(process.cwd(), "test-data");
  if (!fs.existsSync(outputDir)) fs.mkdirSync(outputDir, { recursive: true });
  const outputPath = path.join(outputDir, "k9-penalty-tuning-results.json");
  fs.writeFileSync(outputPath, JSON.stringify({
    penaltyResults: results,
    parityResults,
    best: { insPen: best.insPen, delPen: best.delPen, maxDrift: best.maxDrift },
    nativeViterbi: !!nativeAddon,
    timestamp: new Date().toISOString(),
  }, null, 2));
  console.log(`\n  Results saved to ${outputPath}`);
}

main().catch(err => {
  console.error("Penalty tuning failed:", err);
  process.exit(1);
});
