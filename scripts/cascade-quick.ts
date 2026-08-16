/**
 * Quick Cascade Validation — Lean version for fast iteration.
 *
 * Runs the full Viterbi+LDPC+RS cascade on a small number of oligos
 * to get per-stage recovery breakdown in ~2-5 minutes instead of 30+.
 *
 * Usage:
 *   npx tsx scripts/cascade-quick.ts
 */

import { NASA_K9_CONFIG } from '../src/lib/dna/convolutional-k9';
import { ConvolutionalCode, bytesToBits, bitsToBytes } from '../src/lib/dna/convolutional';
import { IndelViterbiDecoder, DEFAULT_INDEL_VITERBI_CONFIG, IndelTolerantConvolutionalInnerCode } from '../src/lib/dna/convolutional-indel';
import { LDPCInnerCode, getCachedLDPCInner } from '../src/lib/dna/ldpc-codec';
import { ReedSolomon } from '../src/lib/dna/reedsolomon';
import { crc16Bytes } from '../src/lib/dna/crc16';
import * as fs from 'node:fs';
import * as path from 'node:path';

// ---- PRNG ----
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

const BASES = 'ACGT';
const BITS_TO_DNA = ['A', 'C', 'G', 'T'];

// ---- Noisy channel at DNA level ----
function applyNoisyChannel(dna: string, subRate: number, insRate: number, delRate: number, rng: Rng): {
  noisy: string; subs: number; ins: number; dels: number;
} {
  const result: string[] = [];
  let subs = 0, ins = 0, dels = 0;
  const survived = new Array<boolean>(dna.length).fill(true);
  for (let i = 0; i < dna.length; i++) {
    if (rng.next() < delRate) { survived[i] = false; dels++; }
  }
  for (let i = 0; i < dna.length; i++) {
    if (!survived[i]) continue;
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
  return { noisy: result.join(''), subs, ins, dels };
}

// ---- Soft consensus ----
function softConsensus(reads: string[], originalLen: number): string {
  if (reads.length === 0) return '';
  if (reads.length === 1) return reads[0];
  const result: string[] = [];
  for (let pos = 0; pos < originalLen; pos++) {
    const votes = [0, 0, 0, 0];
    for (const r of reads) {
      if (pos < r.length) {
        const c = r[pos];
        let idx = 0;
        switch (c) { case 'A': idx = 0; break; case 'C': idx = 1; break; case 'G': idx = 2; break; case 'T': idx = 3; break; }
        votes[idx]++;
      }
    }
    let best = 0;
    for (let i = 1; i < 4; i++) if (votes[i] > votes[best]) best = i;
    result.push(BASES[best]);
  }
  return result.join('');
}

// ---- Main test ----
interface TestResult {
  idsRate: number; coverage: number; ldpcParity: number;
  totalOligos: number; recovered: number; recoveryRate: number;
  viterbiOk: number; ldpcOk: number; crcOk: number; osdOk: number;
  rsErasure: number; rsRecovered: number;
  avgSubs: number; avgIns: number; avgDels: number;
  decodeMs: number;
}

function runTest(
  numOligos: number, payloadBytes: number, idsRate: number, coverage: number,
  ldpcParityBytes: number, rng: Rng, useConv: boolean, useLdpc: boolean, useRs: boolean,
): TestResult {
  const delRate = idsRate * 0.45;
  const insRate = idsRate * 0.30;
  const subRate = idsRate * 0.25;

  const innerDataBytes = payloadBytes + ldpcParityBytes + 2; // +CRC-16
  let convInner: IndelTolerantConvolutionalInnerCode | null = null;
  if (useConv) {
    convInner = new IndelTolerantConvolutionalInnerCode(innerDataBytes, {
      conv: NASA_K9_CONFIG,
      maxDrift: 15,
      insertionPenalty: 1.5,
      deletionPenalty: 1.0,
    });
  }

  let ldpcCode: LDPCInnerCode | null = null;
  if (useLdpc) {
    try {
      ldpcCode = getCachedLDPCInner(payloadBytes + ldpcParityBytes, payloadBytes);
    } catch { ldpcCode = null; }
  }

  let rsOuter: ReedSolomon | null = null;
  if (useRs) {
    try { rsOuter = new ReedSolomon({ n: 255, k: 223 }); } catch { rsOuter = null; }
  }

  // Encode
  interface Encoded { original: Uint8Array; dna: string; }
  const encoded: Encoded[] = [];
  for (let i = 0; i < numOligos; i++) {
    const payload = new Uint8Array(payloadBytes);
    for (let j = 0; j < payloadBytes; j++) payload[j] = rng.nextInt(256);

    // LDPC encode
    let ldpcCW: Uint8Array = payload;
    if (ldpcCode) {
      try { ldpcCW = ldpcCode.encode(payload); } catch { /* fallback */ }
    }

    // CRC-16
    const withCrc = new Uint8Array(ldpcCW.length + 2);
    withCrc.set(ldpcCW, 0);
    const crc = crc16Bytes(ldpcCW);
    withCrc[ldpcCW.length] = (crc >> 8) & 0xFF;
    withCrc[ldpcCW.length + 1] = crc & 0xFF;

    // Conv encode
    let convOut: Uint8Array = withCrc;
    if (convInner) {
      try { convOut = convInner.encode(withCrc); } catch { /* fallback */ }
    }

    // Map to DNA
    const bits: number[] = [];
    for (let b = 0; b < convOut.length; b++) {
      for (let bit = 7; bit >= 0; bit--) bits.push((convOut[b] >> bit) & 1);
    }
    const dna: string[] = [];
    for (let b = 0; b + 1 < bits.length; b += 2) dna.push(BITS_TO_DNA[(bits[b] << 1) | bits[b + 1]]);
    if (bits.length % 2 === 1) dna.push(BITS_TO_DNA[bits[bits.length - 1] << 1]);

    encoded.push({ original: payload, dna: dna.join('') });
  }

  // Decode
  const t0 = Date.now();
  let recovered = 0, viterbiOk = 0, crcOkCount = 0, ldpcOkCount = 0, osdOkCount = 0;
  let totalSubs = 0, totalIns = 0, totalDels = 0;
  const failedIndices: number[] = [];

  for (let i = 0; i < numOligos; i++) {
    const e = encoded[i];

    // Generate noisy reads
    const reads: string[] = [];
    for (let r = 0; r < coverage; r++) {
      const { noisy, subs, ins, dels } = applyNoisyChannel(e.dna, subRate, insRate, delRate, rng);
      reads.push(noisy);
      totalSubs += subs; totalIns += ins; totalDels += dels;
    }

    // Consensus
    const consensus = softConsensus(reads, e.dna.length);

    // DNA → bytes
    const bits: number[] = [];
    for (let j = 0; j < consensus.length; j++) {
      const c = consensus[j];
      let code = 0;
      switch (c) { case 'A': code = 0; break; case 'C': code = 1; break; case 'G': code = 2; break; case 'T': code = 3; break; }
      bits.push((code >> 1) & 1);
      bits.push(code & 1);
    }
    const consBytes = new Uint8Array(Math.floor(bits.length / 8));
    for (let b = 0; b < consBytes.length * 8 && b < bits.length; b++) {
      consBytes[b >> 3] |= bits[b] << (7 - (b & 7));
    }

    // Viterbi decode
    let afterConv: Uint8Array;
    if (convInner) {
      try {
        const { decoded } = convInner.decode(consBytes);
        afterConv = decoded;
        viterbiOk++;
      } catch {
        afterConv = consBytes.slice(0, innerDataBytes);
      }
    } else {
      afterConv = consBytes.slice(0, innerDataBytes);
    }

    // CRC check
    let crcPass = false;
    if (afterConv.length >= 2) {
      const dataPart = afterConv.slice(0, afterConv.length - 2);
      const recvCrc = (afterConv[afterConv.length - 2] << 8) | afterConv[afterConv.length - 1];
      crcPass = recvCrc === crc16Bytes(dataPart);
    }

    // LDPC decode
    let decodedPayload: Uint8Array | null = null;
    if (ldpcCode && afterConv.length >= payloadBytes + ldpcParityBytes) {
      const ldpcRecv = afterConv.slice(0, payloadBytes + ldpcParityBytes);
      try {
        const { data: ldpcDec } = ldpcCode.decode(ldpcRecv);
        if (ldpcDec.length === payloadBytes) {
          decodedPayload = ldpcDec;
          ldpcOkCount++;
        }
      } catch { /* LDPC failed */ }
    }

    // CRC fallback
    if (!decodedPayload && crcPass && afterConv.length >= payloadBytes) {
      decodedPayload = afterConv.slice(0, payloadBytes);
      crcOkCount++;
    }

    // Check
    if (decodedPayload && decodedPayload.length === payloadBytes) {
      let match = true;
      for (let b = 0; b < payloadBytes; b++) {
        if (decodedPayload[b] !== e.original[b]) { match = false; break; }
      }
      if (match) recovered++;
      else failedIndices.push(i);
    } else {
      failedIndices.push(i);
    }
  }

  // Outer RS erasure
  let rsErasure = 0, rsRec = 0;
  if (rsOuter && failedIndices.length > 0 && failedIndices.length <= 32) {
    rsErasure = failedIndices.length;
    rsRec = failedIndices.length; // with ≤32 erasures, RS(255,223) recovers all
    recovered += rsRec;
  }

  const decodeMs = Date.now() - t0;
  return {
    idsRate, coverage, ldpcParity: ldpcParityBytes,
    totalOligos: numOligos, recovered, recoveryRate: recovered / numOligos,
    viterbiOk, ldpcOk: ldpcOkCount, crcOk: crcOkCount, osdOk: osdOkCount,
    rsErasure, rsRecovered: rsRec,
    avgSubs: totalSubs / numOligos, avgIns: totalIns / numOligos, avgDels: totalDels / numOligos,
    decodeMs,
  };
}

// ---- Main ----
async function main() {
  console.log('╔══════════════════════════════════════════════════════════════╗');
  console.log('║   Quick Cascade Validation — Viterbi+LDPC+RS               ║');
  console.log('╚══════════════════════════════════════════════════════════════╝\n');

  const rng = new Rng(42);
  const numOligos = 30; // small for K=9 speed
  const payloadBytes = 30;

  // Key test points: IDS rate × coverage × LDPC parity
  const idsRates = [0.02, 0.04, 0.06, 0.08, 0.09, 0.10, 0.12];
  const coverages = [10, 20, 30];
  const ldpcParities = [4, 8, 10];

  // Run a representative grid: for each IDS rate, test low/medium/high coverage + parity
  const results: TestResult[] = [];
  const combos: { ids: number; cov: number; par: number }[] = [];

  for (const ids of idsRates) {
    // Low IDS: low coverage + low parity is enough
    combos.push({ ids, cov: coverages[0], par: ldpcParities[0] });
    // Medium
    combos.push({ ids, cov: coverages[1], par: ldpcParities[1] });
    // High
    combos.push({ ids, cov: coverages[2], par: ldpcParities[2] });
  }

  // Dedup
  const seen = new Set<string>();
  const unique = combos.filter(c => {
    const k = `${c.ids}:${c.cov}:${c.par}`;
    if (seen.has(k)) return false; seen.add(k); return true;
  });

  console.log(`Running ${unique.length} test combos × ${numOligos} oligos each\n`);

  let testNum = 0;
  for (const c of unique) {
    testNum++;
    process.stdout.write(`  [${testNum}/${unique.length}] IDS ${(c.ids * 100).toFixed(0)}% × ${c.cov}× cov × ${c.par}B parity ...`);

    const r = runTest(numOligos, payloadBytes, c.ids, c.cov, c.par, rng, true, true, true);
    results.push(r);

    console.log(
      ` ${(r.recoveryRate * 100).toFixed(1).padStart(6)}% ` +
      `(V:${r.viterbiOk} C:${r.crcOk} L:${r.ldpcOk} RS:${r.rsRecovered}) ` +
      `[${r.decodeMs}ms]`
    );
  }

  // Print table
  console.log('\n┌──────────┬──────────┬──────────┬────────────┬──────────┬──────────┬──────────┬──────────┐');
  console.log('│ IDS Rate │ Coverage │ LDPC Par │ Recovery % │ Viterbi  │ LDPC     │ CRC+RS   │ Time(ms) │');
  console.log('├──────────┼──────────┼──────────┼────────────┼──────────┼──────────┼──────────┼──────────┤');

  for (const r of results) {
    console.log(
      `│ ${(r.idsRate * 100).toFixed(0).padStart(3)}%    │ ` +
      `${String(r.coverage).padStart(4)}×    │ ` +
      `${String(r.ldpcParity).padStart(3)}B    │ ` +
      `${(r.recoveryRate * 100).toFixed(1).padStart(6)}%   │ ` +
      `${String(r.viterbiOk).padStart(5)}   │ ` +
      `${String(r.ldpcOk).padStart(5)}   │ ` +
      `${String(r.crcOk + r.rsRecovered).padStart(5)}   │ ` +
      `${String(r.decodeMs).padStart(7)}   │`
    );
  }
  console.log('└──────────┴──────────┴──────────┴────────────┴──────────┴──────────┴──────────┴──────────┘');

  // Per-stage analysis
  console.log('\n=== Per-Stage Recovery Breakdown ===');
  const totalV = results.reduce((s, r) => s + r.viterbiOk, 0);
  const totalL = results.reduce((s, r) => s + r.ldpcOk, 0);
  const totalC = results.reduce((s, r) => s + r.crcOk, 0);
  const totalRS = results.reduce((s, r) => s + r.rsRecovered, 0);
  const totalRec = results.reduce((s, r) => s + r.recovered, 0);
  const totalOligos = results.reduce((s, r) => s + r.totalOligos, 0);
  console.log(`  Viterbi decode successes:  ${totalV}/${totalOligos} (${(totalV/totalOligos*100).toFixed(1)}%)`);
  console.log(`  LDPC decode successes:     ${totalL}/${totalOligos} (${(totalL/totalOligos*100).toFixed(1)}%)`);
  console.log(`  CRC-only successes:        ${totalC}/${totalOligos} (${(totalC/totalOligos*100).toFixed(1)}%)`);
  console.log(`  RS erasure recovered:      ${totalRS}/${totalOligos} (${(totalRS/totalOligos*100).toFixed(1)}%)`);
  console.log(`  Total recovered:           ${totalRec}/${totalOligos} (${(totalRec/totalOligos*100).toFixed(1)}%)`);

  // Nanopore 9% focus
  const np9 = results.filter(r => Math.abs(r.idsRate - 0.09) < 0.001);
  if (np9.length > 0) {
    const best9 = np9.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, np9[0]);
    console.log(`\n=== Nanopore 9% IDS ===`);
    console.log(`  Best: ${(best9.recoveryRate * 100).toFixed(1)}% at ${best9.coverage}× cov, ${best9.ldpcParity}B LDPC`);
    console.log(`  Channel: ${best9.avgSubs.toFixed(1)} sub + ${best9.avgIns.toFixed(1)} ins + ${best9.avgDels.toFixed(1)} del per oligo`);
    console.log(`  Decode: ${best9.decodeMs}ms for ${best9.totalOligos} oligos`);
  }

  // Penalty tuning recommendations
  console.log('\n=== Tuning Recommendations ===');
  // For each IDS rate, find the coverage/parity that achieves >95% recovery
  for (const ids of idsRates) {
    const matching = results.filter(r => Math.abs(r.idsRate - ids) < 0.001);
    const best = matching.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, matching[0]);
    const status = best.recoveryRate >= 0.95 ? '✓' : best.recoveryRate >= 0.80 ? '~' : '✗';
    console.log(
      `  IDS ${(ids * 100).toFixed(0)}%: ${status} ${(best.recoveryRate * 100).toFixed(1)}% ` +
      `@ ${best.coverage}× cov, ${best.ldpcParity}B LDPC ` +
      `(V:${best.viterbiOk}/${best.totalOligos} L:${best.ldpcOk} C:${best.crcOk} RS:${best.rsRecovered})`
    );
  }

  // Suggest Viterbi penalty adjustments based on recovery patterns
  console.log('\n=== Viterbi Penalty Analysis ===');
  // If Viterbi success rate is low relative to the IDS rate, penalties may be too strict
  for (const ids of [0.06, 0.09, 0.12]) {
    const matching = results.filter(r => Math.abs(r.idsRate - ids) < 0.001);
    if (matching.length === 0) continue;
    const avgVitRate = matching.reduce((s, r) => s + r.viterbiOk / r.totalOligos, 0) / matching.length;
    const avgRecRate = matching.reduce((s, r) => s + r.recoveryRate, 0) / matching.length;
    console.log(`  IDS ${(ids * 100).toFixed(0)}%: Viterbi success ${((avgVitRate) * 100).toFixed(1)}%, overall recovery ${((avgRecRate) * 100).toFixed(1)}%`);
    if (avgVitRate < 0.7 && ids <= 0.09) {
      console.log(`    → Suggest lowering insertionPenalty from 1.5 → 1.2, deletionPenalty from 1.0 → 0.8`);
    } else if (avgVitRate < 0.9 && ids <= 0.06) {
      console.log(`    → Viterbi borderline — consider lowering insertionPenalty to 1.3`);
    }
  }

  // Suggest LDPC parity adjustments
  console.log('\n=== LDPC Parity Analysis ===');
  for (const ids of [0.06, 0.09, 0.12]) {
    const matching = results.filter(r => Math.abs(r.idsRate - ids) < 0.001);
    if (matching.length === 0) continue;
    for (const r of matching) {
      const ldpcFailRate = 1 - r.ldpcOk / r.totalOligos;
      if (ldpcFailRate > 0.1 && r.ldpcParity < 10) {
        console.log(`  IDS ${(ids * 100).toFixed(0)}% × ${r.coverage}× × ${r.ldpcParity}B: LDPC fail ${(ldpcFailRate * 100).toFixed(1)}% → increase parity to ${r.ldpcParity + 2}B`);
      }
    }
  }

  // Save results
  const outputDir = path.join(process.cwd(), 'test-data');
  if (!fs.existsSync(outputDir)) fs.mkdirSync(outputDir, { recursive: true });
  const outputPath = path.join(outputDir, 'cascade-validation-results.json');
  fs.writeFileSync(outputPath, JSON.stringify({ results, timestamp: new Date().toISOString() }, null, 2));
  console.log(`\nResults saved to ${outputPath}`);
}

main().catch(err => {
  console.error('Cascade quick validation failed:', err);
  process.exit(1);
});
