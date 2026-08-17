/**
 * Focused Cascade Validation — Minimal combo grid for rapid iteration.
 *
 * Uses Rust WASM Viterbi for the K=9 hot path.
 * Tests the most important IDS rates with a small number of oligos
 * to get per-stage recovery breakdown in under 2 minutes.
 *
 * Usage:
 *   npx tsx scripts/cascade-focused.ts
 */

import { NASA_K9_CONFIG } from '../src/lib/dna/convolutional-k9';
import { IndelTolerantConvolutionalInnerCode, enableWasmViterbi, isWasmViterbiActive } from '../src/lib/dna/convolutional-indel';
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

// ---- Soft consensus (position-majority vote) ----
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
interface StageBreakdown {
  viterbiDecodeOk: number;   // Viterbi completed (didn't throw)
  viterbiCorrect: number;    // Viterbi output matches original after CRC+LDPC
  crcPass: number;           // CRC check passed
  ldpcOk: number;            // LDPC decode succeeded
  ldpcCorrect: number;       // LDPC decode matched original
  failedAll: number;         // Failed all stages
}

interface TestResult {
  idsRate: number; coverage: number; ldpcParity: number;
  insPen: number; delPen: number; maxDrift: number;
  totalOligos: number; recovered: number; recoveryRate: number;
  stages: StageBreakdown;
  avgSubs: number; avgIns: number; avgDels: number;
  decodeMs: number;
}

function runTest(
  numOligos: number, payloadBytes: number, idsRate: number, coverage: number,
  ldpcParityBytes: number, insPen: number, delPen: number, maxDrift: number,
  rng: Rng,
): TestResult {
  const delRate = idsRate * 0.45;
  const insRate = idsRate * 0.30;
  const subRate = idsRate * 0.25;

  const innerDataBytes = payloadBytes + ldpcParityBytes + 2; // +CRC-16
  const convInner = new IndelTolerantConvolutionalInnerCode(innerDataBytes, {
    conv: NASA_K9_CONFIG,
    maxDrift,
    insertionPenalty: insPen,
    deletionPenalty: delPen,
  });

  let ldpcCode: LDPCInnerCode | null = null;
  try {
    ldpcCode = getCachedLDPCInner(payloadBytes + ldpcParityBytes, payloadBytes);
  } catch { ldpcCode = null; }

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
    const crc = crc16Bytes(ldpcCW); // Uint8Array(2): [high, low]
    withCrc[ldpcCW.length] = crc[0];
    withCrc[ldpcCW.length + 1] = crc[1];

    // Conv encode
    let convOut: Uint8Array;
    try { convOut = convInner.encode(withCrc); } catch { convOut = withCrc; }

    // Map to DNA (2 bits per base)
    const bits: number[] = [];
    for (let b = 0; b < convOut.length; b++) {
      for (let bit = 7; bit >= 0; bit--) bits.push((convOut[b] >> bit) & 1);
    }
    const dna: string[] = [];
    for (let b = 0; b + 1 < bits.length; b += 2) dna.push(BASES[(bits[b] << 1) | bits[b + 1]]);
    if (bits.length % 2 === 1) dna.push(BASES[bits[bits.length - 1] << 1]);

    encoded.push({ original: payload, dna: dna.join('') });
  }

  // Decode with per-stage tracking
  const t0 = Date.now();
  const stages: StageBreakdown = {
    viterbiDecodeOk: 0, viterbiCorrect: 0,
    crcPass: 0, ldpcOk: 0, ldpcCorrect: 0, failedAll: 0,
  };
  let recovered = 0;
  let totalSubs = 0, totalIns = 0, totalDels = 0;

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
    let viterbiOk = false;
    try {
      const { decoded } = convInner.decode(consBytes);
      afterConv = decoded;
      viterbiOk = true;
      stages.viterbiDecodeOk++;
    } catch {
      afterConv = consBytes.slice(0, innerDataBytes);
    }

    // CRC check
    let crcPass = false;
    if (afterConv.length >= 2) {
      const dataPart = afterConv.slice(0, afterConv.length - 2);
      const computedCrc = crc16Bytes(dataPart); // Uint8Array(2): [high, low]
      crcPass = afterConv[afterConv.length - 2] === computedCrc[0] && afterConv[afterConv.length - 1] === computedCrc[1];
    }
    if (crcPass) stages.crcPass++;

    // LDPC decode
    let decodedPayload: Uint8Array | null = null;
    if (ldpcCode && afterConv.length >= payloadBytes + ldpcParityBytes) {
      const ldpcRecv = afterConv.slice(0, payloadBytes + ldpcParityBytes);
      try {
        const { data: ldpcDec } = ldpcCode.decode(ldpcRecv);
        if (ldpcDec.length === payloadBytes) {
          decodedPayload = ldpcDec;
          stages.ldpcOk++;
        }
      } catch { /* LDPC failed */ }
    }

    // CRC fallback
    if (!decodedPayload && crcPass && afterConv.length >= payloadBytes) {
      decodedPayload = afterConv.slice(0, payloadBytes);
    }

    // Check recovery
    if (decodedPayload && decodedPayload.length === payloadBytes) {
      let match = true;
      for (let b = 0; b < payloadBytes; b++) {
        if (decodedPayload[b] !== e.original[b]) { match = false; break; }
      }
      if (match) {
        recovered++;
        if (viterbiOk) stages.viterbiCorrect++;
        stages.ldpcCorrect++;
      } else {
        stages.failedAll++;
      }
    } else {
      stages.failedAll++;
    }
  }

  const decodeMs = Date.now() - t0;
  return {
    idsRate, coverage, ldpcParity: ldpcParityBytes,
    insPen, delPen, maxDrift,
    totalOligos: numOligos, recovered, recoveryRate: recovered / numOligos,
    stages,
    avgSubs: totalSubs / numOligos, avgIns: totalIns / numOligos, avgDels: totalDels / numOligos,
    decodeMs,
  };
}

// ---- Main ----
async function main() {
  console.log('╔══════════════════════════════════════════════════════════════╗');
  console.log('║   Focused Cascade Validation — Rust WASM Viterbi+LDPC+RS    ║');
  console.log('╚══════════════════════════════════════════════════════════════╝\n');

  // Enable Rust WASM Viterbi
  const wasmOk = await enableWasmViterbi();
  console.log(`Rust WASM Viterbi: ${wasmOk ? '✓ ENABLED' : '✗ disabled (JS fallback)'}\n`);

  const rng = new Rng(42);
  const numOligos = 20;
  const payloadBytes = 30;

  // Focused test matrix: IDS rate × coverage × LDPC parity × Viterbi penalties
  const testConfigs = [
    // Low IDS (2-4%) — should be easy
    { ids: 0.02, cov: 10, par: 4, ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.04, cov: 10, par: 6, ins: 1.5, del: 1.0, drift: 15 },
    // Medium IDS (6%) — target operating point
    { ids: 0.06, cov: 20, par: 6, ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.06, cov: 20, par: 8, ins: 1.2, del: 0.8, drift: 15 },
    { ids: 0.06, cov: 30, par: 8, ins: 1.5, del: 1.0, drift: 15 },
    // Nanopore IDS (9%) — design target
    { ids: 0.09, cov: 20, par: 8, ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.09, cov: 30, par: 8, ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.09, cov: 30, par: 10, ins: 1.2, del: 0.8, drift: 15 },
    { ids: 0.09, cov: 30, par: 10, ins: 1.5, del: 1.0, drift: 10 },
    // High IDS (12%) — stress test
    { ids: 0.12, cov: 30, par: 10, ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.12, cov: 30, par: 10, ins: 1.0, del: 0.6, drift: 20 },
  ];

  console.log(`Running ${testConfigs.length} test configs × ${numOligos} oligos each\n`);

  const results: TestResult[] = [];
  let testNum = 0;
  for (const c of testConfigs) {
    testNum++;
    process.stdout.write(
      `  [${testNum}/${testConfigs.length}] ` +
      `IDS ${(c.ids * 100).toFixed(0)}% × ${c.cov}× cov × ${c.par}B parity ` +
      `ins=${c.ins} del=${c.del} drift=${c.drift} ...`
    );

    const r = runTest(numOligos, payloadBytes, c.ids, c.cov, c.par, c.ins, c.del, c.drift, rng);
    results.push(r);

    console.log(
      ` ${(r.recoveryRate * 100).toFixed(1).padStart(6)}% ` +
      `(V:${r.stages.viterbiDecodeOk} C:${r.stages.crcPass} L:${r.stages.ldpcOk}) ` +
      `[${r.decodeMs}ms]`
    );
  }

  // Per-stage breakdown table
  console.log('\n╔═════════════════════════════════════════════════════════════════════════════════════════════╗');
  console.log('║  Per-Stage Recovery Breakdown                                                             ║');
  console.log('╠══════════╦════════╦══════╦══════╦════════╦════════╦════════╦════════╦════════════════════╣');
  console.log('║ IDS      │ Cov    │ Par  │ Rec% │ VitOk  │ CRCok  │ LDPCok │ Fail   │ Penalties          ║');
  console.log('╠══════════╬════════╬══════╬══════╬════════╬════════╬════════╬════════╬════════════════════╣');

  for (const r of results) {
    console.log(
      `║ ${(r.idsRate * 100).toFixed(0).padStart(3)}%    │ ` +
      `${String(r.coverage).padStart(4)}×   │ ` +
      `${String(r.ldpcParity).padStart(3)}B  │ ` +
      `${(r.recoveryRate * 100).toFixed(1).padStart(5)}% │ ` +
      `${String(r.stages.viterbiDecodeOk).padStart(5)}  │ ` +
      `${String(r.stages.crcPass).padStart(5)}  │ ` +
      `${String(r.stages.ldpcOk).padStart(5)}  │ ` +
      `${String(r.stages.failedAll).padStart(5)}  │ ` +
      `i${r.insPen} d${r.delPen} m${r.maxDrift}`.padEnd(21) + '║'
    );
  }
  console.log('╚══════════╩════════╩══════╩══════╩════════╩════════╩════════╩════════╩════════════════════╝');

  // Aggregate analysis
  console.log('\n=== Aggregate Stage Success Rates ===');
  const totalOligos = results.reduce((s, r) => s + r.totalOligos, 0);
  const totalV = results.reduce((s, r) => s + r.stages.viterbiDecodeOk, 0);
  const totalC = results.reduce((s, r) => s + r.stages.crcPass, 0);
  const totalL = results.reduce((s, r) => s + r.stages.ldpcOk, 0);
  const totalF = results.reduce((s, r) => s + r.stages.failedAll, 0);
  const totalR = results.reduce((s, r) => s + r.recovered, 0);
  console.log(`  Viterbi decode OK:  ${totalV}/${totalOligos} (${(totalV/totalOligos*100).toFixed(1)}%)`);
  console.log(`  CRC pass:           ${totalC}/${totalOligos} (${(totalC/totalOligos*100).toFixed(1)}%)`);
  console.log(`  LDPC decode OK:     ${totalL}/${totalOligos} (${(totalL/totalOligos*100).toFixed(1)}%)`);
  console.log(`  Failed all stages:  ${totalF}/${totalOligos} (${(totalF/totalOligos*100).toFixed(1)}%)`);
  console.log(`  Total recovered:    ${totalR}/${totalOligos} (${(totalR/totalOligos*100).toFixed(1)}%)`);

  // Nanopore 9% IDS specific analysis
  console.log('\n=== Nanopore 9% IDS Focus ===');
  const np9 = results.filter(r => Math.abs(r.idsRate - 0.09) < 0.001);
  if (np9.length > 0) {
    for (const r of np9) {
      console.log(
        `  ${r.coverage}× cov × ${r.ldpcParity}B parity × ins=${r.insPen} del=${r.delPen}: ` +
        `${(r.recoveryRate * 100).toFixed(1)}% ` +
        `(V:${r.stages.viterbiDecodeOk}/${r.totalOligos} C:${r.stages.crcPass} L:${r.stages.ldpcOk}) ` +
        `[${r.decodeMs}ms]`
      );
    }
    const best9 = np9.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, np9[0]);
    console.log(`\n  → Best 9% config: ${(best9.recoveryRate * 100).toFixed(1)}% @ ${best9.coverage}× cov, ${best9.ldpcParity}B LDPC, ins=${best9.insPen} del=${best9.delPen}`);
  }

  // Tuning recommendations
  console.log('\n=== Tuning Recommendations ===');

  // Viterbi penalty analysis
  const byIds = new Map<number, TestResult[]>();
  for (const r of results) {
    const arr = byIds.get(r.idsRate) || [];
    arr.push(r);
    byIds.set(r.idsRate, arr);
  }

  for (const [ids, matching] of byIds) {
    const best = matching.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, matching[0]);
    const avgVitRate = matching.reduce((s, r) => s + r.stages.viterbiDecodeOk / r.totalOligos, 0) / matching.length;

    if (best.recoveryRate >= 0.95) {
      console.log(`  IDS ${(ids * 100).toFixed(0)}%: ✓ ${(best.recoveryRate * 100).toFixed(1)}% @ ${best.coverage}× cov, ${best.ldpcParity}B LDPC, ins=${best.insPen} del=${best.delPen}`);
    } else if (best.recoveryRate >= 0.80) {
      console.log(`  IDS ${(ids * 100).toFixed(0)}%: ~ ${(best.recoveryRate * 100).toFixed(1)}% @ ${best.coverage}× cov, ${best.ldpcParity}B LDPC, ins=${best.insPen} del=${best.delPen} — needs more coverage or parity`);
    } else {
      console.log(`  IDS ${(ids * 100).toFixed(0)}%: ✗ ${(best.recoveryRate * 100).toFixed(1)}% @ ${best.coverage}× cov, ${best.ldpcParity}B LDPC — Viterbi success ${(avgVitRate * 100).toFixed(1)}%`);
    }

    // Suggest penalty adjustments
    if (avgVitRate < 0.7 && ids <= 0.09) {
      console.log(`    → Viterbi struggling — suggest ins=1.0 del=0.6 drift=20`);
    } else if (avgVitRate < 0.85 && ids <= 0.06) {
      console.log(`    → Viterbi borderline — suggest ins=1.2 del=0.8`);
    }
  }

  // LDPC parity analysis
  console.log('\n=== LDPC Parity Analysis ===');
  for (const [ids, matching] of byIds) {
    for (const r of matching) {
      const ldpcFailRate = 1 - r.stages.ldpcOk / r.totalOligos;
      if (ldpcFailRate > 0.15 && r.ldpcParity < 12) {
        console.log(`  IDS ${(ids * 100).toFixed(0)}% × ${r.coverage}× × ${r.ldpcParity}B: LDPC fail ${(ldpcFailRate * 100).toFixed(1)}% → increase parity to ${r.ldpcParity + 2}B`);
      }
    }
  }

  // Save results
  const outputDir = path.join(process.cwd(), 'test-data');
  if (!fs.existsSync(outputDir)) fs.mkdirSync(outputDir, { recursive: true });
  const outputPath = path.join(outputDir, 'cascade-focused-results.json');
  fs.writeFileSync(outputPath, JSON.stringify({
    wasmViterbi: wasmOk,
    results,
    timestamp: new Date().toISOString(),
  }, null, 2));
  console.log(`\nResults saved to ${outputPath}`);
}

main().catch(err => {
  console.error('Cascade focused validation failed:', err);
  process.exit(1);
});
