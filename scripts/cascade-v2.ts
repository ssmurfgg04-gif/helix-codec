/**
 * Focused Cascade Validation v2 — Per-read Viterbi decode + majority vote.
 *
 * Key insight: The indel-tolerant Viterbi decoder IS the indel corrector.
 * We should NOT do position-majority consensus before Viterbi (it misaligns
 * reads with indels). Instead, decode each read through Viterbi individually,
 * then majority-vote the decoded byte sequences.
 *
 * Pipeline per oligo:
 *   1. Generate coverage × noisy reads
 *   2. For each read: DNA→bytes → Viterbi decode → decoded bytes
 *   3. Majority vote the decoded byte sequences (per-position)
 *   4. CRC check → LDPC decode → verify
 *
 * Uses Rust WASM Viterbi for the K=9 hot path.
 */

import { NASA_K9_CONFIG } from '../src/lib/dna/convolutional-k9';
import { IndelTolerantConvolutionalInnerCode, enableWasmViterbi } from '../src/lib/dna/convolutional-indel';
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

function applyNoisyChannel(dna: string, subRate: number, insRate: number, delRate: number, rng: Rng): string {
  const result: string[] = [];
  for (let i = 0; i < dna.length; i++) {
    // Deletion
    if (rng.next() < delRate) continue;
    // Substitution
    let base = dna[i];
    if (rng.next() < subRate) {
      let nb: string;
      do { nb = BASES[rng.nextInt(4)]; } while (nb === base);
      base = nb;
    }
    result.push(base);
    // Insertion
    if (rng.next() < insRate) {
      result.push(BASES[rng.nextInt(4)]);
    }
  }
  return result.join('');
}

/** Convert DNA string to byte array (2 bits per base: A=00 C=01 G=10 T=11) */
function dnaToBytes(dna: string): Uint8Array {
  const bits: number[] = [];
  for (let i = 0; i < dna.length; i++) {
    let code = 0;
    switch (dna[i]) { case 'A': code = 0; break; case 'C': code = 1; break; case 'G': code = 2; break; case 'T': code = 3; break; }
    bits.push((code >> 1) & 1);
    bits.push(code & 1);
  }
  const bytes = new Uint8Array(Math.floor(bits.length / 8));
  for (let b = 0; b < bytes.length * 8 && b < bits.length; b++) {
    bytes[b >> 3] |= bits[b] << (7 - (b & 7));
  }
  return bytes;
}

/** Convert byte array to DNA string */
function bytesToDna(bytes: Uint8Array): string {
  const bits: number[] = [];
  for (let b = 0; b < bytes.length; b++) {
    for (let bit = 7; bit >= 0; bit--) bits.push((bytes[b] >> bit) & 1);
  }
  const dna: string[] = [];
  for (let b = 0; b + 1 < bits.length; b += 2) dna.push(BASES[(bits[b] << 1) | bits[b + 1]]);
  return dna.join('');
}

/** Majority vote across multiple decoded byte arrays (per-position) */
function majorityVoteDecoded(decodedArrays: Uint8Array[], expectedLen: number): Uint8Array {
  if (decodedArrays.length === 0) return new Uint8Array(expectedLen);
  if (decodedArrays.length === 1) return decodedArrays[0];
  const result = new Uint8Array(expectedLen);
  for (let pos = 0; pos < expectedLen; pos++) {
    const votes = new Array(256).fill(0);
    for (const arr of decodedArrays) {
      if (pos < arr.length) votes[arr[pos]]++;
    }
    let best = 0;
    for (let v = 1; v < 256; v++) if (votes[v] > votes[best]) best = v;
    result[pos] = best;
  }
  return result;
}

interface StageBreakdown {
  viterbiDecodeOk: number;
  viterbiPerReadOk: number;
  totalReads: number;
  crcPass: number;
  ldpcOk: number;
  failedAll: number;
}

interface TestResult {
  idsRate: number; coverage: number; ldpcParity: number;
  insPen: number; delPen: number; maxDrift: number;
  totalOligos: number; recovered: number; recoveryRate: number;
  stages: StageBreakdown;
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

  // Encode all oligos
  interface Encoded { original: Uint8Array; convOut: Uint8Array; dna: string; }
  const encoded: Encoded[] = [];
  for (let i = 0; i < numOligos; i++) {
    const payload = new Uint8Array(payloadBytes);
    for (let j = 0; j < payloadBytes; j++) payload[j] = rng.nextInt(256);

    let ldpcCW: Uint8Array = payload;
    if (ldpcCode) {
      try { ldpcCW = ldpcCode.encode(payload); } catch { /* fallback */ }
    }

    const withCrc = new Uint8Array(ldpcCW.length + 2);
    withCrc.set(ldpcCW, 0);
    const crc = crc16Bytes(ldpcCW); // Uint8Array(2): [high, low]
    withCrc[ldpcCW.length] = crc[0];
    withCrc[ldpcCW.length + 1] = crc[1];

    let convOut: Uint8Array;
    try { convOut = convInner.encode(withCrc); } catch { convOut = withCrc; }

    const dna = bytesToDna(convOut);
    encoded.push({ original: payload, convOut, dna });
  }

  // Decode with per-read Viterbi
  const t0 = Date.now();
  const stages: StageBreakdown = {
    viterbiDecodeOk: 0, viterbiPerReadOk: 0, totalReads: 0,
    crcPass: 0, ldpcOk: 0, failedAll: 0,
  };
  let recovered = 0;

  for (let i = 0; i < numOligos; i++) {
    const e = encoded[i];

    // Generate noisy reads and decode each through Viterbi
    const decodedReads: Uint8Array[] = [];
    for (let r = 0; r < coverage; r++) {
      const noisyDna = applyNoisyChannel(e.dna, subRate, insRate, delRate, rng);
      stages.totalReads++;

      // DNA → bytes
      const recvBytes = dnaToBytes(noisyDna);

      // Viterbi decode this read
      try {
        const { decoded } = convInner.decode(recvBytes);
        if (decoded.length === innerDataBytes) {
          decodedReads.push(decoded);
          stages.viterbiPerReadOk++;
        }
      } catch { /* Viterbi failed for this read */ }
    }

    if (decodedReads.length > 0) {
      stages.viterbiDecodeOk++;
    }

    // Majority vote across decoded reads
    const afterConv = majorityVoteDecoded(decodedReads, innerDataBytes);

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
      if (match) recovered++;
      else stages.failedAll++;
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
    decodeMs,
  };
}

async function main() {
  console.log('╔══════════════════════════════════════════════════════════════════╗');
  console.log('║   Cascade Validation v2 — Per-read Viterbi + Majority Vote      ║');
  console.log('╚══════════════════════════════════════════════════════════════════╝\n');

  const wasmOk = await enableWasmViterbi();
  console.log(`Rust WASM Viterbi: ${wasmOk ? '✓ ENABLED' : '✗ disabled'}\n`);

  const rng = new Rng(42);
  const numOligos = 15;
  const payloadBytes = 30;

  // Test matrix focused on key operating points
  const testConfigs = [
    // Zero noise sanity check
    { ids: 0.00, cov: 1,  par: 4,  ins: 1.5, del: 1.0, drift: 15 },
    // Low IDS
    { ids: 0.02, cov: 5,  par: 4,  ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.04, cov: 10, par: 6,  ins: 1.5, del: 1.0, drift: 15 },
    // Medium IDS (6%)
    { ids: 0.06, cov: 10, par: 6,  ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.06, cov: 10, par: 8,  ins: 1.2, del: 0.8, drift: 15 },
    // Nanopore 9% IDS — design target
    { ids: 0.09, cov: 10, par: 8,  ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.09, cov: 15, par: 8,  ins: 1.5, del: 1.0, drift: 15 },
    { ids: 0.09, cov: 15, par: 10, ins: 1.2, del: 0.8, drift: 15 },
    { ids: 0.09, cov: 20, par: 10, ins: 1.5, del: 1.0, drift: 10 },
    // High IDS (12%)
    { ids: 0.12, cov: 20, par: 10, ins: 1.0, del: 0.6, drift: 20 },
    { ids: 0.12, cov: 20, par: 10, ins: 1.5, del: 1.0, drift: 15 },
  ];

  console.log(`Running ${testConfigs.length} configs × ${numOligos} oligos each\n`);

  const results: TestResult[] = [];
  let testNum = 0;
  for (const c of testConfigs) {
    testNum++;
    process.stdout.write(
      `  [${testNum}/${testConfigs.length}] ` +
      `IDS ${(c.ids * 100).toFixed(0)}% × ${c.cov}× × ${c.par}B ` +
      `ins=${c.ins} del=${c.del} d=${c.drift} ...`
    );

    const r = runTest(numOligos, payloadBytes, c.ids, c.cov, c.par, c.ins, c.del, c.drift, rng);
    results.push(r);

    const readOk = r.stages.totalReads > 0
      ? `${(r.stages.viterbiPerReadOk / r.stages.totalReads * 100).toFixed(0)}%`
      : '-';
    console.log(
      ` ${(r.recoveryRate * 100).toFixed(1).padStart(6)}% ` +
      `(Vread:${readOk} C:${r.stages.crcPass} L:${r.stages.ldpcOk}) ` +
      `[${r.decodeMs}ms]`
    );
  }

  // Summary table
  console.log('\n┌──────────┬──────┬──────┬────────┬──────────┬──────┬──────┬──────┬─────────────┐');
  console.log('│ IDS      │ Cov  │ Par  │ Rec%   │ Vread%   │ CRC  │ LDPC │ Fail │ Penalties   │');
  console.log('├──────────┼──────┼──────┼────────┼──────────┼──────┼──────┼──────┼─────────────┤');

  for (const r of results) {
    const vreadPct = r.stages.totalReads > 0
      ? (r.stages.viterbiPerReadOk / r.stages.totalReads * 100).toFixed(0)
      : '-';
    console.log(
      `│ ${(r.idsRate * 100).toFixed(0).padStart(3)}%    │ ` +
      `${String(r.coverage).padStart(4)} │ ` +
      `${String(r.ldpcParity).padStart(3)}B │ ` +
      `${(r.recoveryRate * 100).toFixed(1).padStart(5)}% │ ` +
      `${vreadPct.padStart(7)}% │ ` +
      `${String(r.stages.crcPass).padStart(4)} │ ` +
      `${String(r.stages.ldpcOk).padStart(4)} │ ` +
      `${String(r.stages.failedAll).padStart(4)} │ ` +
      `i${r.insPen} d${r.delPen} m${r.maxDrift}`.padEnd(14) + '│'
    );
  }
  console.log('└──────────┴──────┴──────┴────────┴──────────┴──────┴──────┴──────┴─────────────┘');

  // Per-stage analysis
  console.log('\n=== Per-Stage Recovery ===');
  const totalOligos = results.reduce((s, r) => s + r.totalOligos, 0);
  const totalReads = results.reduce((s, r) => s + r.stages.totalReads, 0);
  const totalVReadOk = results.reduce((s, r) => s + r.stages.viterbiPerReadOk, 0);
  const totalVOk = results.reduce((s, r) => s + r.stages.viterbiDecodeOk, 0);
  const totalC = results.reduce((s, r) => s + r.stages.crcPass, 0);
  const totalL = results.reduce((s, r) => s + r.stages.ldpcOk, 0);
  const totalR = results.reduce((s, r) => s + r.recovered, 0);
  console.log(`  Per-read Viterbi OK: ${totalVReadOk}/${totalReads} (${(totalVReadOk/totalReads*100).toFixed(1)}%)`);
  console.log(`  Oligo Viterbi OK:    ${totalVOk}/${totalOligos} (${(totalVOk/totalOligos*100).toFixed(1)}%)`);
  console.log(`  CRC pass:            ${totalC}/${totalOligos} (${(totalC/totalOligos*100).toFixed(1)}%)`);
  console.log(`  LDPC OK:             ${totalL}/${totalOligos} (${(totalL/totalOligos*100).toFixed(1)}%)`);
  console.log(`  Total recovered:     ${totalR}/${totalOligos} (${(totalR/totalOligos*100).toFixed(1)}%)`);

  // Nanopore 9% focus
  console.log('\n=== Nanopore 9% IDS ===');
  const np9 = results.filter(r => Math.abs(r.idsRate - 0.09) < 0.001);
  for (const r of np9) {
    const vreadPct = r.stages.totalReads > 0
      ? (r.stages.viterbiPerReadOk / r.stages.totalReads * 100).toFixed(1)
      : '-';
    console.log(
      `  ${r.coverage}× × ${r.ldpcParity}B × i${r.insPen} d${r.delPen}: ` +
      `${(r.recoveryRate * 100).toFixed(1)}% ` +
      `(Vread:${vreadPct}% C:${r.stages.crcPass} L:${r.stages.ldpcOk})`
    );
  }

  // Tuning
  console.log('\n=== Tuning Recommendations ===');
  const byIds = new Map<number, TestResult[]>();
  for (const r of results) { const arr = byIds.get(r.idsRate) || []; arr.push(r); byIds.set(r.idsRate, arr); }
  for (const [ids, matching] of byIds) {
    const best = matching.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, matching[0]);
    const avgVread = matching.reduce((s, r) => s + (r.stages.totalReads > 0 ? r.stages.viterbiPerReadOk / r.stages.totalReads : 0), 0) / matching.length;
    const mark = best.recoveryRate >= 0.95 ? '✓' : best.recoveryRate >= 0.80 ? '~' : '✗';
    console.log(`  IDS ${(ids * 100).toFixed(0)}%: ${mark} ${(best.recoveryRate * 100).toFixed(1)}% @ ${best.coverage}× ${best.ldpcParity}B i${best.insPen} d${best.delPen} (Vread:${(avgVread*100).toFixed(0)}%)`);
    if (avgVread < 0.5) {
      console.log(`    → Viterbi per-read failing — try lower ins/del penalties or increase coverage`);
    }
  }

  // Save
  const outputDir = path.join(process.cwd(), 'test-data');
  if (!fs.existsSync(outputDir)) fs.mkdirSync(outputDir, { recursive: true });
  const outputPath = path.join(outputDir, 'cascade-v2-results.json');
  fs.writeFileSync(outputPath, JSON.stringify({ wasmViterbi: wasmOk, results, timestamp: new Date().toISOString() }, null, 2));
  console.log(`\nResults saved to ${outputPath}`);
}

main().catch(err => { console.error('Failed:', err); process.exit(1); });
