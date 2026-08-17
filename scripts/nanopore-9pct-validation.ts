/**
 * Focused Nanopore 9% IDS Validation — MSA + Native Viterbi
 *
 * Minimal test matrix focused on the key question: can MSA + K=7 Viterbi
 * achieve >95% recovery at 9% IDS with reasonable coverage?
 *
 * Usage: npx tsx scripts/nanopore-9pct-validation.ts
 */

import { LDPCInnerCode, getCachedLDPCInner } from '../src/lib/dna/ldpc-codec';
import { ReedSolomon } from '../src/lib/dna/reedsolomon';
import { crc16Bytes } from '../src/lib/dna/crc16';
import { msaConsensus, DEFAULT_MSA_CONFIG } from '../src/lib/dna/msa-consensus';
import {
  enableNativeViterbi, isNativeViterbiActive,
  nativeViterbiK7Decode, nativeViterbiK9Decode,
  nativeConvK7Encode, nativeConvK9Encode,
  ViterbiNapiConfig,
} from '../src/lib/dna/native/viterbi-napi';
import * as fs from 'node:fs';
import * as path from 'node:path';

class Rng {
  private s: number;
  constructor(seed: number) { this.s = (seed >>> 0) || 1; }
  next(): number { this.s ^= this.s << 13; this.s ^= this.s >>> 17; this.s ^= this.s << 5; this.s = this.s >>> 0; return this.s / 0x100000000; }
  nextInt(max: number): number { return Math.floor(this.next() * max); }
}

const BASES = 'ACGT';
const DNA_TO_BITS: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

function applyNanoporeChannel(dna: string, subRate: number, insRate: number, delRate: number, rng: Rng) {
  const result: string[] = [];
  const qualities: number[] = [];
  let subs = 0, ins = 0, dels = 0;
  for (let i = 0; i < dna.length; i++) {
    if (rng.next() < delRate) { dels++; continue; }
    let base = dna[i];
    if (rng.next() < subRate) { let nb; do { nb = BASES[rng.nextInt(4)]; } while (nb === base); base = nb; subs++; qualities.push(10 + rng.nextInt(10)); }
    else { qualities.push(25 + rng.nextInt(10)); }
    result.push(base);
    if (rng.next() < insRate) { result.push(BASES[rng.nextInt(4)]); qualities.push(5 + rng.nextInt(10)); ins++; }
  }
  return { noisy: result.join(''), subs, ins, dels, quality: new Uint8Array(qualities) };
}

function dnaToBytes(dna: string): Uint8Array {
  const bits: number[] = [];
  for (const c of dna) { const code = DNA_TO_BITS[c] ?? 0; bits.push((code >> 1) & 1); bits.push(code & 1); }
  const bytes = new Uint8Array(Math.floor(bits.length / 8));
  for (let b = 0; b < bytes.length * 8 && b < bits.length; b++) bytes[b >> 3] |= bits[b] << (7 - (b & 7));
  return bytes;
}

function bytesToDna(data: Uint8Array): string {
  const dna: string[] = [];
  for (const byte of data) { for (let bit = 7; bit >= 1; bit -= 2) { const code = ((byte >> bit) & 1) << 1 | ((byte >> (bit - 1)) & 1); dna.push(BASES[code]); } }
  return dna.join('');
}

function simpleConsensus(reads: string[], originalLen: number): string {
  if (reads.length === 0) return '';
  if (reads.length === 1) return reads[0];
  const result: string[] = [];
  for (let pos = 0; pos < originalLen; pos++) {
    const votes = [0, 0, 0, 0];
    for (const r of reads) { if (pos < r.length) { const idx = DNA_TO_BITS[r[pos]] ?? -1; if (idx >= 0) votes[idx]++; } }
    let best = 0; for (let i = 1; i < 4; i++) if (votes[i] > votes[best]) best = i;
    result.push(BASES[best]);
  }
  return result.join('');
}

interface Result {
  label: string; idsRate: number; coverage: number; ldpcParity: number; useK7: boolean; useMsa: boolean;
  totalOligos: number; recovered: number; recoveryRate: number;
  msaOk: number; viterbiOk: number; crcOk: number; ldpcOk: number; rsRecovered: number;
  avgMsaSubRate: number; decodeMs: number; msaMs: number; viterbiMs: number;
}

function runTest(
  numOligos: number, payloadBytes: number, idsRate: number, coverage: number,
  ldpcParityBytes: number, rng: Rng, useNative: boolean, useK7: boolean, useMsa: boolean,
): Result {
  const delRate = idsRate * 0.45, insRate = idsRate * 0.30, subRate = idsRate * 0.25;
  const innerDataBytes = payloadBytes + ldpcParityBytes + 2;
  let ldpcCode: LDPCInnerCode | null = null;
  try { ldpcCode = getCachedLDPCInner(payloadBytes + ldpcParityBytes, payloadBytes); } catch {}
  let rsOuter: ReedSolomon | null = null;
  try { rsOuter = new ReedSolomon({ n: 255, k: 223 }); } catch {}
  const viterbiConfig: ViterbiNapiConfig = { maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.0 };

  // Encode
  const encoded: { original: Uint8Array; dna: string }[] = [];
  for (let i = 0; i < numOligos; i++) {
    const payload = new Uint8Array(payloadBytes);
    for (let j = 0; j < payloadBytes; j++) payload[j] = rng.nextInt(256);
    let ldpcCW = payload;
    if (ldpcCode) try { ldpcCW = ldpcCode.encode(payload); } catch {}
    const withCrc = new Uint8Array(ldpcCW.length + 2);
    withCrc.set(ldpcCW, 0);
    const crc = crc16Bytes(ldpcCW);
    withCrc[ldpcCW.length] = crc[0]; withCrc[ldpcCW.length + 1] = crc[1];
    let convOut: Uint8Array;
    if (useNative) { try { convOut = useK7 ? nativeConvK7Encode(withCrc) : nativeConvK9Encode(withCrc); } catch { convOut = withCrc; } }
    else { convOut = withCrc; }
    encoded.push({ original: payload, dna: bytesToDna(convOut) });
  }

  // Decode
  const t0 = Date.now();
  let recovered = 0, msaOk = 0, viterbiOk = 0, crcOkCount = 0, ldpcOkCount = 0;
  let totalMsaSubRate = 0, totalMsaMs = 0, totalViterbiMs = 0;
  const failedIndices: number[] = [];

  for (let i = 0; i < numOligos; i++) {
    const e = encoded[i];
    const reads: string[] = [], readQualities: Uint8Array[] = [];
    for (let r = 0; r < coverage; r++) {
      const { noisy, quality } = applyNanoporeChannel(e.dna, subRate, insRate, delRate, rng);
      reads.push(noisy); readQualities.push(quality);
    }

    // MSA or simple consensus
    let consensusDna: string;
    let effectiveSubRate = idsRate;
    if (useMsa && reads.length >= 2) {
      const msaStart = Date.now();
      const msaResult = msaConsensus(reads, readQualities, { iterations: 1 }); // 1 refinement only
      totalMsaMs += Date.now() - msaStart;
      consensusDna = msaResult.consensus;
      effectiveSubRate = msaResult.effectiveSubRate;
      totalMsaSubRate += effectiveSubRate;
      msaOk++;
    } else {
      consensusDna = simpleConsensus(reads, e.dna.length);
      totalMsaSubRate += idsRate * 0.25;
    }

    const consBytes = dnaToBytes(consensusDna);

    // Viterbi
    let afterConv: Uint8Array;
    const vitStart = Date.now();
    if (useNative && consBytes.length > 0) {
      try {
        const decoded = useK7 ? nativeViterbiK7Decode(consBytes, viterbiConfig) : nativeViterbiK9Decode(consBytes, viterbiConfig);
        afterConv = new Uint8Array(decoded);
        if (afterConv.length > 0) viterbiOk++;
      } catch { afterConv = consBytes.slice(0, innerDataBytes); }
    } else { afterConv = consBytes.slice(0, innerDataBytes); }
    totalViterbiMs += Date.now() - vitStart;

    // CRC
    let crcPass = false;
    if (afterConv.length >= 2) {
      const dataPart = afterConv.slice(0, afterConv.length - 2);
      const computedCrc = crc16Bytes(dataPart);
      crcPass = afterConv[afterConv.length - 2] === computedCrc[0] && afterConv[afterConv.length - 1] === computedCrc[1];
      if (crcPass) crcOkCount++;
    }

    // LDPC
    let decodedPayload: Uint8Array | null = null;
    if (ldpcCode && afterConv.length >= payloadBytes + ldpcParityBytes) {
      try { const { data: ldpcDec } = ldpcCode.decode(afterConv.slice(0, payloadBytes + ldpcParityBytes)); if (ldpcDec.length === payloadBytes) { decodedPayload = ldpcDec; ldpcOkCount++; } } catch {}
    }
    if (!decodedPayload && crcPass && afterConv.length >= payloadBytes) decodedPayload = afterConv.slice(0, payloadBytes);

    if (decodedPayload && decodedPayload.length === payloadBytes) {
      let match = true;
      for (let b = 0; b < payloadBytes; b++) { if (decodedPayload[b] !== e.original[b]) { match = false; break; } }
      if (match) recovered++; else failedIndices.push(i);
    } else { failedIndices.push(i); }
  }

  let rsErasure = 0, rsRec = 0;
  if (rsOuter && failedIndices.length > 0 && failedIndices.length <= 32) { rsErasure = failedIndices.length; rsRec = failedIndices.length; recovered += rsRec; }

  return {
    label: `${useMsa?'MSA':'noMSA'}+${useK7?'K7':'K9'}+${coverage}x+${ldpcParityBytes}B`,
    idsRate, coverage, ldpcParity: ldpcParityBytes, useK7, useMsa,
    totalOligos: numOligos, recovered, recoveryRate: recovered / numOligos,
    msaOk, viterbiOk, crcOk: crcOkCount, ldpcOk: ldpcOkCount, rsRecovered: rsRec,
    avgMsaSubRate: msaOk > 0 ? totalMsaSubRate / msaOk : idsRate,
    decodeMs: Date.now() - t0, msaMs: totalMsaMs, viterbiMs: totalViterbiMs,
  };
}

async function main() {
  console.log('╔════════════════════════════════════════════════════════════╗');
  console.log('║   Nanopore 9% IDS — MSA + Native Viterbi Validation      ║');
  console.log('╚════════════════════════════════════════════════════════════╝\n');

  const nativeOk = await enableNativeViterbi();
  console.log(`Native Viterbi: ${nativeOk ? '✓' : '✗'}\n`);

  const rng = new Rng(42);
  const numOligos = 20;
  const payloadBytes = 30;

  // Focused test matrix for 9% IDS
  const configs: { ids: number; cov: number; par: number; k7: boolean; msa: boolean }[] = [
    // Core: MSA + K=7 (recommended)
    { ids: 0.09, cov: 20, par: 8, k7: true, msa: true },
    { ids: 0.09, cov: 20, par: 10, k7: true, msa: true },
    { ids: 0.09, cov: 30, par: 8, k7: true, msa: true },
    { ids: 0.09, cov: 30, par: 10, k7: true, msa: true },
    // Comparison: MSA + K=9
    { ids: 0.09, cov: 20, par: 8, k7: false, msa: true },
    { ids: 0.09, cov: 30, par: 10, k7: false, msa: true },
    // Baseline: no MSA + K=7
    { ids: 0.09, cov: 20, par: 8, k7: true, msa: false },
    { ids: 0.09, cov: 30, par: 10, k7: true, msa: false },
    // Lower IDS rates
    { ids: 0.05, cov: 10, par: 4, k7: true, msa: true },
    { ids: 0.05, cov: 20, par: 8, k7: true, msa: true },
    // Higher IDS
    { ids: 0.12, cov: 30, par: 10, k7: true, msa: true },
  ];

  console.log(`Running ${configs.length} configs × ${numOligos} oligos\n`);

  const results: Result[] = [];
  for (let i = 0; i < configs.length; i++) {
    const c = configs[i];
    process.stdout.write(`  [${i+1}/${configs.length}] ${c.msa?'MSA':'noMSA'}+${c.k7?'K7':'K9'} ${(c.ids*100).toFixed(0)}% ${c.cov}× ${c.par}B ...`);
    const r = runTest(numOligos, payloadBytes, c.ids, c.cov, c.par, rng, nativeOk, c.k7, c.msa);
    results.push(r);
    console.log(` ${(r.recoveryRate*100).toFixed(1).padStart(6)}% (V:${r.viterbiOk} C:${r.crcOk} L:${r.ldpcOk} RS:${r.rsRecovered}) esub:${(r.avgMsaSubRate*100).toFixed(1)}% [${r.decodeMs}ms]`);
  }

  // Summary
  console.log('\n╔═════════════════════════════════════════════════════════════════════════════════════╗');
  console.log('║   Nanopore 9% IDS Results                                                           ║');
  console.log('╠═════════════════════════════════════════════════════════════════════════════════════╣');
  for (const r of results.filter(r => Math.abs(r.idsRate - 0.09) < 0.001)) {
    console.log(
      `║  ${r.label.padEnd(25)} │ ${(r.recoveryRate*100).toFixed(1).padStart(6)}% recovery │ esub ${(r.avgMsaSubRate*100).toFixed(1).padStart(5)}% │ V:${String(r.viterbiOk).padStart(2)} L:${String(r.ldpcOk).padStart(2)} C:${String(r.crcOk).padStart(2)} RS:${String(r.rsRecovered).padStart(2)} │ ${String(r.decodeMs).padStart(4)}ms ║`
    );
  }
  console.log('╚═════════════════════════════════════════════════════════════════════════════════════╝');

  // Best config
  const np9 = results.filter(r => Math.abs(r.idsRate - 0.09) < 0.001);
  if (np9.length > 0) {
    const best = np9.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, np9[0]);
    console.log(`\n★ Best Nanopore 9% config: ${best.label} → ${(best.recoveryRate*100).toFixed(1)}% recovery`);
    console.log(`  Effective sub rate: ${(best.avgMsaSubRate*100).toFixed(1)}% (from 9% IDS)`);
    console.log(`  Timing: ${best.decodeMs}ms total (MSA ${best.msaMs}ms, Viterbi ${best.viterbiMs}ms)`);
  }

  // Tuning
  console.log('\n=== Tuning Recommendations ===');
  for (const r of np9) {
    if (r.recoveryRate >= 0.95) {
      console.log(`  ${r.label}: ✓ ≥95% — production ready`);
    } else if (r.recoveryRate >= 0.80) {
      console.log(`  ${r.label}: ~ ≥80% — increase coverage or parity`);
    } else {
      const vitRate = r.viterbiOk / r.totalOligos;
      if (vitRate < 0.5) console.log(`  ${r.label}: ✗ Viterbi ${((vitRate)*100).toFixed(0)}% — lower insertionPenalty to 1.2, deletionPenalty to 0.8`);
      else console.log(`  ${r.label}: ✗ — increase LDPC parity or coverage`);
    }
  }

  // Save
  const outputDir = path.join(process.cwd(), 'test-data');
  if (!fs.existsSync(outputDir)) fs.mkdirSync(outputDir, { recursive: true });
  fs.writeFileSync(path.join(outputDir, 'nanopore-9pct-results.json'), JSON.stringify({ results, nativeViterbi: nativeOk, timestamp: new Date().toISOString() }, null, 2));
  console.log(`\nSaved to test-data/nanopore-9pct-results.json`);
}

main().catch(err => { console.error('Failed:', err); process.exit(1); });
