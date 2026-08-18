/**
 * K=9 Viterbi Penalty Tuning — CORRECTED: bit-level noise simulation
 *
 * The Viterbi decoder operates on conv-encoded BIT STREAMS, not DNA.
 * We simulate noise at the bit level (substitutions + insertions + deletions)
 * which is what the indel-tolerant Viterbi is designed to correct.
 */
import * as path from "node:path";
import * as fs from "node:fs";

const addonPath = path.resolve(process.cwd(), "rust/helix-dna-napi/target/release/libhelix_dna_napi.so");
let nativeAddon: any = null;
try {
  const mod = { exports: {} };
  (process as any).dlopen(mod, addonPath);
  nativeAddon = mod.exports;
  console.log(`[tune] Native Viterbi: ${nativeAddon.napiVersion()}`);
} catch {
  console.log(`[tune] Native Viterbi NOT available`);
  process.exit(1);
}

import { LDPCInnerCode, getCachedLDPCInner } from "../src/lib/dna/ldpc-codec";
import { crc16Bytes } from "../src/lib/dna/crc16";

class Rng {
  private s: number;
  constructor(seed: number) { this.s = (seed >>> 0) || 1; }
  next(): number { this.s ^= this.s << 13; this.s ^= this.s >>> 17; this.s ^= this.s << 5; this.s = this.s >>> 0; return this.s / 0x100000000; }
  nextInt(max: number): number { return Math.floor(this.next() * max); }
}

/** Apply bit-level IDS noise to a byte stream */
function applyBitNoise(data: Uint8Array, subRate: number, insRate: number, delRate: number, rng: Rng): Uint8Array {
  const bits: number[] = [];
  // Unpack to bits
  for (let i = 0; i < data.length * 8; i++) bits.push((data[i >> 3] >> (7 - (i & 7))) & 1);

  const result: number[] = [];
  let pos = 0;
  while (pos < bits.length) {
    // Deletion: skip this bit
    if (rng.next() < delRate) { pos++; continue; }

    // Substitution: flip this bit
    if (rng.next() < subRate) {
      result.push(bits[pos] ^ 1);
    } else {
      result.push(bits[pos]);
    }
    pos++;

    // Insertion: insert a random bit after this one
    if (rng.next() < insRate) {
      result.push(rng.nextInt(2));
    }
  }

  // Pack back to bytes
  const out = new Uint8Array(Math.ceil(result.length / 8));
  for (let i = 0; i < result.length; i++) out[i >> 3] |= result[i] << (7 - (i & 7));
  return out;
}

interface TuneResult {
  insPen: number; delPen: number; maxDrift: number;
  ldpcParity: number; idsRate: number; coverage: number;
  useK7: boolean; recovered: number; totalOligos: number;
  recoveryRate: number; avgVitMs: number;
}

function runTest(
  n: number, payloadBytes: number, idsRate: number, coverage: number,
  ldpcParity: number, insPen: number, delPen: number, maxDrift: number,
  useK7: boolean, rng: Rng,
): TuneResult {
  // IDS composition: 45% deletion, 30% insertion, 25% substitution (ONT R10.4.1)
  const delR = idsRate * 0.45;
  const insR = idsRate * 0.30;
  const subR = idsRate * 0.25;
  const innerBytes = payloadBytes + ldpcParity + 2; // +CRC-16

  let ldpc: LDPCInnerCode | null = null;
  try { ldpc = getCachedLDPCInner(payloadBytes + ldpcParity, payloadBytes); } catch {}

  let recovered = 0, totalVitMs = 0;

  for (let i = 0; i < n; i++) {
    const payload = new Uint8Array(payloadBytes);
    for (let j = 0; j < payloadBytes; j++) payload[j] = rng.nextInt(256);

    // LDPC encode
    let ldpcCW = payload;
    if (ldpc) try { ldpcCW = ldpc.encode(payload); } catch {}

    // CRC-16
    const withCrc = new Uint8Array(ldpcCW.length + 2);
    withCrc.set(ldpcCW, 0);
    const crc = crc16Bytes(ldpcCW);
    withCrc[ldpcCW.length] = crc[0]; withCrc[ldpcCW.length + 1] = crc[1];

    // Convolutional encode
    const convOut = useK7
      ? new Uint8Array(nativeAddon.convK7Encode(withCrc))
      : new Uint8Array(nativeAddon.convK9Encode(withCrc));

    // Simulate noisy reads + consensus (bit-level)
    // Each read gets independent noise; consensus = majority vote per bit
    const readBits: number[][] = [];
    for (let r = 0; r < coverage; r++) {
      const noisy = applyBitNoise(convOut, subR, insR, delR, rng);
      const bits: number[] = [];
      for (let b = 0; b < noisy.length * 8; b++) bits.push((noisy[b >> 3] >> (7 - (b & 7))) & 1);
      readBits.push(bits);
    }

    // Majority-vote consensus at bit level
    const consensusLen = Math.round(readBits.reduce((s, r) => s + r.length, 0) / readBits.length);
    const consBits: number[] = [];
    for (let pos = 0; pos < consensusLen; pos++) {
      let ones = 0, zeros = 0;
      for (const r of readBits) {
        if (pos < r.length) {
          if (r[pos] === 1) ones++; else zeros++;
        }
      }
      consBits.push(ones >= zeros ? 1 : 0);
    }
    const consBytes = new Uint8Array(Math.ceil(consBits.length / 8));
    for (let b = 0; b < consBits.length; b++) consBytes[b >> 3] |= consBits[b] << (7 - (b & 7));

    // Viterbi decode
    const t0 = Date.now();
    let afterConv: Uint8Array;
    try {
      const cfg = { maxDrift, insertionPenalty: insPen, deletionPenalty: delPen, numInfoBits: innerBytes * 8 };
      const dec = useK7 ? nativeAddon.viterbiK7Decode(Buffer.from(consBytes), cfg) : nativeAddon.viterbiK9Decode(Buffer.from(consBytes), cfg);
      afterConv = new Uint8Array(dec);
    } catch { afterConv = consBytes.slice(0, innerBytes); }
    totalVitMs += Date.now() - t0;

    // LDPC decode
    let decoded: Uint8Array | null = null;
    if (ldpc && afterConv.length >= payloadBytes + ldpcParity) {
      try { const { data: d } = ldpc.decode(afterConv.slice(0, payloadBytes + ldpcParity)); if (d.length === payloadBytes) decoded = d; } catch {}
    }
    // CRC fallback
    if (!decoded && afterConv.length >= innerBytes) {
      const dp = afterConv.slice(0, afterConv.length - 2);
      const c = crc16Bytes(dp);
      if (afterConv[afterConv.length - 2] === c[0] && afterConv[afterConv.length - 1] === c[1]) decoded = afterConv.slice(0, payloadBytes);
    }
    if (decoded) { let ok = true; for (let b = 0; b < payloadBytes; b++) if (decoded[b] !== payload[b]) { ok = false; break; } if (ok) recovered++; }
  }
  return { insPen, delPen, maxDrift, ldpcParity, idsRate, coverage, useK7, recovered, totalOligos: n, recoveryRate: recovered / n, avgVitMs: totalVitMs / n };
}

async function main() {
  const rng = new Rng(42);
  const N = 20, payloadBytes = 30;
  const allResults: TuneResult[] = [];

  console.log("╔════════════════════════════════════════════════════════════╗");
  console.log("║   K=9 Viterbi Penalty Tuning (bit-level noise)            ║");
  console.log("╚════════════════════════════════════════════════════════════╝\n");

  // Phase 1: Penalty sweep at 9% IDS
  console.log("=== Phase 1: Penalty Sweep at 9% IDS, 20× coverage, 8B LDPC ===\n");
  const penalties = [
    { ins: 1.0, del: 1.0 }, { ins: 1.2, del: 1.2 }, { ins: 1.5, del: 1.5 },
    { ins: 1.8, del: 1.8 }, { ins: 2.0, del: 2.0 }, { ins: 1.2, del: 1.5 },
  ];

  let bestRate = 0, bestPen = { ins: 1.5, del: 1.5 }, bestMd = 10, bestK7 = false;

  for (const p of penalties) {
    for (const md of [5, 8, 10]) {
      for (const k7 of [false, true]) {
        const r = runTest(N, payloadBytes, 0.09, 20, 8, p.ins, p.del, md, k7, rng);
        const label = `${k7 ? 'K7' : 'K9'} ins=${p.ins} del=${p.del} md=${md}`;
        console.log(`  ${label}: ${(r.recoveryRate*100).toFixed(0).padStart(3)}% vit=${r.avgVitMs.toFixed(1)}ms`);
        allResults.push(r);
        if (r.recoveryRate > bestRate || (r.recoveryRate === bestRate && r.avgVitMs < (allResults.find(x => x.recoveryRate === bestRate)?.avgVitMs ?? Infinity))) {
          bestRate = r.recoveryRate; bestPen = p; bestMd = md; bestK7 = k7;
        }
      }
    }
  }

  console.log(`\n  Best: ${bestK7 ? 'K7' : 'K9'} ins=${bestPen.ins} del=${bestPen.del} md=${bestMd} → ${(bestRate*100).toFixed(0)}%`);

  // Phase 2: Coverage × LDPC sweep with best penalties
  console.log("\n=== Phase 2: IDS × Coverage × LDPC Parity ===\n");
  for (const ids of [0.04, 0.06, 0.09, 0.12]) {
    for (const cov of [10, 20, 30]) {
      for (const par of [4, 8, 10]) {
        const r = runTest(N, payloadBytes, ids, cov, par, bestPen.ins, bestPen.del, bestMd, bestK7, rng);
        const label = `IDS${(ids*100).toFixed(0)}% ${cov}× ${par}B`;
        const status = r.recoveryRate >= 0.95 ? "✓" : r.recoveryRate >= 0.80 ? "~" : "✗";
        console.log(`  ${status} ${label}: ${(r.recoveryRate*100).toFixed(0).padStart(3)}% vit=${r.avgVitMs.toFixed(1)}ms`);
        allResults.push(r);
      }
    }
  }

  // Summary
  console.log("\n" + "=".repeat(60));
  console.log("  TUNING SUMMARY");
  console.log("=".repeat(60));
  console.log(`  Best config: ${bestK7 ? 'K7' : 'K9'} ins=${bestPen.ins} del=${bestPen.del} maxDrift=${bestMd}`);
  console.log(`  At 9% IDS, 20× cov, 8B LDPC: ${(bestRate*100).toFixed(0)}% recovery`);

  // Per-IDS summary
  console.log("\n  Recommended LDPC parity per IDS rate:");
  for (const ids of [0.04, 0.06, 0.09, 0.12]) {
    const matching = allResults.filter(r => Math.abs(r.idsRate - ids) < 0.001 && r.coverage === 20);
    const best = matching.reduce((a, b) => a.recoveryRate > b.recoveryRate ? a : b, matching[0]);
    if (best) {
      const status = best.recoveryRate >= 0.95 ? "✓" : best.recoveryRate >= 0.80 ? "~" : "✗";
      console.log(`    ${(ids*100).toFixed(0)}%: ${status} ${(best.recoveryRate*100).toFixed(0)}% @ ${best.ldpcParity}B LDPC, ${best.coverage}× cov`);
    }
  }

  // Save
  const outDir = path.join(process.cwd(), "test-data");
  if (!fs.existsSync(outDir)) fs.mkdirSync(outDir, { recursive: true });
  fs.writeFileSync(path.join(outDir, "k9-penalty-tuning-results.json"), JSON.stringify({
    results: allResults,
    best: { ...bestPen, maxDrift: bestMd, useK7: bestK7, rate: bestRate },
    timestamp: new Date().toISOString(),
  }, null, 2));
  console.log(`\n  Results saved to test-data/k9-penalty-tuning-results.json`);
}

main().catch(e => { console.error(e); process.exit(1); });
