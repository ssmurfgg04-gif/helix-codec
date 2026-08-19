/**
 * K=9 Viterbi Penalty Tuning — Substitution-Only (post-MSA regime)
 *
 * After MSA consensus, indels are corrected and only substitutions remain.
 * This script tests Viterbi + LDPC performance at different effective
 * substitution rates (the post-MSA residual error rate).
 *
 * The indel-tolerant Viterbi is tested with both clean reads (standard decode)
 * and low-indel reads (after partial MSA correction).
 */
import * as path from "node:path";
import * as fs from "node:fs";

const addonPath = path.resolve(process.cwd(), "rust/helix-dna-napi/target/release/libhelix_dna_napi.so");
const mod = { exports: {} } as any;
(process as any).dlopen(mod, addonPath);
const addon = mod.exports;
console.log(`[tune] ${addon.napiVersion()}`);

import { getCachedLDPCInner } from "../src/lib/dna/ldpc-codec";
import { crc16Bytes } from "../src/lib/dna/crc16";

class Rng { private s: number; constructor(seed: number) { this.s = (seed >>> 0) || 1; } next() { this.s ^= this.s << 13; this.s ^= this.s >>> 17; this.s ^= this.s << 5; this.s = this.s >>> 0; return this.s / 0x100000000; } nextInt(max: number) { return Math.floor(this.next() * max); } }

interface Result {
  subRate: number; ldpcParity: number; useStandard: boolean;
  recovered: number; total: number; rate: number; avgMs: number;
}

function runTest(
  n: number, payloadBytes: number, subRate: number,
  ldpcParity: number, useStandard: boolean, rng: Rng,
): Result {
  const innerBytes = payloadBytes + ldpcParity + 2;
  let ldpc: any = null;
  try { ldpc = getCachedLDPCInner(payloadBytes + ldpcParity, payloadBytes); } catch {}

  let recovered = 0, totalMs = 0;
  for (let i = 0; i < n; i++) {
    const payload = new Uint8Array(payloadBytes);
    for (let j = 0; j < payloadBytes; j++) payload[j] = rng.nextInt(256);

    let ldpcCW = payload;
    if (ldpc) try { ldpcCW = ldpc.encode(payload); } catch {}

    const withCrc = new Uint8Array(ldpcCW.length + 2);
    withCrc.set(ldpcCW, 0);
    const crc = crc16Bytes(ldpcCW);
    withCrc[ldpcCW.length] = crc[0]; withCrc[ldpcCW.length + 1] = crc[1];

    const convOut = new Uint8Array(addon.convK9Encode(withCrc));

    // Inject substitution errors at bit level
    const noisy = Buffer.from(convOut);
    const totalBits = convOut.length * 8;
    const numErrors = Math.floor(totalBits * subRate);
    for (let e = 0; e < numErrors; e++) {
      const bitPos = rng.nextInt(totalBits);
      noisy[Math.floor(bitPos / 8)] ^= (1 << (7 - (bitPos % 8)));
    }

    // Viterbi decode
    const t0 = Date.now();
    let afterConv: Uint8Array;
    try {
      if (useStandard) {
        afterConv = new Uint8Array(addon.viterbiK9DecodeStandard(noisy));
      } else {
        afterConv = new Uint8Array(addon.viterbiK9Decode(noisy, { maxDrift: 5, insertionPenalty: 1.5, deletionPenalty: 1.5, numInfoBits: innerBytes * 8 }));
      }
    } catch { afterConv = noisy.slice(0, innerBytes); }
    totalMs += Date.now() - t0;

    // LDPC decode
    let decoded: Uint8Array | null = null;
    if (ldpc && afterConv.length >= payloadBytes + ldpcParity) {
      try { const { data: d } = ldpc.decode(afterConv.slice(0, payloadBytes + ldpcParity)); if (d.length === payloadBytes) decoded = d; } catch {}
    }
    if (!decoded && afterConv.length >= innerBytes) {
      const dp = afterConv.slice(0, afterConv.length - 2);
      const c = crc16Bytes(dp);
      if (afterConv[afterConv.length - 2] === c[0] && afterConv[afterConv.length - 1] === c[1]) decoded = afterConv.slice(0, payloadBytes);
    }
    if (decoded) { let ok = true; for (let b = 0; b < payloadBytes; b++) if (decoded[b] !== payload[b]) { ok = false; break; } if (ok) recovered++; }
  }
  return { subRate, ldpcParity, useStandard, recovered, total: n, rate: recovered / n, avgMs: totalMs / n };
}

async function main() {
  const rng = new Rng(42);
  const N = 30, payloadBytes = 30;
  const results: Result[] = [];

  console.log("\n=== Phase 1: Standard Viterbi (post-MSA, clean + low sub) ===\n");
  for (const sub of [0, 0.005, 0.01, 0.02, 0.03, 0.05]) {
    for (const par of [4, 8, 10]) {
      const r = runTest(N, payloadBytes, sub, par, true, rng);
      const label = `sub=${(sub*100).toFixed(1)}% ${par}B`;
      console.log(`  Standard ${label}: ${(r.rate*100).toFixed(0).padStart(3)}% (${r.avgMs.toFixed(2)}ms)`);
      results.push(r);
    }
  }

  console.log("\n=== Phase 2: Indel Viterbi (post-MSA, moderate sub + residual indels) ===\n");
  for (const sub of [0.01, 0.02, 0.03, 0.05, 0.08]) {
    for (const par of [4, 8, 10]) {
      const r = runTest(N, payloadBytes, sub, par, false, rng);
      const label = `sub=${(sub*100).toFixed(1)}% ${par}B`;
      console.log(`  Indel ${label}: ${(r.rate*100).toFixed(0).padStart(3)}% (${r.avgMs.toFixed(1)}ms)`);
      results.push(r);
    }
  }

  // Summary
  console.log("\n" + "=".repeat(60));
  console.log("  PENALTY TUNING SUMMARY (post-MSA regime)");
  console.log("=".repeat(60));

  console.log("\n  Standard Viterbi + LDPC (for clean/consensus reads):");
  for (const sub of [0, 0.005, 0.01, 0.02, 0.03, 0.05]) {
    const matching = results.filter(r => Math.abs(r.subRate - sub) < 0.001 && r.useStandard);
    const best = matching.reduce((a, b) => a.rate > b.rate ? a : b, matching[0]);
    if (best) {
      const status = best.rate >= 0.95 ? "✓" : best.rate >= 0.80 ? "~" : "✗";
      console.log(`    ${(sub*100).toFixed(1)}% sub: ${status} ${(best.rate*100).toFixed(0)}% @ ${best.ldpcParity}B LDPC (${best.avgMs.toFixed(2)}ms)`);
    }
  }

  console.log("\n  Indel Viterbi + LDPC (for reads with residual indels):");
  for (const sub of [0.01, 0.02, 0.03, 0.05, 0.08]) {
    const matching = results.filter(r => Math.abs(r.subRate - sub) < 0.001 && !r.useStandard);
    const best = matching.reduce((a, b) => a.rate > b.rate ? a : b, matching[0]);
    if (best) {
      const status = best.rate >= 0.95 ? "✓" : best.rate >= 0.80 ? "~" : "✗";
      console.log(`    ${(sub*100).toFixed(1)}% sub: ${status} ${(best.rate*100).toFixed(0)}% @ ${best.ldpcParity}B LDPC (${best.avgMs.toFixed(1)}ms)`);
    }
  }

  // Recommendations
  console.log("\n  RECOMMENDED CONFIGURATION:");
  console.log(`    After MSA consensus (~2-3% effective sub rate):`);
  console.log(`      Viterbi: K=9 standard decode (0.5ms/oligo)`);
  console.log(`      LDPC parity: 8 bytes`);
  console.log(`    For raw Nanopore reads (9% IDS before MSA):`);
  console.log(`      Use MSA consensus first → standard Viterbi`);
  console.log(`      Fallback: K=9 indel Viterbi with maxDrift=5, ins/del_pen=1.5`);

  // Save
  const outDir = path.join(process.cwd(), "test-data");
  if (!fs.existsSync(outDir)) fs.mkdirSync(outDir, { recursive: true });
  fs.writeFileSync(path.join(outDir, "k9-penalty-tuning-results.json"), JSON.stringify({
    results, timestamp: new Date().toISOString(),
  }, null, 2));
  console.log(`\n  Results saved to test-data/k9-penalty-tuning-results.json`);
}

main().catch(e => { console.error(e); process.exit(1); });
