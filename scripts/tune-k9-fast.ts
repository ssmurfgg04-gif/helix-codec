/**
 * Fast K=9 Penalty Tuning — uses native Viterbi, 10 oligos, focused grid.
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

import { ConvolutionalCode, bytesToBits } from "../src/lib/dna/convolutional";
import { LDPCInnerCode, getCachedLDPCInner } from "../src/lib/dna/ldpc-codec";
import { crc16Bytes } from "../src/lib/dna/crc16";

class Rng {
  private s: number;
  constructor(seed: number) { this.s = (seed >>> 0) || 1; }
  next(): number { this.s ^= this.s << 13; this.s ^= this.s >>> 17; this.s ^= this.s << 5; this.s = this.s >>> 0; return this.s / 0x100000000; }
  nextInt(max: number): number { return Math.floor(this.next() * max); }
}

const BASES = "ACGT";
const DNA_TO_BITS: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

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
  return dna.join("");
}

function applyChannel(dna: string, sub: number, ins: number, del: number, rng: Rng): string {
  const r: string[] = [];
  for (let i = 0; i < dna.length; i++) {
    if (rng.next() < del) continue;
    let base = dna[i];
    if (rng.next() < sub) { let nb; do { nb = BASES[rng.nextInt(4)]; } while (nb === base); base = nb; }
    r.push(base);
    if (rng.next() < ins) r.push(BASES[rng.nextInt(4)]);
  }
  return r.join("");
}

function consensus(reads: string[], len: number): string {
  if (reads.length <= 1) return reads[0]?.slice(0, len) ?? "";
  const r: string[] = [];
  for (let p = 0; p < len; p++) {
    const v = [0, 0, 0, 0];
    for (const rd of reads) { if (p < rd.length) { const i = DNA_TO_BITS[rd[p]] ?? -1; if (i >= 0) v[i]++; } }
    let b = 0; for (let i = 1; i < 4; i++) if (v[i] > v[b]) b = i;
    r.push(BASES[b]);
  }
  return r.join("");
}

function runTest(
  n: number, payloadBytes: number, idsRate: number, coverage: number,
  ldpcParity: number, insPen: number, delPen: number, maxDrift: number,
  useK7: boolean, rng: Rng,
): { recovered: number; rate: number; avgMs: number } {
  const delR = idsRate * 0.45, insR = idsRate * 0.30, subR = idsRate * 0.25;
  const innerBytes = payloadBytes + ldpcParity + 2;
  let ldpc: LDPCInnerCode | null = null;
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

    const convOut = useK7
      ? new Uint8Array(nativeAddon.convK7Encode(withCrc))
      : new Uint8Array(nativeAddon.convK9Encode(withCrc));

    const dna = bytesToDna(convOut);
    const reads: string[] = [];
    for (let r = 0; r < coverage; r++) reads.push(applyChannel(dna, subR, insR, delR, rng));
    const cons = consensus(reads, dna.length);
    const consBytes = dnaToBytes(cons);

    const t0 = Date.now();
    let afterConv: Uint8Array;
    try {
      const cfg = { maxDrift, insertionPenalty: insPen, deletionPenalty: delPen, numInfoBits: innerBytes * 8 };
      const dec = useK7 ? nativeAddon.viterbiK7Decode(Buffer.from(consBytes), cfg) : nativeAddon.viterbiK9Decode(Buffer.from(consBytes), cfg);
      afterConv = new Uint8Array(dec);
    } catch { afterConv = consBytes.slice(0, innerBytes); }
    totalMs += Date.now() - t0;

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
  return { recovered, rate: recovered / n, avgMs: totalMs / n };
}

async function main() {
  const rng = new Rng(42);
  const N = 10, payloadBytes = 30;
  const results: any[] = [];

  console.log("=== Phase 1: Penalty Sweep at 9% IDS, 20× coverage ===\n");
  const penalties = [
    { ins: 1.0, del: 1.0 }, { ins: 1.2, del: 1.2 }, { ins: 1.5, del: 1.5 },
    { ins: 1.8, del: 1.8 }, { ins: 2.0, del: 2.0 }, { ins: 1.2, del: 1.5 },
  ];

  let bestRate = 0, bestPen = { ins: 1.5, del: 1.5 }, bestMd = 10;

  for (const p of penalties) {
    for (const md of [5, 8, 10]) {
      for (const k7 of [false, true]) {
        const r = runTest(N, payloadBytes, 0.09, 20, 8, p.ins, p.del, md, k7, rng);
        const label = `${k7 ? 'K7' : 'K9'} ins=${p.ins} del=${p.del} md=${md}`;
        console.log(`  ${label}: ${(r.rate*100).toFixed(0).padStart(3)}% (${r.avgMs.toFixed(1)}ms)`);
        results.push({ ...r, ins: p.ins, del: p.del, md, k7, label });
        if (r.rate > bestRate) { bestRate = r.rate; bestPen = p; bestMd = md; }
      }
    }
  }

  console.log(`\n  Best: ins=${bestPen.ins} del=${bestPen.del} md=${bestMd} → ${(bestRate*100).toFixed(0)}%`);

  console.log("\n=== Phase 2: IDS × Coverage × LDPC Parity with Best Penalties ===\n");
  for (const ids of [0.04, 0.06, 0.09, 0.12]) {
    for (const cov of [10, 20]) {
      for (const par of [4, 8, 10]) {
        const r = runTest(N, payloadBytes, ids, cov, par, bestPen.ins, bestPen.del, bestMd, false, rng);
        const label = `IDS${(ids*100).toFixed(0)}% ${cov}× ${par}B`;
        console.log(`  ${label}: ${(r.rate*100).toFixed(0).padStart(3)}% (${r.avgMs.toFixed(1)}ms)`);
        results.push({ ...r, ids, cov, par, label });
      }
    }
  }

  // Save
  const outDir = path.join(process.cwd(), "test-data");
  if (!fs.existsSync(outDir)) fs.mkdirSync(outDir, { recursive: true });
  fs.writeFileSync(path.join(outDir, "k9-penalty-tuning-results.json"), JSON.stringify({
    results, best: { ...bestPen, maxDrift: bestMd, rate: bestRate }, timestamp: new Date().toISOString(),
  }, null, 2));
  console.log("\n  Results saved to test-data/k9-penalty-tuning-results.json");
}

main().catch(e => { console.error(e); process.exit(1); });
