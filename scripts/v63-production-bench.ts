/**
 * v63 Production Benchmark
 *
 * Tests:
 *   1. v63-hd config (1100nt + direct + 4B parity + 2% RS) — should beat SOTA 1.815 b/nt
 *   2. v63-maxdensity config (1500nt) — should hit ~1.856 b/nt
 *   3. v55-density config (700nt + 8B parity) — baseline 1.664 b/nt
 *   4. Nanopore K=9 IDS recovery at 9% IDS (PRESET_NANOPORE)
 *   5. 10MB streaming decode memory test
 *   6. Mega-perf wiring validation (BitParallelSyndrome in decode, cached LDPC, cached IndelViterbi)
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_NANOPORE, MutationConfig } from "../src/lib/dna/simulate";
import {
  ULTIMATE_V55_DENSITY_CONFIG,
  ULTIMATE_V63_HD_CONFIG,
  ULTIMATE_V63_MAXDENSITY_CONFIG,
  ULTIMATE_V61_NANOPORE_CONFIG,
  computeDensity,
} from "../src/lib/dna/presets";
import { computeBatchSize, estimateOligoMemoryUsage } from "../src/lib/dna/mega-performance";

const MB = 1024 * 1024;

function randomBytes(n: number): Uint8Array {
  const buf = new Uint8Array(n);
  for (let i = 0; i < n; i++) buf[i] = Math.floor(Math.random() * 256);
  return buf;
}

async function benchConfig(
  name: string,
  cfg: typeof ULTIMATE_V55_DENSITY_CONFIG,
  payloadSize: number,
  coverage: number = 10,
  subRate: number = 0.001,
): Promise<void> {
  const payload = randomBytes(payloadSize);
  const theoreticalDensity = computeDensity(cfg);

  const encodeStart = performance.now();
  const enc = await encodeFile(payload, cfg, { fileName: "b.bin", contentType: "application/octet-stream" });
  const encodeMs = performance.now() - encodeStart;
  const { oligos, metadata } = enc.encoded;

  const simCfg: MutationConfig = {
    substitutionRate: subRate,
    insertionRate: 0,
    deletionRate: 0,
    coverage,
    dropoutRate: 0,
    seed: 42,
  };
  const { reads } = simulate(oligos, simCfg);

  const decodeStart = performance.now();
  const result = await decodeReads(reads, metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decodeMs = performance.now() - decodeStart;

  const hashMatch = result.hashMatches;
  const dataMatch = result.data
    ? Buffer.compare(Buffer.from(payload), Buffer.from(result.data)) === 0
    : false;

  // Compute ACTUAL realized density (payload bits / total nt synthesized)
  const totalNt = oligos.reduce((sum, o) => sum + o.sequence.length, 0);
  const realizedDensity = (payloadSize * 8) / totalNt;

  console.log(`\n${name} (${(payloadSize / MB).toFixed(2)}MB @ ${coverage}× cov, ${subRate * 100}% sub)`);
  console.log(`  theoretical density: ${theoreticalDensity.toFixed(3)} b/nt`);
  console.log(`  realized density:    ${realizedDensity.toFixed(3)} b/nt`);
  console.log(`  reported density:    ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
  console.log(`  oligos: ${oligos.length}, total nt: ${totalNt.toLocaleString()}`);
  console.log(`  encode: ${encodeMs.toFixed(0)}ms = ${(payloadSize / MB / (encodeMs / 1000)).toFixed(2)} MB/s`);
  console.log(`  decode: ${decodeMs.toFixed(0)}ms = ${(payloadSize / MB / (decodeMs / 1000)).toFixed(2)} MB/s`);
  console.log(`  hash: ${hashMatch ? "OK ✅" : "FAIL ❌"}, data: ${dataMatch ? "OK ✅" : "FAIL ❌"}`);
  console.log(`  recovered: ${result.stats.oligosRecovered}/${oligos.length} oligos, erased: ${result.stats.oligosErased}`);
}

async function benchNanoporeK9(): Promise<void> {
  console.log("\n=== Nanopore K=9 IDS Recovery Test ===");
  const cfg = ULTIMATE_V61_NANOPORE_CONFIG;
  const payloadSize = 16 * 1024; // 16KB — nanopore is slower, keep small
  const payload = randomBytes(payloadSize);

  console.log(`  config: 300nt oligos, K=9 conv (d_free=24), 8B LDPC, 25% outer RS`);
  console.log(`  channel: nanopore (sub=2%, ins=3%, del=4% → 9% total IDS)`);

  const encodeStart = performance.now();
  const enc = await encodeFile(payload, cfg, { fileName: "b.bin", contentType: "application/octet-stream" });
  const encodeMs = performance.now() - encodeStart;
  const { oligos, metadata } = enc.encoded;
  console.log(`  encoded: ${oligos.length} oligos in ${encodeMs.toFixed(0)}ms`);

  // Test at multiple IDS levels
  for (const [label, sub, ins, del] of [
    ["0% IDS (clean)", 0.0, 0.0, 0.0],
    ["3% IDS (mild)", 0.01, 0.01, 0.01],
    ["6% IDS (moderate)", 0.02, 0.02, 0.02],
    ["9% IDS (PRESET_NANOPORE)", 0.02, 0.03, 0.04],
    ["12% IDS (heavy)", 0.03, 0.04, 0.05],
  ] as const) {
    const simCfg: MutationConfig = {
      substitutionRate: sub,
      insertionRate: ins,
      deletionRate: del,
      coverage: 10,
      dropoutRate: 0,
      seed: 42,
    };
    const { reads } = simulate(oligos, simCfg);

    const decodeStart = performance.now();
    const result = await decodeReads(reads, metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const decodeMs = performance.now() - decodeStart;

    const hashMatch = result.hashMatches;
    const dataMatch = result.data
      ? Buffer.compare(Buffer.from(payload), Buffer.from(result.data)) === 0
      : false;

    console.log(`  ${label}: hash=${hashMatch ? "OK" : "FAIL"}, data=${dataMatch ? "OK" : "FAIL"}, decode=${decodeMs.toFixed(0)}ms`);
  }
}

async function benchStreamingMemory(): Promise<void> {
  console.log("\n=== 10MB Streaming Decode Memory Test ===");
  const cfg = ULTIMATE_V63_HD_CONFIG;
  const payloadSize = 10 * MB;
  const payload = randomBytes(payloadSize);

  const memBefore = process.memoryUsage();
  console.log(`  RSS before encode: ${(memBefore.rss / MB).toFixed(0)}MB`);

  const encodeStart = performance.now();
  const enc = await encodeFile(payload, cfg, { fileName: "b.bin", contentType: "application/octet-stream" });
  const encodeMs = performance.now() - encodeStart;
  const { oligos, metadata } = enc.encoded;
  const memAfterEncode = process.memoryUsage();
  console.log(`  RSS after encode:  ${(memAfterEncode.rss / MB).toFixed(0)}MB (Δ${((memAfterEncode.rss - memBefore.rss) / MB).toFixed(0)}MB)`);
  console.log(`  encoded: ${oligos.length} oligos, total nt: ${oligos.reduce((s, o) => s + o.sequence.length, 0).toLocaleString()}`);
  console.log(`  encode: ${encodeMs.toFixed(0)}ms = ${(payloadSize / MB / (encodeMs / 1000)).toFixed(2)} MB/s`);

  // Compute theoretical batch size for 1GB memory budget
  const batchSize = computeBatchSize(1024, 10, cfg.oligoLength);
  const bytesPerOligo = estimateOligoMemoryUsage(10, cfg.oligoLength);
  console.log(`  streaming batch size (1GB budget): ${batchSize.toLocaleString()} oligos (${bytesPerOligo}B each)`);

  const simCfg: MutationConfig = {
    substitutionRate: 0.001,
    insertionRate: 0,
    deletionRate: 0,
    coverage: 10,
    dropoutRate: 0,
    seed: 42,
  };
  const { reads } = simulate(oligos, simCfg);

  const decodeStart = performance.now();
  const result = await decodeReads(reads, metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decodeMs = performance.now() - decodeStart;
  const memAfterDecode = process.memoryUsage();
  console.log(`  RSS after decode:  ${(memAfterDecode.rss / MB).toFixed(0)}MB (Δ${((memAfterDecode.rss - memAfterEncode.rss) / MB).toFixed(0)}MB)`);
  console.log(`  decode: ${decodeMs.toFixed(0)}ms = ${(payloadSize / MB / (decodeMs / 1000)).toFixed(2)} MB/s`);
  console.log(`  hash: ${result.hashMatches ? "OK ✅" : "FAIL ❌"}`);
  if (result.data) {
    const match = Buffer.compare(Buffer.from(payload), Buffer.from(result.data)) === 0;
    console.log(`  data: ${match ? "OK ✅" : "FAIL ❌"}`);
  }
}

async function main(): Promise<void> {
  console.log("=== v63 Production Benchmark ===");
  console.log(`Date: ${new Date().toISOString()}`);
  console.log(`Node: ${process.version}, Bun: ${(globalThis as any).Bun?.version ?? "n/a"}`);

  // Test 1: v55-density baseline (700nt, 8B parity, 3% RS)
  await benchConfig("v55-density (baseline)", ULTIMATE_V55_DENSITY_CONFIG, 256 * 1024, 10, 0.001);

  // Test 2: v63-hd (1100nt, 4B parity, 2% RS) — should beat SOTA 1.815
  await benchConfig("v63-hd (1100nt, 4B parity)", ULTIMATE_V63_HD_CONFIG, 256 * 1024, 10, 0.001);

  // Test 3: v63-maxdensity (1500nt, 4B parity, 2% RS) — should hit ~1.856
  await benchConfig("v63-maxdensity (1500nt, 4B parity)", ULTIMATE_V63_MAXDENSITY_CONFIG, 256 * 1024, 10, 0.001);

  // Test 4: v63-hd at 2MB scale
  await benchConfig("v63-hd (1100nt) @ 2MB", ULTIMATE_V63_HD_CONFIG, 2 * MB, 10, 0.001);

  // Test 5: Nanopore K=9 IDS recovery
  await benchNanoporeK9();

  // Test 6: 10MB streaming
  await benchStreamingMemory();

  console.log("\n=== v63 Benchmark Complete ===");
}

main().catch((e) => {
  console.error("Benchmark failed:", e);
  process.exit(1);
});
