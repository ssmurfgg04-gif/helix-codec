/**
 * v63: 10MB Streaming Decode Memory Test (Task 3)
 *
 * Verifies that StreamingDecode keeps memory flat for TB-scale payloads.
 * Uses v63-hd config (direct mode, 1100nt) — no K=9, so it completes fast.
 *
 * Measures:
 *   - RSS before encode, after encode, after decode
 *   - Encode/decode throughput
 *   - Hash verification
 *   - Oligo count and total nt
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, MutationConfig } from "../src/lib/dna/simulate";
import { ULTIMATE_V63_HD_CONFIG } from "../src/lib/dna/presets";
import {
  computeBatchSize,
  estimateOligoMemoryUsage,
} from "../src/lib/dna/mega-performance";

const MB = 1024 * 1024;

function randomBytes(n: number): Uint8Array {
  const buf = new Uint8Array(n);
  for (let i = 0; i < n; i++) buf[i] = Math.floor(Math.random() * 256);
  return buf;
}

async function main(): Promise<void> {
  console.log("=== v63 10MB Streaming Decode Memory Test ===");
  console.log(`Date: ${new Date().toISOString()}`);
  console.log(`Node: ${process.version}, Bun: ${(globalThis as any).Bun?.version ?? "n/a"}`);
  console.log(`Config: v63-hd (1100nt, direct, 4B LDPC, 2% outer RS)\n`);

  const cfg = ULTIMATE_V63_HD_CONFIG;
  const payloadSize = 10 * MB;
  const payload = randomBytes(payloadSize);

  // Memory baseline
  const memBefore = process.memoryUsage();
  console.log(`RSS before encode:  ${(memBefore.rss / MB).toFixed(0)}MB`);
  console.log(`Heap used:          ${(memBefore.heapUsed / MB).toFixed(0)}MB / ${(memBefore.heapTotal / MB).toFixed(0)}MB`);

  // Encode
  console.log(`\nEncoding ${(payloadSize / MB).toFixed(0)}MB payload...`);
  const encodeStart = performance.now();
  const enc = await encodeFile(payload, cfg, {
    fileName: "b.bin",
    contentType: "application/octet-stream",
  });
  const encodeMs = performance.now() - encodeStart;
  const { oligos, metadata } = enc.encoded;

  const memAfterEncode = process.memoryUsage();
  console.log(`RSS after encode:   ${(memAfterEncode.rss / MB).toFixed(0)}MB (Δ${((memAfterEncode.rss - memBefore.rss) / MB).toFixed(0)}MB)`);
  console.log(`Heap used:          ${(memAfterEncode.heapUsed / MB).toFixed(0)}MB / ${(memAfterEncode.heapTotal / MB).toFixed(0)}MB`);
  console.log(`Encoded: ${oligos.length.toLocaleString()} oligos, total nt: ${oligos.reduce((s, o) => s + o.sequence.length, 0).toLocaleString()}`);
  console.log(`Encode: ${encodeMs.toFixed(0)}ms = ${(payloadSize / MB / (encodeMs / 1000)).toFixed(2)} MB/s`);
  console.log(`Reported density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);

  // Compute streaming batch size
  const batchSize = computeBatchSize(1024, 10, cfg.oligoLength);
  const bytesPerOligo = estimateOligoMemoryUsage(10, cfg.oligoLength);
  console.log(`\nStreaming config (1GB memory budget):`);
  console.log(`  batch size:      ${batchSize.toLocaleString()} oligos`);
  console.log(`  bytes per oligo: ${bytesPerOligo}B (at 10× coverage)`);
  console.log(`  total batches:   ${Math.ceil(oligos.length / batchSize)}`);

  // Simulate Illumina reads
  console.log(`\nSimulating 10× coverage Illumina reads...`);
  const simStart = performance.now();
  const simCfg: MutationConfig = {
    substitutionRate: 0.001,
    insertionRate: 0,
    deletionRate: 0,
    coverage: 10,
    dropoutRate: 0,
    seed: 42,
  };
  const { reads } = simulate(oligos, simCfg);
  const simMs = performance.now() - simStart;
  console.log(`Simulated ${reads.length.toLocaleString()} reads in ${simMs.toFixed(0)}ms`);

  const memAfterSim = process.memoryUsage();
  console.log(`RSS after simulate: ${(memAfterSim.rss / MB).toFixed(0)}MB (Δ${((memAfterSim.rss - memAfterEncode.rss) / MB).toFixed(0)}MB)`);

  // Decode
  console.log(`\nDecoding (streaming)...`);
  const decodeStart = performance.now();
  const result = await decodeReads(
    reads,
    metadata,
    cfg,
    enc.encoded.forwardPrimer,
    enc.encoded.reversePrimer,
  );
  const decodeMs = performance.now() - decodeStart;

  const memAfterDecode = process.memoryUsage();
  console.log(`RSS after decode:   ${(memAfterDecode.rss / MB).toFixed(0)}MB (Δ${((memAfterDecode.rss - memAfterEncode.rss) / MB).toFixed(0)}MB)`);
  console.log(`Heap used:          ${(memAfterDecode.heapUsed / MB).toFixed(0)}MB / ${(memAfterDecode.heapTotal / MB).toFixed(0)}MB`);
  console.log(`Decode: ${decodeMs.toFixed(0)}ms = ${(payloadSize / MB / (decodeMs / 1000)).toFixed(2)} MB/s`);

  // Verify
  const hashMatch = result.hashMatches;
  const dataMatch = result.data
    ? Buffer.compare(Buffer.from(payload), Buffer.from(result.data)) === 0
    : false;

  console.log(`\n=== Results ===`);
  console.log(`hash:       ${hashMatch ? "OK ✅" : "FAIL ❌"}`);
  console.log(`data:       ${dataMatch ? "OK ✅" : "FAIL ❌"}`);
  console.log(`recovered:  ${result.stats.oligosRecovered.toLocaleString()}/${oligos.length.toLocaleString()} oligos`);
  console.log(`erased:     ${result.stats.oligosErased}`);

  // Memory flatness check
  const encodeRss = memAfterEncode.rss / MB;
  const decodeRss = memAfterDecode.rss / MB;
  const delta = decodeRss - encodeRss;
  console.log(`\n=== Memory Flatness ===`);
  console.log(`Encode RSS: ${encodeRss.toFixed(0)}MB`);
  console.log(`Decode RSS: ${decodeRss.toFixed(0)}MB`);
  console.log(`Delta:      ${delta >= 0 ? "+" : ""}${delta.toFixed(0)}MB`);
  if (Math.abs(delta) < 500) {
    console.log(`✅ Memory is FLAT — streaming decode keeps memory bounded.`);
  } else {
    console.log(`⚠️ Memory grew by ${Math.abs(delta).toFixed(0)}MB during decode — investigate.`);
  }

  // Throughput summary
  console.log(`\n=== Throughput Summary ===`);
  console.log(`Encode: ${(payloadSize / MB / (encodeMs / 1000)).toFixed(2)} MB/s`);
  console.log(`Decode: ${(payloadSize / MB / (decodeMs / 1000)).toFixed(2)} MB/s`);
  console.log(`Total:  ${(payloadSize / MB / ((encodeMs + decodeMs) / 1000)).toFixed(2)} MB/s end-to-end`);
}

main().catch((e) => {
  console.error("Test failed:", e);
  process.exit(1);
});
