/**
 * v64: Quick K=9 decode timing test
 *
 * Measures actual decode speed of the IndelTolerantConvolutionalInnerCode
 * with the v64 optimizations (buffer pool + maxDrift=15 + PrecomputedTransitionLUT).
 *
 * Encodes a SINGLE oligo's worth of data, simulates 9% IDS, and decodes it.
 * Repeats 10 times to get average timing.
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, MutationConfig } from "../src/lib/dna/simulate";
import { ULTIMATE_V61_NANOPORE_CONFIG } from "../src/lib/dna/presets";
import { getIndelViterbiCacheStats } from "../src/lib/dna/convolutional-indel";

function randomBytes(n: number): Uint8Array {
  const buf = new Uint8Array(n);
  for (let i = 0; i < n; i++) buf[i] = Math.floor(Math.random() * 256);
  return buf;
}

async function main(): Promise<void> {
  console.log("=== v64 K=9 Decode Timing Test ===");
  console.log(`Date: ${new Date().toISOString()}`);
  console.log(`Config: maxDrift=15, buffer pool, PrecomputedTransitionLUT\n`);

  const cfg = ULTIMATE_V61_NANOPORE_CONFIG;

  // Use 1KB payload, 5× coverage
  const payloadSize = 1024;
  const payload = randomBytes(payloadSize);

  console.log(`Encoding ${payloadSize}B payload...`);
  const enc = await encodeFile(payload, cfg, {
    fileName: "b.bin",
    contentType: "application/octet-stream",
  });
  const { oligos, metadata } = enc.encoded;
  console.log(`Encoded: ${oligos.length} oligos\n`);

  // Simulate 9% IDS at 5× coverage
  const simCfg: MutationConfig = {
    substitutionRate: 0.02,
    insertionRate: 0.03,
    deletionRate: 0.04,
    coverage: 5,
    dropoutRate: 0,
    seed: 42,
  };
  const { reads } = simulate(oligos, simCfg);
  console.log(`Simulated ${reads.length} reads at 9% IDS, 5× coverage\n`);

  // Time the decode
  console.log(`Decoding...`);
  const decodeStart = performance.now();
  const result = await decodeReads(
    reads,
    metadata,
    cfg,
    enc.encoded.forwardPrimer,
    enc.encoded.reversePrimer,
  );
  const decodeMs = performance.now() - decodeStart;

  console.log(`\nDecode took ${(decodeMs / 1000).toFixed(1)}s`);
  console.log(`Per-read: ${(decodeMs / reads.length).toFixed(0)}ms/read`);
  console.log(`Hash: ${result.hashMatches ? "OK ✅" : "FAIL ❌"}`);
  console.log(`Recovered: ${result.stats.oligosRecovered}/${oligos.length}`);

  const stats = getIndelViterbiCacheStats();
  console.log(`\nCache stats:`);
  console.log(`  transition cache: ${stats.transitionCacheSize} entries`);
  console.log(`  inner cache: ${stats.innerCacheSize} entries`);
  console.log(`  buffer pool: ${stats.bufferPoolSize} buffers (${(stats.bufferPoolBytes / 1024 / 1024).toFixed(0)}MB)`);
}

main().catch((e) => {
  console.error("Test failed:", e);
  process.exit(1);
});
