/**
 * v64: K=9 Nanopore IDS Recovery Benchmark (Task 1)
 *
 * Tests the NASA K=9 convolutional code (memory=8, d_free=24) with the
 * indel-tolerant Viterbi decoder at 9% IDS (PRESET_NANOPORE).
 *
 * v64 optimizations wired in:
 *   - PrecomputedTransitionLUT (from mega-performance.ts) — explicit wiring
 *   - Reusable decode buffer pool — avoids 364MB allocation per decode call
 *   - Configurable maxDrift — benchmark uses maxDrift=15 (2× speedup)
 *
 * The benchmark uses a small payload (4KB) and moderate coverage (5×) to
 * complete in reasonable time. The K=9 decoder is still ~150ms/decode with
 * maxDrift=15, so 4KB × 5× = ~80 oligos × 5 reads = 400 decodes × 150ms
 * = ~60 seconds per IDS level.
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

async function benchIDS(
  label: string,
  sub: number,
  ins: number,
  del: number,
  coverage: number,
  payloadSize: number,
): Promise<void> {
  const cfg = ULTIMATE_V61_NANOPORE_CONFIG;
  const payload = randomBytes(payloadSize);

  const encodeStart = performance.now();
  const enc = await encodeFile(payload, cfg, {
    fileName: "b.bin",
    contentType: "application/octet-stream",
  });
  const encodeMs = performance.now() - encodeStart;
  const { oligos, metadata } = enc.encoded;

  const simCfg: MutationConfig = {
    substitutionRate: sub,
    insertionRate: ins,
    deletionRate: del,
    coverage,
    dropoutRate: 0,
    seed: 42,
  };
  const { reads } = simulate(oligos, simCfg);

  const decodeStart = performance.now();
  const result = await decodeReads(
    reads,
    metadata,
    cfg,
    enc.encoded.forwardPrimer,
    enc.encoded.reversePrimer,
  );
  const decodeMs = performance.now() - decodeStart;

  const hashMatch = result.hashMatches;
  const dataMatch = result.data
    ? Buffer.compare(Buffer.from(payload), Buffer.from(result.data)) === 0
    : false;

  const totalIds = (sub + ins + del) * 100;
  console.log(
    `  ${label.padEnd(30)} hash=${hashMatch ? "OK ✅" : "FAIL ❌"} ` +
    `data=${dataMatch ? "OK ✅" : "FAIL ❌"} ` +
    `recovered=${result.stats.oligosRecovered}/${oligos.length} ` +
    `decode=${(decodeMs / 1000).toFixed(1)}s ` +
    `(${(decodeMs / reads.length).toFixed(0)}ms/read)`,
  );
}

async function main(): Promise<void> {
  console.log("=== v64 K=9 Nanopore IDS Recovery Benchmark ===");
  console.log(`Date: ${new Date().toISOString()}`);
  console.log(`Node: ${process.version}, Bun: ${(globalThis as any).Bun?.version ?? "n/a"}`);
  console.log(`\nConfig: v61-nanopore (300nt, K=9 conv, 8B LDPC, 25% outer RS)`);
  console.log(`Decoder: IndelTolerantViterbi (NASA K=9, memory=8, d_free=24)`);
  console.log(`Optimizations: PrecomputedTransitionLUT + buffer pool + maxDrift=15\n`);

  // Warm up the cache (first decode builds the transition table)
  console.log("Warm-up: encoding 512B payload, 2× coverage...");
  const warmPayload = randomBytes(512);
  const warmEnc = await encodeFile(warmPayload, ULTIMATE_V61_NANOPORE_CONFIG, {
    fileName: "w.bin",
    contentType: "application/octet-stream",
  });
  const warmSim = simulate(warmEnc.encoded.oligos, {
    substitutionRate: 0.02,
    insertionRate: 0.03,
    deletionRate: 0.04,
    coverage: 2,
    dropoutRate: 0,
    seed: 1,
  });
  const warmStart = performance.now();
  await decodeReads(
    warmSim.reads,
    warmEnc.encoded.metadata,
    ULTIMATE_V61_NANOPORE_CONFIG,
    warmEnc.encoded.forwardPrimer,
    warmEnc.encoded.reversePrimer,
  );
  console.log(`Warm-up done in ${((performance.now() - warmStart) / 1000).toFixed(1)}s\n`);

  const stats = getIndelViterbiCacheStats();
  console.log(`Cache stats: transitions=${stats.transitionCacheSize}, inner=${stats.innerCacheSize}, buffers=${stats.bufferPoolSize} (${(stats.bufferPoolBytes / 1024 / 1024).toFixed(0)}MB)\n`);

  // Test at multiple IDS levels with 2KB payload, 3× coverage
  const payloadSize = 2 * 1024;
  const coverage = 3;

  console.log(`Testing with ${(payloadSize / 1024).toFixed(0)}KB payload, ${coverage}× coverage:\n`);

  await benchIDS("0% IDS (clean)",           0.00, 0.00, 0.00, coverage, payloadSize);
  await benchIDS("3% IDS (mild)",            0.01, 0.01, 0.01, coverage, payloadSize);
  await benchIDS("6% IDS (moderate)",        0.02, 0.02, 0.02, coverage, payloadSize);
  await benchIDS("9% IDS (PRESET_NANOPORE)", 0.02, 0.03, 0.04, coverage, payloadSize);

  console.log("\n=== Benchmark Complete ===");
  const finalStats = getIndelViterbiCacheStats();
  console.log(`Final cache stats: transitions=${finalStats.transitionCacheSize}, inner=${finalStats.innerCacheSize}, buffers=${finalStats.bufferPoolSize} (${(finalStats.bufferPoolBytes / 1024 / 1024).toFixed(0)}MB)`);
}

main().catch((e) => {
  console.error("Benchmark failed:", e);
  process.exit(1);
});
