/**
 * v59 K-mer Clustering Test — verify address recovery at 9% IDS.
 *
 * Tests that k-mer clustering recovers addresses that exact matching misses,
 * unlocking both nanopore 9% IDS tolerance and arithmetic mode.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { ULTIMATE_NANOPORE_V52_CONFIG, ULTIMATE_V55_DENSITY_CONFIG, computeDensity } from "../src/lib/dna/presets";
import * as fs from "fs";

const TAG = "[v59-kmer]";

async function testNanoporeIds() {
  console.log(`${TAG} === Test 1: Nanopore 9% IDS with k-mer clustering ===\n`);

  // Small payload for fast iteration
  const payload = new Uint8Array(16 * 1024); // 16KB
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;

  const cfg = ULTIMATE_NANOPORE_V52_CONFIG;
  console.log(`${TAG} Config: ${cfg.oligoLength}nt oligos, channel=${cfg.channel}, conv=${cfg.useConvolutionalInner}`);
  console.log(`${TAG} Inner parity: ${cfg.innerParityBytes}B, outer RS: ${(cfg.outerParityRatio * 100).toFixed(0)}%`);

  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  console.log(`${TAG} Encoded: ${enc.encoded.oligos.length} oligos, ${(encMs / 1000).toFixed(2)}s`);

  // Simulate nanopore reads at 9% IDS, 10× coverage
  const coverage = 10;
  const sim = simulate(enc.encoded.oligos, { ...PRESET_NANOPORE, coverage, seed: 42 });
  console.log(`${TAG} Simulated ${sim.reads.length} reads at ${coverage}× coverage, 9% IDS\n`);

  // Decode with k-mer clustering (auto-enabled for nanopore channel)
  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;

  console.log(`${TAG} Result:`);
  console.log(`${TAG}   decode: ${(decMs / 1000).toFixed(2)}s = ${((payload.length / 1024 / 1024) / (decMs / 1000)).toFixed(2)} MB/s`);
  console.log(`${TAG}   hash: ${dec.hashMatches ? "OK ✓" : "FAIL ✗"}`);
  console.log(`${TAG}   oligos recovered: ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
  console.log();
}

async function testArithmeticMode() {
  console.log(`${TAG} === Test 2: Arithmetic mode at 1.85+ b/nt ===\n`);

  // Use the v55 omni config (arithmetic mode)
  const cfg = {
    ...ULTIMATE_V55_DENSITY_CONFIG,
    mappingMode: "arithmetic" as const,
    innerParityBytes: 8, // m=64 for erasure correction
    outerParityRatio: 0.05,
  };

  const payload = new Uint8Array(16 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;

  console.log(`${TAG} Config: ${cfg.oligoLength}nt, mapping=${cfg.mappingMode}, parity=${cfg.innerParityBytes}B`);
  const density = computeDensity(cfg);
  console.log(`${TAG} Theoretical density: ${density.toFixed(3)} b/nt`);

  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  console.log(`${TAG} Encoded: ${enc.encoded.oligos.length} oligos, ${(encMs / 1000).toFixed(2)}s`);
  console.log(`${TAG} Realized density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);

  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });

  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;

  console.log(`${TAG} Result:`);
  console.log(`${TAG}   decode: ${(decMs / 1000).toFixed(2)}s`);
  console.log(`${TAG}   hash: ${dec.hashMatches ? "OK ✓" : "FAIL ✗"}`);
  console.log();
}

async function main() {
  console.log(`${TAG} v59 K-mer Clustering Validation\n`);
  console.log(`${TAG} This test verifies that k-mer based address recovery`);
  console.log(`${TAG} unblocks both nanopore 9% IDS tolerance and arithmetic mode.\n`);

  await testNanoporeIds();
  await testArithmeticMode();

  console.log(`${TAG} === Summary ===`);
  console.log(`${TAG} v59 improvements:`);
  console.log(`${TAG}   - K-mer clustering wired into decode path (clusterReadsWithKmer)`);
  console.log(`${TAG}   - Pre-computed reference addresses for all N oligo indices`);
  console.log(`${TAG}   - K-mer inverted index for O(R·k) matching`);
  console.log(`${TAG}   - Arithmetic mode + nanopore channel auto-route to JS path`);
}
main().catch(e => { console.error(e); process.exit(1); });
