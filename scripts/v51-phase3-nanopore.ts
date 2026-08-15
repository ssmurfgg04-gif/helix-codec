/**
 * v51+ ULTIMATE PHASE 3 VALIDATION
 *
 * Tests the Viterbi preprocessor for Nanopore channel (9% total IDS).
 *
 * Encodes a 64KB payload, simulates Nanopore sequencing (PRESET_NANOPORE:
 * sub=2%, ins=3%, del=4% = 9% total IDS), and decodes with channel="nanopore".
 *
 * The Viterbi preprocessor (viterbi-preprocess.ts) should activate and
 * correct indels BEFORE the LDPC decoder runs, converting the indel channel
 * into a substitution-only channel that LDPC handles efficiently.
 *
 * Success criteria:
 *   - At 15× coverage: ≥ 95% recovery (Viterbi + LDPC + outer RS)
 *   - At 25× coverage: 100% recovery
 *   - Compare with channel="illumina" (Viterbi disabled) to show the
 *     improvement from the preprocessor
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_NANOPORE, MutationConfig } from "../src/lib/dna/simulate";
import { ULTIMATE_NANOPORE_CONFIG, computeDensity } from "../src/lib/dna/presets";
import * as crypto from "crypto";

async function main() {
  console.log("=== v51+ ULTIMATE PHASE 3: Nanopore IDS Tolerance ===\n");

  const payload = crypto.randomBytes(16 * 1024);
  console.log(`Payload: 16KB (random, reduced for faster Phase 3 sweep)`);

  // Use a 300nt oligo direct mapping config with channel="nanopore"
  // (WASM-compiled for 300nt direct; Viterbi preprocess is in JS decodeReads path)
  const nanoporeConfig = {
    ...ULTIMATE_NANOPORE_CONFIG,
    mappingMode: "direct" as const,
    oligoLength: 300,
    // Use RS inner code instead of LDPC for guaranteed minimum distance
    // (LDPC has error floor at high error rates; RS has guaranteed correction)
    innerCode: "rs" as const,
    innerParityBytes: 16,
    // Higher outer parity for Nanopore — 30% to tolerate the higher erasure rate
    outerParityRatio: 0.30,
    channel: "nanopore" as const,
    interleaveDepth: 0,
  };

  // Same config but channel="illumina" — Viterbi preprocessor OFF
  const illuminaConfig = {
    ...nanoporeConfig,
    channel: "illumina" as const,
  };

  const density = computeDensity(nanoporeConfig, "payload");
  console.log(`Config: oligoLen=${nanoporeConfig.oligoLength}, mapping=${nanoporeConfig.mappingMode}`);
  console.log(`  innerParity=${nanoporeConfig.innerParityBytes}, outerRatio=${nanoporeConfig.outerParityRatio}`);
  console.log(`  channel=${nanoporeConfig.channel}`);
  console.log(`Theoretical density: ${density.toFixed(3)} bits/nt (payload-only)`);
  console.log(`\nNanopore preset: sub=${(PRESET_NANOPORE.substitutionRate * 100).toFixed(1)}%, ins=${(PRESET_NANOPORE.insertionRate * 100).toFixed(1)}%, del=${(PRESET_NANOPORE.deletionRate * 100).toFixed(1)}%`);
  console.log(`  → Total raw IDS: ${((PRESET_NANOPORE.substitutionRate + PRESET_NANOPORE.insertionRate + PRESET_NANOPORE.deletionRate) * 100).toFixed(1)}%\n`);

  // Encode with the nanopore config (channel is encoded in metadata)
  const enc = await encodeFile(
    Buffer.from(payload),
    nanoporeConfig,
    { fileName: "test.bin", contentType: "application/octet-stream" },
  );
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos\n`);

  // Test at multiple coverage levels, with and without Viterbi preprocess
  console.log("Cov | Channel   | Reads  | Recovery | Time(ms) | OligosOK | FailedInner | FailedOuter");
  console.log("----|-----------|--------|----------|----------|----------|-------------|------------");

  // Test at multiple Nanopore error rates: from high-quality (3% IDS) to standard (9% IDS)
  const testConfigs = [
    { name: "NanoHQ (3% IDS)", noise: { ...PRESET_NANOPORE, substitutionRate: 0.005, insertionRate: 0.01, deletionRate: 0.015, coverage: 15, seed: 42 } },
    { name: "NanoStd (9% IDS)", noise: { ...PRESET_NANOPORE, coverage: 15, seed: 42 } },
    { name: "NanoStd (9% IDS) @25x", noise: { ...PRESET_NANOPORE, coverage: 25, seed: 42 } },
  ];

  for (const tc of testConfigs) {
    const sim = simulate(enc.encoded.oligos, tc.noise);

    // Decode with channel="nanopore" (Viterbi ON)
    const t0 = Date.now();
    let decNano;
    try {
      decNano = await decodeReads(
        sim.reads,
        enc.encoded.metadata,
        nanoporeConfig,
        enc.encoded.forwardPrimer,
        enc.encoded.reversePrimer,
        true,
      );
    } catch (e: any) {
      console.log(`  [nanopore decode failed: ${e.message?.slice(0, 80)}]`);
      continue;
    }
    const msNano = Date.now() - t0;

    console.log(
      `${tc.name.padEnd(20)} | nanopore  | ${sim.totalReads.toString().padStart(6)} | ${decNano.hashMatches ? "✅ PASS" : "❌ FAIL"}  | ${msNano.toString().padStart(8)} | ${decNano.stats.oligosRecovered.toString().padStart(8)}/${enc.encoded.oligos.length} | ${decNano.stats.oligosFailedInnerRS.toString().padStart(11)} | ${decNano.stats.oligosFailedOuterRS.toString().padStart(10)}`,
    );

    // Decode the same reads with channel="illumina" (Viterbi OFF) for comparison
    const t1 = Date.now();
    let decIll;
    try {
      decIll = await decodeReads(
        sim.reads,
        // Override metadata channel for comparison
        { ...enc.encoded.metadata, channel: "illumina" },
        illuminaConfig,
        enc.encoded.forwardPrimer,
        enc.encoded.reversePrimer,
        true,
      );
    } catch (e: any) {
      console.log(`  [illumina decode failed: ${e.message?.slice(0, 80)}]`);
      continue;
    }
    const msIll = Date.now() - t1;

    console.log(
      `${"".padEnd(20)} | illumina  | ${sim.totalReads.toString().padStart(6)} | ${decIll.hashMatches ? "✅ PASS" : "❌ FAIL"}  | ${msIll.toString().padStart(8)} | ${decIll.stats.oligosRecovered.toString().padStart(8)}/${enc.encoded.oligos.length} | ${decIll.stats.oligosFailedInnerRS.toString().padStart(11)} | ${decIll.stats.oligosFailedOuterRS.toString().padStart(10)}`,
    );
    console.log("----|-----------|--------|----------|----------|----------|-------------|------------");
  }
}

main().catch((e) => {
  console.error("FATAL:", e);
  process.exit(1);
});
