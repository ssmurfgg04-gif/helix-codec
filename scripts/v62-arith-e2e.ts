/**
 * v62: End-to-end test for arithmetic-v2 mode (address outside arithmetic stream).
 *
 * Verifies that the production encode→decode pipeline works for arithmetic mode,
 * which was previously broken (v57-v61) because the address was inside the
 * arithmetic stream.
 *
 * Run: bun scripts/v62-arith-e2e.ts
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { ULTIMATE_V61_ARITHMETIC_CONFIG } from "../src/lib/dna/presets";
import { simulate } from "../src/lib/dna/simulate";

async function main() {
  console.log("=== v62 Arithmetic-v2 End-to-End Test ===\n");

  // Test 1: Small payload (1KB)
  console.log("Test 1: 1KB payload, arithmetic-v2 mode, 10× coverage, 0% noise");
  const data1 = new Uint8Array(1024);
  for (let i = 0; i < data1.length; i++) data1[i] = (i * 37 + 13) & 0xff;

  const cfg1 = { ...ULTIMATE_V61_ARITHMETIC_CONFIG };
  try {
    const encodeResult = await encodeFile(data1, cfg1, {
      fileName: "test1k.bin",
      contentType: "application/octet-stream",
    });
    console.log(`  Encoded: ${encodeResult.encoded.oligos.length} oligos`);
    console.log(`  Density: ${encodeResult.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
    console.log(`  Encode time: ${encodeResult.stats.encodeTimeMs}ms`);

    // Simulate reads
    const simResult = simulate(encodeResult.encoded.oligos, {
      coverage: 10,
      substitutionRate: 0,
      insertionRate: 0,
      deletionRate: 0,
      dropoutRate: 0,
      seed: 42,
    });
    const reads = simResult.reads;
    console.log(`  Simulated ${reads.length} reads (10× coverage)`);

    // Decode
    process.env.HELIX_DEBUG = "1";
    const decodeResult = await decodeReads(
      reads,
      encodeResult.encoded.metadata,
      cfg1,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );
    console.log(`  Decoded: ${decodeResult.data.length} bytes`);
    console.log(`  Hash match: ${decodeResult.hashMatches ? "✅ YES" : "❌ NO"}`);
    console.log(`  Oligo count: ${encodeResult.encoded.oligos.length}`);
    console.log(`  Layout: payloadBytes=${encodeResult.encoded.metadata.payloadBytesPerOligo}, innerRS=${JSON.stringify(encodeResult.encoded.metadata.innerRS)}`);
    if (decodeResult.data.length === data1.length) {
      let match = true;
      for (let i = 0; i < data1.length; i++) {
        if (decodeResult.data[i] !== data1[i]) { match = false; break; }
      }
      console.log(`  Data match: ${match ? "✅ YES" : "❌ NO"}`);
    }
  } catch (e) {
    console.log(`  ❌ FAILED: ${e}`);
  }

  // Test 2: Medium payload (64KB) with low noise
  console.log("\nTest 2: 64KB payload, arithmetic-v2 mode, 10× coverage, 0.1% sub");
  const data2 = new Uint8Array(65536);
  for (let i = 0; i < data2.length; i++) data2[i] = (i * 73 + 17) & 0xff;

  const cfg2 = { ...ULTIMATE_V61_ARITHMETIC_CONFIG };
  try {
    const encodeResult = await encodeFile(data2, cfg2, {
      fileName: "test64k.bin",
      contentType: "application/octet-stream",
    });
    console.log(`  Encoded: ${encodeResult.encoded.oligos.length} oligos`);
    console.log(`  Density: ${encodeResult.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);

    const simResult2 = simulate(encodeResult.encoded.oligos, {
      coverage: 10,
      substitutionRate: 0.001,
      insertionRate: 0,
      deletionRate: 0,
      dropoutRate: 0,
      seed: 123,
    });
    const reads = simResult2.reads;
    console.log(`  Simulated ${reads.length} reads (10× coverage, 0.1% sub)`);

    const decodeResult = await decodeReads(
      reads,
      encodeResult.encoded.metadata,
      cfg2,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );
    console.log(`  Hash match: ${decodeResult.hashMatches ? "✅ YES" : "❌ NO"}`);
  } catch (e) {
    console.log(`  ❌ FAILED: ${e}`);
  }

  // Test 3: Compare density vs direct mode
  console.log("\nTest 3: Density comparison (arithmetic-v2 vs direct)");
  const ULTIMATE_V55_DENSITY_CONFIG = (await import("../src/lib/dna/presets")).ULTIMATE_V55_DENSITY_CONFIG;
  const cfgDirect = { ...ULTIMATE_V55_DENSITY_CONFIG };
  const cfgArith = { ...ULTIMATE_V61_ARITHMETIC_CONFIG };
  try {
    const data = new Uint8Array(1024);
    for (let i = 0; i < data.length; i++) data[i] = i & 0xff;

    const r1 = await encodeFile(data, cfgDirect, { fileName: "t", contentType: "x" });
    const r2 = await encodeFile(data, cfgArith, { fileName: "t", contentType: "x" });
    console.log(`  Direct mode:       ${r1.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
    console.log(`  Arithmetic-v2:     ${r2.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
    console.log(`  SOTA (Yi Ding):    1.815 b/nt`);
  } catch (e) {
    console.log(`  ❌ FAILED: ${e}`);
  }

  console.log("\n=== v62 E2E Test Complete ===");
}

main().catch(console.error);
