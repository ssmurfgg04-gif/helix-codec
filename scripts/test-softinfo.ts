// Test soft-information passthrough vs. standard decoding.
// Compares recovery success at Real 2024 error rates with and without soft-info.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA, PRESET_REAL_2024 } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";

async function main() {
  console.log("=== Soft-Information Passthrough Test ===\n");

  // Use a larger payload to get meaningful oligo counts
  const text = `The Synthetic DNA Data Storage Codec is a system for encoding digital information into synthetic DNA sequences. DNA has an information density of approximately 455 exabytes per gram, making it an attractive medium for long-term archival storage. Unlike magnetic tape or optical disks, DNA can last for thousands of years if stored in cool, dry conditions.`;
  const payload = new TextEncoder().encode(text.repeat(3));
  console.log(`Payload: ${payload.length} bytes`);

  const encodeResult = await encodeFile(payload, DEFAULT_CONFIG, {
    fileName: "softinfo-test.txt",
    contentType: "text/plain",
  });
  console.log(`Encoded: ${encodeResult.stats.oligoCount} oligos\n`);

  // Test 1: Illumina (low error) — both should pass
  console.log("--- Illumina (0.25% total, 20x) ---");
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_ILLUMINA, seed: 42 });
    console.log(`  Reads: ${sim.totalReads}, Errors: ${sim.totalErrors}`);

    const decStandard = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, false,
    );
    console.log(`  Standard:  ${decStandard.hashMatches ? "PASS" : "FAIL"} (${decStandard.stats.oligosRecovered} recovered, ${decStandard.stats.oligosFailedInnerRS} inner RS fails)`);

    const decSoftInfo = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, true,
    );
    console.log(`  Soft-info: ${decSoftInfo.hashMatches ? "PASS" : "FAIL"} (${decSoftInfo.stats.oligosRecovered} recovered, ${decSoftInfo.stats.oligosFailedInnerRS} inner RS fails)`);
  }

  // Test 2: Real 2024 (high error) — soft-info should help significantly
  console.log("\n--- Real 2024 (12.3% total, 25x) ---");
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_REAL_2024, seed: 42 });
    console.log(`  Reads: ${sim.totalReads}, Errors: ${sim.totalErrors}`);

    const decStandard = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, false,
    );
    console.log(`  Standard:  ${decStandard.hashMatches ? "PASS" : "FAIL"} (${decStandard.stats.oligosRecovered} recovered, ${decStandard.stats.oligosFailedInnerRS} inner RS fails)`);

    const decSoftInfo = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, true,
    );
    console.log(`  Soft-info: ${decSoftInfo.hashMatches ? "PASS" : "FAIL"} (${decSoftInfo.stats.oligosRecovered} recovered, ${decSoftInfo.stats.oligosFailedInnerRS} inner RS fails)`);
  }

  // Test 3: Real 2024 with higher coverage (40x) — soft-info should help more
  console.log("\n--- Real 2024 (12.3% total, 40x coverage) ---");
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_REAL_2024, coverage: 40, seed: 42 });
    console.log(`  Reads: ${sim.totalReads}, Errors: ${sim.totalErrors}`);

    const decStandard = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, false,
    );
    console.log(`  Standard:  ${decStandard.hashMatches ? "PASS" : "FAIL"} (${decStandard.stats.oligosRecovered} recovered, ${decStandard.stats.oligosFailedInnerRS} inner RS fails)`);

    const decSoftInfo = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, true,
    );
    console.log(`  Soft-info: ${decSoftInfo.hashMatches ? "PASS" : "FAIL"} (${decSoftInfo.stats.oligosRecovered} recovered, ${decSoftInfo.stats.oligosFailedInnerRS} inner RS fails)`);
  }

  // Test 4: Moderate sub-only (1% sub, 0 indel, 10x) — soft-info sweet spot
  console.log("\n--- Sub-only (1% sub, 0 indel, 10x) — soft-info erasure sweet spot ---");
  {
    const sim = simulate(encodeResult.encoded.oligos, {
      substitutionRate: 0.01,
      insertionRate: 0,
      deletionRate: 0,
      coverage: 10,
      dropoutRate: 0,
      seed: 42,
    });
    console.log(`  Reads: ${sim.totalReads}, Errors: ${sim.totalErrors}`);

    const decStandard = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, false,
    );
    console.log(`  Standard:  ${decStandard.hashMatches ? "PASS" : "FAIL"} (${decStandard.stats.oligosRecovered} recovered, ${decStandard.stats.oligosFailedInnerRS} inner RS fails)`);

    const decSoftInfo = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, true,
    );
    console.log(`  Soft-info: ${decSoftInfo.hashMatches ? "PASS" : "FAIL"} (${decSoftInfo.stats.oligosRecovered} recovered, ${decSoftInfo.stats.oligosFailedInnerRS} inner RS fails)`);
  }

  // Test 6: Heavy sub (5% sub, 0 indel, 8x) — few reads, many errors per read
  // This is where soft-info should help: consensus has residual errors, erasure hints
  // tell RS where they are (2x capacity)
  console.log("\n--- Heavy sub (5% sub, 0 indel, 8x) — soft-info should help ---");
  {
    const sim = simulate(encodeResult.encoded.oligos, {
      substitutionRate: 0.05,
      insertionRate: 0,
      deletionRate: 0,
      coverage: 8,
      dropoutRate: 0,
      seed: 42,
    });
    console.log(`  Reads: ${sim.totalReads}, Errors: ${sim.totalErrors}`);

    const decStandard = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, false,
    );
    console.log(`  Standard:  ${decStandard.hashMatches ? "PASS" : "FAIL"} (${decStandard.stats.oligosRecovered} recovered, ${decStandard.stats.oligosFailedInnerRS} inner RS fails, ${decStandard.stats.oligosErased} erased)`);

    const decSoftInfo = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, true,
    );
    console.log(`  Soft-info: ${decSoftInfo.hashMatches ? "PASS" : "FAIL"} (${decSoftInfo.stats.oligosRecovered} recovered, ${decSoftInfo.stats.oligosFailedInnerRS} inner RS fails, ${decSoftInfo.stats.oligosErased} erased)`);
  }

  // Test 7: Very heavy sub (8% sub, 0 indel, 5x) — extreme case
  console.log("\n--- Very heavy sub (8% sub, 0 indel, 5x) ---");
  {
    const sim = simulate(encodeResult.encoded.oligos, {
      substitutionRate: 0.08,
      insertionRate: 0,
      deletionRate: 0,
      coverage: 5,
      dropoutRate: 0,
      seed: 42,
    });
    console.log(`  Reads: ${sim.totalReads}, Errors: ${sim.totalErrors}`);

    const decStandard = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, false,
    );
    console.log(`  Standard:  ${decStandard.hashMatches ? "PASS" : "FAIL"} (${decStandard.stats.oligosRecovered} recovered, ${decStandard.stats.oligosFailedInnerRS} inner RS fails)`);

    const decSoftInfo = await decodeReads(
      sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer, true,
    );
    console.log(`  Soft-info: ${decSoftInfo.hashMatches ? "PASS" : "FAIL"} (${decSoftInfo.stats.oligosRecovered} recovered, ${decSoftInfo.stats.oligosFailedInnerRS} inner RS fails)`);
  }
}

main().catch(console.error);
