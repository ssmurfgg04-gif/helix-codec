// End-to-end test for the DNA storage codec.
// Run with: npx tsx /home/z/my-project/scripts/test-codec.ts

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG, CodecConfig } from "../src/lib/dna/types";
import { bytesToDna, dnaToBytes, satisfiesConstraints, gcContent, maxHomopolymerRun } from "../src/lib/dna/mapping";
import { crc16, crc16Bytes, verifyCrc16 } from "../src/lib/dna/crc16";

function assert(cond: boolean, msg: string) {
  if (!cond) {
    console.error("FAIL:", msg);
    process.exit(1);
  } else {
    console.log("PASS:", msg);
  }
}

function equalArrays(a: Uint8Array, b: Uint8Array): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

async function main() {
  console.log("\n=== Unit tests ===\n");

  // Mapping round-trip
  {
    const data = new Uint8Array([0x12, 0x34, 0x56, 0x78, 0xab, 0xcd, 0xef, 0x00, 0xff]);
    const dna = bytesToDna(data);
    const recovered = dnaToBytes(dna);
    assert(equalArrays(recovered, data), "DNA mapping round-trip");
    assert(dna.length === data.length * 4, "DNA length is 4x byte length");
  }

  // CRC round-trip
  {
    const data = new Uint8Array([0x01, 0x02, 0x03, 0x04, 0x05]);
    const crc = crc16(data);
    assert(crc > 0, "CRC-16 produces nonzero");
    const withCrc = new Uint8Array(data.length + 2);
    withCrc.set(data, 0);
    const crcB = crc16Bytes(data);
    withCrc.set(crcB, data.length);
    assert(verifyCrc16(withCrc), "CRC-16 verifies");
    withCrc[0] ^= 0xff;
    assert(!verifyCrc16(withCrc), "CRC-16 detects corruption");
  }

  console.log("\n=== Encoding tests ===\n");

  // Encode a small file (use random data for realistic entropy)
  const testPayload = new Uint8Array(4096);
  let seed = 42;
  for (let i = 0; i < testPayload.length; i++) {
    seed ^= seed << 13; seed ^= seed >>> 17; seed ^= seed << 5; seed >>>= 0;
    testPayload[i] = seed & 0xff;
  }
  console.log(`Test payload: ${testPayload.length} bytes`);

  // Use Goldman mapping (guaranteed homopolymer-free, no screening needed)
  const testConfig = { ...DEFAULT_CONFIG, mappingMode: "goldman" as const, goldmanMode: "fast" as const, oligoLength: 208, primerLength: 20, innerParityBytes: 4, outerParityRatio: 0.1 };
  const encodeResult = await encodeFile(testPayload, testConfig, {
    fileName: "hello.txt",
    contentType: "text/plain",
  });

  console.log(`Encoded: ${encodeResult.stats.oligoCount} oligos`);
  console.log(`Net density: ${encodeResult.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`Overhead: ${encodeResult.stats.overheadPercent.toFixed(1)}%`);
  console.log(`Screening retries: ${encodeResult.stats.screeningRetries}`);

  // All oligos must satisfy constraints
  let allSatisfy = true;
  for (const oligo of encodeResult.encoded.oligos) {
    const inner = oligo.sequence.slice(
      testConfig.primerLength,
      oligo.sequence.length - testConfig.primerLength,
    );
    if (!satisfiesConstraints(inner, {
      gcMin: testConfig.constraints.gcMin,
      gcMax: testConfig.constraints.gcMax,
      maxHomopolymer: testConfig.constraints.maxHomopolymer,
    })) {
      allSatisfy = false;
      console.log(`  Oligo ${oligo.index} FAILS: GC=${gcContent(inner).toFixed(2)} maxHp=${maxHomopolymerRun(inner)}`);
    }
  }
  assert(allSatisfy, "All oligos satisfy constraints (GC + homopolymer)");

  console.log("\n=== Decoding tests ===\n");

  // Test 1: Clean decode (no errors, coverage=1)
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
    const decodeResult = await decodeReads(
      sim.reads,
      encodeResult.encoded.metadata,
      testConfig,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );
    assert(decodeResult.data !== null, "Clean decode produces data");
    assert(decodeResult.hashMatches, "Clean decode hash matches");
    assert(equalArrays(decodeResult.data!, testPayload), "Clean decode data matches original");
    console.log(`  Clean: ${sim.totalReads} reads, decode ${decodeResult.stats.decodeTimeMs}ms`);
  }

  // Test 2: Illumina errors, 20x coverage
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_ILLUMINA, seed: 42 });
    console.log(`  Illumina: ${sim.totalReads} reads, ${sim.totalErrors} total errors`);
    const decodeResult = await decodeReads(
      sim.reads,
      encodeResult.encoded.metadata,
      testConfig,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );
    assert(decodeResult.hashMatches, "Illumina 20x decode hash matches");
    console.log(`  Recovered: ${decodeResult.stats.oligosRecovered}/${encodeResult.encoded.metadata.oligoCount} oligos`);
    console.log(`  Erased: ${decodeResult.stats.oligosErased}, Failed inner RS: ${decodeResult.stats.oligosFailedInnerRS}`);
  }

  // Test 3: High error rate (5% total)
  {
    const sim = simulate(encodeResult.encoded.oligos, {
      substitutionRate: 0.02,
      insertionRate: 0.015,
      deletionRate: 0.015,
      coverage: 30,
      dropoutRate: 0.0,
      seed: 42,
    });
    console.log(`  High-error 5%: ${sim.totalReads} reads, ${sim.totalErrors} errors`);
    const decodeResult = await decodeReads(
      sim.reads,
      encodeResult.encoded.metadata,
      testConfig,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );
    console.log(`  Recovered: ${decodeResult.stats.oligosRecovered}, hash match: ${decodeResult.hashMatches}`);
    // At 5% errors, we may or may not recover — just report.
  }

  console.log("\n=== Spec verification ===\n");

  // Verify oligo length is exactly as configured
  for (const oligo of encodeResult.encoded.oligos) {
    assert(oligo.length === testConfig.oligoLength, `Oligo ${oligo.index} length correct`);
    break; // just check first
  }

  // Verify primer presence
  const firstOligo = encodeResult.encoded.oligos[0];
  assert(
    firstOligo.sequence.startsWith(encodeResult.encoded.forwardPrimer),
    "Forward primer present",
  );
  assert(
    firstOligo.sequence.endsWith(encodeResult.encoded.reversePrimer),
    "Reverse primer present",
  );

  // Verify address structure
  console.log(`  Oligo count: ${encodeResult.encoded.metadata.oligoCount}`);
  console.log(`  Payload bytes/oligo: ${encodeResult.encoded.metadata.payloadBytesPerOligo}`);
  console.log(`  Inner RS: (${encodeResult.encoded.metadata.innerRS.n}, ${encodeResult.encoded.metadata.innerRS.k})`);
  console.log(`  Outer RS: (${encodeResult.encoded.metadata.outerRS.n}, ${encodeResult.encoded.metadata.outerRS.k})`);
  console.log(`  Parity oligos: ${encodeResult.encoded.metadata.parityOligos}`);
  console.log(`  File hash: ${encodeResult.encoded.metadata.fileHash.slice(0, 16)}...`);

  console.log("\nAll codec tests passed.");
}

main().catch((e) => {
  console.error("FATAL:", e);
  process.exit(1);
});
