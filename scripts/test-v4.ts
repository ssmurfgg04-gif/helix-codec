// Test v4.0 modules: CRC-32, GF(2^16), Fountain codes, Goldman mapping, Profile HMM.
import { crc32, crc32Bytes, verifyCrc32 } from "../src/lib/dna/crc32";
import { fountainEncode, fountainDecode, robustSolitonCDF, DEFAULT_FOUNTAIN_CONFIG } from "../src/lib/dna/fountain";
import { bytesToGoldmanDna, goldmanDnaToBytes, hasHomopolymer } from "../src/lib/dna/goldman";
import { forwardBackward, fusePosteriors, posteriorsToLLRs, isConfident } from "../src/lib/dna/profileHmm";

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

console.log("=== v4.0 Module Tests ===\n");

// --- CRC-32 ---
console.log("--- CRC-32 ---\n");
{
  const data = new Uint8Array([0x01, 0x02, 0x03, 0x04, 0x05]);
  const crc = crc32(data);
  assert(crc > 0, "CRC-32 produces nonzero");
  // Known test vector: CRC-32 of "123456789" = 0xCBF43926
  const testVec = new TextEncoder().encode("123456789");
  assert(crc32(testVec) === 0xcbf43926, `CRC-32 of "123456789" = 0xcbf43926 (got 0x${crc32(testVec).toString(16)})`);

  const withCrc = new Uint8Array(data.length + 4);
  withCrc.set(data, 0);
  withCrc.set(crc32Bytes(data), data.length);
  assert(verifyCrc32(withCrc), "CRC-32 verifies");

  withCrc[0] ^= 0xff;
  assert(!verifyCrc32(withCrc), "CRC-32 detects corruption");
}

// --- Fountain codes ---
console.log("\n--- LT/Fountain Codes ---\n");
{
  const data = new TextEncoder().encode("Fountain code test data for DNA storage!");
  const encoding = fountainEncode(data, { ...DEFAULT_FOUNTAIN_CONFIG, seed: 42 }, Math.ceil(data.length / 32) * 1.5);
  console.log(`  Encoded ${data.length} bytes → ${encoding.numChunks} chunks, ${encoding.droplets.length} droplets`);

  const recovered = fountainDecode(encoding);
  assert(recovered !== null, "Fountain decode succeeds");
  if (recovered) {
    assert(equalArrays(recovered, data), "Fountain round-trip matches");
  }

  // Test with fewer droplets (near threshold)
  const minDroplets = Math.ceil(encoding.numChunks * 1.05);
  const encoding2 = fountainEncode(data, { ...DEFAULT_FOUNTAIN_CONFIG, seed: 42 }, minDroplets);
  const recovered2 = fountainDecode(encoding2);
  // Should usually work with 5% overhead
  if (recovered2) {
    assert(equalArrays(recovered2, data), "Fountain with 5% overhead succeeds");
  } else {
    console.log("  (5% overhead decode failed — expected occasionally with fountain codes)");
  }

  // RSD
  const { cdf, S } = robustSolitonCDF(100, 0.1, 0.5);
  assert(cdf[100] === 1.0, "RSD CDF ends at 1.0");
  assert(S > 0, `RSD S parameter > 0 (got ${S.toFixed(2)})`);
}

// --- Goldman rotational mapping ---
console.log("\n--- Goldman Rotational Mapping ---\n");
{
  const data = new Uint8Array([0x48, 0x65, 0x6c, 0x6c, 0x6f]); // "Hello"
  const dna = bytesToGoldmanDna(data);
  console.log(`  "Hello" → ${dna}`);

  assert(!hasHomopolymer(dna), "Goldman DNA has no homopolymers");
  assert(dna.length === data.length * 6, `DNA length = ${data.length * 6} (6 trits/byte)`);

  const recovered = goldmanDnaToBytes(dna);
  assert(equalArrays(recovered, data), "Goldman round-trip matches");

  // Larger test
  const bigData = new Uint8Array(100);
  for (let i = 0; i < 100; i++) bigData[i] = (i * 31 + 17) & 0xff;
  const bigDna = bytesToGoldmanDna(bigData);
  assert(!hasHomopolymer(bigDna), "Larger Goldman DNA has no homopolymers");
  assert(equalArrays(goldmanDnaToBytes(bigDna), bigData), "Larger Goldman round-trip");

  // All same byte (stress test for homopolymer avoidance)
  const sameData = new Uint8Array(20).fill(0xaa);
  const sameDna = bytesToGoldmanDna(sameData);
  assert(!hasHomopolymer(sameDna), "All-0xAA still no homopolymers");
}

// --- Profile HMM ---
console.log("\n--- Profile HMM ---\n");
{
  const ref = "ACGTACGTACGT";
  const read = "ACGTACGTACGT"; // perfect match
  const quality = new Uint8Array(read.length).fill(30);

  const result = forwardBackward(read, ref, quality);
  assert(result.logLikelihood > -100, "Perfect match has reasonable log-likelihood");
  assert(result.path.length === ref.length, "Path length = ref length");
  assert(result.posteriors.length === ref.length * 4, "Posteriors length = ref * 4");

  // Fused posteriors
  const fused = fusePosteriors([result.posteriors, result.posteriors]);
  assert(fused.length === result.posteriors.length, "Fused posteriors same length");

  // LLRs
  const llrs = posteriorsToLLRs(fused);
  assert(llrs.length === ref.length * 2, "LLRs length = ref * 2");

  // Confidence check
  assert(isConfident(fused, 0.99), "Two perfect reads → confident at 0.99");
  assert(isConfident(fused, 0.999), "Two perfect reads → confident at 0.999");

  // Read with substitution
  const subRead = "ACGTACGTACGA"; // last base substituted
  const subQuality = new Uint8Array(subRead.length).fill(30);
  subQuality[subQuality.length - 1] = 5; // low Q at substitution
  const subResult = forwardBackward(subRead, ref, subQuality);
  assert(subResult.logLikelihood < result.logLikelihood, "Substituted read has lower log-likelihood");

  // Fuse a perfect read with a substituted read — should still be confident
  const fused2 = fusePosteriors([result.posteriors, subResult.posteriors]);
  assert(isConfident(fused2, 0.9), "Perfect + substituted → confident at 0.9");
}

console.log("\n=== All v4.0 module tests passed ===");
