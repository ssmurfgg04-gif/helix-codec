// Comprehensive edge-case and fuzzing tests.
// Run: npx tsx scripts/test-comprehensive.ts

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import { bytesToDna, dnaToBytes, gcContent, maxHomopolymerRun, satisfiesConstraints } from "../src/lib/dna/mapping";
import { crc16, crc16Bytes, verifyCrc16 } from "../src/lib/dna/crc16";
import { detectHairpins, reverseComplement, checkStructureConstraints } from "../src/lib/dna/structure";
import { kmerToBits, extractKmers, clusterByAddress } from "../src/lib/dna/kmer";
import { generateKeyPair, signArchive, verifyArchive } from "../src/lib/dna/signing";
import { holographicEncode, holographicDecode } from "../src/lib/dna/holographic";
import { insertCRCMarkers, stripCRCMarkers, crc8 } from "../src/lib/dna/crcmarker";

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
  console.log("=== Comprehensive Test Suite ===\n");

  // --- Edge cases ---
  console.log("--- Edge Cases ---\n");

  // Empty payload
  {
    const data = new Uint8Array(0);
    try {
      const enc = await encodeFile(data, DEFAULT_CONFIG, { fileName: "empty.bin", contentType: "application/octet-stream" });
      assert(enc.stats.oligoCount > 0, "Empty payload encodes (at least 1 oligo)");
    } catch (e) {
      assert(false, `Empty payload should not crash: ${(e as Error).message}`);
    }
  }

  // 1 byte
  {
    const data = new Uint8Array([42]);
    const enc = await encodeFile(data, DEFAULT_CONFIG, { fileName: "one.bin", contentType: "application/octet-stream" });
    assert(enc.stats.oligoCount > 0, "1 byte encodes");
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, DEFAULT_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    assert(dec.hashMatches, "1 byte round-trips");
  }

  // All zeros
  {
    const data = new Uint8Array(100);
    const enc = await encodeFile(data, DEFAULT_CONFIG, { fileName: "zeros.bin", contentType: "application/octet-stream" });
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, DEFAULT_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    assert(dec.hashMatches, "All-zeros round-trips");
  }

  // All 0xFF
  {
    const data = new Uint8Array(100).fill(0xff);
    const enc = await encodeFile(data, DEFAULT_CONFIG, { fileName: "ff.bin", contentType: "application/octet-stream" });
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, DEFAULT_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    assert(dec.hashMatches, "All-0xFF round-trips");
  }

  // Random binary
  {
    const data = new Uint8Array(500);
    for (let i = 0; i < 500; i++) data[i] = Math.floor(Math.random() * 256);
    const enc = await encodeFile(data, DEFAULT_CONFIG, { fileName: "random.bin", contentType: "application/octet-stream" });
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, DEFAULT_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    assert(dec.hashMatches, "Random binary round-trips");
  }

  // --- DNA mapping ---
  console.log("\n--- DNA Mapping ---\n");

  // All 256 byte values
  {
    const data = new Uint8Array(256);
    for (let i = 0; i < 256; i++) data[i] = i;
    const dna = bytesToDna(data);
    const recovered = dnaToBytes(dna);
    assert(equalArrays(recovered, data), "All 256 byte values round-trip");
    assert(dna.length === 1024, "256 bytes → 1024 nt");
  }

  // Only ACGT characters
  {
    const data = new Uint8Array([0x00, 0x55, 0xaa, 0xff]);
    const dna = bytesToDna(data);
    assert(/^[ACGT]+$/.test(dna), "DNA string only contains ACGT");
  }

  // --- CRC ---
  console.log("\n--- CRC-16 ---\n");

  // CRC detects single-bit errors
  {
    const data = new Uint8Array([0x01, 0x02, 0x03, 0x04, 0x05]);
    const withCrc = new Uint8Array(data.length + 2);
    withCrc.set(data, 0);
    withCrc.set(crc16Bytes(data), data.length);
    assert(verifyCrc16(withCrc), "CRC verifies correct data");

    // Flip one bit
    const corrupted = withCrc.slice();
    corrupted[0] ^= 0x01;
    assert(!verifyCrc16(corrupted), "CRC detects single-bit error");

    // Flip a bit in the CRC itself
    const corrupted2 = withCrc.slice();
    corrupted2[data.length] ^= 0x01;
    assert(!verifyCrc16(corrupted2), "CRC detects corrupted CRC bytes");
  }

  // --- Structure detection ---
  console.log("\n--- Structure Detection ---\n");

  // Hairpin detection
  {
    // ACGT is its own reverse complement (palindrome)
    const seq = "ACGTAAAACGT"; // arm1=ACGT, loop=AAA, arm2=ACGT (revcomp of ACGT)
    const hairpins = detectHairpins(seq);
    assert(hairpins.length > 0, "Hairpin detected in palindromic sequence");
  }

  // Reverse complement
  {
    assert(reverseComplement("ACGT") === "ACGT", "ACGT is its own reverse complement");
    assert(reverseComplement("AAAA") === "TTTT", "AAAA revcomp = TTTT");
    assert(reverseComplement("ACGTACGT") === "ACGTACGT", "ACGTACGT is a palindrome");
  }

  // Structure constraints
  {
    const badSeq = "AAAAAAAAGGGGGGGG"; // homopolymer + GC bias
    const issues = checkStructureConstraints(badSeq);
    assert(issues.length > 0, "Bad sequence has structure issues");
    assert(issues.some((i) => i.type === "homopolymer"), "Homopolymer detected");
  }

  // --- K-mer indexing ---
  console.log("\n--- K-mer Indexing ---\n");

  {
    const kmer = "ACGTACGT";
    const bits = kmerToBits(kmer);
    assert(bits >= 0, "K-mer encodes to non-negative bits");

    const kmers = extractKmers("ACGTACGTACGT", 4);
    assert(kmers.length === 9, `Extract 4-mers from 12nt: got ${kmers.length} (expected 9)`);
  }

  // Clustering by address
  {
    const reads = [
      { sequence: "ACGTACGTACGTAAAA" }, // address = ACGTACGT
      { sequence: "ACGTACGTACGTAAAC" }, // same address
      { sequence: "TTTTTTTTACGTAAAG" }, // different address
    ];
    const clusters = clusterByAddress(reads, 8);
    assert(clusters.size === 2, `Clustered into 2 groups: got ${clusters.size}`);
  }

  // --- Signing ---
  console.log("\n--- Ed25519 Signing ---\n");

  {
    const keypair = generateKeyPair();
    assert(keypair.publicKey.length === 32, "Public key is 32 bytes");
    assert(keypair.privateKey.length === 32, "Private key is 32 bytes");

    const manifest = { format: "bioarchive/v1", archiveId: "test" };
    const signed = signArchive(manifest, keypair.privateKey);
    assert(signed.signature.length === 128, "Signature is 64 bytes (128 hex chars)");
    assert(verifyArchive(signed), "Valid signature verifies");

    // Tamper
    const tampered = { ...signed, archiveId: "hacked" };
    assert(!verifyArchive(tampered), "Tampered manifest fails verification");
  }

  // --- Holographic ---
  console.log("\n--- Holographic Sharding ---\n");

  // Round-trip
  {
    const data = new Uint8Array(100);
    for (let i = 0; i < 100; i++) data[i] = i;
    const enc = holographicEncode(data, { dataShards: 10, totalShards: 15, blockSize: 10 });
    const dec = holographicDecode(enc.shards, enc);
    assert(equalArrays(dec, data), "Holographic round-trip");
  }

  // Boundary: exactly K shards
  {
    const data = new Uint8Array(50);
    for (let i = 0; i < 50; i++) data[i] = i;
    const enc = holographicEncode(data, { dataShards: 10, totalShards: 15, blockSize: 10 });
    const dec = holographicDecode(enc.shards.slice(0, 10), enc); // exactly K
    assert(equalArrays(dec, data), "Exactly K shards recovers");
  }

  // --- CRC markers ---
  console.log("\n--- CRC Markers ---\n");

  {
    const payload = new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8]);
    const marked = insertCRCMarkers(payload, { segmentSize: 4, markerSize: 1, checkpoint: 3 });
    assert(marked.length === payload.length + 2, `Marked length: ${marked.length} (expected ${payload.length + 2})`);

    const stripped = stripCRCMarkers(marked, { segmentSize: 4, markerSize: 1, checkpoint: 3 });
    assert(equalArrays(stripped.payload, payload), "CRC marker round-trip");
    assert(stripped.failedSegments.length === 0, "No failed segments (clean data)");
  }

  // CRC-8
  {
    const data = new Uint8Array([0x01, 0x02, 0x03]);
    const crc = crc8(data);
    assert(crc >= 0 && crc <= 255, "CRC-8 produces valid byte");
  }

  // --- Determinism ---
  console.log("\n--- Determinism ---\n");

  {
    const data = new TextEncoder().encode("Determinism test");
    const enc1 = await encodeFile(data, DEFAULT_CONFIG, { fileName: "det.txt", contentType: "text/plain" });
    const enc2 = await encodeFile(data, DEFAULT_CONFIG, { fileName: "det.txt", contentType: "text/plain" });
    // Primers are deterministic
    assert(enc1.encoded.forwardPrimer === enc2.encoded.forwardPrimer, "Primers are deterministic");
    // Oligo sequences depend on whitening (deterministic) and screening (seed-based, deterministic)
    assert(
      enc1.encoded.oligos[0].sequence === enc2.encoded.oligos[0].sequence,
      "Oligo sequences are deterministic",
    );
  }

  // --- Mutation simulation determinism ---
  {
    const data = new TextEncoder().encode("Sim determinism");
    const enc = await encodeFile(data, DEFAULT_CONFIG, { fileName: "sim.txt", contentType: "text/plain" });
    const sim1 = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, seed: 42 });
    const sim2 = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, seed: 42 });
    assert(sim1.totalReads === sim2.totalReads, "Simulation read count is deterministic");
    assert(sim1.totalErrors === sim2.totalErrors, "Simulation error count is deterministic");
  }

  console.log("\n=== All comprehensive tests passed ===");
}

main().catch((e) => {
  console.error("FATAL:", e);
  process.exit(1);
});
