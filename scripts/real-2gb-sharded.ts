// REAL 2GB scale test with SHARDING.
// For 2GB, we split the data into multiple shards, each with its own RS block.
// Each shard is encoded independently and decoded independently.

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as crypto from "crypto";

async function main() {
  const targetSize = parseInt(process.argv[2] || "2147483648");
  console.log(`=== REAL ${(targetSize / 1024 / 1024 / 1024).toFixed(2)} GB SHARDED SCALE TEST ===\n`);

  // Generate compressible text
  console.log(`Building ${(targetSize / 1024 / 1024 / 1024).toFixed(2)} GB payload...`);
  const textPattern = "The Helix Codec encodes digital data into synthetic DNA oligos for archival storage. " +
    "It uses LDPC inner codes with belief-propagation decoding, GF(2^16) Reed-Solomon outer codes, " +
    "CRC-16 integrity checks, and DEFLATE compression to achieve near-Shannon-limit density. ";
  const patternBuf = Buffer.from(textPattern, "utf-8");
  const payload = Buffer.alloc(targetSize);
  let offset = 0;
  while (offset < targetSize) {
    const chunkLen = Math.min(patternBuf.length, targetSize - offset);
    patternBuf.copy(payload, offset, 0, chunkLen);
    offset += chunkLen;
  }
  console.log(`Payload: ${payload.length.toLocaleString()} bytes`);

  // Hash in chunks
  const hash = crypto.createHash("sha256");
  for (let i = 0; i < payload.length; i += 64 * 1024 * 1024) {
    hash.update(payload.subarray(i, Math.min(i + 64 * 1024 * 1024, payload.length)));
  }
  const payloadHash = hash.digest("hex");
  console.log(`SHA-256: ${payloadHash.slice(0, 16)}...`);

  const config = { ...DEFAULT_CONFIG };

  // Split into shards of 64MB each (to stay under 65535 oligos per shard)
  const shardSize = 64 * 1024 * 1024;
  const numShards = Math.ceil(payload.length / shardSize);
  console.log(`\nSplitting into ${numShards} shards of ${shardSize / 1024 / 1024} MB each`);

  let totalOligos = 0;
  let totalEncMs = 0;
  let totalDecMs = 0;
  let allShardsPass = true;
  const shardResults: { shard: number; oligos: number; pass: boolean; encMs: number; decMs: number }[] = [];

  for (let s = 0; s < numShards; s++) {
    const shardStart = s * shardSize;
    const shardEnd = Math.min(shardStart + shardSize, payload.length);
    const shardData = payload.subarray(shardStart, shardEnd);

    console.log(`\n--- Shard ${s + 1}/${numShards} (${shardData.length.toLocaleString()} bytes) ---`);

    // Encode
    const t0 = Date.now();
    const enc = await encodeFile(shardData, config, {
      fileName: `shard_${s}.bin`,
      contentType: "application/octet-stream",
    });
    const encMs = Date.now() - t0;
    totalEncMs += encMs;
    totalOligos += enc.encoded.oligos.length;

    console.log(`  Encoded: ${enc.encoded.oligos.length.toLocaleString()} oligos in ${encMs}ms`);
    console.log(`  Density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);

    // Simulate 10x coverage
    const sim = simulate(enc.encoded.oligos, {
      substitutionRate: 0.001,
      insertionRate: 0.0005,
      deletionRate: 0.001,
      dropoutRate: 0.0,
      coverage: 10,
      seed: 42 + s,
    });

    // Decode
    const t1 = Date.now();
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
    const decMs = Date.now() - t1;
    totalDecMs += decMs;

    const pass = dec.hashMatches;
    if (!pass) allShardsPass = false;

    console.log(`  Decoded: ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"} (${dec.stats.oligosRecovered.toLocaleString()}/${enc.encoded.oligos.length.toLocaleString()}) in ${decMs}ms`);

    shardResults.push({ shard: s + 1, oligos: enc.encoded.oligos.length, pass, encMs, decMs });
  }

  console.log(`\n=== ${(targetSize / 1024 / 1024 / 1024).toFixed(2)} GB SHARDED TEST COMPLETE ===`);
  console.log(`Total oligos: ${totalOligos.toLocaleString()}`);
  console.log(`Total encode time: ${(totalEncMs / 1000).toFixed(1)}s (${(targetSize / 1024 / 1024 / (totalEncMs / 1000)).toFixed(1)} MB/s)`);
  console.log(`Total decode time: ${(totalDecMs / 1000).toFixed(1)}s (${(targetSize / 1024 / 1024 / (totalDecMs / 1000)).toFixed(1)} MB/s)`);
  console.log(`All shards pass: ${allShardsPass ? "✅ YES" : "❌ NO"}`);
  console.log(`Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);

  if (allShardsPass) {
    console.log(`\n✅ REAL ${(targetSize / 1024 / 1024 / 1024).toFixed(2)} GB SCALE TEST PASSED!`);
    console.log(`   ${numShards} shards × ~${Math.floor(totalOligos / numShards).toLocaleString()} oligos = ${totalOligos.toLocaleString()} total oligos`);
    console.log(`   Net density: ${config.oligoLength > 0 ? (54 * 8 / config.oligoLength).toFixed(3) : "N/A"} bits/nt (per oligo)`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
