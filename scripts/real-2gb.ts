// REAL 2GB scale test — encode and decode actual 2GB of data.
// This is NOT a 64MB sample extrapolated — it's a real 2GB encode+decode.

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as crypto from "crypto";

async function main() {
  const targetSize = parseInt(process.argv[2] || "2147483648"); // 2GB
  console.log(`=== REAL ${(targetSize / 1024 / 1024 / 1024).toFixed(2)} GB SCALE TEST ===\n`);

  // Generate compressible text data (real archives are text-heavy)
  console.log(`Generating ${(targetSize / 1024 / 1024).toFixed(0)} MB of compressible text...`);
  const textPattern = "The Helix Codec encodes digital data into synthetic DNA oligos for archival storage. " +
    "It uses LDPC inner codes with belief-propagation decoding, GF(2^16) Reed-Solomon outer codes, " +
    "CRC-16 integrity checks, and DEFLATE compression to achieve near-Shannon-limit density. " +
    "The codec runs entirely in the browser with zero dependencies on external servers. ";

  // For 2GB, we can't use String.repeat (JS string length limit ~512MB).
  // Instead, build the payload in chunks using Buffer.
  console.log(`Building ${(targetSize / 1024 / 1024 / 1024).toFixed(2)} GB payload...`);
  const patternBuf = Buffer.from(textPattern, "utf-8");
  const patternLen = patternBuf.length;
  const payload = Buffer.alloc(targetSize);
  let offset = 0;
  while (offset < targetSize) {
    const chunkLen = Math.min(patternLen, targetSize - offset);
    patternBuf.copy(payload, offset, 0, chunkLen);
    offset += chunkLen;
  }
  console.log(`Payload: ${payload.length.toLocaleString()} bytes (${(payload.length / 1024 / 1024 / 1024).toFixed(2)} GB)`);

  const hash = crypto.createHash("sha256");
  // Hash in chunks to avoid "data is too long" error
  for (let i = 0; i < payload.length; i += 64 * 1024 * 1024) {
    hash.update(payload.subarray(i, Math.min(i + 64 * 1024 * 1024, payload.length)));
  }
  const hashHex = hash.digest("hex");
  console.log(`SHA-256: ${hashHex.slice(0, 16)}...`);
  console.log(`Memory after payload: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);

  const config = { ...DEFAULT_CONFIG };
  console.log(`\nConfig: oligoLength=${config.oligoLength}, mappingMode=${config.mappingMode}, parity=${config.outerParityRatio}`);

  // Encode
  console.log(`\nEncoding ${ (targetSize / 1024 / 1024 / 1024).toFixed(2)} GB...`);
  const t0 = Date.now();
  const enc = await encodeFile(payload, config, { fileName: "2gb.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  console.log(`  Encoded: ${enc.encoded.oligos.length.toLocaleString()} oligos in ${(encMs / 1000).toFixed(1)}s`);
  console.log(`  Encode throughput: ${(payload.length / 1024 / 1024 / (encMs / 1000)).toFixed(1)} MB/s`);
  console.log(`  Compressed: ${(enc.stats.compressedSize / 1024 / 1024).toFixed(2)} MB (ratio: ${(enc.stats.compressedSize / payload.length).toFixed(4)})`);
  console.log(`  Net density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`  Payload density: ${(enc.stats.payloadBytesPerOligo * 8 / (config.oligoLength - 2 * config.primerLength)).toFixed(3)} bits/nt`);
  console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);

  // Simulate at 10x coverage
  const coverage = 10;
  console.log(`\nSimulating ${coverage}x coverage...`);
  const simT0 = Date.now();
  const sim = simulate(enc.encoded.oligos, {
    substitutionRate: 0.001,
    insertionRate: 0.0005,
    deletionRate: 0.001,
    dropoutRate: 0.0,
    coverage,
    seed: 42,
  });
  const simMs = Date.now() - simT0;
  console.log(`  Simulated ${sim.totalReads.toLocaleString()} reads in ${(simMs / 1000).toFixed(1)}s`);

  // Decode
  console.log(`\nDecoding ${sim.totalReads.toLocaleString()} reads...`);
  const decT0 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
  const decMs = Date.now() - decT0;
  console.log(`  Decoded in ${(decMs / 1000).toFixed(1)}s`);
  console.log(`  Decode throughput: ${(payload.length / 1024 / 1024 / (decMs / 1000)).toFixed(1)} MB/s`);
  console.log(`  Recovery: ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"} (${dec.stats.oligosRecovered.toLocaleString()}/${enc.encoded.oligos.length.toLocaleString()})`);
  console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);

  if (dec.hashMatches) {
    console.log(`\n✅ REAL ${(targetSize / 1024 / 1024 / 1024).toFixed(2)} GB SCALE TEST PASSED!`);
    console.log(`   Encoded ${payload.length.toLocaleString()} bytes → ${enc.encoded.oligos.length.toLocaleString()} oligos`);
    console.log(`   Simulated ${sim.totalReads.toLocaleString()} reads at ${coverage}x coverage`);
    console.log(`   Decoded with 100% recovery`);
    console.log(`   Net density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
    console.log(`   Total time: ${((encMs + simMs + decMs) / 1000).toFixed(1)}s`);
  } else {
    console.log(`\n❌ FAILED`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
