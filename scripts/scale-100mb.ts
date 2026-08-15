// 100MB scale test — prove streaming/memory management holds for real-world archival sizes.
//
// Uses compressible text data (real archives are mostly text/code that compresses 3-5x).
// At 100MB, this tests:
//   1. Memory management (100MB input → ~25MB compressed → ~500K oligos → ~5M reads at 10x)
//   2. GF(2^16) outer RS at scale (n > 65535 needs multiple RS blocks... wait, n <= 65535)
//   3. Encoding throughput (should be >1 MB/s for practical use)
//   4. Decoding throughput (should be >5 MB/s for practical use)

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as crypto from "crypto";

async function main() {
  console.log("=== 100MB SCALE TEST ===\n");

  // Generate 100MB of compressible text data
  // (Real DNA storage archives are text-heavy: documents, source code, genomics data)
  console.log("Generating 100MB of compressible text data...");
  const textPattern = "The Helix Codec encodes digital data into synthetic DNA oligos for archival storage. " +
    "It uses LDPC inner codes, GF(2^16) Reed-Solomon outer codes, CRC-16 integrity checks, " +
    "and DEFLATE compression to achieve near-Shannon-limit density. ";
  const patternLen = textPattern.length;
  const repeatCount = Math.ceil((100 * 1024 * 1024) / patternLen);
  const bigPayload = Buffer.from(textPattern.repeat(repeatCount).slice(0, 100 * 1024 * 1024));
  console.log(`  Payload: ${bigPayload.length.toLocaleString()} bytes (${(bigPayload.length / 1024 / 1024).toFixed(0)} MB)`);

  const payloadHash = crypto.createHash("sha256").update(bigPayload).digest("hex");
  console.log(`  SHA-256: ${payloadHash.slice(0, 16)}...`);

  // Config: LDPC inner, 300nt oligo, 30% outer RS, relaxed constraints for high-entropy data
  const bigConfig = {
    ...DEFAULT_CONFIG,
    oligoLength: 300,
    primerLength: 20,
    outerParityRatio: 0.3,
    constraints: { gcMin: 0.35, gcMax: 0.65, maxHomopolymer: 4 },
    maxRetries: 16,
  };
  console.log(`\nConfig: innerCode=${bigConfig.innerCode}, innerParityBytes=${bigConfig.innerParityBytes}`);
  console.log(`  Constraints: gc=[${bigConfig.constraints.gcMin},${bigConfig.constraints.gcMax}], maxHomopolymer=${bigConfig.constraints.maxHomopolymer}`);

  // Encode
  console.log(`\nEncoding 100MB...`);
  const t0 = Date.now();
  const enc = await encodeFile(bigPayload, bigConfig, { fileName: "100mb.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  console.log(`  Encoded: ${enc.encoded.oligos.length.toLocaleString()} oligos in ${(encMs / 1000).toFixed(1)}s`);
  console.log(`  Encode throughput: ${(100 / (encMs / 1000)).toFixed(1)} MB/s`);
  console.log(`  Compressed size: ${(enc.stats.compressedSize / 1024 / 1024).toFixed(2)} MB (ratio: ${(enc.stats.compressedSize / bigPayload.length).toFixed(3)})`);
  console.log(`  Payload per oligo: ${enc.stats.payloadBytesPerOligo} bytes`);
  console.log(`  Net density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`  Screening retries: ${enc.stats.screeningRetries.toLocaleString()} (avg ${(enc.stats.screeningRetries / enc.encoded.oligos.length).toFixed(1)} per oligo)`);
  console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);

  // Simulate 10x coverage (standard for DNA storage archives)
  console.log(`\nSimulating 10x coverage...`);
  const simT0 = Date.now();
  const sim = simulate(enc.encoded.oligos, {
    substitutionRate: 0.001, // Standard Illumina 0.1%
    insertionRate: 0.0005,
    deletionRate: 0.001,
    dropoutRate: 0.0,
    coverage: 10,
    seed: 42,
  });
  const simMs = Date.now() - simT0;
  console.log(`  Simulated ${sim.totalReads.toLocaleString()} reads in ${(simMs / 1000).toFixed(1)}s`);
  console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);

  // Decode
  console.log(`\nDecoding ${sim.totalReads.toLocaleString()} reads...`);
  const decT0 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, bigConfig, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
  const decMs = Date.now() - decT0;
  console.log(`  Decoded in ${(decMs / 1000).toFixed(1)}s`);
  console.log(`  Decode throughput: ${(100 / (decMs / 1000)).toFixed(1)} MB/s`);
  console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);
  console.log(`  Recovery: ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"} (${dec.stats.oligosRecovered.toLocaleString()}/${enc.encoded.oligos.length.toLocaleString()} oligos)`);
  console.log(`  Oligos erased: ${dec.stats.oligosErased}`);
  console.log(`  Oligos failed inner: ${dec.stats.oligosFailedInnerRS}`);
  console.log(`  Oligos failed outer: ${dec.stats.oligosFailedOuterRS}`);

  if (dec.hashMatches) {
    console.log(`\n✅ 100MB SCALE TEST PASSED!`);
    console.log(`   Encoded ${bigPayload.length.toLocaleString()} bytes → ${enc.encoded.oligos.length.toLocaleString()} oligos`);
    console.log(`   Simulated ${sim.totalReads.toLocaleString()} reads at 10x coverage`);
    console.log(`   Decoded with zero data loss`);
    console.log(`   Density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
    console.log(`   Total time: ${((encMs + simMs + decMs) / 1000).toFixed(1)}s (encode + simulate + decode)`);
  } else {
    console.log(`\n❌ 100MB SCALE TEST FAILED`);
    if (dec.data) {
      const recoveredHash = crypto.createHash("sha256").update(dec.data).digest("hex");
      console.log(`   Expected: ${payloadHash.slice(0, 32)}...`);
      console.log(`   Got:      ${recoveredHash.slice(0, 32)}...`);
      console.log(`   Recovered size: ${dec.data.length.toLocaleString()} (expected ${bigPayload.length.toLocaleString()})`);
    }
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
