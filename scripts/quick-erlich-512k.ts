// Quick full-pipeline test with 512KB subset of Erlich payload.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";
import * as crypto from "crypto";

async function main() {
  console.log("=== Quick Full-Pipeline Test (512KB Erlich subset) ===\n");

  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  const subset = payload.slice(0, 512 * 1024); // 512KB
  console.log(`Payload: ${subset.length.toLocaleString()} bytes (512KB subset of 2.1MB Erlich)`);

  const payloadHash = crypto.createHash("sha256").update(subset).digest("hex");
  console.log(`SHA-256: ${payloadHash.slice(0, 16)}...`);

  const testConfig = {
    ...DEFAULT_CONFIG,
    oligoLength: 300,
    primerLength: 20,
    outerParityRatio: 0.3,
    constraints: { gcMin: 0.35, gcMax: 0.65, maxHomopolymer: 4 },
    maxRetries: 16,
  };

  console.log(`\nEncoding...`);
  const t0 = Date.now();
  const enc = await encodeFile(subset, testConfig, { fileName: "erlich_512k.bin", contentType: "application/octet-stream" });
  console.log(`  Encoded: ${enc.encoded.oligos.length.toLocaleString()} oligos in ${Date.now() - t0}ms`);
  console.log(`  Screening retries: ${enc.stats.screeningRetries.toLocaleString()} (avg ${(enc.stats.screeningRetries / enc.encoded.oligos.length).toFixed(1)} per oligo)`);
  console.log(`  Density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);

  // Simulate 20x coverage with empirical Erlich noise
  console.log(`\nSimulating 20x coverage...`);
  const simT0 = Date.now();
  const sim = simulate(enc.encoded.oligos, {
    substitutionRate: 0.000167,
    insertionRate: 0.0001,
    deletionRate: 0.0001,
    dropoutRate: 0.0,
    coverage: 20,
    seed: 42,
  });
  console.log(`  Simulated ${sim.totalReads.toLocaleString()} reads in ${Date.now() - simT0}ms`);

  // Decode
  console.log(`\nDecoding...`);
  const decT0 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, testConfig, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
  const decMs = Date.now() - decT0;
  console.log(`  Decoded in ${decMs}ms`);
  console.log(`  Hash matches: ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"}`);
  console.log(`  Oligos recovered: ${dec.stats.oligosRecovered.toLocaleString()}/${enc.encoded.oligos.length.toLocaleString()}`);
  console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);

  if (dec.hashMatches) {
    console.log(`\n✅ 512KB Erlich subset RECOVERED at 20x coverage with LDPC!`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
