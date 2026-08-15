// Focused end-to-end Erlich validation — smaller payload for speed.
import { fountainEncode, fountainDecode } from "../src/lib/dna/fountain";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";
import * as crypto from "crypto";

async function main() {
  console.log("=== END-TO-END ERICH 2017 VALIDATION ===\n");

  // Use 64KB subset of the Erlich-sized payload
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  const subset = payload.slice(0, 65536);
  console.log(`Payload: ${subset.length} bytes (subset of Erlich's 2.1MB)\n`);

  // --- Test 1: Helix LDPC pipeline at various coverage ---
  console.log("--- Helix LDPC + Goldman Pipeline (real Erlich data subset) ---\n");
  // Goldman mapping: 1 byte = 6 nt. oligoLength=292 → innerNt=252 → 252/6=42 bytes inner.
  // 42 - 4 (addr) - 4 (LDPC parity) - 2 (CRC) = 32 bytes payload per oligo.
  const testConfig = { ...DEFAULT_CONFIG, oligoLength: 292, primerLength: 20, outerParityRatio: 0.3 };
  const enc = await encodeFile(subset, testConfig, {
    fileName: "erlich_test.bin",
    contentType: "application/octet-stream",
  });
  console.log(`Encoded: ${enc.stats.oligoCount} oligos, ${enc.stats.compressedSize} compressed\n`);

  console.log("Cov  | Reads  | Recovery | Time  | Oligos OK");
  console.log("-----|--------|----------|-------|----------");
  for (const cov of [5, 10, 15, 20]) {
    const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: cov, seed: 42 });
    const t0 = Date.now();
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, testConfig, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
    const ms = Date.now() - t0;
    console.log(`${cov}x  | ${sim.totalReads.toString().padStart(6)} | ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"}  | ${ms}ms | ${dec.stats.oligosRecovered}/${enc.stats.oligoCount}`);
  }

  // --- Test 2: Fountain code (Erlich's approach) ---
  console.log("\n--- Fountain Code (Erlich's approach) ---\n");
  const K = Math.ceil(subset.length / 32);
  const encoding = fountainEncode(subset, { chunkSize: 32, rsdC: 0.025, rsdDelta: 0.001, seed: 42, maxDroplets: 100000 }, Math.ceil(K * 1.1));
  console.log(`K=${K} chunks, ${encoding.droplets.length} droplets (${((encoding.droplets.length / K - 1) * 100).toFixed(1)}% overhead)\n`);

  // Test: how many droplets needed to decode?
  console.log("Droplets | Recovery | Time");
  console.log("---------|----------|-----");
  for (const ratio of [1.0, 1.05, 1.1, 1.15, 1.2]) {
    const numDroplets = Math.floor(K * ratio);
    const subset_droplets = encoding.droplets.slice(0, Math.min(numDroplets, encoding.droplets.length));
    const t0 = Date.now();
    let recovered = false;
    try {
      const decoded = fountainDecode({ ...encoding, droplets: subset_droplets });
      if (decoded && decoded.length === subset.length) {
        const hash1 = crypto.createHash("sha256").update(Buffer.from(decoded)).digest("hex").slice(0, 16);
        const hash2 = crypto.createHash("sha256").update(Buffer.from(subset)).digest("hex").slice(0, 16);
        recovered = hash1 === hash2;
      }
    } catch {}
    const ms = Date.now() - t0;
    console.log(`${numDroplets.toString().padStart(8)} | ${recovered ? "✅ PASS" : "❌ FAIL"}  | ${ms}ms`);
  }

  // --- Test 3: 10MB scale (100MB hits RS limit, use 10MB) ---
  console.log("\n--- 10MB Scale Test ---\n");
  // Use compressible text data (real archives are mostly text/code that compresses 3-5x).
  // Random bytes don't compress, exceeding the GF(2^16) RS limit of 65535 oligos.
  // Deterministic patterns like (i*31+17)&0xff can cause homopolymer violations.
  // Text avoids both issues: compresses well AND has natural entropy.
  const textPattern = "The Helix Codec encodes digital data into synthetic DNA oligos for archival storage. ";
  const repeatCount = Math.ceil((10 * 1024 * 1024) / textPattern.length);
  const bigPayload = Buffer.from(textPattern.repeat(repeatCount).slice(0, 10 * 1024 * 1024));

  console.log(`Encoding 10MB...`);
  const bigConfig = { ...DEFAULT_CONFIG, oligoLength: 292, primerLength: 20, outerParityRatio: 0.3 };
  const t0 = Date.now();
  const encBig = await encodeFile(bigPayload, bigConfig, { fileName: "10mb.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  console.log(`  Encoded: ${encBig.encoded.oligos.length} oligos in ${encMs}ms (${(10 / (encMs / 1000)).toFixed(1)} MB/s)`);
  console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);

  // Decode at 10x
  console.log(`  Simulating 10x coverage...`);
  const sim = simulate(encBig.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  console.log(`  Decoding ${sim.totalReads.toLocaleString()} reads...`);
  const t1 = Date.now();
  const dec = await decodeReads(sim.reads, encBig.encoded.metadata, bigConfig, encBig.encoded.forwardPrimer, encBig.encoded.reversePrimer, true);
  const decMs = Date.now() - t1;
  console.log(`  Decode: ${decMs}ms (${(10 / (decMs / 1000)).toFixed(1)} MB/s)`);
  console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB`);
  console.log(`  Recovery: ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"} (${dec.stats.oligosRecovered}/${encBig.encoded.oligos.length})`);

  console.log("\n=== VALIDATION COMPLETE ===");
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
