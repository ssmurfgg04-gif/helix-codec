// Comprehensive benchmark for paper — all real numbers
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  const fullPayload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  console.log("=== HELIX CODEC — COMPREHENSIVE BENCHMARK ===\n");
  console.log(`Payload sizes: 64KB, 256KB, 1MB, 2.1MB`);
  console.log(`Config: oligoLength=300, primerLength=20, LDPC(62,58), outerRS(10%), direct mapping\n`);

  // === Throughput benchmark ===
  console.log("=== THROUGHPUT BENCHMARK ===\n");
  console.log("Payload | Encode | Decode | Enc MB/s | Dec MB/s | Oligos | Density | Result");
  console.log("--------|--------|--------|----------|----------|--------|---------|-------");

  for (const size of [65536, 262144, 1048576, 2116608]) {
    const payload = fullPayload.slice(0, size);
    const sizeLabel = size >= 1048576 ? `${(size / 1048576).toFixed(1)}MB` : `${(size / 1024).toFixed(0)}KB`;
    const config = { ...DEFAULT_CONFIG };

    const t0 = Date.now();
    const enc = await encodeFile(payload, config, { fileName: "test.bin", contentType: "application/octet-stream" });
    const encMs = Date.now() - t0;

    const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });

    const t1 = Date.now();
    const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const decMs = Date.now() - t1;

    const encThroughput = (size / 1024 / 1024) / (encMs / 1000);
    const decThroughput = (size / 1024 / 1024) / (decMs / 1000);

    console.log(`${sizeLabel.padEnd(7)} | ${encMs.toString().padStart(5)}ms | ${decMs.toString().padStart(5)}ms | ${encThroughput.toFixed(2).padStart(7)}  | ${decThroughput.toFixed(2).padStart(7)}  | ${enc.encoded.oligos.length.toString().padStart(6)} | ${enc.stats.netDensityBitsPerNt.toFixed(3)}   | ${dec.hashMatches ? "PASS" : "FAIL"}`);
  }

  // === Coverage sweep (Erlich 2.1MB) ===
  console.log("\n=== COVERAGE SWEEP (Erlich 2.1MB, empirical noise) ===\n");
  const erlichConfig = { ...DEFAULT_CONFIG, oligoLength: 300, primerLength: 20, outerParityRatio: 0.1, maxRetries: 1 };
  const enc = await encodeFile(fullPayload, erlichConfig, { fileName: "erlich.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos, density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt\n`);

  console.log("Coverage | Reads    | Decode  | MB/s   | Result");
  console.log("---------|----------|---------|--------|-------");

  for (const coverage of [3, 5, 10, 20]) {
    const sim = simulate(enc.encoded.oligos, {
      substitutionRate: 0.000176, insertionRate: 0.0001, deletionRate: 0.0001,
      dropoutRate: 0,
      coverage, seed: 42,
    });
    const t0 = Date.now();
    const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, erlichConfig, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const decMs = Date.now() - t0;
    const throughput = (fullPayload.length / 1024 / 1024) / (decMs / 1000);
    console.log(`${coverage.toString().padStart(8)}x | ${sim.totalReads.toString().padStart(8)} | ${decMs.toString().padStart(5)}ms | ${throughput.toFixed(2).padStart(6)} | ${dec.hashMatches ? "PASS" : "FAIL"}`);
  }

  // === Interleaving benchmark ===
  console.log("\n=== INTERLEAVING BENCHMARK (256KB, 10x) ===\n");
  const payload256 = fullPayload.slice(0, 262144);
  console.log("Depth | Perfect | 0.0176% sub | 0.1% sub | 0.5% sub");
  console.log("------|---------|-------------|----------|----------");

  for (const depth of [0, 2, 4]) {
    const config = { ...DEFAULT_CONFIG, interleaveDepth: depth };
    const enc = await encodeFile(payload256, config, { fileName: "test.bin", contentType: "application/octet-stream" });

    // Perfect reads
    const perfectReads = enc.encoded.oligos.map(o => ({ sequence: o.sequence, quality: [] as number[], oligoIndex: 0, substitutions: 0, insertions: 0, deletions: 0 }));
    const dec1 = await decodeReadsUltra(perfectReads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);

    // Erlich noise
    const sim1 = simulate(enc.encoded.oligos, { substitutionRate: 0.000176, insertionRate: 0.0001, deletionRate: 0.0001, dropoutRate: 0, coverage: 10, seed: 42 });
    const dec2 = await decodeReadsUltra(sim1.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);

    // 0.1% sub
    const sim2 = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
    const dec3 = await decodeReadsUltra(sim2.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);

    // 0.5% sub
    const sim3 = simulate(enc.encoded.oligos, { substitutionRate: 0.005, insertionRate: 0.001, deletionRate: 0.001, dropoutRate: 0, coverage: 10, seed: 42 });
    const dec4 = await decodeReadsUltra(sim3.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);

    console.log(`${depth.toString().padStart(5)} | ${dec1.hashMatches ? "PASS   " : "FAIL   "} | ${dec2.hashMatches ? "PASS       " : "FAIL       "} | ${dec3.hashMatches ? "PASS    " : "FAIL    "} | ${dec4.hashMatches ? "PASS    " : "FAIL    "}`);
  }

  // === 8-byte LDPC parity benchmark ===
  console.log("\n=== LDPC PARITY BENCHMARK (256KB, 10x) ===\n");
  console.log("Parity | Density | 0.1% sub | 0.5% sub");
  console.log("-------|---------|----------|----------");
  for (const parity of [4, 8]) {
    const config = { ...DEFAULT_CONFIG, innerParityBytes: parity };
    const enc = await encodeFile(payload256, config, { fileName: "test.bin", contentType: "application/octet-stream" });
    const sim1 = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
    const dec1 = await decodeReadsUltra(sim1.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const sim2 = simulate(enc.encoded.oligos, { substitutionRate: 0.005, insertionRate: 0.001, deletionRate: 0.001, dropoutRate: 0, coverage: 10, seed: 42 });
    const dec2 = await decodeReadsUltra(sim2.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    console.log(`${parity.toString().padStart(6)}B | ${enc.stats.netDensityBitsPerNt.toFixed(3)}   | ${dec1.hashMatches ? "PASS    " : "FAIL    "} | ${dec2.hashMatches ? "PASS    " : "FAIL    "}`);
  }

  // === Encryption benchmark ===
  console.log("\n=== ENCRYPTION BENCHMARK (256KB, 10x) ===\n");
  const encConfig = { ...DEFAULT_CONFIG, encryptPassword: "benchmark-test" };
  const encEnc = await encodeFile(payload256, encConfig, { fileName: "test.bin", contentType: "application/octet-stream" });
  const simEnc = simulate(encEnc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const t0 = Date.now();
  const decEnc = await decodeReadsUltra(simEnc.reads, encEnc.encoded.metadata, encConfig, encEnc.encoded.forwardPrimer, encEnc.encoded.reversePrimer);
  const encMs = Date.now() - t0;
  console.log(`Encrypted decode: ${encMs}ms, ${decEnc.hashMatches ? "PASS" : "FAIL"}`);

  console.log("\n=== BENCHMARK COMPLETE ===");
}

main().catch((e) => { console.error(e); process.exit(1); });
