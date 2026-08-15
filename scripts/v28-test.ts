// Test the triage + neural LDPC + attention consensus integration.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  console.log("=== v28 Integration Test ===\n");

  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 65536);
  const config = { ...DEFAULT_CONFIG };
  console.log(`Config: ${config.oligoLength}nt, ${config.mappingMode}, ${config.innerCode}/${config.ldpcDecoder}, ${(config.outerParityRatio * 100).toFixed(0)}% parity\n`);

  const enc = await encodeFile(payload, config, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos`);
  console.log(`Net density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`Payload density: ${(enc.stats.payloadBytesPerOligo * 8 / (config.oligoLength - 2 * config.primerLength)).toFixed(3)} bits/nt\n`);

  console.log("Cov  | Reads  | Recovery | Time  | Oligos OK | Failed");
  console.log("-----|--------|----------|-------|----------|------");
  for (const cov of [5, 10, 20]) {
    const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: cov, seed: 42 });
    const t0 = Date.now();
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
    const ms = Date.now() - t0;
    console.log(`${cov}x  | ${sim.totalReads.toString().padStart(6)} | ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"}  | ${ms}ms | ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}    | ${dec.stats.oligosFailedInnerRS}`);
  }

  // Test full 2.1MB Erlich
  console.log("\n--- Full 2.1MB Erlich ---");
  const fullPayload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  const t0 = Date.now();
  const fullEnc = await encodeFile(fullPayload, config, { fileName: "erlich.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${fullEnc.encoded.oligos.length.toLocaleString()} oligos in ${((Date.now()-t0)/1000).toFixed(1)}s`);

  const sim = simulate(fullEnc.encoded.oligos, {
    substitutionRate: 0.000167, insertionRate: 0.0001, deletionRate: 0.0001,
    dropoutRate: 0.0, coverage: 10, seed: 42,
  });
  const t1 = Date.now();
  const dec = await decodeReads(sim.reads, fullEnc.encoded.metadata, config, fullEnc.encoded.forwardPrimer, fullEnc.encoded.reversePrimer, true);
  console.log(`Decoded in ${((Date.now()-t1)/1000).toFixed(1)}s`);
  console.log(`Hash: ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"}`);
  console.log(`Recovered: ${dec.stats.oligosRecovered.toLocaleString()}/${fullEnc.encoded.oligos.length.toLocaleString()}`);
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
