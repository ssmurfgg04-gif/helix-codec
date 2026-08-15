// Benchmark: Ultra WASM decode vs original JS decode
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 262144);
  const config = { ...DEFAULT_CONFIG };

  const t0 = Date.now();
  const enc = await encodeFile(payload, config, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  console.log(`Encode: ${encMs}ms`);
  console.log(`Oligos: ${enc.encoded.oligos.length}, Reads: ${enc.encoded.oligos.length * 10}`);

  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });

  // Original JS decode
  const t1 = Date.now();
  const dec1 = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
  const dec1Ms = Date.now() - t1;
  console.log(`\nJS decode:   ${dec1Ms}ms (${(payload.length / 1024 / 1024 / (dec1Ms / 1000)).toFixed(1)} MB/s)`);
  console.log(`  Result: ${dec1.hashMatches ? "PASS" : "FAIL"}`);
  console.log(`  Per-read: ${(dec1Ms / sim.totalReads).toFixed(4)}ms`);

  // Ultra WASM decode
  const t2 = Date.now();
  const dec2 = await decodeReadsUltra(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const dec2Ms = Date.now() - t2;
  console.log(`\nWASM decode: ${dec2Ms}ms (${(payload.length / 1024 / 1024 / (dec2Ms / 1000)).toFixed(1)} MB/s)`);
  console.log(`  Result: ${dec2.hashMatches ? "PASS" : "FAIL"}`);
  console.log(`  Per-read: ${(dec2Ms / sim.totalReads).toFixed(4)}ms`);

  // Speedup
  console.log(`\n=== SPEEDUP ===`);
  console.log(`Decode: ${dec1Ms}ms → ${dec2Ms}ms = ${(dec1Ms / dec2Ms).toFixed(1)}x faster`);
  console.log(`Total:  ${encMs + dec1Ms}ms → ${encMs + dec2Ms}ms = ${((encMs + dec1Ms) / (encMs + dec2Ms)).toFixed(1)}x faster`);

  // vs v34 baseline
  const v34enc = 5290, v34dec = 2765;
  console.log(`\n=== vs v34 BASELINE ===`);
  console.log(`Encode: ${v34enc}ms → ${encMs}ms = ${(v34enc / encMs).toFixed(1)}x`);
  console.log(`Decode: ${v34dec}ms → ${dec2Ms}ms = ${(v34dec / dec2Ms).toFixed(1)}x`);
  console.log(`Total:  ${v34enc + v34dec}ms → ${encMs + dec2Ms}ms = ${((v34enc + v34dec) / (encMs + dec2Ms)).toFixed(1)}x`);
}

main().catch((e) => { console.error(e); process.exit(1); });
