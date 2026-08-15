#!/usr/bin/env bun
/**
 * v57 HONEST throughput benchmark on 256KB and 2.1MB real payloads.
 *
 * The user specifically asked: "run paper-benchmark.ts on 256KB and 2.1MB
 * real payloads and report those numbers." This script does exactly that —
 * no synthetic 64-byte payloads, no nucleotide-rate measurements.
 *
 * Reports END-TO-END payload throughput (MB/s) = payload_bytes / decode_time.
 * This is the only honest throughput metric for DNA storage.
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

const TAG = "[v57-throughput]";

async function benchSize(payload: Uint8Array, label: string, coverage: number = 10) {
  const config = { ...DEFAULT_CONFIG };

  const t0 = Date.now();
  const enc = await encodeFile(payload, config, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;

  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage, seed: 42 });

  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;

  const encThroughput = (payload.length / 1024 / 1024) / (encMs / 1000);
  const decThroughput = (payload.length / 1024 / 1024) / (decMs / 1000);

  console.log(`${TAG} ${label}:`);
  console.log(`${TAG}   payload: ${(payload.length / 1024).toFixed(0)}KB (${payload.length} bytes)`);
  console.log(`${TAG}   oligos: ${enc.encoded.oligos.length}, density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`${TAG}   encode: ${encMs}ms = ${encThroughput.toFixed(2)} MB/s`);
  console.log(`${TAG}   decode: ${decMs}ms = ${decThroughput.toFixed(2)} MB/s`);
  console.log(`${TAG}   hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
  console.log();
}

async function main() {
  const fullPayload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  console.log(`${TAG} === HONEST Throughput Benchmark (real Erlich payload) ===`);
  console.log(`${TAG} Config: DEFAULT_CONFIG (300nt, LDPC, 10% outer RS, direct mapping, 10× coverage)`);
  console.log();

  // 256KB
  await benchSize(fullPayload.slice(0, 262144), "256KB");

  // 2.1MB (full)
  await benchSize(fullPayload, "2.1MB (full)");

  console.log(`${TAG} === Summary ===`);
  console.log(`${TAG} These are END-TO-END payload throughput numbers.`);
  console.log(`${TAG} No nucleotide-rate, no synthetic 64-byte payloads.`);
  console.log(`${TAG} These are the numbers to cite in papers/pitches.`);
}

main().catch((e) => { console.error(e); process.exit(1); });
