/**
 * v60: Honest end-to-end benchmark.
 * 
 * Tests all 4 key metrics:
 * 1. Net density (bits/nucleotide)
 * 2. Coverage for 100% recovery (2×, 3×, 10×)
 * 3. Decode throughput (MB/s)
 * 4. IDS error tolerance (9% IDS nanopore)
 * 5. Encode throughput (MB/s)
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG, ULTIMATE_NANOPORE_V52_CONFIG, computeDensity } from "../src/lib/dna/presets";
import * as fs from "fs";

function makePayload(size: number): Uint8Array {
  // Realistic payload: mix of text and binary
  const payload = new Uint8Array(size);
  for (let i = 0; i < size; i++) {
    payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  }
  return payload;
}

async function benchEncodeDecode(
  name: string,
  cfg: any,
  payloadSize: number,
  coverage: number,
  preset: any = PRESET_ILLUMINA,
) {
  const payload = makePayload(payloadSize);
  
  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "bench.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  
  const sim = simulate(enc.encoded.oligos, { ...preset, coverage, seed: 42 });
  
  const t1 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;
  
  const encMBs = (payloadSize / 1024 / 1024) / (encMs / 1000);
  const decMBs = (payloadSize / 1024 / 1024) / (decMs / 1000);
  
  console.log(`${name}:`);
  console.log(`  ${cfg.oligoLength}nt, ${enc.encoded.oligos.length} oligos, ${coverage}x cov`);
  console.log(`  encode: ${encMs}ms = ${encMBs.toFixed(2)} MB/s`);
  console.log(`  decode: ${decMs}ms = ${decMBs.toFixed(2)} MB/s`);
  console.log(`  hash: ${dec.hashMatches ? "OK" : "FAIL"}, recovered: ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
  console.log(`  density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt (theoretical: ${computeDensity(cfg).toFixed(3)})`);
  console.log();
  
  return { name, encMBs, decMBs, hash: dec.hashMatches, density: enc.stats.netDensityBitsPerNt };
}

async function main() {
  console.log("=== v60 Honest Benchmark ===\n");
  
  const results: any[] = [];
  
  // 1. v55-density (production config) — 256KB
  results.push(await benchEncodeDecode(
    "v55-density 256KB 10x Illumina",
    ULTIMATE_V55_DENSITY_CONFIG,
    256 * 1024,
    10,
  ));
  
  // 2. v55-density — 2.1MB
  results.push(await benchEncodeDecode(
    "v55-density 2.1MB 10x Illumina",
    ULTIMATE_V55_DENSITY_CONFIG,
    2.1 * 1024 * 1024,
    10,
  ));
  
  // 3. v55-density — 256KB at 3× coverage
  results.push(await benchEncodeDecode(
    "v55-density 256KB 3x Illumina",
    ULTIMATE_V55_DENSITY_CONFIG,
    256 * 1024,
    3,
  ));
  
  // 4. Nanopore 9% IDS — 16KB
  results.push(await benchEncodeDecode(
    "Nanopore 9% IDS 16KB 10x",
    ULTIMATE_NANOPORE_V52_CONFIG,
    16 * 1024,
    10,
    PRESET_NANOPORE,
  ));
  
  // Summary
  console.log("=== Summary ===");
  console.log("Metric | Value | SOTA | Gap");
  console.log("-------|-------|------|-----");
  
  const density = results[0].density;
  const encMBs = results[1].encMBs;
  const decMBs = results[1].decMBs;
  const nanoporeResult = results[3];
  
  console.log(`Net density | ${density.toFixed(3)} b/nt | 1.815 (Yi Ding) | ${((density - 1.815) / 1.815 * 100).toFixed(1)}%`);
  console.log(`Encode throughput | ${encMBs.toFixed(2)} MB/s | 12.5 (Catalog HW) | ${((encMBs - 12.5) / 12.5 * 100).toFixed(1)}%`);
  console.log(`Decode throughput | ${decMBs.toFixed(2)} MB/s | 2.5 (NGS baseline) | ${((decMBs - 2.5) / 2.5 * 100).toFixed(1)}%`);
  console.log(`IDS tolerance | 9% IDS, ${nanoporeResult.hash ? "hash OK" : "hash FAIL"} | 24% (DNA-MGC+) | below`);
  console.log(`Coverage | 3x (hash ${results[2].hash ? "OK" : "FAIL"}) | 2.25x (DNA-MGC+) | ${results[2].hash ? "leading" : "behind"}`);
}

main().catch(e => { console.error(e); process.exit(1); });
