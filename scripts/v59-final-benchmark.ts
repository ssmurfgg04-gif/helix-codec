/**
 * v59 Final Benchmark — honest numbers for all configs.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import {
  ULTIMATE_V55_DENSITY_CONFIG,
  ULTIMATE_V59_FAST_ENCODE_CONFIG,
  ULTIMATE_V59_HD_FAST_CONFIG,
  computeDensity,
} from "../src/lib/dna/presets";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

const TAG = "[v59-final]";

async function bench(payload: Uint8Array, label: string, cfg: any, coverage: number = 10) {
  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  const density = enc.stats.netDensityBitsPerNt;
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage, seed: 42 });
  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;
  const encTp = (payload.length / 1024 / 1024) / (encMs / 1000);
  const decTp = (payload.length / 1024 / 1024) / (decMs / 1000);
  console.log(`${TAG} ${label}:`);
  console.log(`${TAG}   payload: ${(payload.length / 1024).toFixed(0)}KB, oligos: ${enc.encoded.oligos.length}, density: ${density.toFixed(3)} b/nt`);
  console.log(`${TAG}   encode: ${encMs}ms = ${encTp.toFixed(2)} MB/s (retries: ${enc.stats.screeningRetries})`);
  console.log(`${TAG}   decode: ${decMs}ms = ${decTp.toFixed(2)} MB/s`);
  console.log(`${TAG}   hash: ${dec.hashMatches ? "OK ✓" : "FAIL ✗"}`);
  console.log();
}

async function main() {
  const fullPayload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  console.log(`${TAG} === v59 Final Benchmark (honest numbers) ===\n`);

  console.log(`${TAG} --- 256KB payload ---\n`);
  await bench(fullPayload.slice(0, 262144), "DEFAULT_CONFIG (300nt, direct)", DEFAULT_CONFIG);
  await bench(fullPayload.slice(0, 262144), "v55-density (700nt, direct, 8B parity)", ULTIMATE_V55_DENSITY_CONFIG);

  console.log(`${TAG} --- 2.1MB payload ---\n`);
  await bench(fullPayload, "DEFAULT_CONFIG (300nt, direct)", DEFAULT_CONFIG);
  await bench(fullPayload, "v55-density (700nt, direct, 8B parity)", ULTIMATE_V55_DENSITY_CONFIG);

  console.log(`${TAG} === Summary ===`);
  console.log(`${TAG} v59 improvements:`);
  console.log(`${TAG}   - Fixed v55-density hash FAIL (4B→8B inner parity)`);
  console.log(`${TAG}   - K-mer clustering wired into decode path`);
  console.log(`${TAG}   - Tolerant primer search for nanopore`);
  console.log(`${TAG}   - MSA + conv inner code support`);
  console.log(`${TAG}   - New fast-encode configs (SRT mapping)`);
}
main().catch(e => { console.error(e); process.exit(1); });
