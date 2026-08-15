// Benchmark: Arithmetic coding density + parallel decode
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { decodeReadsParallel } from "../src/lib/dna/ultra-decode-parallel";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";
import * as os from "os";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 262144);

  // === Test 1: Direct mapping (baseline) ===
  console.log("=== Test 1: Direct mapping (baseline) ===");
  const directConfig = { ...DEFAULT_CONFIG, mappingMode: "direct" as const };
  const enc1 = await encodeFile(payload, directConfig, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`  Oligos: ${enc1.encoded.oligos.length}, Density: ${enc1.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);

  const sim1 = simulate(enc1.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const t1 = Date.now();
  const dec1 = await decodeReadsUltra(sim1.reads, enc1.encoded.metadata, directConfig, enc1.encoded.forwardPrimer, enc1.encoded.reversePrimer);
  console.log(`  Decode: ${Date.now() - t1}ms, Result: ${dec1.hashMatches ? "PASS" : "FAIL"}`);

  // === Test 2: Arithmetic mapping (higher density) ===
  console.log("\n=== Test 2: Arithmetic mapping (higher density) ===");
  const arithConfig = { ...DEFAULT_CONFIG, mappingMode: "arithmetic" as const };
  let enc2;
  try {
    enc2 = await encodeFile(payload, arithConfig, { fileName: "test.bin", contentType: "application/octet-stream" });
    console.log(`  Oligos: ${enc2.encoded.oligos.length}, Density: ${enc2.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
    console.log(`  Density improvement: ${((enc2.stats.netDensityBitsPerNt / enc1.stats.netDensityBitsPerNt - 1) * 100).toFixed(1)}%`);
  } catch (e: any) {
    console.log(`  Encode FAILED: ${e.message}`);
    return;
  }

  const sim2 = simulate(enc2.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const t2 = Date.now();
  const dec2 = await decodeReadsUltra(sim2.reads, enc2.encoded.metadata, arithConfig, enc2.encoded.forwardPrimer, enc2.encoded.reversePrimer);
  console.log(`  Decode: ${Date.now() - t2}ms, Result: ${dec2.hashMatches ? "PASS" : "FAIL"}`);

  // === Test 3: Parallel decode (direct mapping, multi-core) ===
  console.log("\n=== Test 3: Parallel decode (direct mapping) ===");
  const numWorkers = Math.min(os.cpus().length, 8);
  console.log(`  Workers: ${numWorkers} (CPU cores: ${os.cpus().length})`);

  const t3 = Date.now();
  const dec3 = await decodeReadsParallel(sim1.reads, enc1.encoded.metadata, directConfig, enc1.encoded.forwardPrimer, enc1.encoded.reversePrimer, numWorkers);
  console.log(`  Decode: ${Date.now() - t3}ms, Result: ${dec3.hashMatches ? "PASS" : "FAIL"}`);

  // === Summary ===
  console.log("\n=== SUMMARY ===");
  console.log(`  Direct decode (single-thread):  ${dec1.stats.decodeTimeMs}ms, ${dec1.hashMatches ? "PASS" : "FAIL"}, ${enc1.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`  Arithmetic decode (single):     ${dec2.stats.decodeTimeMs}ms, ${dec2.hashMatches ? "PASS" : "FAIL"}, ${enc2.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`  Parallel decode (${numWorkers} workers): ${dec3.stats.decodeTimeMs}ms, ${dec3.hashMatches ? "PASS" : "FAIL"}, ${enc1.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);

  if (dec3.stats.decodeTimeMs > 0) {
    console.log(`  Parallel speedup: ${(dec1.stats.decodeTimeMs / dec3.stats.decodeTimeMs).toFixed(1)}x`);
  }
}

main().catch((e) => { console.error(e); process.exit(1); });
