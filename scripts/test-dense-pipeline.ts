// Test the full pipeline with dense Goldman mapping.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  console.log("=== Dense Goldman Mapping Pipeline Test ===\n");

  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  const subset = payload.slice(0, 65536); // 64KB
  console.log(`Payload: ${subset.length} bytes\n`);

  // Test with dense Goldman mode (default)
  const config = { ...DEFAULT_CONFIG };
  console.log(`Config: mappingMode=${config.mappingMode}, goldmanMode=${config.goldmanMode}, oligoLength=${config.oligoLength}`);

  const enc = await encodeFile(subset, config, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos`);
  console.log(`  Payload per oligo: ${enc.stats.payloadBytesPerOligo} bytes`);
  console.log(`  Density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`  Screening retries: ${enc.stats.screeningRetries}\n`);

  console.log("Cov  | Reads  | Recovery | Time  | Oligos OK");
  console.log("-----|--------|----------|-------|----------");
  for (const cov of [5, 10, 15, 20]) {
    const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: cov, seed: 42 });
    const t0 = Date.now();
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
    const ms = Date.now() - t0;
    console.log(`${cov}x  | ${sim.totalReads.toString().padStart(6)} | ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"}  | ${ms}ms | ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
