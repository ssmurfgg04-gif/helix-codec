// Nanopore validation: test codec at high indel rates typical of Nanopore sequencing.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  console.log("=== Nanopore Validation ===\n");

  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  const subset = payload.slice(0, 65536);

  const config = { ...DEFAULT_CONFIG };
  console.log(`Config: innerCode=${config.innerCode}, ldpcDecoder=${config.ldpcDecoder}`);
  console.log(`Nanopore preset: sub=${(PRESET_NANOPORE.substitutionRate * 100).toFixed(1)}%, ins=${(PRESET_NANOPORE.insertionRate * 100).toFixed(1)}%, del=${(PRESET_NANOPORE.deletionRate * 100).toFixed(1)}%\n`);

  const enc = await encodeFile(subset, config, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos\n`);

  console.log("Cov  | Reads  | Recovery | Time  | Oligos OK");
  console.log("-----|--------|----------|-------|----------");
  for (const cov of [10, 20, 30]) {
    const sim = simulate(enc.encoded.oligos, { ...PRESET_NANOPORE, coverage: cov, seed: 42 });
    const t0 = Date.now();
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
    const ms = Date.now() - t0;
    console.log(`${cov}x  | ${sim.totalReads.toString().padStart(6)} | ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"}  | ${ms}ms | ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
  }

  // Also test with moderate indel rates (between Illumina and Nanopore)
  console.log("\n--- Moderate Indel Rates ---\n");
  for (const indelRate of [0.005, 0.01, 0.02]) {
    const sim = simulate(enc.encoded.oligos, {
      substitutionRate: 0.001,
      insertionRate: indelRate,
      deletionRate: indelRate,
      dropoutRate: 0.0,
      coverage: 20,
      seed: 42,
    });
    const t0 = Date.now();
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
    const ms = Date.now() - t0;
    console.log(`indel=${(indelRate * 100).toFixed(1)}%: ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"} | ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length} | ${ms}ms`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
