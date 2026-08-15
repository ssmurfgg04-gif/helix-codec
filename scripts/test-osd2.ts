// Test OSD-2 soft-decision decoding at elevated error rates.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  console.log("=== OSD-2 Soft-Decision Decoding Test ===\n");

  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  const subset = payload.slice(0, 65536); // 64KB

  const testConfig = { ...DEFAULT_CONFIG, oligoLength: 292, primerLength: 20, outerParityRatio: 0.3 };

  console.log(`Config: innerCode=${testConfig.innerCode}, mappingMode=${testConfig.mappingMode}, oligoLength=${testConfig.oligoLength}`);
  console.log(`  (OSD-2 is automatically triggered when hard-decision LDPC fails)\n`);

  const enc = await encodeFile(subset, testConfig, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos\n`);

  // Test at various error rates — OSD-2 should help at higher rates
  console.log("Sub Rate | Coverage | Recovery | Time  | Oligos OK");
  console.log("---------|----------|----------|-------|----------");

  for (const subRate of [0.001, 0.005, 0.01, 0.02, 0.05]) {
    for (const cov of [10, 20]) {
      const sim = simulate(enc.encoded.oligos, {
        substitutionRate: subRate,
        insertionRate: subRate * 0.5,
        deletionRate: subRate,
        dropoutRate: 0.0,
        coverage: cov,
        seed: 42,
      });
      const t0 = Date.now();
      const dec = await decodeReads(sim.reads, enc.encoded.metadata, testConfig, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
      const ms = Date.now() - t0;
      const subPct = (subRate * 100).toFixed(1) + "%";
      console.log(`${subPct.padStart(8)} | ${cov}x       | ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"}  | ${ms}ms | ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
    }
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
