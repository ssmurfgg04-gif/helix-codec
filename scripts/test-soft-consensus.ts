// Test soft-consensus at various error rates.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 16384);
  const config = { ...DEFAULT_CONFIG };
  const enc = await encodeFile(payload, config, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log("Encoded:", enc.encoded.oligos.length, "oligos");

  for (const subRate of [0.001, 0.005, 0.01]) {
    const sim = simulate(enc.encoded.oligos, {
      substitutionRate: subRate,
      insertionRate: subRate * 0.5,
      deletionRate: subRate,
      dropoutRate: 0.0,
      coverage: 20,
      seed: 42,
    });
    const t0 = Date.now();
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
    console.log(`sub=${(subRate * 100).toFixed(1)}%: ${dec.hashMatches ? "PASS" : "FAIL"} ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length} ${Date.now() - t0}ms failed=${dec.stats.oligosFailedInnerRS}`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
