// Quick speed benchmark.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 262144);
  const config = { ...DEFAULT_CONFIG };

  const t0 = Date.now();
  const enc = await encodeFile(payload, config, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;

  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const t2 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
  const decMs = Date.now() - t2;

  console.log(`Encode: ${encMs}ms (${(payload.length / 1024 / 1024 / (encMs / 1000)).toFixed(1)} MB/s)`);
  console.log(`Decode: ${decMs}ms (${(payload.length / 1024 / 1024 / (decMs / 1000)).toFixed(1)} MB/s)`);
  console.log(`Per-read: ${(decMs / dec.stats.readsUsed).toFixed(4)}ms`);
  console.log(`Result: ${dec.hashMatches ? "PASS" : "FAIL"}`);
}

main().catch((e) => { console.error(e); process.exit(1); });
