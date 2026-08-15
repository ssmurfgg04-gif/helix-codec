/**
 * v59 JS path test — force JS decode for 2.1MB v55-density.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";

const TAG = "[v59-js]";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  console.log(`${TAG} Payload: ${(payload.length / 1024 / 1024).toFixed(2)}MB`);

  const cfg = { ...ULTIMATE_V55_DENSITY_CONFIG, lowCoverageTrigger: 999 }; // force JS path
  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  console.log(`${TAG} Encode: ${encMs}ms = ${((payload.length / 1024 / 1024) / (encMs / 1000)).toFixed(2)} MB/s`);

  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const t1 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
  const decMs = Date.now() - t1;
  console.log(`${TAG} Decode (JS): ${decMs}ms = ${((payload.length / 1024 / 1024) / (decMs / 1000)).toFixed(2)} MB/s`);
  console.log(`${TAG} Hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
  console.log(`${TAG} Oligos recovered: ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
}
main().catch(e => { console.error(e); process.exit(1); });
