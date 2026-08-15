/**
 * v59 v55-density 256KB test.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";

const TAG = "[v59-256]";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 262144);
  console.log(`${TAG} Payload: ${(payload.length / 1024).toFixed(0)}KB`);

  const cfg = ULTIMATE_V55_DENSITY_CONFIG;
  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  console.log(`${TAG} Encode: ${encMs}ms = ${((payload.length / 1024 / 1024) / (encMs / 1000)).toFixed(2)} MB/s`);
  console.log(`${TAG} Oligos: ${enc.encoded.oligos.length}, density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);

  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;
  console.log(`${TAG} Decode: ${decMs}ms = ${((payload.length / 1024 / 1024) / (decMs / 1000)).toFixed(2)} MB/s`);
  console.log(`${TAG} Hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
}
main().catch(e => { console.error(e); process.exit(1); });
