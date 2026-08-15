import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import * as crypto from "crypto";

async function test(label: string, payload: Uint8Array, cfg: any, coverage: number) {
  const t0 = Date.now();
  const enc = await encodeFile(Buffer.from(payload), cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const density = enc.stats.netDensityBitsPerNt;
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage, seed: 42 });
  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;
  console.log(`${label.padEnd(35)} | hash=${dec.hashMatches?"OK":"FAIL"} | density=${density.toFixed(3)} b/nt | oligos=${enc.encoded.oligos.length} | dec=${decMs}ms`);
}

const payload = crypto.randomBytes(64 * 1024);
console.log("=== v58 Arithmetic — larger payload (more oligos for outer RS) ===\n");

await test("64KB @ 10x", payload, ULTIMATE_V55_DENSITY_CONFIG, 10);
await test("64KB @ 15x", payload, ULTIMATE_V55_DENSITY_CONFIG, 15);
await test("64KB @ 20x", payload, ULTIMATE_V55_DENSITY_CONFIG, 20);
