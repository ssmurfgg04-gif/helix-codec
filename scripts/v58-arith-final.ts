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
  const tp = (payload.length / 1024 / 1024) / (decMs / 1000);
  console.log(`${label.padEnd(30)} | hash=${dec.hashMatches?"OK":"FAIL"} | density=${density.toFixed(3)} b/nt | dec=${decMs}ms (${tp.toFixed(2)} MB/s) | oligos=${enc.encoded.oligos.length}`);
}

const payload = crypto.randomBytes(8 * 1024);

console.log("=== v58 Arithmetic Mode (innerParityBytes=8) ===\n");
await test("arith 8KB @ 10x", payload, ULTIMATE_V55_DENSITY_CONFIG, 10);
await test("arith 8KB @ 5x", payload, ULTIMATE_V55_DENSITY_CONFIG, 5);
await test("arith 8KB @ 3x", payload, ULTIMATE_V55_DENSITY_CONFIG, 3);

const bigPayload = crypto.randomBytes(256 * 1024);
await test("arith 256KB @ 10x", bigPayload, ULTIMATE_V55_DENSITY_CONFIG, 10);
