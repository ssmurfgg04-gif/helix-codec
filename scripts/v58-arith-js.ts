import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import * as crypto from "crypto";

async function test(label: string, payload: Uint8Array, cfg: any, coverage: number) {
  const enc = await encodeFile(Buffer.from(payload), cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage, seed: 42 });
  const t0 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
  const ms = Date.now() - t0;
  console.log(`${label.padEnd(30)} | hash=${dec.hashMatches?"OK":"FAIL"} | density=${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt | dec=${ms}ms | oligos=${dec.stats.oligosRecovered}`);
}

const payload = crypto.randomBytes(16 * 1024);
const cfg = {
  oligoLength: 700, primerLength: 20,
  innerCode: "ldpc", ldpcDecoder: "auto", mappingMode: "arithmetic",
  innerParityBytes: 8, outerParityRatio: 0.05,
  constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
  compress: true, maxRetries: 1, interleaveDepth: 0,
  channel: "illumina", lowCoverageTrigger: 999, useConvolutionalInner: false,
};

console.log("=== v58 Arithmetic — JS path (lowCoverageTrigger=999) ===\n");
await test("16KB arith JS @ 10x", payload, cfg, 10);
await test("16KB arith JS @ 20x", payload, cfg, 20);
