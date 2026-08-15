import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import * as crypto from "crypto";

async function test(label: string, payload: Uint8Array, cfg: any, coverage: number = 10) {
  const t0 = Date.now();
  const enc = await encodeFile(Buffer.from(payload), cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  const density = enc.stats.netDensityBitsPerNt;
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage, seed: 42 });
  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;
  console.log(`${label.padEnd(30)} | hash=${dec.hashMatches?"OK":"FAIL"} | density=${density.toFixed(3)} b/nt | enc=${encMs}ms dec=${decMs}ms | oligos=${enc.encoded.oligos.length}`);
}

const payload = crypto.randomBytes(8 * 1024);

// 700nt + arithmetic + 4B LDPC + 3% outer RS (matches ULTIMATE_V55 but with arithmetic)
const cfg = {
  oligoLength: 700, primerLength: 20,
  innerCode: "ldpc", ldpcDecoder: "auto", mappingMode: "arithmetic",
  innerParityBytes: 4, outerParityRatio: 0.03,
  constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
  compress: true, maxRetries: 1, interleaveDepth: 0,
  channel: "illumina", lowCoverageTrigger: 999, useConvolutionalInner: false,
};

await test("700nt arith @ 10x", payload, cfg, 10);
await test("700nt arith @ 5x", payload, cfg, 5);
await test("700nt arith @ 3x", payload, cfg, 3);
