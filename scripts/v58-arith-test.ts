// Test arithmetic mode end-to-end with hash verification
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import * as crypto from "crypto";
import { computeDensity } from "./presets";

async function test(label: string, payload: Uint8Array, cfg: any, coverage: number = 10) {
  const t0 = Date.now();
  const enc = await encodeFile(Buffer.from(payload), cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  const density = enc.stats.netDensityBitsPerNt;

  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage, seed: 42 });
  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;
  const tp = (payload.length / 1024 / 1024) / (decMs / 1000);
  console.log(`${label.padEnd(30)} | hash=${dec.hashMatches?"OK":"FAIL"} | density=${density.toFixed(3)} b/nt | enc=${encMs}ms dec=${decMs}ms (${tp.toFixed(2)} MB/s) | oligos=${dec.stats.oligosRecovered}`);
}

async function main() {
  console.log("=== v58 Arithmetic Mode E2E ===\n");
  const payload = crypto.randomBytes(8 * 1024); // 8KB

  // Arithmetic mode with low coverage trigger (forces JS path with LDPC erasure)
  const arithCfg = {
    oligoLength: 300, primerLength: 20,
    innerCode: "ldpc", ldpcDecoder: "auto", mappingMode: "arithmetic",
    innerParityBytes: 4, outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true, maxRetries: 1, interleaveDepth: 0,
    channel: "illumina", lowCoverageTrigger: 999, useConvolutionalInner: false,
  };

  // Test at multiple coverages
  await test("arith @ 10x", payload, arithCfg, 10);
  await test("arith @ 5x", payload, arithCfg, 5);
  await test("arith @ 3x", payload, arithCfg, 3);

  // Test with larger payload
  const bigPayload = crypto.randomBytes(64 * 1024);
  await test("arith 64KB @ 10x", bigPayload, arithCfg, 10);
}
main().catch(e => { console.error(e); process.exit(1); });
