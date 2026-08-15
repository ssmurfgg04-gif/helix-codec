// Trace arithmetic mode decode to find the failure point
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { computeLayoutAuto } from "../src/lib/dna/types";
import * as crypto from "crypto";

const payload = crypto.randomBytes(1024);
const cfg = {
  oligoLength: 700, primerLength: 20,
  innerCode: "ldpc", ldpcDecoder: "auto", mappingMode: "arithmetic",
  innerParityBytes: 4, outerParityRatio: 0.03,
  constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
  compress: true, maxRetries: 1, interleaveDepth: 0,
  channel: "illumina", lowCoverageTrigger: 999, useConvolutionalInner: false,
};

const enc = await encodeFile(Buffer.from(payload), cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
console.log("oligos:", enc.encoded.oligos.length, "density:", enc.stats.netDensityBitsPerNt.toFixed(3));
const layout = computeLayoutAuto(cfg);
console.log("layout:", JSON.stringify(layout));
console.log("innerN:", layout.addressBytes + layout.payloadBytes + layout.innerParityBytes, "totalInnerBytes:", layout.totalInnerBytes);

const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
console.log("reads:", sim.reads.length);

const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
console.log("hash:", dec.hashMatches ? "OK" : "FAIL");
console.log("stats:", JSON.stringify(dec.stats));
console.log("per-oligo (first 5):");
for (let i = 0; i < Math.min(5, dec.perOligo.length); i++) {
  const p = dec.perOligo[i];
  console.log(`  oligo[${i}]: reads=${p.readCount} crcPassed=${p.crcPassed} corrected=${p.innerRS?.corrected} success=${p.innerRS?.success}`);
}
console.log("first 16 of decoded:", Array.from(dec.data?.slice(0,16) ?? []).map(b=>b.toString(16).padStart(2,'0')).join(' '));
console.log("first 16 expected:  ", Array.from(payload.slice(0,16)).map(b=>b.toString(16).padStart(2,'0')).join(' '));
