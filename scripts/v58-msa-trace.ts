// Trace which strategies are firing and where the failures are
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import * as crypto from "crypto";

const payload = crypto.randomBytes(4 * 1024); // 4KB

const cfg = {
  oligoLength: 300, primerLength: 20,
  innerCode: "ldpc", ldpcDecoder: "auto", mappingMode: "direct",
  innerParityBytes: 8, outerParityRatio: 0.30,
  constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
  compress: true, maxRetries: 1, interleaveDepth: 0,
  channel: "nanopore", lowCoverageTrigger: 999, useConvolutionalInner: false,
};

const enc = await encodeFile(Buffer.from(payload), cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
console.log(`oligos: ${enc.encoded.oligos.length}, outerN=${enc.encoded.metadata.outerRS.n}, outerK=${enc.encoded.metadata.outerRS.k}, nsym=${enc.encoded.metadata.outerRS.n - enc.encoded.metadata.outerRS.k}`);

// 5% IDS at 10x — easier case
const noise = { ...PRESET_NANOPORE, substitutionRate: 0.01, insertionRate: 0.015, deletionRate: 0.025, coverage: 10, seed: 42 };
const sim = simulate(enc.encoded.oligos, noise);
console.log(`reads: ${sim.reads.length}, avg cluster: ${(sim.reads.length/enc.encoded.oligos.length).toFixed(2)}`);

// Per-oligo: how many reads, and what's the avg length variation?
const oligoReads = new Map<number, any[]>();
for (const r of sim.reads) {
  // Just group by index 0..N-1 using r.address or position
}
// Actually the reads don't carry oligo index — let me check
console.log("first read:", { seq: sim.reads[0].sequence.slice(0,50)+"...", len: sim.reads[0].sequence.length, hasQ: !!sim.reads[0].quality });

const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
console.log(`recovered: ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}, failInner=${dec.stats.oligosFailedInnerRS}, failOuter=${dec.stats.oligosFailedOuterRS}, hash=${dec.hashMatches?"OK":"FAIL"}`);

// Check per-oligo details — how many used which strategy?
const strategyCounts: Record<string, number> = {};
let totalCorrected = 0;
for (const p of dec.perOligo) {
  if (p.crcPassed) totalCorrected += p.innerRS?.corrected ?? 0;
}
console.log(`totalCorrected bits: ${totalCorrected}`);
console.log(`per-oligo sample (first 10):`);
for (let i = 0; i < Math.min(10, dec.perOligo.length); i++) {
  const p = dec.perOligo[i];
  console.log(`  oligo[${i}]: reads=${p.readCount} crcPassed=${p.crcPassed} corrected=${p.innerRS?.corrected} success=${p.innerRS?.success}`);
}
