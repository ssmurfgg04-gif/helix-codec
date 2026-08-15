// Profile 2.1MB to find scaling bottleneck
import { encodeFile } from "../src/lib/dna/codec";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG, computeLayoutAuto } from "../src/lib/dna/types";
import { fullDecode } from "../src/lib/dna/wasm-batch-decode";
import * as fs from "fs";

const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
console.log("payload:", payload.length, "bytes");

const t0 = Date.now();
const enc = await encodeFile(payload, { ...DEFAULT_CONFIG }, { fileName: "test.bin", contentType: "application/octet-stream" });
console.log(`encode: ${Date.now()-t0}ms`);

const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
console.log(`reads: ${sim.reads.length}`);

const t1 = Date.now();
let totalLen = 0;
for (const r of sim.reads) totalLen += r.sequence.length;
const allReads = new Uint8Array(totalLen);
const readOffsets = new Uint8Array(sim.reads.length * 4);
const readLengths = new Uint8Array(sim.reads.length * 4);
let offset = 0;
for (let i = 0; i < sim.reads.length; i++) {
  const seq = sim.reads[i].sequence;
  new DataView(readOffsets.buffer, i * 4, 4).setUint32(0, offset, true);
  new DataView(readLengths.buffer, i * 4, 4).setUint32(0, seq.length, true);
  for (let j = 0; j < seq.length; j++) allReads[offset + j] = seq.charCodeAt(j);
  offset += seq.length;
}
console.log(`flatten: ${Date.now()-t1}ms, allReads: ${totalLen} bytes`);

const layout = computeLayoutAuto({ ...DEFAULT_CONFIG });
const innerK = layout.addressBytes + layout.payloadBytes;
const innerN = innerK + layout.innerParityBytes;

const t2 = Date.now();
const result = fullDecode(
  sim.reads, enc.encoded.forwardPrimer, enc.encoded.reversePrimer,
  enc.encoded.metadata.oligoCount, innerN, innerK, layout.totalInnerBytes,
  enc.encoded.metadata.outerRS.n, enc.encoded.metadata.outerRS.k, layout.payloadBytes,
  enc.encoded.metadata.fileSize, enc.encoded.metadata.compression === "deflate",
);
console.log(`WASM fullDecode: ${Date.now()-t2}ms = ${(payload.length/1024/1024)/((Date.now()-t2)/1000)} MB/s`);
console.log(`result len: ${result.length}, hash match: ${result[0]===payload[0] && result[result.length-1]===payload[payload.length-1]}`);
