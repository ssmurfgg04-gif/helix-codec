// Profile each phase of encode + decode for 2.1MB to find bottlenecks
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG, computeLayoutAuto } from "../src/lib/dna/types";
import { fullDecode } from "../src/lib/dna/wasm-batch-decode";
import * as fs from "fs";

const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 262144); // 256KB

const t0 = Date.now();
const enc = await encodeFile(payload, { ...DEFAULT_CONFIG }, { fileName: "test.bin", contentType: "application/octet-stream" });
const encMs = Date.now() - t0;
console.log(`encode: ${encMs}ms = ${(payload.length/1024/1024)/(encMs/1000)} MB/s`);

// Profile read flattening (JS overhead)
const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
console.log(`simulate: ${sim.reads.length} reads`);

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
const flattenMs = Date.now() - t1;
console.log(`flatten: ${flattenMs}ms (allReads: ${totalLen} bytes)`);

// WASM call only
const layout = computeLayoutAuto({ ...DEFAULT_CONFIG });
const innerK = layout.addressBytes + layout.payloadBytes;
const innerN = innerK + layout.innerParityBytes;
const fwd = new Uint8Array(enc.encoded.forwardPrimer.length);
for (let i=0; i<fwd.length; i++) fwd[i] = enc.encoded.forwardPrimer.charCodeAt(i);
const rev = new Uint8Array(enc.encoded.reversePrimer.length);
for (let i=0; i<rev.length; i++) rev[i] = enc.encoded.reversePrimer.charCodeAt(i);

const t2 = Date.now();
const result = fullDecode(
  sim.reads, enc.encoded.forwardPrimer, enc.encoded.reversePrimer,
  enc.encoded.metadata.oligoCount, innerN, innerK, layout.totalInnerBytes,
  enc.encoded.metadata.outerRS.n, enc.encoded.metadata.outerRS.k, layout.payloadBytes,
  enc.encoded.metadata.fileSize, enc.encoded.metadata.compression === "deflate",
);
const wasmMs = Date.now() - t2;
console.log(`WASM fullDecode: ${wasmMs}ms = ${(payload.length/1024/1024)/(wasmMs/1000)} MB/s`);
console.log(`result len: ${result.length}, hash: ${Array.from(result.slice(0,16)).map(b=>b.toString(16).padStart(2,'0')).join(' ')}`);
console.log(`expected:     ${Array.from(payload.slice(0,16)).map(b=>b.toString(16).padStart(2,'0')).join(' ')}`);

const t3 = Date.now();
const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, { ...DEFAULT_CONFIG }, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
const decMs = Date.now() - t3;
console.log(`\ndecodeReadsUltra (end-to-end): ${decMs}ms = ${(payload.length/1024/1024)/(decMs/1000)} MB/s`);
console.log(`hash: ${dec.hashMatches?"OK":"FAIL"}`);
