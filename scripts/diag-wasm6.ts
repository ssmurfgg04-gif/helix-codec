// Use full_decode_pre_deflate to see what comes out of LDPC+RS without inflate
import { encodeFile } from "../src/lib/dna/codec";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG, computeLayoutAuto } from "../src/lib/dna/types";
import * as fs from "fs";
const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");

const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 1024);
const cfg = { ...DEFAULT_CONFIG, compress: false };  // disable deflate for clearer diagnostic
const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
const layout = computeLayoutAuto(cfg);
const innerK = layout.addressBytes + layout.payloadBytes;
const innerN = innerK + layout.innerParityBytes;

const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });

// Flatten reads
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
const fwd = new Uint8Array(enc.encoded.forwardPrimer.length);
for (let i=0; i<fwd.length; i++) fwd[i] = enc.encoded.forwardPrimer.charCodeAt(i);
const rev = new Uint8Array(enc.encoded.reversePrimer.length);
for (let i=0; i<rev.length; i++) rev[i] = enc.encoded.reversePrimer.charCodeAt(i);

// Pre-deflate = post-RS concatenated payloads
const preDeflate = wasm.full_decode_pre_deflate(
  allReads, readOffsets, readLengths,
  fwd, rev,
  enc.encoded.metadata.oligoCount,
  innerN, innerK, layout.totalInnerBytes,
  enc.encoded.metadata.outerRS.n, enc.encoded.metadata.outerRS.k,
  layout.payloadBytes,
);
console.log("pre-deflate bytes:", preDeflate.length, "expected:", enc.encoded.metadata.outerRS.k * layout.payloadBytes);
console.log("first 32 bytes:", Array.from(preDeflate.slice(0,32)).map(b=>b.toString(16).padStart(2,'0')).join(' '));

// Compare to expected (raw payload bytes, since compress=false)
console.log("expected first 32:", Array.from(payload.slice(0,32)).map(b=>b.toString(16).padStart(2,'0')).join(' '));
console.log("all zero?:", Array.from(preDeflate).every(b=>b===0));

// Count zero blocks (per oligo = 54 bytes payload)
const blockCount = Math.floor(preDeflate.length / layout.payloadBytes);
let zeroBlocks = 0;
for (let i=0; i<blockCount; i++) {
  const off = i * layout.payloadBytes;
  let isZero = true;
  for (let j=0; j<layout.payloadBytes; j++) {
    if (preDeflate[off+j] !== 0) { isZero = false; break; }
  }
  if (isZero) zeroBlocks++;
}
console.log(`zero blocks: ${zeroBlocks}/${blockCount}`);
