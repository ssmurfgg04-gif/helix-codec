// Trace the full arithmetic decode pipeline step by step
import { encodeFile } from "../src/lib/dna/codec";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { computeLayoutAuto } from "../src/lib/dna/types";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { arithmeticDnaToBytesCrc } from "../src/lib/dna/markov-arithmetic";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import * as crypto from "crypto";

const payload = crypto.randomBytes(1024);
const cfg = ULTIMATE_V55_DENSITY_CONFIG;
const enc = await encodeFile(Buffer.from(payload), cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
console.log("oligos:", enc.encoded.oligos.length, "density:", enc.stats.netDensityBitsPerNt.toFixed(3));

const layout = computeLayoutAuto(cfg);
const innerN = layout.addressBytes + layout.payloadBytes + layout.innerParityBytes;
console.log("innerN:", innerN, "payloadBytes:", layout.payloadBytes);

const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
console.log("reads:", sim.reads.length);

// Take the first read, trim primers, arithmetic decode, LDPC decode
const read = sim.reads[0];
const fwd = enc.encoded.forwardPrimer;
const rev = enc.encoded.reversePrimer;
console.log("read length:", read.sequence.length, "fwd primer:", fwd, "rev primer:", rev);

// Trim primers
const inner = read.sequence.slice(fwd.length, read.sequence.length - rev.length);
console.log("inner DNA length:", inner.length, "expected:", layout.totalInnerBytes * 4);

// Arithmetic decode
const arithBlockSize = Math.floor((layout.totalInnerBytes * 4) / 2);
const result = arithmeticDnaToBytesCrc(inner, 3, innerN, arithBlockSize);
console.log("arith decode: len=", result.data.length, "erasures:", result.erasures.filter(e=>e).length);

// Check erasure positions
for (let i = 0; i < result.erasures.length; i++) {
  if (result.erasures[i]) console.log(`  erased byte ${i}`);
}

// LDPC erasure decode
const ldpc = new LDPCInnerCode({ n: innerN, k: innerN - layout.innerParityBytes });
const rsCodeword = result.data.slice(0, innerN);
const erasePos: number[] = [];
for (let i = 0; i < innerN * 8; i++) {
  if (result.erasures[Math.floor(i/8)]) erasePos.push(i);
}
console.log("erase bit positions:", erasePos.length);

try {
  const r = ldpc.decodeWithErasures(rsCodeword, erasePos);
  console.log("LDPC erasure decode: SUCCESS, erased=", r.erased);
  // Check address
  const uw0 = r.data[0] ^ 0x1b;
  const uw1 = r.data[1] ^ 0x4b;
  const uw2 = r.data[2] ^ 0x24;
  const oi = (uw0 << 16) | (uw1 << 8) | uw2;
  console.log("decoded oligo index:", oi, "(expected 0-", enc.encoded.oligos.length-1, ")");
} catch (e: any) {
  console.log("LDPC erasure decode FAILED:", e.message);
  // Try without erasures
  try {
    const r = ldpc.decode(rsCodeword);
    console.log("LDPC hard decode: SUCCESS, corrected=", r.corrected);
  } catch (e2: any) {
    console.log("LDPC hard decode FAILED:", e2.message);
  }
}
