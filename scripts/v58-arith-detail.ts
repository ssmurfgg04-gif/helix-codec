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
const layout = computeLayoutAuto(cfg);
const innerN = layout.addressBytes + layout.payloadBytes + layout.innerParityBytes;

const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });

// Try multiple reads to find one that decodes
const ldpc = new LDPCInnerCode({ n: innerN, k: innerN - layout.innerParityBytes });
const fwd = enc.encoded.forwardPrimer;
const rev = enc.encoded.reversePrimer;
const arithBlockSize = Math.floor((layout.totalInnerBytes * 4) / 2);

let successCount = 0;
let failCount = 0;
let noErasureCount = 0;
for (let r = 0; r < sim.reads.length; r++) {
  const read = sim.reads[r];
  if (read.sequence.length < fwd.length + rev.length) continue;
  const inner = read.sequence.slice(fwd.length, read.sequence.length - rev.length);
  
  try {
    const result = arithmeticDnaToBytesCrc(inner, 3, innerN, arithBlockSize);
    const rsCodeword = result.data.slice(0, innerN);
    const erasePos: number[] = [];
    for (let i = 0; i < innerN; i++) {
      if (result.erasures[i]) {
        for (let bit = 0; bit < 8; bit++) erasePos.push(i * 8 + bit);
      }
    }
    
    if (erasePos.length === 0) {
      // No erasures — hard decode
      const decoded = ldpc.decode(rsCodeword);
      const uw0 = decoded.data[0] ^ 0x1b;
      const uw1 = decoded.data[1] ^ 0x4b;
      const uw2 = decoded.data[2] ^ 0x24;
      const oi = (uw0 << 16) | (uw1 << 8) | uw2;
      if (oi < enc.encoded.oligos.length) {
        successCount++;
        noErasureCount++;
      }
      continue;
    }
    
    try {
      const decoded = ldpc.decodeWithErasures(rsCodeword, erasePos);
      const uw0 = decoded.data[0] ^ 0x1b;
      const uw1 = decoded.data[1] ^ 0x4b;
      const uw2 = decoded.data[2] ^ 0x24;
      const oi = (uw0 << 16) | (uw1 << 8) | uw2;
      if (oi < enc.encoded.oligos.length) {
        successCount++;
      } else {
        failCount++;
      }
    } catch (e: any) {
      failCount++;
      if (failCount <= 3) console.log(`read ${r}: erasure decode failed: ${e.message.slice(0,80)}`);
    }
  } catch {
    failCount++;
  }
}
console.log(`success: ${successCount}, fail: ${failCount}, noErasure: ${noErasureCount}`);
