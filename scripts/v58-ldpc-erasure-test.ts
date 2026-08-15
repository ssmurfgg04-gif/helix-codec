// Test LDPC erasure decoder directly
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";

const innerN = 156;
const innerK = 152;
const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });

// Encode random data
const data = new Uint8Array(innerK);
for (let i = 0; i < innerK; i++) data[i] = (i * 37 + 13) & 0xff;
const codeword = ldpc.encode(data);
console.log(`encoded: ${codeword.length} bytes (expected ${innerN})`);

// Test 1: no erasures — should decode cleanly
const r1 = ldpc.decode(codeword);
console.log(`no erasures: match=${r1.data.every((b,i)=>b===data[i])}, corrected=${r1.corrected}`);

// Test 2: erase 1 byte (8 bits) at position 12
const cw2 = new Uint8Array(codeword);
cw2[12] = 0; // zero out byte 12
const erasePos: number[] = [];
for (let bit = 0; bit < 8; bit++) erasePos.push(12 * 8 + bit);
try {
  const r2 = ldpc.decodeWithErasures(cw2, erasePos);
  console.log(`1 byte erasure @12: match=${r2.data.every((b,i)=>b===data[i])}, erased=${r2.erased}`);
} catch (e: any) {
  console.log(`1 byte erasure @12: FAILED — ${e.message}`);
}

// Test 3: erase 2 bytes at positions 12 and 91 (typical for 2-block arithmetic)
const cw3 = new Uint8Array(codeword);
cw3[12] = 0; cw3[91] = 0;
const erasePos3: number[] = [];
for (let bit = 0; bit < 8; bit++) {
  erasePos3.push(12 * 8 + bit);
  erasePos3.push(91 * 8 + bit);
}
try {
  const r3 = ldpc.decodeWithErasures(cw3, erasePos3);
  console.log(`2 byte erasures @12,91: match=${r3.data.every((b,i)=>b===data[i])}, erased=${r3.erased}`);
} catch (e: any) {
  console.log(`2 byte erasures @12,91: FAILED — ${e.message}`);
}

// Test 4: erase 3 bytes (worst case for 3-block arithmetic)
const cw4 = new Uint8Array(codeword);
cw4[12] = 0; cw4[60] = 0; cw4[120] = 0;
const erasePos4: number[] = [];
for (let bit = 0; bit < 8; bit++) {
  erasePos4.push(12 * 8 + bit);
  erasePos4.push(60 * 8 + bit);
  erasePos4.push(120 * 8 + bit);
}
try {
  const r4 = ldpc.decodeWithErasures(cw4, erasePos4);
  console.log(`3 byte erasures @12,60,120: match=${r4.data.every((b,i)=>b===data[i])}, erased=${r4.erased}`);
} catch (e: any) {
  console.log(`3 byte erasures @12,60,120: FAILED — ${e.message}`);
}
