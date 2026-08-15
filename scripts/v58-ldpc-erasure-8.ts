import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";

// Test with innerParityBytes=8 (m=64)
const innerN = 156;
const innerK = 148; // 156 - 8
const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });

const data = new Uint8Array(innerK);
for (let i = 0; i < innerK; i++) data[i] = (i * 37 + 13) & 0xff;
const codeword = ldpc.encode(data);

// Test 2 byte erasures
const cw3 = new Uint8Array(codeword);
cw3[12] = 0; cw3[91] = 0;
const erasePos3: number[] = [];
for (let bit = 0; bit < 8; bit++) {
  erasePos3.push(12 * 8 + bit);
  erasePos3.push(91 * 8 + bit);
}
try {
  const r3 = ldpc.decodeWithErasures(cw3, erasePos3);
  console.log(`m=64, 2 byte erasures @12,91: match=${r3.data.every((b,i)=>b===data[i])}, erased=${r3.erased}`);
} catch (e: any) {
  console.log(`m=64, 2 byte erasures @12,91: FAILED — ${e.message}`);
}

// Test 3 byte erasures
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
  console.log(`m=64, 3 byte erasures @12,60,120: match=${r4.data.every((b,i)=>b===data[i])}, erased=${r4.erased}`);
} catch (e: any) {
  console.log(`m=64, 3 byte erasures @12,60,120: FAILED — ${e.message}`);
}

// Test 4 byte erasures (worst case for 4-block arithmetic)
const cw5 = new Uint8Array(codeword);
cw5[12] = 0; cw5[50] = 0; cw5[100] = 0; cw5[140] = 0;
const erasePos5: number[] = [];
for (let bit = 0; bit < 8; bit++) {
  erasePos5.push(12 * 8 + bit);
  erasePos5.push(50 * 8 + bit);
  erasePos5.push(100 * 8 + bit);
  erasePos5.push(140 * 8 + bit);
}
try {
  const r5 = ldpc.decodeWithErasures(cw5, erasePos5);
  console.log(`m=64, 4 byte erasures: match=${r5.data.every((b,i)=>b===data[i])}, erased=${r5.erased}`);
} catch (e: any) {
  console.log(`m=64, 4 byte erasures: FAILED — ${e.message}`);
}
