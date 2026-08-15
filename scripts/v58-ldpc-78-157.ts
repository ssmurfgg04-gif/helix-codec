import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";

const innerN = 158;
const innerK = 150;
const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });

const data = new Uint8Array(innerK);
for (let i = 0; i < innerK; i++) data[i] = (i * 37 + 13) & 0xff;
const codeword = ldpc.encode(data);

// Test erasure at positions 78 and 157 (the actual arithmetic termination positions)
const cw = new Uint8Array(codeword);
cw[78] = 0; cw[157] = 0;
const erasePos: number[] = [];
for (let bit = 0; bit < 8; bit++) {
  erasePos.push(78 * 8 + bit);
  erasePos.push(157 * 8 + bit);
}
try {
  const r = ldpc.decodeWithErasures(cw, erasePos);
  let match = true;
  for (let i = 0; i < innerK; i++) {
    if (r.data[i] !== data[i]) { match = false; break; }
  }
  console.log(`erasure @78,157: match=${match}, erased=${r.erased}`);
} catch (e: any) {
  console.log(`erasure @78,157: FAILED — ${e.message}`);
}

// Test with positions 78 and 157 set to WRONG values (not zero)
const cw2 = new Uint8Array(codeword);
cw2[78] = codeword[78] - 1; // off by 1 (matches arithmetic termination corruption)
cw2[157] = codeword[157] - 1;
try {
  const r = ldpc.decodeWithErasures(cw2, erasePos);
  let match = true;
  for (let i = 0; i < innerK; i++) {
    if (r.data[i] !== data[i]) { match = false; break; }
  }
  console.log(`wrong-val @78,157: match=${match}, erased=${r.erased}`);
} catch (e: any) {
  console.log(`wrong-val @78,157: FAILED — ${e.message}`);
}

// Check syndrome of the wrong-val codeword
let syndromeNonZero = 0;
for (let i = 0; i < ldpc.mBits; i++) {
  let s = 0;
  const cols = ldpc.rowCols[i];
  for (let idx = 0; idx < cols.length; idx++) {
    const byteIdx = cols[idx] >> 3;
    const bitIdx = 7 - (cols[idx] & 7);
    s ^= (cw2[byteIdx] >> bitIdx) & 1;
  }
  if (s) syndromeNonZero++;
}
console.log(`wrong-val syndrome: ${syndromeNonZero} non-zero bits out of ${ldpc.mBits}`);
