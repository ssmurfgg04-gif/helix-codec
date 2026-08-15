/**
 * v60: Quick test for indel-tolerant Viterbi decoder.
 *
 * Encode a random byte sequence, inject various error patterns
 * (subs only, ins only, del only, mixed IDS), and verify the decoder
 * recovers the original.
 */

import { IndelTolerantConvolutionalInnerCode } from "../src/lib/dna/convolutional-indel";
import { ConvolutionalInnerCode } from "../src/lib/dna/convolutional";

function randomBytes(n: number): Uint8Array {
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) out[i] = Math.floor(Math.random() * 256);
  return out;
}

function injectErrors(bits: Uint8Array, subRate: number, insRate: number, delRate: number): Uint8Array {
  const out: number[] = [];
  for (let i = 0; i < bits.length; i++) {
    const r = Math.random();
    if (r < delRate) {
      // Deletion: skip this bit
      continue;
    }
    if (r < delRate + insRate) {
      // Insertion: insert a random bit BEFORE this bit
      out.push(Math.random() < 0.5 ? 0 : 1);
    }
    // Output the bit (possibly with substitution)
    let b = bits[i];
    if (Math.random() < subRate) b ^= 1;
    out.push(b);
  }
  return new Uint8Array(out);
}

function bytesToBits(data: Uint8Array): Uint8Array {
  const bits = new Uint8Array(data.length * 8);
  for (let i = 0; i < data.length; i++) {
    for (let b = 0; b < 8; b++) {
      bits[i * 8 + b] = (data[i] >> (7 - b)) & 1;
    }
  }
  return bits;
}

function bitsToBytes(bits: Uint8Array): Uint8Array {
  const out = new Uint8Array(Math.floor(bits.length / 8));
  for (let i = 0; i < out.length; i++) {
    let byte = 0;
    for (let b = 0; b < 8; b++) {
      byte |= bits[i * 8 + b] << (7 - b);
    }
    out[i] = byte;
  }
  return out;
}

const inputBytes = 20;
const input = randomBytes(inputBytes);
console.log("Input:", Buffer.from(input).toString("hex"));

// Standard conv code (substitution-only Viterbi)
const standardConv = new ConvolutionalInnerCode(inputBytes);
// Indel-tolerant conv code (augmented Viterbi with drift)
const indelConv = new IndelTolerantConvolutionalInnerCode(inputBytes);

const encoded = standardConv.encode(input);
const encodedBits = bytesToBits(encoded);
console.log(`Encoded: ${encoded.length} bytes (${encodedBits.length} bits)`);

// Test 1: No errors
{
  const received = bitsToBytes(encodedBits);
  const decoded = standardConv.decode(received);
  const match = Buffer.from(input).equals(Buffer.from(decoded));
  console.log(`Test 1 (no errors, standard Viterbi): ${match ? "PASS" : "FAIL"}`);
}

// Test 2: 1% substitution only (standard Viterbi should handle)
{
  const noisyBits = injectErrors(encodedBits, 0.01, 0, 0);
  const received = bitsToBytes(noisyBits);
  const decodedStandard = standardConv.decode(received);
  const decodedIndel = indelConv.decode(received);
  const matchStandard = Buffer.from(input).equals(Buffer.from(decodedStandard));
  const matchIndel = Buffer.from(input).equals(Buffer.from(decodedIndel));
  console.log(`Test 2 (1% sub, standard Viterbi): ${matchStandard ? "PASS" : "FAIL"}`);
  console.log(`Test 2 (1% sub, indel Viterbi):    ${matchIndel ? "PASS" : "FAIL"}`);
}

// Test 3: 3% substitution + 3% insertion + 3% deletion
{
  const noisyBits = injectErrors(encodedBits, 0.03, 0.03, 0.03);
  const received = bitsToBytes(noisyBits);
  const decodedStandard = standardConv.decode(received);
  const decodedIndel = indelConv.decode(received);
  const matchStandard = Buffer.from(input).equals(Buffer.from(decodedStandard));
  const matchIndel = Buffer.from(input).equals(Buffer.from(decodedIndel));
  console.log(`Test 3 (3%+3%+3% IDS, standard Viterbi): ${matchStandard ? "PASS" : "FAIL"}`);
  console.log(`Test 3 (3%+3%+3% IDS, indel Viterbi):    ${matchIndel ? "PASS" : "FAIL"}`);
}

// Test 4: 5% substitution + 5% insertion + 5% deletion (9% total IDS)
{
  const noisyBits = injectErrors(encodedBits, 0.05, 0.05, 0.05);
  const received = bitsToBytes(noisyBits);
  const decodedStandard = standardConv.decode(received);
  const decodedIndel = indelConv.decode(received);
  const matchStandard = Buffer.from(input).equals(Buffer.from(decodedStandard));
  const matchIndel = Buffer.from(input).equals(Buffer.from(decodedIndel));
  console.log(`Test 4 (9% IDS, standard Viterbi): ${matchStandard ? "PASS" : "FAIL"}`);
  console.log(`Test 4 (9% IDS, indel Viterbi):    ${matchIndel ? "PASS" : "FAIL"}`);
}

// Test 5: Heavy insertions (5% sub + 10% ins + 5% del)
{
  const noisyBits = injectErrors(encodedBits, 0.05, 0.10, 0.05);
  const received = bitsToBytes(noisyBits);
  const decodedStandard = standardConv.decode(received);
  const decodedIndel = indelConv.decode(received);
  const matchStandard = Buffer.from(input).equals(Buffer.from(decodedStandard));
  const matchIndel = Buffer.from(input).equals(Buffer.from(decodedIndel));
  console.log(`Test 5 (5%+10%+5%, standard Viterbi): ${matchStandard ? "PASS" : "FAIL"}`);
  console.log(`Test 5 (5%+10%+5%, indel Viterbi):    ${matchIndel ? "PASS" : "FAIL"}`);
}

// Test 6: Pure indels (no subs) — 5% ins + 5% del
{
  const noisyBits = injectErrors(encodedBits, 0, 0.05, 0.05);
  const received = bitsToBytes(noisyBits);
  const decodedStandard = standardConv.decode(received);
  const decodedIndel = indelConv.decode(received);
  const matchStandard = Buffer.from(input).equals(Buffer.from(decodedStandard));
  const matchIndel = Buffer.from(input).equals(Buffer.from(decodedIndel));
  console.log(`Test 6 (5%ins+5%del, standard Viterbi): ${matchStandard ? "PASS" : "FAIL"}`);
  console.log(`Test 6 (5%ins+5%del, indel Viterbi):    ${matchIndel ? "PASS" : "FAIL"}`);
}

// Run tests 3-6 multiple times to get statistics
console.log("\n=== Statistical tests (20 trials each) ===");
for (const [name, sub, ins, del] of [
  ["3%+3%+3% IDS", 0.03, 0.03, 0.03],
  ["5%+5%+5% IDS", 0.05, 0.05, 0.05],
  ["3%+6%+3% IDS", 0.03, 0.06, 0.03],
] as [string, number, number, number][]) {
  let standardPass = 0;
  let indelPass = 0;
  const trials = 20;
  for (let trial = 0; trial < trials; trial++) {
    const testInput = randomBytes(inputBytes);
    const enc = standardConv.encode(testInput);
    const encBits = bytesToBits(enc);
    const noisyBits = injectErrors(encBits, sub, ins, del);
    // Pad bits to byte boundary
    const paddedBits = new Uint8Array(Math.ceil(noisyBits.length / 8) * 8);
    paddedBits.set(noisyBits);
    const received = bitsToBytes(paddedBits);
    try {
      const decStandard = standardConv.decode(received);
      if (Buffer.from(testInput).equals(Buffer.from(decStandard))) standardPass++;
    } catch {}
    try {
      const decIndel = indelConv.decode(received);
      if (Buffer.from(testInput).equals(Buffer.from(decIndel))) indelPass++;
    } catch {}
  }
  console.log(`${name}: standard ${standardPass}/${trials} (${(standardPass/trials*100).toFixed(0)}%), indel ${indelPass}/${trials} (${(indelPass/trials*100).toFixed(0)}%)`);
}
