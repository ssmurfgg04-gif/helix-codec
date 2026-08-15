/**
 * v60: Granular tests for indel-tolerant Viterbi.
 * Test each error type in isolation to find the bug.
 */

import { IndelTolerantConvolutionalInnerCode } from "../src/lib/dna/convolutional-indel";
import { ConvolutionalInnerCode } from "../src/lib/dna/convolutional";

function randomBytes(n: number): Uint8Array {
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) out[i] = Math.floor(Math.random() * 256);
  return out;
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

function injectSubsOnly(bits: Uint8Array, rate: number): Uint8Array {
  const out = new Uint8Array(bits.length);
  for (let i = 0; i < bits.length; i++) {
    let b = bits[i];
    if (Math.random() < rate) b ^= 1;
    out[i] = b;
  }
  return out;
}

function injectInsOnly(bits: Uint8Array, rate: number): Uint8Array {
  const out: number[] = [];
  for (let i = 0; i < bits.length; i++) {
    if (Math.random() < rate) {
      out.push(Math.random() < 0.5 ? 0 : 1); // inserted bit
    }
    out.push(bits[i]);
  }
  return new Uint8Array(out);
}

function injectDelOnly(bits: Uint8Array, rate: number): Uint8Array {
  const out: number[] = [];
  for (let i = 0; i < bits.length; i++) {
    if (Math.random() < rate) continue; // deleted bit
    out.push(bits[i]);
  }
  return new Uint8Array(out);
}

const inputBytes = 20;
const input = randomBytes(inputBytes);
console.log("Input:", Buffer.from(input).toString("hex"));

const standardConv = new ConvolutionalInnerCode(inputBytes);
const indelConv = new IndelTolerantConvolutionalInnerCode(inputBytes);

const encoded = standardConv.encode(input);
const encodedBits = bytesToBits(encoded);
console.log(`Encoded: ${encoded.length} bytes (${encodedBits.length} bits)`);

// Test pure substitutions at various rates
console.log("\n=== Pure substitutions ===");
for (const rate of [0.01, 0.02, 0.03, 0.05, 0.10]) {
  let standardPass = 0, indelPass = 0;
  const trials = 20;
  for (let t = 0; t < trials; t++) {
    const testInput = randomBytes(inputBytes);
    const enc = standardConv.encode(testInput);
    const encBits = bytesToBits(enc);
    const noisyBits = injectSubsOnly(encBits, rate);
    const received = bitsToBytes(noisyBits);
    try {
      const dec1 = standardConv.decode(received);
      if (Buffer.from(testInput).equals(Buffer.from(dec1))) standardPass++;
    } catch {}
    try {
      const dec2 = indelConv.decode(received);
      if (Buffer.from(testInput).equals(Buffer.from(dec2))) indelPass++;
    } catch {}
  }
  console.log(`  ${rate*100}% sub: standard ${standardPass}/${trials}, indel ${indelPass}/${trials}`);
}

// Test pure insertions
console.log("\n=== Pure insertions ===");
for (const rate of [0.01, 0.02, 0.03, 0.05, 0.10]) {
  let indelPass = 0;
  const trials = 20;
  for (let t = 0; t < trials; t++) {
    const testInput = randomBytes(inputBytes);
    const enc = standardConv.encode(testInput);
    const encBits = bytesToBits(enc);
    const noisyBits = injectInsOnly(encBits, rate);
    const padded = new Uint8Array(Math.ceil(noisyBits.length / 8) * 8);
    padded.set(noisyBits);
    const received = bitsToBytes(padded);
    try {
      const dec = indelConv.decode(received);
      if (Buffer.from(testInput).equals(Buffer.from(dec))) indelPass++;
    } catch {}
  }
  console.log(`  ${rate*100}% ins: indel ${indelPass}/${trials}`);
}

// Test pure deletions
console.log("\n=== Pure deletions ===");
for (const rate of [0.01, 0.02, 0.03, 0.05, 0.10]) {
  let indelPass = 0;
  const trials = 20;
  for (let t = 0; t < trials; t++) {
    const testInput = randomBytes(inputBytes);
    const enc = standardConv.encode(testInput);
    const encBits = bytesToBits(enc);
    const noisyBits = injectDelOnly(encBits, rate);
    const padded = new Uint8Array(Math.ceil(noisyBits.length / 8) * 8);
    padded.set(noisyBits);
    const received = bitsToBytes(padded);
    try {
      const dec = indelConv.decode(received);
      if (Buffer.from(testInput).equals(Buffer.from(dec))) indelPass++;
    } catch {}
  }
  console.log(`  ${rate*100}% del: indel ${indelPass}/${trials}`);
}

// Test combined IDS at low rates
console.log("\n=== Combined IDS (low rates) ===");
for (const [sub, ins, del] of [[0.005, 0.005, 0.005], [0.01, 0.01, 0.01], [0.015, 0.015, 0.015]] as number[][]) {
  let indelPass = 0;
  const trials = 20;
  for (let t = 0; t < trials; t++) {
    const testInput = randomBytes(inputBytes);
    const enc = standardConv.encode(testInput);
    const encBits = bytesToBits(enc);
    let noisy = injectSubsOnly(encBits, sub);
    noisy = injectInsOnly(noisy, ins);
    noisy = injectDelOnly(noisy, del);
    const padded = new Uint8Array(Math.ceil(noisy.length / 8) * 8);
    padded.set(noisy);
    const received = bitsToBytes(padded);
    try {
      const dec = indelConv.decode(received);
      if (Buffer.from(testInput).equals(Buffer.from(dec))) indelPass++;
    } catch {}
  }
  console.log(`  ${(sub+ins+del)*100}% IDS (${sub*100}+${ins*100}+${del*100}): indel ${indelPass}/${trials}`);
}
