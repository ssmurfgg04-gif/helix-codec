/**
 * v60: Debug test — single insertion, all-zeros input.
 */

import { IndelViterbiDecoder } from "../src/lib/dna/convolutional-indel";
import { ConvolutionalCode } from "../src/lib/dna/convolutional";

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

const conv = new ConvolutionalCode();
const decoder = new IndelViterbiDecoder({});

// Encode all-zeros
const inputBytes = 20;
const input = new Uint8Array(inputBytes); // all zeros
const inputBits = bytesToBits(input);
const encodedBits = conv.encode(inputBits);
console.log(`Encoded: ${encodedBits.length} bits`);
console.log(`First 20 encoded bits: ${Array.from(encodedBits.slice(0, 20)).join("")}`);

// Test 1: No errors
{
  const decoded = decoder.decode(encodedBits, inputBytes * 8);
  const decodedBytes = bitsToBytes(decoded);
  const match = Buffer.from(input).equals(Buffer.from(decodedBytes));
  console.log(`Test 1 (no errors): ${match ? "PASS" : "FAIL"}`);
  if (!match) {
    console.log(`  Decoded bits: ${Array.from(decoded.slice(0, 20)).join("")}`);
  }
}

// Test 2: 1 insertion at position 50
{
  const noisy = new Uint8Array(encodedBits.length + 1);
  for (let i = 0; i < 50; i++) noisy[i] = encodedBits[i];
  noisy[50] = 1; // inserted bit
  for (let i = 50; i < encodedBits.length; i++) noisy[i + 1] = encodedBits[i];
  const decoded = decoder.decode(noisy, inputBytes * 8);
  const decodedBytes = bitsToBytes(decoded);
  const match = Buffer.from(input).equals(Buffer.from(decodedBytes));
  console.log(`Test 2 (1 ins @ 50): ${match ? "PASS" : "FAIL"}`);
  if (!match) {
    console.log(`  Decoded bits: ${Array.from(decoded.slice(0, 20)).join("")}`);
    // Count errors
    let errors = 0;
    for (let i = 0; i < decoded.length; i++) {
      if (decoded[i] !== inputBits[i]) errors++;
    }
    console.log(`  Bit errors: ${errors}/${decoded.length}`);
  }
}

// Test 3: 1 deletion at position 50
{
  const noisy = new Uint8Array(encodedBits.length - 1);
  for (let i = 0; i < 50; i++) noisy[i] = encodedBits[i];
  for (let i = 51; i < encodedBits.length; i++) noisy[i - 1] = encodedBits[i];
  const decoded = decoder.decode(noisy, inputBytes * 8);
  const decodedBytes = bitsToBytes(decoded);
  const match = Buffer.from(input).equals(Buffer.from(decodedBytes));
  console.log(`Test 3 (1 del @ 50): ${match ? "PASS" : "FAIL"}`);
  if (!match) {
    console.log(`  Decoded bits: ${Array.from(decoded.slice(0, 20)).join("")}`);
    let errors = 0;
    for (let i = 0; i < decoded.length; i++) {
      if (decoded[i] !== inputBits[i]) errors++;
    }
    console.log(`  Bit errors: ${errors}/${decoded.length}`);
  }
}

// Test 4: 1 substitution at position 50
{
  const noisy = new Uint8Array(encodedBits.length);
  noisy.set(encodedBits);
  noisy[50] ^= 1;
  const decoded = decoder.decode(noisy, inputBytes * 8);
  const decodedBytes = bitsToBytes(decoded);
  const match = Buffer.from(input).equals(Buffer.from(decodedBytes));
  console.log(`Test 4 (1 sub @ 50): ${match ? "PASS" : "FAIL"}`);
  if (!match) {
    console.log(`  Decoded bits: ${Array.from(decoded.slice(0, 20)).join("")}`);
    let errors = 0;
    for (let i = 0; i < decoded.length; i++) {
      if (decoded[i] !== inputBits[i]) errors++;
    }
    console.log(`  Bit errors: ${errors}/${decoded.length}`);
  }
}

// Test 5: 3 insertions spread out
{
  const insPositions = [50, 100, 200];
  const noisy: number[] = [];
  let j = 0;
  for (let i = 0; i < encodedBits.length; i++) {
    if (j < insPositions.length && i === insPositions[j]) {
      noisy.push(Math.random() < 0.5 ? 0 : 1); // inserted bit
      j++;
    }
    noisy.push(encodedBits[i]);
  }
  const decoded = decoder.decode(new Uint8Array(noisy), inputBytes * 8);
  const decodedBytes = bitsToBytes(decoded);
  const match = Buffer.from(input).equals(Buffer.from(decodedBytes));
  console.log(`Test 5 (3 ins spread): ${match ? "PASS" : "FAIL"}`);
  if (!match) {
    let errors = 0;
    for (let i = 0; i < decoded.length; i++) {
      if (decoded[i] !== inputBits[i]) errors++;
    }
    console.log(`  Bit errors: ${errors}/${decoded.length}`);
  }
}
