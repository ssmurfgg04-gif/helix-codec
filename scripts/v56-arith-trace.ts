#!/usr/bin/env bun
/**
 * v56 — Trace encoder/decoder state divergence for all-FF data.
 * Compare bit-by-bit to find where the encoder input and decoder output diverge.
 */

import {
  bytesToArithmeticDna,
  arithmeticDnaToBytes,
} from "../src/lib/dna/markov-arithmetic";

const TAG = "[v56-trace]";

// Test: encode 1 byte of 0xFF, decode, compare bit-by-bit
function traceByte(value: number) {
  const data = new Uint8Array([value]);
  const numBases = 8; // generous: 8 nt for 1 byte
  const dna = bytesToArithmeticDna(data, 3, numBases);
  const decoded = arithmeticDnaToBytes(dna, 3, 1);

  const expectedBits: number[] = [];
  for (let b = 7; b >= 0; b--) expectedBits.push((value >> b) & 1);

  // Get the raw decoder output bits by encoding a longer expectedBytes
  const decodedRaw = arithmeticDnaToBytes(dna, 3, 4); // 4 bytes = 32 bits
  const actualBits: number[] = [];
  for (let i = 0; i < 4; i++) {
    for (let b = 7; b >= 0; b--) actualBits.push((decodedRaw[i] >> b) & 1);
  }

  const match = decoded[0] === value;
  console.log(`${TAG} byte=0x${value.toString(16).padStart(2, '0')} dna=${dna} match=${match}`);
  console.log(`  expected bits: ${expectedBits.join("")} (first 8)`);
  console.log(`  actual bits:   ${actualBits.slice(0, 8).join("")} | ${actualBits.slice(8, 16).join("")} | ${actualBits.slice(16, 24).join("")} | ${actualBits.slice(24, 32).join("")}`);

  // Find first bit mismatch
  for (let i = 0; i < 8; i++) {
    if (expectedBits[i] !== actualBits[i]) {
      console.log(`  first mismatch at bit ${i} (byte ${Math.floor(i/8)}, bit ${7 - (i % 8)}): expected ${expectedBits[i]}, got ${actualBits[i]}`);
      break;
    }
  }
  console.log();
}

// Trace several byte values
for (const v of [0x00, 0xFF, 0x55, 0xAA, 0x42, 0x01, 0x80, 0x7F, 0xFE, 0x12]) {
  traceByte(v);
}

// Trace 2-byte values
function traceTwoBytes(v1: number, v2: number) {
  const data = new Uint8Array([v1, v2]);
  const numBases = 12;
  const dna = bytesToArithmeticDna(data, 3, numBases);
  const decoded = arithmeticDnaToBytes(dna, 3, 2);
  const match = decoded[0] === v1 && decoded[1] === v2;
  console.log(`${TAG} bytes=0x${v1.toString(16).padStart(2, '0')} 0x${v2.toString(16).padStart(2, '0')} match=${match}`);
  if (!match) {
    console.log(`  expected: ${v1.toString(16).padStart(2, '0')} ${v2.toString(16).padStart(2, '0')}`);
    console.log(`  decoded:  ${decoded[0].toString(16).padStart(2, '0')} ${decoded[1].toString(16).padStart(2, '0')}`);
  }
}

console.log(`\n--- 2-byte traces ---`);
traceTwoBytes(0xFF, 0xFF);
traceTwoBytes(0x00, 0x00);
traceTwoBytes(0xFF, 0x00);
traceTwoBytes(0x00, 0xFF);
traceTwoBytes(0x42, 0x42);
