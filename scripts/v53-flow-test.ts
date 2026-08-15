#!/usr/bin/env bun
/**
 * v53 Codec-Flow Arithmetic Round-Trip Test
 *
 * Mimics the exact codec.ts encode flow + Rust decode flow to isolate
 * the encoder/decoder mismatch.
 */

import { bytesToArithmeticDnaCrc } from "../src/lib/dna/markov-arithmetic";

const TAG = "[v53-flow]";

const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");

function crc8ts(data: Uint8Array): number {
  let crc = 0;
  for (const b of data) {
    crc ^= b;
    for (let i = 0; i < 8; i++) {
      if (crc & 0x80) crc = ((crc << 1) ^ 0x07) & 0xff;
      else crc = (crc << 1) & 0xff;
    }
  }
  return crc & 0xff;
}

function dnaToBytes(dna: string): Uint8Array {
  const out = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) out[i] = dna.charCodeAt(i);
  return out;
}

function test(name: string, innerN: number, totalInnerBytes: number, blockSize: number) {
  console.log(`\n--- ${TAG} ${name} ---`);
  console.log(`${TAG} innerN=${innerN}, totalInnerBytes=${totalInnerBytes}, blockSize=${blockSize}`);

  const innerDnaLen = totalInnerBytes * 4;
  const bytesPerBlockTotal = Math.floor(blockSize / 4);
  const bytesPerBlockData = bytesPerBlockTotal - 3;
  const numBlocks = Math.floor(innerDnaLen / blockSize);
  const dataCapacity = numBlocks * bytesPerBlockData;
  console.log(`${TAG} bytesPerBlockTotal=${bytesPerBlockTotal}, bytesPerBlockData=${bytesPerBlockData}, numBlocks=${numBlocks}, dataCapacity=${dataCapacity}`);

  if (dataCapacity < innerN) {
    console.log(`${TAG} SKIP: dataCapacity ${dataCapacity} < innerN ${innerN}`);
    return false;
  }

  // Create a fake LDPC codeword (random bytes)
  const rsCodeword = new Uint8Array(innerN);
  for (let i = 0; i < innerN; i++) rsCodeword[i] = (i * 37 + 11) & 0xff;

  // Mimic codec.ts encode flow
  const innerBlock = new Uint8Array(totalInnerBytes);
  innerBlock.set(rsCodeword, 0);
  // CRC at position innerN..innerN+1 (not used in arithmetic mode, but set for completeness)
  innerBlock[innerN] = 0xab;
  innerBlock[innerN + 1] = 0xcd;

  const truncatedBlock = innerBlock.slice(0, innerN);  // = rsCodeword
  const padded = new Uint8Array(dataCapacity);  // = dataCapacity bytes
  padded.set(truncatedBlock, 0);  // copy rsCodeword, rest is 0

  console.log(`${TAG} padded.length=${padded.length}, first 16: ${Array.from(padded).slice(0, 16).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);

  // Encode
  const dna = bytesToArithmeticDnaCrc(padded, 3, numBlocks * blockSize, blockSize);
  console.log(`${TAG} encoded: ${dna.length} nt`);

  // Decode via WASM test_arithmetic_decode_crc
  const dnaBytes = dnaToBytes(dna);
  const decoded = wasm.test_arithmetic_decode_crc(dnaBytes, 3, innerN, blockSize);

  // Split: first innerN bytes = data, next innerN bytes = erasures
  const dataBytes = decoded.slice(0, innerN);
  const erasureFlags = Array.from(decoded.slice(innerN));
  const erasedCount = erasureFlags.filter((f: number) => f === 1).length;
  console.log(`${TAG} decoded: ${decoded.length} bytes, erasures: ${erasedCount}/${innerN}`);

  // Compare with original rsCodeword (NOT padded — padded has zeros after innerN)
  const match = dataBytes.length === rsCodeword.length && Array.from(dataBytes).every((b: number, i: number) => b === rsCodeword[i]);
  console.log(`${TAG} match: ${match}`);
  if (!match) {
    // Find first mismatch
    let firstMismatch = -1;
    for (let i = 0; i < innerN; i++) {
      if (dataBytes[i] !== rsCodeword[i]) { firstMismatch = i; break; }
    }
    console.log(`${TAG} first mismatch at byte: ${firstMismatch}`);
    if (firstMismatch >= 0) {
      console.log(`${TAG} expected[${firstMismatch}]: 0x${rsCodeword[firstMismatch].toString(16).padStart(2, '0')}`);
      console.log(`${TAG} decoded[${firstMismatch}]:  0x${dataBytes[firstMismatch].toString(16).padStart(2, '0')}`);
    }
    console.log(`${TAG} expected[0..16]: ${Array.from(rsCodeword).slice(0, 16).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    console.log(`${TAG} decoded[0..16]:  ${Array.from(dataBytes).slice(0, 16).map((b: number) => b.toString(16).padStart(2, '0')).join(' ')}`);
  }
  return match;
}

function main() {
  console.log(`${TAG} === v53 Codec-Flow Arithmetic Round-Trip Test ===`);

  // Mimic ULTIMATE_DENSITY_CONFIG: 500nt, 20nt primer, 8 parity, 5% outer
  // After v53 -3 fix: payloadBytes=96, innerN=108, totalInnerBytes=115
  const r1 = test("Test 1: 500nt + arithmetic (ULTIMATE_DENSITY)", 108, 115, 230);

  // Mimic ULTIMATE_LOW_COVERAGE: 500nt, 10 parity, 10% outer
  // After v53 -3 fix: payloadBytes=94, innerN=108, totalInnerBytes=115
  const r2 = test("Test 2: 500nt + arithmetic + low-cov (ULTIMATE_LOW_COV)", 108, 115, 230);

  // Test with smaller block size
  const r3 = test("Test 3: 500nt + arithmetic, blockSize=153", 108, 115, 153);

  // Test with very small block size
  const r4 = test("Test 4: 500nt + arithmetic, blockSize=92", 108, 115, 92);

  console.log(`\n=== ${TAG} Summary: ${[r1, r2, r3, r4].filter(Boolean).length}/4 PASSED ===`);
}

main();
