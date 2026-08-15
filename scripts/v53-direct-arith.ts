#!/usr/bin/env bun
/**
 * Direct test of bytesToArithmeticDnaCrc + Rust decoder round-trip.
 * Bypasses the full codec to isolate encoder/decoder mismatch.
 */

import { bytesToArithmeticDnaCrc, arithmeticDnaToBytes } from "../src/lib/dna/markov-arithmetic";

const TAG = "[v53-direct]";

// Test 1: TS encoder + TS decoder round-trip (should always work)
function testTsRoundTrip() {
  console.log(`\n--- ${TAG} Test 1: TS encode + TS decode round-trip ---`);
  const data = new Uint8Array(110);
  for (let i = 0; i < 110; i++) data[i] = (i * 37 + 11) & 0xff;

  const blockSize = 230;
  const dna = bytesToArithmeticDnaCrc(data, 3, 460, blockSize);
  console.log(`${TAG} encoded: ${dna.length} nt`);
  console.log(`${TAG} first 60 nt: ${dna.slice(0, 60)}`);
  console.log(`${TAG} last 60 nt: ${dna.slice(-60)}`);

  const decoded = arithmeticDnaToBytes(dna, 3, 110, blockSize);
  console.log(`${TAG} decoded: ${decoded.length} bytes`);

  const match = decoded.length === data.length && Array.from(decoded).every((b, i) => b === data[i]);
  console.log(`${TAG} match: ${match}`);
  if (!match) {
    console.log(`${TAG} data    : ${Array.from(data).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    console.log(`${TAG} decoded : ${Array.from(decoded).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  }
  return match;
}

// Test 2: Different block sizes
function testBlockSizes() {
  console.log(`\n--- ${TAG} Test 2: Block size scan ---`);
  const data = new Uint8Array(110);
  for (let i = 0; i < 110; i++) data[i] = (i * 37 + 11) & 0xff;

  for (const blockSize of [230, 153, 115, 92, 57, 46]) {
    const bytesPerBlockTotal = Math.floor(blockSize / 4);
    const bytesPerBlockData = bytesPerBlockTotal - 2;
    if (bytesPerBlockData <= 0) continue;
    const numBlocks = Math.ceil(data.length / bytesPerBlockData);
    const totalDnaLen = numBlocks * blockSize;
    if (totalDnaLen > 460) {
      console.log(`${TAG} blockSize=${blockSize}: totalDnaLen=${totalDnaLen} > 460, skip`);
      continue;
    }
    const dna = bytesToArithmeticDnaCrc(data, 3, totalDnaLen, blockSize);
    const decoded = arithmeticDnaToBytes(dna, 3, 110, blockSize);
    const match = decoded.length === data.length && Array.from(decoded).every((b, i) => b === data[i]);
    console.log(`${TAG} blockSize=${blockSize}, bpbData=${bytesPerBlockData}, numBlocks=${numBlocks}, dnaLen=${dna.length}, match=${match}`);
  }
}

// Test 3: Encode/decode with targetLen > totalDnaLen (forces padding)
function testPadding() {
  console.log(`\n--- ${TAG} Test 3: Padding behavior ---`);
  const data = new Uint8Array(110);
  for (let i = 0; i < 110; i++) data[i] = (i * 37 + 11) & 0xff;

  // Force targetLen = 460, but only 1 block needed
  const blockSize = 460;  // 1 huge block
  const bytesPerBlockTotal = Math.floor(blockSize / 4);  // 115
  const bytesPerBlockData = bytesPerBlockTotal - 2;  // 113
  console.log(`${TAG} blockSize=${blockSize}, bpbData=${bytesPerBlockData}`);

  const dna = bytesToArithmeticDnaCrc(data, 3, 460, blockSize);
  console.log(`${TAG} encoded: ${dna.length} nt`);

  const decoded = arithmeticDnaToBytes(dna, 3, 110, blockSize);
  console.log(`${TAG} decoded: ${decoded.length} bytes`);
  const match = decoded.length === data.length && Array.from(decoded).every((b, i) => b === data[i]);
  console.log(`${TAG} match: ${match}`);
}

// Test 4: Reproduce what codec.ts does
function testCodecFlow() {
  console.log(`\n--- ${TAG} Test 4: Reproduce codec.ts flow ---`);
  // Simulate the codec.ts flow:
  // 1. innerBlock = [rsCodeword(110) + crc(2) + padding(3)] = 115 bytes
  // 2. truncatedBlock = innerBlock.slice(0, 110) = rsCodeword
  // 3. padded = new Uint8Array(110); padded.set(truncatedBlock, 0)
  // 4. dna = bytesToArithmeticDnaCrc(padded, 3, 460, 230)
  const rsCodeword = new Uint8Array(110);
  for (let i = 0; i < 110; i++) rsCodeword[i] = (i * 37 + 11) & 0xff;
  const crc = new Uint8Array([0xab, 0xcd]);
  const innerBlock = new Uint8Array(115);
  innerBlock.set(rsCodeword, 0);
  innerBlock.set(crc, 110);

  const truncatedBlock = innerBlock.slice(0, 110);
  const padded = new Uint8Array(110);
  padded.set(truncatedBlock, 0);
  console.log(`${TAG} padded.length=${padded.length}, all set: ${Array.from(padded).every((b, i) => b === rsCodeword[i])}`);

  const blockSize = 230;
  const dna = bytesToArithmeticDnaCrc(padded, 3, 460, blockSize);
  console.log(`${TAG} encoded: ${dna.length} nt`);
  console.log(`${TAG} first 60: ${dna.slice(0, 60)}`);
  console.log(`${TAG} nt 230-290: ${dna.slice(230, 290)}`);

  const decoded = arithmeticDnaToBytes(dna, 3, 110, blockSize);
  console.log(`${TAG} decoded: ${decoded.length} bytes`);
  const match = decoded.length === rsCodeword.length && Array.from(decoded).every((b, i) => b === rsCodeword[i]);
  console.log(`${TAG} match: ${match}`);
  if (!match) {
    console.log(`${TAG} expected: ${Array.from(rsCodeword).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    console.log(`${TAG} decoded : ${Array.from(decoded).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  }
}

// Main
const r1 = testTsRoundTrip();
testBlockSizes();
testPadding();
testCodecFlow();

console.log(`\n=== ${TAG} Summary: TS round-trip ${r1 ? "PASS" : "FAIL"} ===`);
process.exit(r1 ? 0 : 1);
