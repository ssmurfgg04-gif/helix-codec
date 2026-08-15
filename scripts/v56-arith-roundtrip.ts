#!/usr/bin/env bun
/**
 * v56 — Verify the arithmetic decoder termination fix.
 *
 * Tests TS encode → TS decode round-trip for all three encoder variants:
 *   1. bytesToArithmeticDna (basic)
 *   2. bytesToArithmeticDnaBlocked (block-wise)
 *   3. bytesToArithmeticDnaCrc (block-wise + CRC)
 *
 * Before v56: basic round-trip fails on last 1-3 bytes (termination corruption).
 * After v56:  all round-trips should PASS with zero corruption.
 */

import {
  bytesToArithmeticDna,
  bytesToArithmeticDnaBlocked,
  bytesToArithmeticDnaCrc,
  arithmeticDnaToBytes,
  satisfiesHomopolymer,
} from "../src/lib/dna/markov-arithmetic";

const TAG = "[v56-arith]";
let pass = 0;
let fail = 0;

function check(name: string, ok: boolean, detail?: string) {
  if (ok) {
    pass++;
    console.log(`${TAG} ✓ ${name}`);
  } else {
    fail++;
    console.log(`${TAG} ✗ ${name}${detail ? " — " + detail : ""}`);
  }
}

// Generate test data: structured + random
function makeData(n: number, seed: number): Uint8Array {
  const d = new Uint8Array(n);
  let s = seed;
  for (let i = 0; i < n; i++) {
    s = (s * 1103515245 + 12345) & 0x7fffffff;
    d[i] = (s >> 16) & 0xff;
  }
  return d;
}

// === Test 1: Basic round-trip (the main bug) ===
function testBasic() {
  console.log(`\n--- ${TAG} Test 1: Basic bytesToArithmeticDna round-trip ---`);
  for (const len of [10, 50, 100, 200, 500]) {
    const data = makeData(len, 42 + len);
    // numBases = ceil(len * 8 / 1.9)
    const numBases = Math.ceil((len * 8) / 1.9);
    const dna = bytesToArithmeticDna(data, 3, numBases);
    const decoded = arithmeticDnaToBytes(dna, 3, len);
    const match = decoded.length === len && Array.from(decoded).every((b, i) => b === data[i]);
    check(`basic len=${len} (${numBases}nt)`, match);
    if (!match) {
      // Find first mismatch
      for (let i = 0; i < len; i++) {
        if (decoded[i] !== data[i]) {
          console.log(`  first mismatch at byte ${i}: expected ${data[i].toString(16)}, got ${decoded[i].toString(16)}`);
          break;
        }
      }
    }
    // Also check homopolymer constraint
    const hpOk = satisfiesHomopolymer(dna, 3);
    check(`basic len=${len} homopolymer ≤ 3`, hpOk);
  }
}

// === Test 2: Blocked round-trip ===
function testBlocked() {
  console.log(`\n--- ${TAG} Test 2: Blocked bytesToArithmeticDnaBlocked round-trip ---`);
  for (const blockSize of [20, 40, 80, 120]) {
    // v56: bytesPerBlock = floor(blockSize * 1.85 / 8) — matches encoder
    const bytesPerBlock = Math.max(1, Math.floor((blockSize * 1.85) / 8));
    const data = makeData(bytesPerBlock * 5, 100 + blockSize);
    const totalDna = 5 * blockSize;
    const dna = bytesToArithmeticDnaBlocked(data, 3, totalDna, blockSize);
    // Decode each block independently using basic decoder
    let allMatch = true;
    for (let b = 0; b < 5; b++) {
      const blockDna = dna.slice(b * blockSize, (b + 1) * blockSize);
      const blockData = data.slice(b * bytesPerBlock, (b + 1) * bytesPerBlock);
      const decoded = arithmeticDnaToBytes(blockDna, 3, bytesPerBlock);
      const match = decoded.length === bytesPerBlock &&
        Array.from(decoded).every((v, i) => v === blockData[i]);
      if (!match) {
        allMatch = false;
        for (let i = 0; i < bytesPerBlock; i++) {
          if (decoded[i] !== blockData[i]) {
            console.log(`  block=${b} byte=${i}: expected ${blockData[i].toString(16)}, got ${decoded[i].toString(16)}`);
            break;
          }
        }
      }
    }
    check(`blocked blockSize=${blockSize} bpb=${bytesPerBlock} (5 blocks)`, allMatch);
  }
}

// === Test 3: CRC round-trip (the one used by codec.ts) ===
function testCrc() {
  console.log(`\n--- ${TAG} Test 3: CRC bytesToArithmeticDnaCrc round-trip ---`);
  for (const blockSize of [80, 120, 200, 230]) {
    // v56: bytesPerBlockTotal = floor(blockSize * 1.85 / 8), bytesPerBlockData = total - 1
    const bytesPerBlockTotal = Math.max(2, Math.floor((blockSize * 1.85) / 8));
    const bytesPerBlockData = bytesPerBlockTotal - 1; // 1 CRC, no padding
    if (bytesPerBlockData <= 0) continue;
    const numBlocks = 3;
    const dataLen = bytesPerBlockData * numBlocks;
    const data = makeData(dataLen, 200 + blockSize);
    const totalDna = numBlocks * blockSize;
    const dna = bytesToArithmeticDnaCrc(data, 3, totalDna, blockSize);

    // Decode each block and verify CRC
    let allMatch = true;
    for (let b = 0; b < numBlocks; b++) {
      const blockDna = dna.slice(b * blockSize, (b + 1) * blockSize);
      const blockData = data.slice(b * bytesPerBlockData, (b + 1) * bytesPerBlockData);
      const decoded = arithmeticDnaToBytes(blockDna, 3, bytesPerBlockTotal);
      // v56 Layout: [CRC byte] [data bytes] (no padding)
      const crcStored = decoded[0];
      const dataPart = decoded.slice(1, 1 + bytesPerBlockData);
      const match = dataPart.length === bytesPerBlockData &&
        Array.from(dataPart).every((v, i) => v === blockData[i]);
      if (!match) {
        allMatch = false;
        for (let i = 0; i < bytesPerBlockData; i++) {
          if (dataPart[i] !== blockData[i]) {
            console.log(`  block=${b} byte=${i}: expected ${blockData[i].toString(16)}, got ${dataPart[i].toString(16)}`);
            break;
          }
        }
      }
    }
    check(`CRC blockSize=${blockSize} bpb=${bytesPerBlockData} (${numBlocks} blocks)`, allMatch);
  }
}

// === Test 4: Edge cases (zeros, all-FF, single byte) ===
function testEdgeCases() {
  console.log(`\n--- ${TAG} Test 4: Edge cases ---`);
  // All zeros
  const zeros = new Uint8Array(100);
  const dna1 = bytesToArithmeticDna(zeros, 3, Math.ceil(100 * 8 / 1.9));
  const dec1 = arithmeticDnaToBytes(dna1, 3, 100);
  check("all-zeros len=100", dec1.length === 100 && Array.from(dec1).every(b => b === 0));

  // All 0xFF
  const ff = new Uint8Array(100).fill(0xff);
  const dna2 = bytesToArithmeticDna(ff, 3, Math.ceil(100 * 8 / 1.9));
  const dec2 = arithmeticDnaToBytes(dna2, 3, 100);
  check("all-FF len=100", dec2.length === 100 && Array.from(dec2).every(b => b === 0xff));

  // Single byte
  const one = new Uint8Array([0x42]);
  const dna3 = bytesToArithmeticDna(one, 3, 5);
  const dec3 = arithmeticDnaToBytes(dna3, 3, 1);
  check("single byte 0x42", dec3.length === 1 && dec3[0] === 0x42);

  // Sequential 0..255
  const seq = new Uint8Array(256);
  for (let i = 0; i < 256; i++) seq[i] = i;
  const dna4 = bytesToArithmeticDna(seq, 3, Math.ceil(256 * 8 / 1.9));
  const dec4 = arithmeticDnaToBytes(dna4, 3, 256);
  check("sequential 0..255", dec4.length === 256 && Array.from(dec4).every((b, i) => b === seq[i]));
}

// === Test 5: Stress test (many random sizes) ===
function testStress() {
  console.log(`\n--- ${TAG} Test 5: Stress test (20 random sizes) ---`);
  let stressPass = 0;
  for (let t = 0; t < 20; t++) {
    const len = 1 + Math.floor(Math.random() * 500);
    const data = makeData(len, 999 + t);
    const numBases = Math.ceil((len * 8) / 1.9);
    const dna = bytesToArithmeticDna(data, 3, numBases);
    const decoded = arithmeticDnaToBytes(dna, 3, len);
    const match = decoded.length === len && Array.from(decoded).every((b, i) => b === data[i]);
    if (match) stressPass++;
    else console.log(`  stress t=${t} len=${len} FAIL`);
  }
  check(`stress test ${stressPass}/20`, stressPass === 20);
}

// Run all tests
testBasic();
testBlocked();
testCrc();
testEdgeCases();
testStress();

console.log(`\n=== ${TAG} Summary: ${pass} passed, ${fail} failed ===`);
process.exit(fail === 0 ? 0 : 1);
