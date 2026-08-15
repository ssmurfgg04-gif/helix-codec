#!/usr/bin/env bun
/**
 * v53 TS vs Rust Arithmetic Decoder Comparison
 *
 * Encode with TS, decode with BOTH TS and Rust, compare byte-by-byte.
 * This isolates whether the Rust arithmetic decoder produces different
 * bytes than the TS arithmetic decoder.
 */

import { bytesToArithmeticDnaCrc, arithmeticDnaToBytes } from "../src/lib/dna/markov-arithmetic";

const TAG = "[v53-compare]";

const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");

function dnaToBytes(dna: string): Uint8Array {
  const out = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) out[i] = dna.charCodeAt(i);
  return out;
}

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

function compare(name: string, data: Uint8Array, blockSize: number, targetLen: number) {
  console.log(`\n--- ${TAG} ${name} ---`);
  console.log(`${TAG} data.length=${data.length}, blockSize=${blockSize}, targetLen=${targetLen}`);

  // TS encode
  const dna = bytesToArithmeticDnaCrc(data, 3, targetLen, blockSize);
  console.log(`${TAG} encoded: ${dna.length} nt`);

  // TS decode (no CRC stripping, raw bytes)
  const tsDecoded = arithmeticDnaToBytes(dna, 3, 57);  // 57 = bytesPerBlockTotal
  console.log(`${TAG} TS decoded: ${tsDecoded.length} bytes`);
  console.log(`${TAG} TS first 16: ${Array.from(tsDecoded).slice(0, 16).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);

  // TS CRC check (mimics Rust)
  const tsDataPart = tsDecoded.slice(1, 1 + data.length);
  const tsStoredCrc = tsDecoded[0];
  const tsComputedCrc = crc8ts(tsDataPart);
  console.log(`${TAG} TS CRC: stored=0x${tsStoredCrc.toString(16).padStart(2, '0')}, computed=0x${tsComputedCrc.toString(16).padStart(2, '0')}, match=${tsStoredCrc === tsComputedCrc}`);

  // Rust decode (uses test_arithmetic_decode, no CRC verification)
  const dnaBytes = dnaToBytes(dna);
  const rustDecoded = wasm.test_arithmetic_decode(dnaBytes, 3, 57);
  console.log(`${TAG} Rust decoded: ${rustDecoded.length} bytes`);
  console.log(`${TAG} Rust first 16: ${Array.from(rustDecoded).slice(0, 16).map((b: number) => b.toString(16).padStart(2, '0')).join(' ')}`);

  // Compare TS vs Rust
  const tsArr = Array.from(tsDecoded);
  const rustArr = Array.from(rustDecoded);
  let mismatches = 0;
  for (let i = 0; i < Math.min(tsArr.length, rustArr.length); i++) {
    if (tsArr[i] !== rustArr[i]) {
      if (mismatches < 5) {
        console.log(`${TAG} mismatch at byte ${i}: TS=0x${tsArr[i].toString(16).padStart(2, '0')} Rust=0x${rustArr[i].toString(16).padStart(2, '0')}`);
      }
      mismatches++;
    }
  }
  console.log(`${TAG} TS vs Rust: ${mismatches} mismatches out of ${Math.min(tsArr.length, rustArr.length)} bytes`);
}

function main() {
  console.log(`${TAG} === v53 TS vs Rust Arithmetic Decoder Comparison ===`);

  // Test A: 36B data, blockSize=153 (works in Test 4)
  const dataA = new Uint8Array(36);
  for (let i = 0; i < 36; i++) dataA[i] = (i * 37 + 11) & 0xff;
  compare("Test A: 36B, blockSize=153 (works)", dataA, 153, 153);

  // Test B: 55B data, blockSize=230 (fails in Test 3)
  const dataB = new Uint8Array(55);
  for (let i = 0; i < 55; i++) dataB[i] = (i * 37 + 11) & 0xff;
  compare("Test B: 55B, blockSize=230 (fails)", dataB, 230, 230);

  // Test C: 53B data, blockSize=230 (fails in Test 1)
  const dataC = new Uint8Array(53);
  for (let i = 0; i < 53; i++) dataC[i] = (i * 37 + 11) & 0xff;
  compare("Test C: 53B, blockSize=230 (fails)", dataC, 230, 230);
}

main();
