#!/usr/bin/env bun
/**
 * v53 Direct WASM Arithmetic Decoder Test
 *
 * Encodes data with TS bytesToArithmeticDnaCrc, then decodes directly with
 * the Rust WASM test_arithmetic_decode_crc function. Bypasses the full codec
 * to isolate the encoder/decoder CRC mismatch.
 */

import { bytesToArithmeticDnaCrc } from "../src/lib/dna/markov-arithmetic";

const TAG = "[v53-wasm-direct]";

// Load WASM module directly
const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");

function dnaToBytes(dna: string): Uint8Array {
  const out = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) out[i] = dna.charCodeAt(i);
  return out;
}

function test(name: string, data: Uint8Array, blockSize: number, targetLen: number) {
  console.log(`\n--- ${TAG} ${name} ---`);
  console.log(`${TAG} data.length=${data.length}, blockSize=${blockSize}, targetLen=${targetLen}`);

  // TS encode
  const dna = bytesToArithmeticDnaCrc(data, 3, targetLen, blockSize);
  console.log(`${TAG} encoded: ${dna.length} nt`);
  console.log(`${TAG} first 40 nt: ${dna.slice(0, 40)}`);

  // Convert DNA string to Uint8Array (ASCII bytes)
  const dnaBytes = dnaToBytes(dna);

  // Rust decode via test_arithmetic_decode_crc
  // Signature: (dna: &[u8], max_homopolymer: usize, expected_bytes: usize, block_size: usize) -> Vec<u8>
  // Returns: [bytes (expected_bytes)] [erasures (expected_bytes)] — packed
  const expectedBytes = data.length;
  try {
    const decoded = wasm.test_arithmetic_decode_crc(dnaBytes, 3, expectedBytes, blockSize);
    console.log(`${TAG} decoded: ${decoded.length} bytes (expected ${expectedBytes * 2} = ${expectedBytes} data + ${expectedBytes} erasures)`);

    // Split: first half = data, second half = erasures
    const dataBytes = decoded.slice(0, expectedBytes);
    const erasureFlags = Array.from(decoded.slice(expectedBytes));
    const erasedCount = erasureFlags.filter((f: number) => f === 1).length;
    console.log(`${TAG} erasures: ${erasedCount}/${expectedBytes}`);

    // Compare
    const match = dataBytes.length === data.length && Array.from(dataBytes).every((b: number, i: number) => b === data[i]);
    console.log(`${TAG} match: ${match}`);
    if (!match) {
      console.log(`${TAG} expected: ${Array.from(data).slice(0, 16).map(b => b.toString(16).padStart(2, '0')).join(' ')}...`);
      console.log(`${TAG} decoded : ${Array.from(dataBytes).slice(0, 16).map((b: number) => b.toString(16).padStart(2, '0')).join(' ')}...`);
    }
    return match;
  } catch (e: any) {
    console.log(`${TAG} decode threw: ${e.message ?? e}`);
    return false;
  }
}

function main() {
  console.log(`${TAG} === v53 Direct WASM Arithmetic Decoder Test ===`);

  // Test 1: Small data, large block (1 block)
  const data1 = new Uint8Array(53);
  for (let i = 0; i < 53; i++) data1[i] = (i * 37 + 11) & 0xff;
  test("Test 1: 53B data, blockSize=230, targetLen=230", data1, 230, 230);

  // Test 2: 110B data, 2 blocks of 230nt each
  const data2 = new Uint8Array(110);
  for (let i = 0; i < 110; i++) data2[i] = (i * 37 + 11) & 0xff;
  test("Test 2: 110B data, blockSize=230, targetLen=460", data2, 230, 460);

  // Test 3: 55B data, 1 block of 230nt
  const data3 = new Uint8Array(55);
  for (let i = 0; i < 55; i++) data3[i] = (i * 37 + 11) & 0xff;
  test("Test 3: 55B data, blockSize=230, targetLen=230", data3, 230, 230);

  // Test 4: smaller block size
  const data4 = new Uint8Array(36);
  for (let i = 0; i < 36; i++) data4[i] = (i * 37 + 11) & 0xff;
  test("Test 4: 36B data, blockSize=153, targetLen=153", data4, 153, 153);
}

main();
