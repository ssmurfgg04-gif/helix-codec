#!/usr/bin/env bun
/**
 * v56 — Find the exact byte where round-trip fails for longer data.
 */

import {
  bytesToArithmeticDna,
  arithmeticDnaToBytes,
} from "../src/lib/dna/markov-arithmetic";

const TAG = "[v56-find]";

function makeData(n: number, seed: number): Uint8Array {
  const d = new Uint8Array(n);
  let s = seed;
  for (let i = 0; i < n; i++) {
    s = (s * 1103515245 + 12345) & 0x7fffffff;
    d[i] = (s >> 16) & 0xff;
  }
  return d;
}

// Test different lengths for each pattern
function findFailLength(name: string, makeFn: (n: number) => Uint8Array) {
  console.log(`\n--- ${TAG} ${name} ---`);
  for (const len of [1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 50, 100]) {
    const data = makeFn(len);
    const numBases = Math.ceil((len * 8) / 1.5) + 8; // very generous
    const dna = bytesToArithmeticDna(data, 3, numBases);
    const decoded = arithmeticDnaToBytes(dna, 3, len);
    const match = decoded.length === len && Array.from(decoded).every((b, i) => b === data[i]);
    if (!match) {
      // Find first mismatch byte
      let firstBad = -1;
      for (let i = 0; i < len; i++) {
        if (decoded[i] !== data[i]) { firstBad = i; break; }
      }
      console.log(`  len=${len} FAIL at byte ${firstBad}: expected 0x${data[firstBad].toString(16)}, got 0x${decoded[firstBad].toString(16)}, numBases=${numBases}`);
      return; // stop at first failure
    } else {
      console.log(`  len=${len} OK (numBases=${numBases})`);
    }
  }
  console.log(`  all lengths pass`);
}

findFailLength("all-FF", (n) => new Uint8Array(n).fill(0xff));
findFailLength("all-zeros", (n) => new Uint8Array(n));
findFailLength("sequential", (n) => { const d = new Uint8Array(n); for (let i = 0; i < n; i++) d[i] = i & 0xff; return d; });
findFailLength("random-42", (n) => makeData(n, 42));
findFailLength("random-999", (n) => makeData(n, 999));
findFailLength("alternating-AA-55", (n) => { const d = new Uint8Array(n); for (let i = 0; i < n; i++) d[i] = i % 2 ? 0xaa : 0x55; return d; });
