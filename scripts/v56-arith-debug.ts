#!/usr/bin/env bun
/**
 * v56 — Debug the arithmetic termination issue.
 * Find the minimum rate that produces correct round-trips.
 */

import {
  bytesToArithmeticDna,
  arithmeticDnaToBytes,
} from "../src/lib/dna/markov-arithmetic";

const TAG = "[v56-debug]";

function makeData(n: number, seed: number): Uint8Array {
  const d = new Uint8Array(n);
  let s = seed;
  for (let i = 0; i < n; i++) {
    s = (s * 1103515245 + 12345) & 0x7fffffff;
    d[i] = (s >> 16) & 0xff;
  }
  return d;
}

// Test basic round-trip at different rates
function testRate(len: number, data: Uint8Array, rate: number): boolean {
  const numBases = Math.ceil((len * 8) / rate);
  const dna = bytesToArithmeticDna(data, 3, numBases);
  const decoded = arithmeticDnaToBytes(dna, 3, len);
  return decoded.length === len && Array.from(decoded).every((b, i) => b === data[i]);
}

// Test with different data patterns
const patterns: { name: string; data: Uint8Array }[] = [
  { name: "random-seed42", data: makeData(100, 42) },
  { name: "all-zeros", data: new Uint8Array(100) },
  { name: "all-FF", data: new Uint8Array(100).fill(0xff) },
  { name: "sequential", data: (() => { const d = new Uint8Array(256); for (let i = 0; i < 256; i++) d[i] = i; return d; })() },
  { name: "alternating", data: (() => { const d = new Uint8Array(100); for (let i = 0; i < 100; i++) d[i] = i % 2 ? 0xaa : 0x55; return d; })() },
  { name: "random-seed999", data: makeData(200, 999) },
];

console.log(`${TAG} Testing different rates to find minimum viable:\n`);
console.log(`${"pattern".padEnd(20)} ${"len".padStart(4)} ${"1.95".padStart(6)} ${"1.9".padStart(6)} ${"1.85".padStart(6)} ${"1.8".padStart(6)} ${"1.75".padStart(6)} ${"1.7".padStart(6)} ${"1.6".padStart(6)} ${"1.5".padStart(6)}`);

for (const { name, data } of patterns) {
  const len = data.length;
  const results = [1.95, 1.9, 1.85, 1.8, 1.75, 1.7, 1.6, 1.5].map(r =>
    testRate(len, data, r) ? "✓" : "✗"
  );
  console.log(`${name.padEnd(20)} ${len.toString().padStart(4)} ${results.map(r => r.padStart(6)).join(" ")}`);
}

// Also test block-wise: encode N bytes into exactly K nt, decode back
console.log(`\n${TAG} Block-wise test (fixed nt, find min viable bytesPerBlock):\n`);
console.log(`${"blockSize".padEnd(12)} ${"bpb=blk/4".padStart(12)} ${"bpb=1.85".padStart(12)} ${"bpb=1.7".padStart(12)} ${"bpb=1.5".padStart(12)}`);

for (const blockSize of [40, 80, 120, 200]) {
  const results: string[] = [];
  for (const bpbCalc of [
    { label: "blk/4", bpb: Math.floor(blockSize / 4) },
    { label: "1.85", bpb: Math.floor((blockSize * 1.85) / 8) },
    { label: "1.7", bpb: Math.floor((blockSize * 1.7) / 8) },
    { label: "1.5", bpb: Math.floor((blockSize * 1.5) / 8) },
  ]) {
    const bpb = Math.max(1, bpbCalc.bpb);
    const data = makeData(bpb, 42 + blockSize);
    const dna = bytesToArithmeticDna(data, 3, blockSize);
    const decoded = arithmeticDnaToBytes(dna, 3, bpb);
    const ok = decoded.length === bpb && Array.from(decoded).every((v, i) => v === data[i]);
    results.push(ok ? "✓" : "✗");
  }
  console.log(`${("blk=" + blockSize).padEnd(12)} ${results.map(r => r.padStart(12)).join(" ")}`);
}
