/**
 * v60: Find the breaking point of the indel Viterbi.
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

function randomBytes(n: number): Uint8Array {
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) out[i] = Math.floor(Math.random() * 256);
  return out;
}

const conv = new ConvolutionalCode();
const decoder = new IndelViterbiDecoder({});

const inputBytes = 20;
const numInfoBits = inputBytes * 8;

// Inject N errors at random positions, all of one type
function testNErrors(n: number, type: "sub" | "ins" | "del", trials: number = 30): number {
  let pass = 0;
  for (let t = 0; t < trials; t++) {
    const input = randomBytes(inputBytes);
    const inputBits = bytesToBits(input);
    const encodedBits = conv.encode(inputBits);
    
    const positions = new Set<number>();
    while (positions.size < n) {
      positions.add(Math.floor(Math.random() * encodedBits.length));
    }
    const posArr = Array.from(positions).sort((a, b) => a - b);
    
    let noisy: number[];
    if (type === "sub") {
      noisy = Array.from(encodedBits);
      for (const p of posArr) noisy[p] ^= 1;
    } else if (type === "ins") {
      noisy = [];
      let j = 0;
      for (let i = 0; i < encodedBits.length; i++) {
        if (j < posArr.length && i === posArr[j]) {
          noisy.push(Math.random() < 0.5 ? 0 : 1);
          j++;
        }
        noisy.push(encodedBits[i]);
      }
    } else {
      noisy = [];
      const delSet = new Set(posArr);
      for (let i = 0; i < encodedBits.length; i++) {
        if (!delSet.has(i)) noisy.push(encodedBits[i]);
      }
    }
    
    try {
      const decoded = decoder.decode(new Uint8Array(noisy), numInfoBits);
      const decodedBytes = bitsToBytes(decoded);
      if (Buffer.from(input).equals(Buffer.from(decodedBytes))) pass++;
    } catch {}
  }
  return pass;
}

console.log("=== Substitution count sweep ===");
for (const n of [1, 2, 3, 5, 8, 10, 15, 20]) {
  const pass = testNErrors(n, "sub");
  console.log(`  ${n} subs: ${pass}/30`);
}

console.log("\n=== Insertion count sweep ===");
for (const n of [1, 2, 3, 5, 8, 10, 15, 20]) {
  const pass = testNErrors(n, "ins");
  console.log(`  ${n} ins: ${pass}/30`);
}

console.log("\n=== Deletion count sweep ===");
for (const n of [1, 2, 3, 5, 8, 10, 15, 20]) {
  const pass = testNErrors(n, "del");
  console.log(`  ${n} del: ${pass}/30`);
}

console.log("\n=== Mixed (sub+ins+del) count sweep ===");
for (const n of [1, 2, 3, 5, 8, 10]) {
  // n of each type, so 3n total errors
  let pass = 0;
  const trials = 30;
  for (let t = 0; t < trials; t++) {
    const input = randomBytes(inputBytes);
    const inputBits = bytesToBits(input);
    let encodedBits = conv.encode(inputBits);
    
    // Apply n subs
    for (let i = 0; i < n; i++) {
      const p = Math.floor(Math.random() * encodedBits.length);
      encodedBits[p] ^= 1;
    }
    // Apply n insertions
    const insArr: number[] = [];
    let noisy: number[] = [];
    for (let i = 0; i < n; i++) insArr.push(Math.floor(Math.random() * encodedBits.length));
    insArr.sort((a, b) => a - b);
    let j = 0;
    for (let i = 0; i < encodedBits.length; i++) {
      while (j < insArr.length && insArr[j] <= i) {
        noisy.push(Math.random() < 0.5 ? 0 : 1);
        j++;
      }
      noisy.push(encodedBits[i]);
    }
    while (j < insArr.length) {
      noisy.push(Math.random() < 0.5 ? 0 : 1);
      j++;
    }
    // Apply n deletions
    const delSet = new Set<number>();
    while (delSet.size < n) delSet.add(Math.floor(Math.random() * noisy.length));
    const final: number[] = [];
    for (let i = 0; i < noisy.length; i++) {
      if (!delSet.has(i)) final.push(noisy[i]);
    }
    
    try {
      const decoded = decoder.decode(new Uint8Array(final), numInfoBits);
      const decodedBytes = bitsToBytes(decoded);
      if (Buffer.from(input).equals(Buffer.from(decodedBytes))) pass++;
    } catch {}
  }
  console.log(`  ${n}+${n}+${n} (3n total): ${pass}/${trials}`);
}
