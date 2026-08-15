const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");
import { bytesToArithmeticDnaCrc, ARITH_CAPACITY_RATE } from "../src/lib/dna/markov-arithmetic";

const innerN = 158;
const blockSize = 330;
const numBlocks = 2;
const data = new Uint8Array(innerN);
for (let i = 0; i < innerN; i++) data[i] = (i * 37 + 13) & 0xff;

const dna = bytesToArithmeticDnaCrc(data, 3, numBlocks * blockSize, blockSize);
const dnaBytes = new Uint8Array(dna.length);
for (let i = 0; i < dna.length; i++) dnaBytes[i] = dna.charCodeAt(i);

// Rust CRC decode
const rustDecoded = wasm.test_arithmetic_decode_crc(dnaBytes, 3, innerN, blockSize);
console.log(`Rust decoded: len=${rustDecoded.length} (expect ${innerN * 2})`);

let mismatches = 0;
let erasures = 0;
for (let i = 0; i < innerN; i++) {
  if (rustDecoded[i] !== data[i]) {
    mismatches++;
    console.log(`  MISMATCH at ${i}: expected ${data[i]}, got ${rustDecoded[i]}, erased=${rustDecoded[innerN+i]===1}`);
  }
  if (rustDecoded[innerN+i] === 1) erasures++;
}
console.log(`mismatches: ${mismatches}, erasures: ${erasures}`);
