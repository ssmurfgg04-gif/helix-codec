// Test the Rust arithmetic decoder via WASM test_arithmetic_decode_crc
const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");
import { bytesToArithmeticDnaCrc, ARITH_CAPACITY_RATE } from "../src/lib/dna/markov-arithmetic";

const data = new Uint8Array(50);
for (let i = 0; i < 50; i++) data[i] = (i * 37) & 0xff;

const blockSize = 60;
const bytesPerBlockTotal = Math.floor(blockSize * ARITH_CAPACITY_RATE / 8);
const bytesPerBlockData = bytesPerBlockTotal - 1;
const numBlocks = Math.ceil(50 / bytesPerBlockData);
const totalDnaLen = numBlocks * blockSize;

const dna = bytesToArithmeticDnaCrc(data, 3, totalDnaLen, blockSize);
console.log(`JS encoded DNA: len=${dna.length}`);

// Convert DNA to bytes for Rust
const dnaBytes = new Uint8Array(dna.length);
for (let i = 0; i < dna.length; i++) dnaBytes[i] = dna.charCodeAt(i);

// Test Rust arithmetic decode (non-blocked)
const rustDecoded = wasm.test_arithmetic_decode(dnaBytes, 3, 50);
console.log(`Rust decoded (non-blocked): len=${rustDecoded.length}`);
let rustMatch = true;
for (let i = 0; i < 50; i++) {
  if (rustDecoded[i] !== data[i]) {
    console.log(`  Rust MISMATCH at ${i}: expected ${data[i]}, got ${rustDecoded[i]}`);
    rustMatch = false;
    if (i > 5) break;
  }
}
console.log(`Rust match: ${rustMatch}`);

// Test Rust blocked arithmetic decode
const rustDecodedBlocked = wasm.test_arithmetic_decode_blocked(dnaBytes, 3, 50, blockSize);
console.log(`\nRust decoded (blocked, bs=${blockSize}): len=${rustDecodedBlocked.length}`);
let rustBlockMatch = true;
for (let i = 0; i < 50; i++) {
  if (rustDecodedBlocked[i] !== data[i]) {
    console.log(`  Rust blocked MISMATCH at ${i}: expected ${data[i]}, got ${rustDecodedBlocked[i]}`);
    rustBlockMatch = false;
    if (i > 5) break;
  }
}
console.log(`Rust blocked match: ${rustBlockMatch}`);

// Test Rust CRC arithmetic decode
const rustDecodedCrc = wasm.test_arithmetic_decode_crc(dnaBytes, 3, 50, blockSize);
console.log(`\nRust decoded (CRC, bs=${blockSize}): len=${rustDecodedCrc.length}`);
// The first 50 bytes are data, rest is erasure bitmap
let rustCrcMatch = true;
for (let i = 0; i < 50; i++) {
  if (rustDecodedCrc[i] !== data[i]) {
    console.log(`  Rust CRC MISMATCH at ${i}: expected ${data[i]}, got ${rustDecodedCrc[i]}, erased=${rustDecodedCrc[50+i]===1}`);
    rustCrcMatch = false;
    if (i > 10) break;
  }
}
console.log(`Rust CRC match: ${rustCrcMatch}`);
let erasedCount = 0;
for (let i = 50; i < 50+50; i++) if (rustDecodedCrc[i] === 1) erasedCount++;
console.log(`Rust CRC erasures: ${erasedCount}`);
