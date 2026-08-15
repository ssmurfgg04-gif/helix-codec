// Unit test: arithmetic encode → decode roundtrip (no LDPC, no channel)
import { bytesToArithmeticDnaCrc, arithmeticDnaToBytesCrc, ARITH_CAPACITY_RATE } from "../src/lib/dna/markov-arithmetic";

// Test 1: encode 50 bytes, decode, verify
const data = new Uint8Array(50);
for (let i = 0; i < 50; i++) data[i] = (i * 37) & 0xff;

const blockSize = 60; // nt per block
const bytesPerBlockTotal = Math.floor(blockSize * ARITH_CAPACITY_RATE / 8);
const bytesPerBlockData = bytesPerBlockTotal - 1;
const numBlocks = Math.ceil(50 / bytesPerBlockData);
const totalDnaLen = numBlocks * blockSize;

console.log(`blockSize=${blockSize}, bytesPerBlockTotal=${bytesPerBlockTotal}, bytesPerBlockData=${bytesPerBlockData}, numBlocks=${numBlocks}, totalDnaLen=${totalDnaLen}`);

const dna = bytesToArithmeticDnaCrc(data, 3, totalDnaLen, blockSize);
console.log(`encoded DNA: len=${dna.length} (expected ${totalDnaLen})`);

const result = arithmeticDnaToBytesCrc(dna, 3, 50, blockSize);
console.log(`decoded: len=${result.data.length}, erasures: ${result.erasures.filter(e=>e).length}`);

let match = true;
for (let i = 0; i < 50; i++) {
  if (result.data[i] !== data[i]) {
    console.log(`MISMATCH at byte ${i}: expected ${data[i]}, got ${result.data[i]}, erased=${result.erasures[i]}`);
    match = false;
    if (i > 10) break;
  }
}
console.log(`roundtrip match: ${match}`);

// Test 2: simulate 1 substitution error and see what happens
const dnaWithErr = dna.split('');
dnaWithErr[30] = dnaWithErr[30] === 'A' ? 'C' : 'A';
const dnaErr = dnaWithErr.join('');
const result2 = arithmeticDnaToBytesCrc(dnaErr, 3, 50, blockSize);
console.log(`\nWith 1 sub at nt 30: erasures=${result2.erasures.filter(e=>e).length}`);
let mismatches = 0;
for (let i = 0; i < 50; i++) {
  if (result2.data[i] !== data[i] && !result2.erasures[i]) mismatches++;
}
console.log(`silent errors (wrong byte, not erased): ${mismatches}`);
