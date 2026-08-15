import { bytesToArithmeticDnaCrc, arithmeticDnaToBytesCrc, ARITH_CAPACITY_RATE } from "../src/lib/dna/markov-arithmetic";

// Match the actual 700nt config
const innerN = 158;
const innerDnaLen = 660;
const blockSize = 330;
const bytesPerBlockTotal = Math.max(2, Math.floor((blockSize * ARITH_CAPACITY_RATE) / 8));
const bytesPerBlockData = bytesPerBlockTotal - 1;
const numBlocks = Math.floor(innerDnaLen / blockSize);
const dataCapacity = numBlocks * bytesPerBlockData;

console.log(`innerN=${innerN}, blockSize=${blockSize}, bTotal=${bytesPerBlockTotal}, bData=${bytesPerBlockData}, numBlocks=${numBlocks}, capacity=${dataCapacity}`);

const data = new Uint8Array(innerN);
for (let i = 0; i < innerN; i++) data[i] = (i * 37 + 13) & 0xff;

const dna = bytesToArithmeticDnaCrc(data, 3, numBlocks * blockSize, blockSize);
console.log(`DNA length: ${dna.length} (expected ${numBlocks * blockSize})`);

const result = arithmeticDnaToBytesCrc(dna, 3, innerN, blockSize);
console.log(`decoded: len=${result.data.length}, erasures: ${result.erasures.filter(e=>e).length}`);

let match = true;
let firstMismatch = -1;
for (let i = 0; i < innerN; i++) {
  if (result.data[i] !== data[i]) {
    if (firstMismatch === -1) firstMismatch = i;
    console.log(`  MISMATCH at ${i}: expected ${data[i]}, got ${result.data[i]}, erased=${result.erasures[i]}`);
    match = false;
    if (i - firstMismatch > 5) break;
  }
}
console.log(`roundtrip match: ${match}`);
