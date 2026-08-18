// Detailed trace of encode/decode for byte 15 at the failing position
const BASES = ['A', 'C', 'G', 'T'] as const;
const BASE_TO_IDX: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };
const MAP_4TO3_IDX = [
  [1, 2, 3, 1], [0, 2, 3, 0], [0, 1, 3, 0], [0, 1, 2, 0],
];
const GC_CODEBOOKS = [[0, 3, 1, 2], [0, 1, 2, 3], [1, 2, 0, 3]];
const INV_GC_CODEBOOKS = GC_CODEBOOKS.map((cb) => {
  const inv: number[] = [0,0,0,0];
  for (let c = 0; c < 4; c++) inv[cb[c]] = c;
  return inv;
});

// Byte 15 = 85 = 01010101
// bit pairs: 01, 01, 01, 01 → codes 1,1,1,1
// Balanced codebook (1): 01→C, so all 4 pairs → C
// Expected output: CCCC
// But encoded output: CCCG

// The last base being G instead of C means the derangement was applied.
// At position 63, what was the run state?

// Let's manually trace from the encoded DNA
const dna = "AGGGCATTCTCAGCGCCGGCGCAGAACAACGTGGCAGTGTAGTCTCCGTGCGAAGTATAACCCG";

// Check run state at position 60 (start of byte 15)
let prev = -1, run = 0;
for (let i = 0; i < 60; i++) {
  const baseIdx = BASE_TO_IDX[dna[i]];
  if (baseIdx === prev) run++; else { run = 1; prev = baseIdx; }
}
console.log(`State at position 60: prev=${BASES[prev]} (${prev}), runLen=${run}`);

// Trace encode of byte 15 (85 = 01010101)
const byte15 = 85;
const cbIdx = 1; // Balanced
const cb = GC_CODEBOOKS[cbIdx];
console.log(`\nEncoding byte 15 (${byte15} = ${byte15.toString(2).padStart(8,'0')}):`);
console.log(`  Codebook: ${cbIdx} (Balanced: ${cb})`);
console.log(`  State: prev=${BASES[prev]} (${prev}), runLen=${run}`);

for (let pair = 0; pair < 4; pair++) {
  const bits = (byte15 >> (6 - pair * 2)) & 0b11;
  let baseIdx = cb[bits];
  const atLimit = run >= 3 && prev >= 0 && baseIdx === prev;
  const originalBaseIdx = baseIdx;
  if (atLimit) {
    baseIdx = MAP_4TO3_IDX[prev][bits];
  }
  console.log(`  pair ${pair}: bits=${bits.toString(2).padStart(2,'0')}, cb[bits]=${BASES[originalBaseIdx]}(${originalBaseIdx})` +
    (atLimit ? ` → DERANGEMENT → ${BASES[baseIdx]}(${baseIdx})` : '') +
    `, runLen=${run}${atLimit ? ' (AT LIMIT!)' : ''}`);
  if (baseIdx === prev) run++; else { run = 1; prev = baseIdx; }
}

// Now trace the decode of byte 15
console.log(`\nDecoding byte 15 from DNA: ${dna.slice(60, 64)}`);
// Reset state
prev = -1; run = 0;
for (let i = 0; i < 60; i++) {
  const baseIdx = BASE_TO_IDX[dna[i]];
  if (baseIdx === prev) run++; else { run = 1; prev = baseIdx; }
}
console.log(`  State: prev=${BASES[prev]} (${prev}), runLen=${run}`);

const invCb = INV_GC_CODEBOOKS[cbIdx];
let decodedByte = 0;
for (let pair = 0; pair < 4; pair++) {
  const base = dna[60 + pair];
  const baseIdx = BASE_TO_IDX[base];
  let code: number;

  if (run >= 3 && prev >= 0) {
    code = invCb[baseIdx];
    const cbWouldProduce = cb[code];
    console.log(`  pair ${pair}: base=${base}(${baseIdx}), invCb[${baseIdx}]=${code}, cb[${code}]=${BASES[cbWouldProduce]}(${cbWouldProduce})`);
    if (cbWouldProduce === prev) {
      console.log(`    → cb[code] === prev (${BASES[prev]}), DERANGEMENT was applied`);
      for (let c = 0; c < 4; c++) {
        if (MAP_4TO3_IDX[prev][c] === baseIdx) {
          code = c;
          console.log(`    → Found derangement code: ${c} (MAP_4TO3_IDX[${prev}][${c}] = ${baseIdx})`);
          break;
        }
      }
    } else {
      console.log(`    → cb[code] !== prev, no derangement`);
    }
  } else {
    code = invCb[baseIdx];
    console.log(`  pair ${pair}: base=${base}(${baseIdx}), code=${code}`);
  }

  decodedByte = (decodedByte << 2) | code;
  if (baseIdx === prev) run++; else { run = 1; prev = baseIdx; }
}
console.log(`  Decoded byte: ${decodedByte} (${decodedByte.toString(2).padStart(8,'0')})`);
console.log(`  Expected:     ${byte15} (${byte15.toString(2).padStart(8,'0')})`);
