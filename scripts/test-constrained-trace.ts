/**
 * Trace encode/decode for a single byte to find the bug.
 */
const BASES = ['A', 'C', 'G', 'T'] as const;
const BASE_TO_IDX: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

const MAP_4TO3_IDX = [
  [1, 2, 3, 1],  // prev=A
  [0, 2, 3, 0],  // prev=C
  [0, 1, 3, 0],  // prev=G
  [0, 1, 2, 0],  // prev=T
];

const GC_CODEBOOKS = [
  [0, 3, 1, 2],  // A-rich
  [0, 1, 2, 3],  // Balanced
  [1, 2, 0, 3],  // C-rich
];

const INV_GC_CODEBOOKS = GC_CODEBOOKS.map((cb) => {
  const inv = new Array(4) as number[];
  for (let code = 0; code < 4; code++) inv[cb[code]] = code;
  return inv;
});

function selectCodebookByPosition(position: number, cycleLen: number = 16): number {
  const phase = Math.floor(position / cycleLen) % 3;
  if (phase === 1) return 0;
  if (phase === 2) return 2;
  return 1;
}

// Encode the full 16-byte sequence to get the state at byte 15
const data = new Uint8Array(16);
for (let i = 0; i < 16; i++) data[i] = (i * 37 + 42) & 0xff;

let prevIdx = -1;
let runLen = 0;
let totalBases = 0;
const encodedBases: string[] = [];

for (let i = 0; i < 16; i++) {
  const byte = data[i];
  const cbIdx = selectCodebookByPosition(totalBases);
  const cb = GC_CODEBOOKS[cbIdx];

  for (let pair = 0; pair < 4; pair++) {
    const bits = (byte >> (6 - pair * 2)) & 0b11;
    let baseIdx = cb[bits];
    if (runLen >= 3 && prevIdx >= 0 && baseIdx === prevIdx) {
      baseIdx = MAP_4TO3_IDX[prevIdx][bits];
    }
    if (baseIdx === prevIdx) { runLen++; } else { runLen = 1; prevIdx = baseIdx; }
    totalBases++;
    encodedBases.push(BASES[baseIdx]);
  }
}

const dna = encodedBases.join('');
console.log("Encoded DNA:", dna);
console.log("Last 4 bases (byte 15):", dna.slice(60, 64));

// Now decode byte 15
// Need to reconstruct the encoder state at position 60
// Let's just decode the full thing
let decPrevIdx = -1;
let decRunLen = 0;
let decTotalBases = 0;

for (let i = 0; i < 16; i++) {
  const cbIdx = selectCodebookByPosition(decTotalBases);
  const cb = GC_CODEBOOKS[cbIdx];
  const invCb = INV_GC_CODEBOOKS[cbIdx];
  let byte = 0;

  for (let pair = 0; pair < 4; pair++) {
    const base = dna[i * 4 + pair];
    const baseIdx = BASE_TO_IDX[base];
    const bitIdx = i * 8 + pair * 2;

    let code: number;
    if (decRunLen >= 3 && decPrevIdx >= 0) {
      code = invCb[baseIdx];
      if (cb[code] === decPrevIdx) {
        // Derangement was applied
        for (let c = 0; c < 4; c++) {
          if (MAP_4TO3_IDX[decPrevIdx][c] === baseIdx) {
            code = c;
            break;
          }
        }
      }
    } else {
      code = invCb[baseIdx];
    }

    byte = (byte << 2) | code;

    if (baseIdx === decPrevIdx) { decRunLen++; } else { decRunLen = 1; decPrevIdx = baseIdx; }
    decTotalBases++;
  }

  if (i === 15 || data[i] !== byte) {
    console.log(`Byte ${i}: expected=${data[i]} (${data[i].toString(2).padStart(8,'0')}), got=${byte} (${byte.toString(2).padStart(8,'0')})`);
    if (data[i] !== byte) {
      console.log("  MISMATCH!");
      // Trace the encode/decode for this byte
      const encByte = data[i];
      const cbIdx2 = selectCodebookByPosition(i * 4);
      console.log(`  Encode codebook: ${cbIdx2} (${cbIdx2===0?'A-rich':cbIdx2===1?'Balanced':'C-rich'})`);
      console.log(`  Decode codebook: ${cbIdx} (${cbIdx===0?'A-rich':cbIdx===1?'Balanced':'C-rich'})`);
    }
  }
}
