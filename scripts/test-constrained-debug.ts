/**
 * Debug the constrained mapping encode/decode mismatch.
 */
import { bytesToConstrainedDna, constrainedDnaToBytesWithErasure } from "../src/lib/dna/constrained-mapping";

const len = 16;
const data = new Uint8Array(len);
for (let i = 0; i < len; i++) data[i] = (i * 37 + 42) & 0xff;

const encoded = bytesToConstrainedDna(data, 3, 0.4, 0.6);
const decoded = constrainedDnaToBytesWithErasure(encoded.dna, 3, len, encoded.codebookSequence);

console.log("Input data:", Array.from(data));
console.log("Encoded DNA:", encoded.dna);
console.log("Codebook seq:", encoded.codebookSequence);
console.log("Decoded data:", Array.from(decoded.data));
console.log("Erasures:", decoded.erasures.filter(e => e).length);

// Find the first mismatch
for (let i = 0; i < len; i++) {
  if (data[i] !== decoded.data[i]) {
    console.log(`\nFirst mismatch at byte ${i}: expected=${data[i]}, got=${decoded.data[i]}`);
    console.log(`  expected binary: ${data[i].toString(2).padStart(8, '0')}`);
    console.log(`  got binary:      ${decoded.data[i].toString(2).padStart(8, '0')}`);
    // Find which bit pair differs
    for (let p = 0; p < 4; p++) {
      const expBits = (data[i] >> (6 - p * 2)) & 0b11;
      const gotBits = (decoded.data[i] >> (6 - p * 2)) & 0b11;
      if (expBits !== gotBits) {
        console.log(`  bit pair ${p}: expected=${expBits.toString(2).padStart(2,'0')}, got=${gotBits.toString(2).padStart(2,'0')}`);
        console.log(`  DNA pos ${i*4+p}: ${encoded.dna[i*4+p]}`);
      }
    }
    break;
  }
}

// Also test WITHOUT codebooks (gcMin=0, gcMax=1 to effectively disable GC steering)
console.log("\n\n=== Without GC codebooks (gcMin=0, gcMax=1) ===");
const encoded2 = bytesToConstrainedDna(data, 3, 0, 1);
const decoded2 = constrainedDnaToBytesWithErasure(encoded2.dna, 3, len, encoded2.codebookSequence);
let match2 = true;
for (let i = 0; i < len; i++) {
  if (data[i] !== decoded2.data[i]) { match2 = false; break; }
}
console.log("Roundtrip:", match2 ? "PASS" : "FAIL");
if (!match2) {
  for (let i = 0; i < len; i++) {
    if (data[i] !== decoded2.data[i]) {
      console.log(`  Mismatch at byte ${i}: expected=${data[i]}, got=${decoded2.data[i]}`);
      break;
    }
  }
}
