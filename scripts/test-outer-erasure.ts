// Direct test of outer RS erasure decoding with parity erasures.
import { ReedSolomon } from "../src/lib/dna/reedsolomon";

const rs = new ReedSolomon({ n: 15, k: 10 });
const data = new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
const encoded = rs.encode(data);
console.log("encoded:", Array.from(encoded));

// Erase 2 data + 2 parity oligos (positions 2, 5, 11, 13)
const erased = [2, 5, 11, 13];
const recv = encoded.slice();
for (const p of erased) recv[p] = 0;
console.log("recv (erased):", Array.from(recv));
console.log("erased positions:", erased);

try {
  const result = rs.decodeWithErasures(recv, erased);
  console.log("decoded data:", Array.from(result.data));
  console.log("erased count:", result.erased);
  console.log("match:", JSON.stringify(Array.from(result.data)) === JSON.stringify(Array.from(data)));
} catch (e) {
  console.error("Failed:", (e as Error).message);
}
