// Debug script for RS decode with full pipeline logging.
import { ReedSolomon } from "../src/lib/dna/reedsolomon";
import { gfPow, gfMul, gfInverse } from "../src/lib/dna/gf256";

// Use a smaller example for traceability: RS(7, 3), nsym=4, can correct 2 errors
const rs = new ReedSolomon({ n: 7, k: 3 });
const data = new Uint8Array([0x10, 0x20, 0x30]);
const encoded = rs.encode(data);
console.log("encoded:", Array.from(encoded));

// Inject 2 errors at positions 0 and 5 (BE)
const recv = encoded.slice();
recv[0] ^= 0xAB;
recv[5] ^= 0xCD;
console.log("recv:   ", Array.from(recv));
console.log("errors at positions: 0 (delta 0xAB), 5 (delta 0xCD)");

// Compute syndromes
const synd = new Uint8Array(rs.nsym);
for (let i = 0; i < rs.nsym; i++) {
  const x = gfPow(rs.alpha, rs.fcr + i);
  let s = 0;
  for (let j = 0; j < recv.length; j++) s = gfMul(s, x) ^ recv[j];
  synd[i] = s;
}
console.log("syndromes:", Array.from(synd));

// Lambda should be (1 + alpha^0 * x)(1 + alpha^5 * x) for positions 0 and 5 in BE
// Wait — let me think. Wikiversity: position p (BE) corresponds to root alpha^p in Chien search.
// Actually Wikiversity: chien returns nmess - 1 - i where alpha^i is the root.
// So if alpha^i is a root, error position is (n-1-i).
// For error at BE position 0: i = n-1-0 = n-1 = 6. So root is alpha^6.
// For error at BE position 5: i = n-1-5 = 1. So root is alpha^1.
// Lambda(x) = (x - alpha^6)(x - alpha^1) = (x + alpha^6)(x + alpha^1) in GF(2)
// In BE: [1, alpha^6 XOR alpha^1, alpha^6 * alpha^1]
const a1 = gfPow(2, 1);
const a6 = gfPow(2, 6);
console.log("\nExpected Lambda: roots alpha^1 =", a1, ", alpha^6 =", a6);
console.log("Lambda(x) = (x + a6)(x + a1) = x^2 + (a1+a6)x + a1*a6");
console.log("  a1+a6 =", a1 ^ a6, " a1*a6 =", gfMul(a1, a6));
console.log("Expected Lambda BE: [1,", a1 ^ a6, ",", gfMul(a1, a6), "]");

// Run BM
const errLoc = (rs as any).berlekampMassey.call(rs, synd);
console.log("\nBM errLoc:", Array.from(errLoc));

// Run Chien
const errPos = (rs as any).chienSearch.call(rs, errLoc);
console.log("Chien positions:", errPos, "(expected [0, 5])");

// Try full decode
try {
  const result = rs.decode(recv);
  console.log("\nDecode result:", result);
} catch (e) {
  console.error("\nDecode failed:", (e as Error).message);
}
