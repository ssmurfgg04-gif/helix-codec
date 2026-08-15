// Test GF(2^16) and Reed-Solomon over GF(2^16).
import { ReedSolomon216 } from "../src/lib/dna/reedsolomon216";
import { gf16Mul, gf16Div, gf16Inverse, gf16Pow, gf16Add } from "../src/lib/dna/gf216";

function assert(cond: boolean, msg: string) {
  if (!cond) {
    console.error("FAIL:", msg);
    process.exit(1);
  } else {
    console.log("PASS:", msg);
  }
}

function equalArrays(a: Uint16Array, b: Uint16Array): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

console.log("\n=== GF(2^16) Tests ===\n");

// Basic field operations
assert(gf16Mul(0, 5) === 0, "gf16Mul(0,5) === 0");
assert(gf16Mul(1, 7) === 7, "gf16Mul(1,7) === 7");
assert(gf16Mul(2, 3) === 6, "gf16Mul(2,3) === 6");
assert(gf16Inverse(2) !== 0 && gf16Mul(2, gf16Inverse(2)) === 1, "2 * 2^-1 === 1");
assert(gf16Pow(2, 0) === 1, "2^0 === 1");
assert(gf16Pow(2, 65535) === 1, "2^65535 === 1 (field identity)");
assert(gf16Add(5, 3) === 6, "gf16Add(5,3) === 6 (XOR)");

// RS(10, 4) — small test
console.log("\n=== RS(10,4) over GF(2^16) ===\n");
const rs = new ReedSolomon216({ n: 10, k: 4 });
const data = new Uint16Array([0x1234, 0x5678, 0x9abc, 0xdef0]);
const encoded = rs.encode(data);
assert(encoded.length === 10, "encoded length 10");

// No errors
{
  const { data: dec, corrected, erased } = rs.decode(encoded);
  assert(equalArrays(dec, data), "no-error decode matches");
  assert(corrected === 0, "no errors corrected");
  assert(erased === 0, "no erasures");
}

// Erasures only: 5 erasures (max = nsym = 6)
{
  const recv = encoded.slice();
  const positions = [0, 2, 4, 6, 8];
  for (const p of positions) recv[p] = 0;
  const { data: dec, erased } = rs.decode(recv, positions);
  assert(equalArrays(dec, data), "5-erasure decode matches");
  assert(erased === 5, `${erased} erasures (expected 5)`);
}

// 6 erasures (boundary: nsym = 6)
{
  const recv = encoded.slice();
  const positions = [0, 1, 2, 3, 4, 5];
  for (const p of positions) recv[p] = 0;
  const { data: dec, erased } = rs.decode(recv, positions);
  assert(equalArrays(dec, data), "6-erasure decode matches (boundary)");
  assert(erased === 6, `${erased} erasures (expected 6)`);
}

// Larger: RS(40, 32) with 8 erasures
console.log("\n=== RS(40,32) over GF(2^16) ===\n");
{
  const rs2 = new ReedSolomon216({ n: 40, k: 32 });
  const d = new Uint16Array(32);
  for (let i = 0; i < 32; i++) d[i] = (i * 7919 + 31) & 0xffff;
  const enc = rs2.encode(d);

  const recv = enc.slice();
  const positions = [0, 5, 10, 15, 20, 25, 30, 35];
  for (const p of positions) recv[p] = 0;
  const { data: dec, erased } = rs2.decode(recv, positions);
  assert(equalArrays(dec, d), "RS(40,32) 8-erasure decode matches");
  assert(erased === 8, `RS(40,32) ${erased} erasures (expected 8)`);
}

// Large archive: RS(1000, 800) — would be impossible with GF(2^8)
console.log("\n=== RS(1000,800) over GF(2^16) — large archive ===\n");
{
  const rs3 = new ReedSolomon216({ n: 1000, k: 800 });
  const d = new Uint16Array(800);
  for (let i = 0; i < 800; i++) d[i] = (i * 7919 + 31) & 0xffff;
  const enc = rs3.encode(d);

  // 100 erasures (10% dropout)
  const recv = enc.slice();
  const positions: number[] = [];
  for (let i = 0; i < 100; i++) positions.push(i * 10);
  for (const p of positions) recv[p] = 0;
  const { data: dec, erased } = rs3.decode(recv, positions);
  assert(equalArrays(dec, d), "RS(1000,800) 100-erasure decode matches");
  assert(erased === 100, `RS(1000,800) ${erased} erasures (expected 100)`);
}

console.log("\nAll GF(2^16) + RS216 tests passed.");
