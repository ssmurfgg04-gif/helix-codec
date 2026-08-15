// Quick test for GF256 and ReedSolomon implementations.
import { ReedSolomon } from "../src/lib/dna/reedsolomon";
import { gfMul, gfDiv, gfInverse, gfPow } from "../src/lib/dna/gf256";

function assert(cond: boolean, msg: string) {
  if (!cond) {
    console.error("FAIL:", msg);
    process.exit(1);
  } else {
    console.log("PASS:", msg);
  }
}

function equalArrays(a: Uint8Array, b: Uint8Array): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

// --- GF256 tests ---
assert(gfMul(0, 5) === 0, "gfMul(0,5) === 0");
assert(gfMul(1, 7) === 7, "gfMul(1,7) === 7");
assert(gfInverse(2) !== 0 && gfMul(2, gfInverse(2)) === 1, "2 * 2^-1 === 1");
assert(gfPow(2, 0) === 1, "2^0 === 1");
assert(gfPow(2, 255) === 1, "2^255 === 1 (field identity)");
assert(gfDiv(6, 3) === gfMul(6, gfInverse(3)), "div consistent with mul*inv");

// --- RS test: RS(10, 4) ---
const rs = new ReedSolomon({ n: 10, k: 4 });
const data = new Uint8Array([0x10, 0x20, 0x30, 0x40]);
const encoded = rs.encode(data);
assert(encoded.length === 10, "encoded length 10");
assert(equalArrays(encoded.slice(0, 4), data), "systematic: data preserved");

// No errors -> decode works
{
  const { data: dec, corrected, erased } = rs.decode(encoded);
  assert(equalArrays(dec, data), "no-error decode matches");
  assert(corrected === 0, "no errors corrected");
  assert(erased === 0, "no erasures");
}

// Inject 3 errors (max correctable = floor((10-4)/2) = 3)
{
  const recv = encoded.slice();
  recv[0] ^= 0xAB;
  recv[3] ^= 0xCD;
  recv[7] ^= 0x12;
  const { data: dec, corrected } = rs.decode(recv);
  assert(equalArrays(dec, data), "3-error decode matches");
  assert(corrected === 3, `3 errors corrected, got ${corrected}`);
}

// Uncorrectable: 4 errors (max correctable = 3)
{
  const recv = encoded.slice();
  recv[0] ^= 0xAB;
  recv[1] ^= 0xCD;
  recv[2] ^= 0x12;
  recv[3] ^= 0x34;
  try {
    rs.decode(recv);
    assert(false, "4 errors should fail");
  } catch (e) {
    assert(true, "4 errors correctly fails");
  }
}

// Erasures: 5 erasures (max all-erasure = nsym = 6)
{
  const recv = encoded.slice();
  recv[0] = 0; // mark as erased (set to 0)
  recv[2] = 0;
  recv[4] = 0;
  recv[6] = 0;
  recv[8] = 0;
  const { data: dec, erased } = rs.decodeWithErasures(recv, [0, 2, 4, 6, 8]);
  assert(equalArrays(dec, data), "5-erasure decode matches");
  assert(erased === 5, `5 erasures, got ${erased}`);
}

// 6 erasures (max all-erasure = nsym = 6) -- boundary
{
  const recv = encoded.slice();
  recv[0] = 0;
  recv[1] = 0;
  recv[2] = 0;
  recv[3] = 0;
  recv[4] = 0;
  recv[5] = 0;
  const { data: dec, erased } = rs.decodeWithErasures(recv, [0, 1, 2, 3, 4, 5]);
  assert(equalArrays(dec, data), "6-erasure decode matches (boundary)");
  assert(erased === 6, `6 erasures, got ${erased}`);
}

// 7 erasures (too many, should fail)
{
  const recv = encoded.slice();
  for (let i = 0; i < 7; i++) recv[i] = 0;
  try {
    rs.decodeWithErasures(recv, [0, 1, 2, 3, 4, 5, 6]);
    assert(false, "7 erasures should fail");
  } catch (e) {
    assert(true, "7 erasures correctly fails");
  }
}

// Larger: RS(40, 32) — close to typical DNA-storage inner code
{
  const rs2 = new ReedSolomon({ n: 40, k: 32 });
  const d = new Uint8Array(32);
  for (let i = 0; i < 32; i++) d[i] = (i * 7 + 13) & 0xff;
  const enc = rs2.encode(d);
  // Inject 4 errors (max correctable = floor(8/2) = 4)
  const recv = enc.slice();
  recv[0] ^= 0xFF;
  recv[10] ^= 0xAA;
  recv[20] ^= 0x55;
  recv[35] ^= 0x33;
  const { data: dec, corrected } = rs2.decode(recv);
  assert(equalArrays(dec, d), "RS(40,32) 4-error decode matches");
  assert(corrected === 4, `RS(40,32) 4 errors corrected, got ${corrected}`);
}

// Erasures in larger code: RS(40, 32), 8 erasures
{
  const rs2 = new ReedSolomon({ n: 40, k: 32 });
  const d = new Uint8Array(32);
  for (let i = 0; i < 32; i++) d[i] = (i * 7 + 13) & 0xff;
  const enc = rs2.encode(d);
  const recv = enc.slice();
  const positions = [0, 5, 10, 15, 20, 25, 30, 35];
  for (const p of positions) recv[p] = 0;
  const { data: dec, erased } = rs2.decodeWithErasures(recv, positions);
  assert(equalArrays(dec, d), "RS(40,32) 8-erasure decode matches");
  assert(erased === 8, `RS(40,32) 8 erasures, got ${erased}`);
}

console.log("\nAll RS tests passed.");
