/**
 * GF(2^8) Galois Field arithmetic for Reed-Solomon error correction.
 *
 * Uses primitive polynomial 0x12D (x^8 + x^5 + x^3 + x^2 + 1), matching the
 * AZTEC_DATA_8 field used by the `reedsolomon` npm package (ZXing port).
 *
 * Performance optimizations:
 *   - EXP/LOG tables (O(1) mul/div via log/exp lookup)
 *   - 256-entry MUL table (gfMulTable[a*256+b]) for batch operations
 *   - Double-size EXP table (avoids modulo in mul)
 *
 * References:
 *   - L. Rizzo, "Effective Erasure Codes for Reliable Computer Communication
 *     Protocols", ACM SIGCOMM CCR, 1997.
 *   - Plank & Luo, "User's Guide to the Reed-Solomon Codec", UTK, 2009.
 */

// GF(2^8) with primitive polynomial 0x12D (x^8 + x^5 + x^3 + x^2 + 1).
const PRIMITIVE_POLY = 0x12d;

// Eagerly initialize at module load to avoid data race in worker threads
// (previously lazy-init could race if two workers call ensureTables concurrently)
let expTable: Uint8Array;
let logTable: Uint8Array;
let mulTable: Uint8Array; // 256x256 multiplication table

function buildTables(): void {
  const exp = new Uint8Array(512); // double-size to avoid modulo in mult
  const log = new Uint8Array(256);
  let x = 1;
  for (let i = 0; i < 255; i++) {
    exp[i] = x;
    log[x] = i;
    x <<= 1;
    if (x & 0x100) x ^= PRIMITIVE_POLY;
  }
  for (let i = 255; i < 512; i++) exp[i] = exp[i - 255];
  expTable = exp;
  logTable = log;

  // Build full 64KB multiplication table for O(1) batch multiply.
  // mulTable[a * 256 + b] = a * b in GF(2^8).
  // This is 64KB but pays for itself when doing many multiplications
  // (e.g., during RS encode/decode where the same `a` multiplies many `b`s).
  const mul = new Uint8Array(256 * 256);
  for (let a = 0; a < 256; a++) {
    const aBase = a * 256;
    if (a === 0) {
      // All zeros (0 * b = 0)
      continue;
    }
    const logA = log[a];
    for (let b = 0; b < 256; b++) {
      if (b === 0) {
        mul[aBase + b] = 0;
      } else {
        mul[aBase + b] = exp[logA + log[b]];
      }
    }
  }
  mulTable = mul;
}

// Initialize immediately at module load — tables are always available,
// no null checks needed. This eliminates the data race that existed
// when two worker threads could call ensureTables() concurrently
// and both see expTable === null simultaneously.
buildTables();

export function gfAdd(a: number, b: number): number {
  return a ^ b;
}
export const gfSub = gfAdd;

export function gfMul(a: number, b: number): number {
  if (a === 0 || b === 0) return 0;
  return mulTable[a * 256 + b];
}

/** Fast: multiply many bytes by a single GF(2^8) constant. Returns new array. */
export function gfMulSlice(constant: number, data: Uint8Array): Uint8Array {
  const out = new Uint8Array(data.length);
  if (constant === 0) return out;
  const base = constant * 256;
  for (let i = 0; i < data.length; i++) {
    out[i] = mulTable[base + data[i]];
  }
  return out;
}

/** Fast: multiply many bytes by a single GF(2^8) constant, XOR into output. */
export function gfMulSliceXor(constant: number, data: Uint8Array, output: Uint8Array, outputOffset = 0): void {
  if (constant === 0) return;
  const base = constant * 256;
  for (let i = 0; i < data.length; i++) {
    output[outputOffset + i] ^= mulTable[base + data[i]];
  }
}

export function gfDiv(a: number, b: number): number {
  if (a === 0) return 0;
  if (b === 0) throw new Error("GF division by zero");
  return expTable[(logTable[a] - logTable[b] + 255) % 255];
}

export function gfPow(a: number, n: number): number {
  if (n === 0) return 1;
  if (a === 0) return 0;
  return expTable[(logTable[a] * n) % 255];
}

export function gfInverse(a: number): number {
  if (a === 0) throw new Error("GF inverse of zero");
  return expTable[255 - logTable[a]];
}

/** Polynomial multiply in GF(256). Coefficients are little-endian (index 0 = lowest power). */
export function gfPolyMul(a: Uint8Array, b: Uint8Array): Uint8Array {
  const result = new Uint8Array(a.length + b.length - 1);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      result[i + j] ^= gfMul(a[i], b[j]);
    }
  }
  return result;
}

/** Polynomial divide in GF(256) (used for syndrome calc / remainder). */
export function gfPolyDiv(
  dividend: Uint8Array,
  divisor: Uint8Array,
): { quotient: Uint8Array; remainder: Uint8Array } {
  if (divisor.length > dividend.length) {
    return { quotient: new Uint8Array(0), remainder: dividend.slice() };
  }
  const out = dividend.slice();
  const divisorLead = divisor[divisor.length - 1];
  const divisorLeadInv = gfInverse(divisorLead);
  const quotientLen = dividend.length - divisor.length + 1;
  const quotient = new Uint8Array(quotientLen);

  for (let i = quotientLen - 1; i >= 0; i--) {
    const coeff = gfMul(out[i + divisor.length - 1], divisorLeadInv);
    quotient[i] = coeff;
    if (coeff !== 0) {
      for (let j = 0; j < divisor.length; j++) {
        out[i + j] ^= gfMul(divisor[j], coeff);
      }
    }
  }
  const remainder = out.slice(0, divisor.length - 1);
  return { quotient, remainder };
}

/** Evaluate polynomial at point x using Horner's method (little-endian coeffs). */
export function gfPolyEval(poly: Uint8Array, x: number): number {
  let y = poly[poly.length - 1];
  for (let i = poly.length - 2; i >= 0; i--) {
    y = gfMul(y, x) ^ poly[i];
  }
  return y;
}
