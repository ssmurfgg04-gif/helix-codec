/**
 * GF(2^16) Galois Field arithmetic for large-alphabet Reed-Solomon.
 *
 * Supports RS(n, k) where n ≤ 65535 symbols (vs. GF(2^8)'s 255-symbol cap).
 * This is essential for large DNA archives where the outer RS code must span
 * more than 255 oligos.
 *
 * Uses primitive polynomial 0x1100B (x^16 + x^12 + x^3 + x + 1), which is the
 * CCSDS-recommended polynomial for GF(2^16).
 *
 * Implementation:
 *   - EXP and LOG tables of 65536 entries each (128 KB total)
 *   - Multiplication via log/exp lookup (O(1))
 *   - Double-size EXP table to avoid modulo
 *
 * Reference:
 *   - CCSDS 131.0-B-3 (TM Synchronization and Channel Coding)
 *   - Reed & Solomon (1960)
 *   - Mahoraga codec (Banal 2026, arXiv:2604.20810) uses GF(2^16) for outer RS
 */

const PRIMITIVE_POLY = 0x1100b; // x^16 + x^12 + x^3 + x + 1
const FIELD_SIZE = 65536;
const MAX_EXP = 65535;

// Eagerly initialize at module load to avoid data race in worker threads
// (previously lazy-init could race if two workers call ensureTables concurrently)
let expTable: Uint16Array;
let logTable: Uint16Array;

function buildTables(): void {
  const exp = new Uint16Array(MAX_EXP * 2); // double-size to avoid modulo
  const log = new Uint16Array(FIELD_SIZE);
  let x = 1;
  for (let i = 0; i < MAX_EXP; i++) {
    exp[i] = x;
    log[x] = i;
    // Multiply by 2 (shift left), reduce mod primitive poly if bit 16 set
    x <<= 1;
    if (x & FIELD_SIZE) {
      x ^= PRIMITIVE_POLY;
    }
    x &= 0xffff;
  }
  // Replicate for wraparound-free lookup
  for (let i = MAX_EXP; i < MAX_EXP * 2; i++) {
    exp[i] = exp[i - MAX_EXP];
  }
  expTable = exp;
  logTable = log;
}

// Initialize immediately at module load — tables are always available,
// no null checks needed. This eliminates the data race that existed
// when two worker threads could call ensureTables() concurrently
// and both see expTable === null simultaneously.
buildTables();

export function gf16Add(a: number, b: number): number {
  return (a ^ b) & 0xffff;
}

export const gf16Sub = gf16Add;

export function gf16Mul(a: number, b: number): number {
  if (a === 0 || b === 0) return 0;
  return expTable[logTable[a] + logTable[b]];
}

export function gf16Div(a: number, b: number): number {
  if (a === 0) return 0;
  if (b === 0) throw new Error("GF(2^16) division by zero");
  return expTable[(logTable[a] - logTable[b] + MAX_EXP) % MAX_EXP];
}

export function gf16Inverse(a: number): number {
  if (a === 0) throw new Error("GF(2^16) inverse of zero");
  return expTable[MAX_EXP - logTable[a]];
}

export function gf16Pow(a: number, n: number): number {
  if (n === 0) return 1;
  if (a === 0) return 0;
  return expTable[(logTable[a] * n) % MAX_EXP];
}

/**
 * Evaluate a polynomial at x using Horner's method.
 * Coefficients are big-endian (poly[0] = leading).
 */
export function gf16PolyEval(poly: Uint16Array, x: number): number {
  let y = 0;
  for (let i = 0; i < poly.length; i++) {
    y = gf16Mul(y, x) ^ poly[i];
  }
  return y;
}

/** Multiply two polynomials over GF(2^16). Big-endian. */
export function gf16PolyMul(a: Uint16Array, b: Uint16Array): Uint16Array {
  const result = new Uint16Array(a.length + b.length - 1);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      result[i + j] ^= gf16Mul(a[i], b[j]);
    }
  }
  return result;
}

/** Check if the tables are built (for testing). Always true after eager init. */
export function isInitialized(): boolean {
  return true;
}

/** Force table initialization. No-op after eager init (kept for API compat). */
export function init(): void {
  // Tables are already initialized at module load. No-op.
}
