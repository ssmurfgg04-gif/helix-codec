/**
 * DNA base mapping with biological constraint enforcement.
 *
 * Uses direct 2-bit -> base mapping (00=A, 01=C, 10=G, 11=T) for high density
 * (2 bits/base), then enforces constraints via per-oligo screening with
 * seed-based re-encoding retries (DNA Fountain approach, Erlich & Zielinski 2017).
 *
 * Constraints enforced at the oligo level:
 *   - GC content: 40-60% (configurable)
 *   - Max homopolymer run: 3 (configurable)
 *
 * If an oligo fails screening, its source bytes are XOR-ed with a pseudo-random
 * keystream derived from a seed, and re-encoded. The seed (1 byte) is appended
 * to the oligo's index field, so the decoder can reverse the XOR after recovery.
 *
 * Reference:
 *   - Erlich & Zielinski (2017), "DNA Fountain enables a robust and efficient
 *     storage architecture", Science 355:6328.
 */

export type Base = "A" | "C" | "G" | "T";

const BASE_TO_BITS: Record<Base, number> = { A: 0b00, C: 0b01, G: 0b10, T: 0b11 };
const BITS_TO_BASE: Base[] = ["A", "C", "G", "T"];

/** Map 2 bits -> base (direct mapping). */
export function bitsToBase(bits: number): Base {
  return BITS_TO_BASE[bits & 0b11];
}

/** Inverse: base -> 2 bits. */
export function baseToBits(base: Base): number {
  return BASE_TO_BITS[base];
}

/**
 * Encode bytes to DNA using direct 2-bit mapping.
 * Each byte -> 4 bases. Output length = 4 * input length.
 */
// Pre-allocated buffer for bytesToDna (reused across calls)
const DNA_BUF = new Uint8Array(1024);

export function bytesToDna(data: Uint8Array): string {
  const len = data.length;
  // Use String.fromCharCode with spread for fast string building
  const chars: string[] = new Array(len * 4);
  for (let i = 0; i < len; i++) {
    const byte = data[i];
    const off = i * 4;
    chars[off] = BITS_TO_BASE[(byte >> 6) & 0b11];
    chars[off + 1] = BITS_TO_BASE[(byte >> 4) & 0b11];
    chars[off + 2] = BITS_TO_BASE[(byte >> 2) & 0b11];
    chars[off + 3] = BITS_TO_BASE[byte & 0b11];
  }
  return chars.join("");
}

/**
 * Decode DNA string to bytes (inverse of bytesToDna).
 * @param dna string of ACGT (length must be multiple of 4)
 */
// Lookup table for fast DNA char → bits conversion (avoids property access)
const DNA_CHAR_LOOKUP: Int8Array = new Int8Array(128).fill(-1);
DNA_CHAR_LOOKUP[65] = 0; // A
DNA_CHAR_LOOKUP[67] = 1; // C
DNA_CHAR_LOOKUP[71] = 2; // G
DNA_CHAR_LOOKUP[84] = 3; // T

export function dnaToBytes(dna: string): Uint8Array {
  if (dna.length % 4 !== 0) {
    throw new Error(`DNA length ${dna.length} is not a multiple of 4`);
  }
  const out = new Uint8Array(dna.length / 4);
  for (let i = 0; i < out.length; i++) {
    const off = i * 4;
    const b0 = DNA_CHAR_LOOKUP[dna.charCodeAt(off)];
    const b1 = DNA_CHAR_LOOKUP[dna.charCodeAt(off + 1)];
    const b2 = DNA_CHAR_LOOKUP[dna.charCodeAt(off + 2)];
    const b3 = DNA_CHAR_LOOKUP[dna.charCodeAt(off + 3)];
    if (b0 < 0 || b1 < 0 || b2 < 0 || b3 < 0) {
      throw new Error(`Invalid DNA base at position ${off}`);
    }
    out[i] = (b0 << 6) | (b1 << 4) | (b2 << 2) | b3;
  }
  return out;
}

/** Compute GC content of a DNA string (fraction, 0..1). Optimized with charCodeAt. */
export function gcContent(dna: string): number {
  if (dna.length === 0) return 0;
  let gc = 0;
  for (let i = 0; i < dna.length; i++) {
    const c = dna.charCodeAt(i);
    // G=71, C=67
    if (c === 71 || c === 67) gc++;
  }
  return gc / dna.length;
}

/** Find the longest homopolymer run in a DNA string. Optimized with charCodeAt. */
export function maxHomopolymerRun(dna: string): number {
  if (dna.length === 0) return 0;
  let max = 1;
  let current = 1;
  for (let i = 1; i < dna.length; i++) {
    if (dna.charCodeAt(i) === dna.charCodeAt(i - 1)) {
      current++;
      if (current > max) max = current;
    } else {
      current = 1;
    }
  }
  return max;
}

/** Count homopolymer runs of length >= threshold. */
export function countHomopolymers(dna: string, threshold = 3): number {
  let count = 0;
  let run = 1;
  for (let i = 1; i < dna.length; i++) {
    if (dna[i] === dna[i - 1]) {
      run++;
    } else {
      if (run >= threshold) count++;
      run = 1;
    }
  }
  if (run >= threshold) count++;
  return count;
}

export interface DnaConstraints {
  gcMin: number; // minimum acceptable GC fraction
  gcMax: number;
  maxHomopolymer: number; // maximum allowed homopolymer run
}

export const DEFAULT_CONSTRAINTS: DnaConstraints = {
  gcMin: 0.4,
  gcMax: 0.6,
  maxHomopolymer: 3,
};

/** Check if a DNA string satisfies the given constraints. Optimized with early exit. */
export function satisfiesConstraints(
  dna: string,
  constraints: DnaConstraints = DEFAULT_CONSTRAINTS,
): boolean {
  // Check homopolymer first (more likely to fail, cheaper to check)
  let current = 1;
  let gc = 0;
  for (let i = 0; i < dna.length; i++) {
    const c = dna.charCodeAt(i);
    if (c === 71 || c === 67) gc++; // G=71, C=67
    if (i > 0) {
      if (c === dna.charCodeAt(i - 1)) {
        current++;
        if (current > constraints.maxHomopolymer) return false; // early exit
      } else {
        current = 1;
      }
    }
  }
  // Check GC content
  const gcFrac = gc / dna.length;
  if (gcFrac < constraints.gcMin || gcFrac > constraints.gcMax) return false;
  return true;
}

/**
 * XOR a byte array with a pseudo-random keystream derived from `seed`.
 * Used to "re-encode" oligos that fail constraint screening: XOR source bytes
 * with the keystream, then re-encode to DNA. The seed (1 byte) is stored in
 * the oligo header so the decoder can reverse this operation after recovery.
 *
 * Uses xorshift32 for deterministic, fast keystream generation.
 */
export function xorWithSeed(data: Uint8Array, seed: number): Uint8Array {
  const out = new Uint8Array(data.length);
  let state = (seed >>> 0) || 1; // xorshift requires nonzero state
  for (let i = 0; i < data.length; i++) {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    state = state >>> 0;
    out[i] = data[i] ^ (state & 0xff);
  }
  return out;
}

/** Reverse the XOR-with-seed operation (same function — XOR is its own inverse). */
export const unxorWithSeed = xorWithSeed;

/**
 * Fixed XOR pattern applied to address bytes to break up zero runs (which
 * would otherwise create long homopolymer sequences for small oligo indices).
 *
 * Each pattern is chosen so that:
 *   1. Consecutive 2-bit groups WITHIN the byte differ (no intra-byte homopolymer).
 *   2. The byte ENDS with a base that DIFFERS from the next byte's STARTING base
 *      (no inter-byte homopolymer at address boundaries).
 *
 * Verification (2-bit groups, A=00, C=01, G=10, T=11):
 *   0x1B = 00 01 10 11 -> A C G T  (ends with T)
 *   0x4B = 01 00 10 11 -> C A G T  (starts with C != T, ends with T)
 *   0x24 = 00 10 01 00 -> A G C A  (starts with A != T, ends with A)
 *   0x6D = 01 10 11 01 -> C G T C  (starts with C != A, ends with C)
 *
 * For oligo 0 (rawAddress = [0,0,0,0]), whitened = [0x1B, 0x4B, 0x24, 0x6D]
 * -> DNA "ACGTCAGTAGCACGTC" with NO homopolymer.
 */
export const ADDRESS_WHITENING = [0x1b, 0x4b, 0x24, 0x6d];

/**
 * Apply fixed whitening to address bytes (3-byte index + 1-byte seed).
 * The whitening breaks up zero patterns that would create homopolymers.
 * Reversible: call again to undo (XOR is its own inverse).
 */
export function whitenAddress(address: Uint8Array): Uint8Array {
  const out = new Uint8Array(address.length);
  for (let i = 0; i < address.length; i++) {
    out[i] = address[i] ^ ADDRESS_WHITENING[i % ADDRESS_WHITENING.length];
  }
  return out;
}

/** Reverse the address whitening (same function). */
export const unwhitenAddress = whitenAddress;
