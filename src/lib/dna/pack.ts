/**
 * 2-bit DNA packing and bit-parallel operations.
 *
 * Pack encoding (same as mapping.ts dnaToBytes convention):
 *   A = 00, C = 01, G = 10, T = 11
 *
 * Four bases per byte, MSB-first:
 *   byte = (b0 << 6) | (b1 << 4) | (b2 << 2) | b3
 *   where b0 is the leftmost (most significant) base.
 *
 * Bit-parallel operations exploit the 2-bit representation for fast
 * DNA pattern matching:
 *   - Hamming distance via popcount(XOR)
 *   - Exact pattern matching via shift-and algorithm
 *   - Complement via bit-flip (A↔T, C↔G) since 00↔11, 01↔10
 *
 * Reference:
 *   - Wu & Manber (1992). "A Fast Algorithm for Multi-Pattern Searching"
 *     (shift-and / BNDM approach for small alphabets).
 *   - Myers (1999). "A Fast Bit-Vector Algorithm for Approximate String
 *     Matching" (bit-parallel approximate matching).
 */

import { simdUnpack, initSimdUnpack, isSimdAvailable } from "./simd-unpack";

// ---------------------------------------------------------------------------
// Lookup tables for fast char ↔ 2-bit conversion
// ---------------------------------------------------------------------------

/** DNA character → 2-bit code (128-entry ASCII lookup, -1 for invalid). */
const CHAR_TO_BITS: Int8Array = new Int8Array(128).fill(-1);
CHAR_TO_BITS[65] = 0b00; // A
CHAR_TO_BITS[67] = 0b01; // C
CHAR_TO_BITS[71] = 0b10; // G
CHAR_TO_BITS[84] = 0b11; // T

/** 2-bit code → DNA character. */
const BITS_TO_CHAR: string[] = ['A', 'C', 'G', 'T'];

// ---------------------------------------------------------------------------
// Pack / Unpack
// ---------------------------------------------------------------------------

/**
 * Pack a DNA string into 2-bit bytes (4 bases per byte, MSB-first).
 *
 * If the length is not a multiple of 4, the final byte is right-padded
 * with zeros in the least-significant bit pairs.
 *
 * @param dna DNA string of A/C/G/T characters
 * @returns Packed Uint8Array, length = ceil(dna.length / 4)
 *
 * @example
 *   packDnaToBits('ACGT')  // Uint8Array [0b00_01_10_11] = [0x1B]
 *   packDnaToBits('AAAA')  // Uint8Array [0x00]
 */
export function packDnaToBits(dna: string): Uint8Array {
  const numBytes = Math.ceil(dna.length / 4);
  const out = new Uint8Array(numBytes);

  for (let i = 0; i < dna.length; i++) {
    const bits = CHAR_TO_BITS[dna.charCodeAt(i)];
    if (bits < 0) {
      throw new Error(`Invalid DNA base '${dna[i]}' at position ${i}`);
    }
    const byteIdx = i >> 2; // i / 4
    const shift = 6 - ((i & 3) << 1); // 6, 4, 2, 0
    out[byteIdx] |= bits << shift;
  }

  return out;
}

/**
 * Unpack 2-bit bytes back into a DNA string.
 *
 * @param bits Packed byte array (from packDnaToBits)
 * @param numBases Number of bases to decode (may be less than bits.length * 4
 *                 if the original DNA length was not a multiple of 4)
 * @returns DNA string of length numBases
 *
 * @example
 *   unpackBitsToDna(new Uint8Array([0x1B]), 4)  // 'ACGT'
 */
export function unpackBitsToDna(bits: Uint8Array, numBases: number): string {
  if (numBases > bits.length * 4) {
    throw new Error(
      `numBases (${numBases}) exceeds packed capacity (${bits.length * 4})`,
    );
  }

  // Try SIMD path first if available
  if (isSimdAvailable()) {
    const asciiBytes = simdUnpack(bits, numBases);
    // Convert ASCII bytes to string
    let result = '';
    for (let i = 0; i < asciiBytes.length; i++) {
      result += String.fromCharCode(asciiBytes[i]);
    }
    return result;
  }

  const chars: string[] = new Array(numBases);
  for (let i = 0; i < numBases; i++) {
    const byteIdx = i >> 2;
    const shift = 6 - ((i & 3) << 1);
    const code = (bits[byteIdx] >> shift) & 0b11;
    chars[i] = BITS_TO_CHAR[code];
  }

  return chars.join('');
}

// ---------------------------------------------------------------------------
// Rolling hash (Rabin-Karp) for k-mer search
// ---------------------------------------------------------------------------

/**
 * Compute Rabin-Karp rolling hashes over a packed DNA bit array.
 *
 * Each element of the returned array is the hash of the k-mer window
 * starting at that base position (0-indexed in base space, not byte space).
 * The hash uses a 2-bit rolling update: subtract the outgoing base's
 * contribution, shift, and add the incoming base.
 *
 * @param bits Packed DNA byte array
 * @param windowSize K-mer length (number of bases)
 * @returns Uint32Array of length (numBases - windowSize + 1), where
 *          numBases = bits.length * 4. Each entry is a 32-bit hash.
 *
 * @example
 *   const bits = packDnaToBits('ACGTACGT');
 *   const hashes = rollingHash(bits, 3); // 6 hashes for 3-mers
 */
export function rollingHash(bits: Uint8Array, windowSize: number): Uint32Array {
  const numBases = bits.length * 4;
  if (windowSize > numBases || windowSize <= 0) {
    return new Uint32Array(0);
  }

  const resultLen = numBases - windowSize + 1;
  const out = new Uint32Array(resultLen);

  // Rolling hash parameters (B = 4 for 2-bit alphabet, large odd modulus)
  const B = 4;
  const MOD = 0xFFFFFFFB; // large prime < 2^32
  // Precompute B^(windowSize-1) mod MOD
  let bPow = 1;
  for (let i = 0; i < windowSize - 1; i++) {
    bPow = Math.imul(bPow, B) % MOD;
    if (bPow < 0) bPow += MOD;
  }

  // Compute initial hash for the first window
  let hash = 0;
  for (let i = 0; i < windowSize; i++) {
    const byteIdx = i >> 2;
    const shift = 6 - ((i & 3) << 1);
    const base = (bits[byteIdx] >> shift) & 0b11;
    hash = (Math.imul(hash, B) + base) % MOD;
    if (hash < 0) hash += MOD;
  }
  out[0] = hash;

  // Roll: for each subsequent position, update hash in O(1)
  for (let i = 1; i < resultLen; i++) {
    // Outgoing base (position i-1)
    const outByteIdx = (i - 1) >> 2;
    const outShift = 6 - (((i - 1) & 3) << 1);
    const outBase = (bits[outByteIdx] >> outShift) & 0b11;

    // Incoming base (position i + windowSize - 1)
    const inPos = i + windowSize - 1;
    const inByteIdx = inPos >> 2;
    const inShift = 6 - ((inPos & 3) << 1);
    const inBase = (bits[inByteIdx] >> inShift) & 0b11;

    // hash = (hash - outgoing * B^(k-1)) * B + incoming
    let subtract = Math.imul(outBase, bPow) % MOD;
    if (subtract < 0) subtract += MOD;
    hash = (hash - subtract) % MOD;
    if (hash < 0) hash += MOD;
    hash = (Math.imul(hash, B) + inBase) % MOD;
    if (hash < 0) hash += MOD;

    out[i] = hash;
  }

  return out;
}

// ---------------------------------------------------------------------------
// Bit-parallel Hamming distance
// ---------------------------------------------------------------------------

/**
 * Compute the Hamming distance between two packed DNA arrays via popcount(XOR).
 *
 * Uses the bit-parallel trick: XOR the two byte arrays, then count set bits.
 * Since each base is 2 bits, the XOR has 2 set bits per mismatching base.
 * We divide the total popcount by 2 to get the Hamming distance.
 *
 * IMPORTANT: Both arrays must have the same length and encode the same number
 * of bases. Padding zeros in the last byte are counted — use matching lengths.
 *
 * @param a First packed DNA array
 * @param b Second packed DNA array (must be same length as a)
 * @returns Hamming distance (number of mismatching bases)
 *
 * @example
 *   // ACGT vs ACAT → 1 mismatch (G↔A at position 2)
 *   bitParallelHamming(packDnaToBits('ACGT'), packDnaToBits('ACAT')) // 1
 */
export function bitParallelHamming(a: Uint8Array, b: Uint8Array): number {
  if (a.length !== b.length) {
    throw new Error(
      `Arrays must have same length: ${a.length} vs ${b.length}`,
    );
  }

  let popcount = 0;
  for (let i = 0; i < a.length; i++) {
    let x = a[i] ^ b[i];
    // Count set bits in byte (Brian Kernighan's method)
    while (x !== 0) {
      popcount++;
      x &= x - 1;
    }
  }

  // Each mismatching base contributes 2 set bits
  return popcount >> 1;
}

// ---------------------------------------------------------------------------
// Bit-parallel exact pattern matching (shift-and)
// ---------------------------------------------------------------------------

/**
 * Find all exact occurrences of a pattern in a text using the shift-and
 * algorithm (bit-parallel exact matching).
 *
 * Works on packed 2-bit DNA arrays. The pattern and text are both in
 * the 2-bit packed representation (A=00, C=01, G=10, T=11).
 *
 * The shift-and algorithm uses a bit-mask D where bit j is 1 if the
 * first j+1 characters of the pattern match the suffix of the text
 * ending at the current position. When the highest bit (bit m-1)
 * becomes 1, we have a match.
 *
 * Time complexity: O(n) where n = text length in bases (each step is O(1)
 * regardless of pattern length, as long as pattern fits in a machine word).
 * Pattern length limit: 32 bases (fits in a 32-bit integer).
 *
 * @param pattern Packed pattern (must fit in 32 bases)
 * @param text Packed text to search in
 * @returns Array of base positions where the pattern starts (0-indexed)
 *
 * @example
 *   const pat = packDnaToBits('ACG');
 *   const txt = packDnaToBits('TACGTACG');
 *   bitParallelMatch(pat, txt) // [1, 5]
 */
export function bitParallelMatch(pattern: Uint8Array, text: Uint8Array): number[] {
  const patLen = pattern.length * 4;
  const textLen = text.length * 4;

  if (patLen > 32) {
    throw new Error(
      `Pattern length ${patLen} exceeds shift-and limit of 32 bases. ` +
      `Use a longer-pattern algorithm (e.g., BNDM or Boyer-Moore).`,
    );
  }
  if (patLen === 0 || textLen === 0 || patLen > textLen) {
    return [];
  }

  // Build character masks: mask[c] has bit j set if pattern[j] == c
  const mask: Uint32Array = new Uint32Array(4);
  for (let j = 0; j < patLen; j++) {
    const byteIdx = j >> 2;
    const shift = 6 - ((j & 3) << 1);
    const code = (pattern[byteIdx] >> shift) & 0b11;
    mask[code] |= (1 << j) >>> 0;
  }

  const matches: number[] = [];
  const matchBit = (1 << (patLen - 1)) >>> 0;
  let D: number = 0;

  for (let i = 0; i < textLen; i++) {
    const byteIdx = i >> 2;
    const shift = 6 - ((i & 3) << 1);
    const code = (text[byteIdx] >> shift) & 0b11;

    // Shift-and: D = ((D << 1) | 1) & mask[currentChar]
    D = (((D << 1) | 1) >>> 0) & mask[code];

    if ((D & matchBit) !== 0) {
      matches.push(i - patLen + 1);
    }
  }

  return matches;
}

// ---------------------------------------------------------------------------
// Complement and reverse complement
// ---------------------------------------------------------------------------

/**
 * Compute the complement of a packed DNA array via bit-flip.
 *
 * In the 2-bit encoding, complement is simply XOR with 0xFF per byte:
 *   A (00) ↔ T (11)  →  XOR 11
 *   C (01) ↔ G (10)  →  XOR 11
 * Since each 2-bit pair is flipped, XOR with 0b11111111 = 0xFF works
 * for entire bytes.
 *
 * @param bits Packed DNA byte array
 * @returns Complemented packed DNA (same length)
 *
 * @example
 *   complement(packDnaToBits('ACGT'))  // packDnaToBits('TGCA')
 */
export function complement(bits: Uint8Array): Uint8Array {
  const out = new Uint8Array(bits.length);
  for (let i = 0; i < bits.length; i++) {
    out[i] = bits[i] ^ 0xFF;
  }
  return out;
}

/**
 * Compute the reverse complement of a packed DNA array.
 *
 * Reverse complement = reverse the base order + complement each base.
 * This is the biologically relevant operation for double-stranded DNA.
 *
 * The operation is performed entirely on the packed representation:
 *   1. Complement (XOR 0xFF per byte)
 *   2. Reverse the byte order
 *   3. Reverse the bit pairs within each byte (swap 2-bit positions)
 *
 * @param bits Packed DNA byte array
 * @param numBases Number of valid bases (may be < bits.length * 4 if
 *                 the original DNA length was not a multiple of 4)
 * @returns Reverse-complemented packed DNA
 *
 * @example
 *   reverseComplement(packDnaToBits('ACGT'), 4)  // packDnaToBits('ACGT') (palindrome!)
 *   reverseComplement(packDnaToBits('AACG'), 4)  // packDnaToBits('CGTT')
 */
export function reverseComplement(bits: Uint8Array, numBases: number): Uint8Array {
  if (numBases === 0) return new Uint8Array(0);

  // Step 1: Complement (XOR each byte)
  const comp = complement(bits);

  // Step 2 & 3: Reverse base order in the packed representation
  // We read bases from comp in reverse order and pack them into out
  const numBytes = Math.ceil(numBases / 4);
  const out = new Uint8Array(numBytes);

  for (let i = 0; i < numBases; i++) {
    // Source position in the complement (reading backwards)
    const srcPos = numBases - 1 - i;
    const srcByteIdx = srcPos >> 2;
    const srcShift = 6 - ((srcPos & 3) << 1);
    const base = (comp[srcByteIdx] >> srcShift) & 0b11;

    // Destination position in output
    const dstByteIdx = i >> 2;
    const dstShift = 6 - ((i & 3) << 1);
    out[dstByteIdx] |= base << dstShift;
  }

  return out;
}

// ---------------------------------------------------------------------------
// SIMD unpack re-exports
// ---------------------------------------------------------------------------

export { initSimdUnpack, isSimdAvailable } from "./simd-unpack";
