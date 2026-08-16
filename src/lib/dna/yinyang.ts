/**
 * Yin-Yang Coding (YYC) — DNA encoding at 2.0 bits/nt with alternating rules
 *
 * Yin-Yang coding maps 2 binary bits to 1 DNA base using two alternating
 * rule tables. Each base position uses either the Yin or Yang rule,
 * alternating by position.
 *
 * Algorithm:
 *   1. Convert input bytes to a bit stream (MSB first)
 *   2. Group bits into pairs (2 bits per DNA base)
 *   3. For even-positioned pairs, use the Yin rule:
 *        00→A, 01→C, 10→G, 11→T
 *   4. For odd-positioned pairs, use the Yang rule:
 *        00→C, 01→G, 10→T, 11→A
 *   5. The alternating rules ensure that the same 2-bit input produces
 *      different bases at consecutive positions (Yin[x] ≠ Yang[x] for all x),
 *      which reduces homopolymer probability vs. fixed mapping.
 *
 * Density: 2.0 bits/nt (each base encodes exactly 2 bits)
 * Speed: O(n) encode and decode
 *
 * Key property: Yin and Yang rules are complementary — they never produce
 * the same base for the same 2-bit input:
 *   Yin[0]=A ≠ Yang[0]=C,  Yin[1]=C ≠ Yang[1]=G,
 *   Yin[2]=G ≠ Yang[2]=T,  Yin[3]=T ≠ Yang[3]=A
 * This means that if two consecutive positions receive the SAME 2-bit input,
 * the output bases will always differ. Homopolymers can only arise from
 * different inputs at consecutive positions that happen to map to the same
 * base under alternating rules.
 *
 * Reference:
 *   - Ping et al. (2022). Nature Computational Science.
 *   - ntpz870817/Chamaeleo — DNA-storage-YYC
 */

// ---------------------------------------------------------------------------
// Rule tables
// ---------------------------------------------------------------------------

/**
 * Yin rule: maps a 2-bit value (0-3) to a DNA base index.
 *   0b00 → A(0), 0b01 → C(1), 0b10 → G(2), 0b11 → T(3)
 */
const YIN_RULE = [0, 1, 2, 3] as const; // A, C, G, T

/**
 * Yang rule: maps a 2-bit value (0-3) to a DNA base index.
 *   0b00 → C(1), 0b01 → G(2), 0b10 → T(3), 0b11 → A(0)
 *
 * Key property: Yin[x] ≠ Yang[x] for all x ∈ {0,1,2,3}.
 * This guarantees that the SAME input at consecutive positions
 * (which use alternating rules) always produces DIFFERENT bases.
 */
const YANG_RULE = [1, 2, 3, 0] as const; // C, G, T, A

/**
 * Inverse Yin rule: base index → 2-bit value.
 *   A(0)→0b00, C(1)→0b01, G(2)→0b10, T(3)→0b11
 */
const INV_YIN = [0, 1, 2, 3];

/**
 * Inverse Yang rule: base index → 2-bit value.
 *   A(0)→0b11, C(1)→0b00, G(2)→0b01, T(3)→0b10
 */
const INV_YANG = [3, 0, 1, 2];

const BASES = ["A", "C", "G", "T"] as const;
const BASE_TO_IDX = { A: 0, C: 1, G: 2, T: 3 } as const;

// ---------------------------------------------------------------------------
// Encode
// ---------------------------------------------------------------------------

/**
 * Encode bytes to DNA using Yin-Yang Coding (YYC).
 *
 * For each pair of bits from the input (MSB-first bit stream):
 *   - Even-positioned pairs use the Yin rule:  00→A, 01→C, 10→G, 11→T
 *   - Odd-positioned pairs use the Yang rule:   00→C, 01→G, 10→T, 11→A
 *
 * @param data Input bytes to encode
 * @returns DNA string encoded with Yin-Yang coding
 */
export function yinyangEncode(data: Uint8Array): string {
  if (data.length === 0) return "";

  const numBases = data.length * 4; // 2 bits/base, 8 bits/byte → 4 bases/byte
  const output: string[] = new Array(numBases);

  let pos = 0; // base position (determines Yin/Yang rule selection)

  for (let byteIdx = 0; byteIdx < data.length; byteIdx++) {
    const byte = data[byteIdx];
    // Process 8 bits as 4 pairs (MSB first)
    for (let pair = 0; pair < 4; pair++) {
      // Extract 2-bit pair: bits [7-2*pair, 6-2*pair] from the byte
      const twoBits = (byte >>> (6 - 2 * pair)) & 0b11;

      // Alternate: even positions → Yin, odd positions → Yang
      const baseIdx = (pos & 1) === 0
        ? YIN_RULE[twoBits]
        : YANG_RULE[twoBits];

      output[pos] = BASES[baseIdx];
      pos++;
    }
  }

  return output.join("");
}

// ---------------------------------------------------------------------------
// Decode
// ---------------------------------------------------------------------------

/**
 * Decode Yin-Yang encoded DNA back to bytes.
 *
 * For each DNA base at position i:
 *   - Even positions were encoded with the Yin rule  → use inverse Yin
 *   - Odd positions were encoded with the Yang rule   → use inverse Yang
 *
 * The recovered 2-bit values are reassembled into bytes (MSB first).
 *
 * @param dna Yin-Yang encoded DNA string
 * @param numBytes Expected number of output bytes
 * @returns Decoded byte array of length numBytes
 */
export function yinyangDecode(dna: string, numBytes: number): Uint8Array {
  if (dna.length === 0 || numBytes === 0) return new Uint8Array(0);

  const result = new Uint8Array(numBytes);

  for (let pos = 0; pos < dna.length; pos++) {
    const baseIdx = BASE_TO_IDX[dna[pos] as keyof typeof BASE_TO_IDX];
    if (baseIdx === undefined) {
      throw new Error(`Invalid DNA base at position ${pos}: "${dna[pos]}"`);
    }

    // Recover 2-bit value using the inverse of the rule used at this position
    const twoBits = (pos & 1) === 0
      ? INV_YIN[baseIdx]
      : INV_YANG[baseIdx];

    // Pack the 2 bits into the correct byte
    const bitOffset = pos * 2;          // absolute bit offset in output
    const byteIdx = bitOffset >>> 3;    // bitOffset / 8
    const shift = 7 - (bitOffset & 7);  // MSB-first: bit 7 is highest

    if (byteIdx < numBytes) {
      if (shift >= 1) {
        // Both bits fit in the same byte
        result[byteIdx] |= twoBits << (shift - 1);
      } else {
        // shift === 0: bits span two bytes
        result[byteIdx] |= (twoBits >>> 1) & 1;
        if (byteIdx + 1 < numBytes) {
          result[byteIdx + 1] |= (twoBits & 1) << 7;
        }
      }
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// Density and length computations
// ---------------------------------------------------------------------------

/**
 * Compute the DNA length needed to encode n bytes.
 *
 * YYC encodes 2 bits per DNA base.
 * For n bytes = 8n bits, we need 8n / 2 = 4n bases.
 *
 * @param numBytes Number of bytes to encode
 * @returns Required DNA length in bases
 */
export function yinyangDnaLength(numBytes: number): number {
  return numBytes * 4;
}

/**
 * Compute the number of bytes that can be stored in a DNA sequence of given length.
 *
 * @param dnaLength DNA length in bases
 * @returns Number of bytes that can be encoded
 */
export function yinyangBytesFromDnaLength(dnaLength: number): number {
  return Math.floor(dnaLength / 4);
}

/**
 * Compute the effective density in bits per nucleotide.
 *
 * For YYC: density = 2.0 bits/nt (constant, each base encodes exactly 2 bits).
 */
export function yinyangDensity(): number {
  return 2.0;
}

// ---------------------------------------------------------------------------
// Constraint verification
// ---------------------------------------------------------------------------

/**
 * Verify that a Yin-Yang encoded DNA sequence satisfies constraints.
 *
 * Checks:
 *   1. No homopolymers > maxHomo (default: 3)
 *   2. GC content within [0.4, 0.6]
 *   3. All bases are valid {A, C, G, T}
 *
 * YYC's alternating rules reduce homopolymer probability (the same input
 * at consecutive positions always produces different bases), but different
 * inputs can still produce the same base, so homopolymers of length 2+
 * are possible though relatively rare for random data.
 *
 * @param dna DNA string to verify
 * @param maxHomo Maximum allowed homopolymer length (default: 3)
 * @returns Verification result
 */
export function yinyangVerifyConstraints(dna: string, maxHomo: number = 3): {
  valid: boolean;
  maxHomopolymer: number;
  gcContent: number;
  invalidBases: number;
} {
  // Check for invalid bases
  let invalidBases = 0;
  for (let i = 0; i < dna.length; i++) {
    if (!(dna[i] in BASE_TO_IDX)) {
      invalidBases++;
    }
  }

  // Compute max homopolymer run
  let maxRun = dna.length > 0 ? 1 : 0;
  let run = 1;
  for (let i = 1; i < dna.length; i++) {
    if (dna[i] === dna[i - 1]) {
      run++;
      if (run > maxRun) maxRun = run;
    } else {
      run = 1;
    }
  }

  // Compute GC content
  let gc = 0;
  for (let i = 0; i < dna.length; i++) {
    if (dna[i] === "G" || dna[i] === "C") gc++;
  }
  const gcFrac = dna.length > 0 ? gc / dna.length : 0;

  return {
    valid: invalidBases === 0 && maxRun <= maxHomo && gcFrac >= 0.4 && gcFrac <= 0.6,
    maxHomopolymer: maxRun,
    gcContent: gcFrac,
    invalidBases,
  };
}
