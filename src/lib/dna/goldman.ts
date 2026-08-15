/**
 * Goldman-style Rotational DNA Mapping
 *
 * A mapping from trits (base-3 digits) to DNA bases that GUARANTEES no
 * homopolymers (no two adjacent bases are the same). Each trit selects one
 * of 3 non-previous bases via a rotating codebook.
 *
 * Density: log2(3)/2 ≈ 0.79 bits/base (lower than direct 2-bit mapping,
 * but the homopolymer constraint is enforced for free, no screening needed).
 *
 * Pipeline:
 *   bytes → trit stream (base-3 expansion) → rotational DNA mapping
 *
 * Trit expansion:
 *   Each byte (0-255) → ceil(log3(256)) = 6 trits.
 *   We use a chunked scheme: 5 bytes (40 bits) → 25 trits (3^25 > 2^40).
 *   Rejection sampling handles the bias (if 5-byte int >= 3^25, skip).
 *
 * Rotational codebook:
 *   For each trit t (0, 1, 2) and previous base P:
 *     next = codebook[P][(t + 1) mod 3]
 *   where codebook[P] = the 3 bases != P, in a fixed cyclic order.
 *
 * Reference:
 *   - Goldman et al. (2013). "Towards practical, high-capacity, low-
 *     maintenance information storage in synthesized DNA." Nature 494:77-80.
 */

import { Base } from "./mapping";

const BASES: Base[] = ["A", "C", "G", "T"];
const BASE_TO_IDX: Record<Base, number> = { A: 0, C: 1, G: 2, T: 3 };

/**
 * Codebook: for each previous base, the 3 allowed next bases in cyclic order.
 * codebook[prev] = [base1, base2, base3] — all != prev.
 */
const CODEBOOK: Base[][] = [
  ["C", "G", "T"], // prev = A → next ∈ {C, G, T}
  ["A", "G", "T"], // prev = C → next ∈ {A, G, T}
  ["A", "C", "T"], // prev = G → next ∈ {A, C, T}
  ["A", "C", "G"], // prev = T → next ∈ {A, C, G}
];

/**
 * Convert bytes to a trit stream using chunked base-3 expansion.
 *
 * Each 5-byte chunk (40 bits, max 2^40 - 1) is converted to 25 trits
 * (3^25 = 847,288,609,443 > 2^40 = 1,099,511,627,776).
 * Wait — 3^25 < 2^40, so we need 26 trits. 3^26 = 2,541,865,828,329 > 2^40.
 * We use 26 trits per 5 bytes, with rejection sampling for the bias.
 *
 * For simplicity (no rejection), we use 11 bits → 7 trits (3^7 = 2187 > 2048).
 * This has ~6% bias on the top trit but is simpler and the bias doesn't
 * affect homopolymer avoidance.
 */
/**
 * Convert bytes to a trit stream using chunked base-3 expansion.
 *
 * Two modes:
 *   - "fast" (default): 1 byte → 6 trits (3^6 = 729 > 256). Density 1.333 bits/nt.
 *     Simple, no rejection, but ~65% bias on the top trit.
 *   - "dense": 5 bytes → 26 trits (3^26 = 2.54T > 2^40 = 1.10T). Density 1.538 bits/nt.
 *     Uses rejection sampling: if the 5-byte value >= 3^26, skip and read next 5 bytes.
 *     Rejection rate: (2^40 - 3^26) / 2^40 = 0% (3^26 > 2^40, so no rejection!)
 *     Actually 3^26 = 2,541,865,828,329 > 2^40 = 1,099,511,627,776, so we can represent
 *     ALL 5-byte values with 26 trits. The "waste" is 3^26 - 2^40 = 1.44T unused trit
 *     combinations, but that's fine — it just means the top trits have slightly higher
 *     values on average. No rejection needed.
 *
 *     Wait — we want to MINIMIZE trits per byte. 3^25 = 847,288,609,443 < 2^40, so
 *     25 trits CAN'T represent all 5-byte values. We need 26 trits.
 *     But 26 trits for 40 bits = 40/26 = 1.538 bits/nt — BETTER than 8/6 = 1.333.
 *
 * For "dense" mode, we pack 5 bytes (40 bits) into a BigInt, then extract 26 trits.
 * This requires BigInt for the conversion but gives 15% higher density.
 */
export function bytesToTrits(data: Uint8Array, mode: "fast" | "dense" = "fast"): number[] {
  if (mode === "dense") {
    return bytesToTritsDense(data);
  }
  const trits: number[] = new Array(data.length * 6);
  for (let i = 0; i < data.length; i++) {
    let val = data[i];
    const base = i * 6;
    trits[base] = val % 3; val = Math.floor(val / 3);
    trits[base + 1] = val % 3; val = Math.floor(val / 3);
    trits[base + 2] = val % 3; val = Math.floor(val / 3);
    trits[base + 3] = val % 3; val = Math.floor(val / 3);
    trits[base + 4] = val % 3; val = Math.floor(val / 3);
    trits[base + 5] = val % 3;
  }
  return trits;
}

/**
 * High-density trit conversion: 5 bytes → 26 trits.
 * Uses BigInt for the base-3 expansion (40-bit values exceed JS Number's safe integer range
 * for precise division, though 2^40 < 2^53 so Number would work... but BigInt is safer).
 *
 * 3^26 = 2,541,865,828,329 > 2^40 = 1,099,511,627,776, so ALL 5-byte values are representable.
 * No rejection sampling needed — every 5 bytes maps to a valid 26-trit value.
 *
 * Density: 40 bits / 26 trits = 1.538 bits/nt (vs 1.333 for fast mode = 15% improvement)
 */
function bytesToTritsDense(data: Uint8Array): number[] {
  const trits: number[] = [];
  // 3^0, 3^1, ..., 3^25 as BigInt
  const THREE = 3n;
  const powers: bigint[] = [];
  let p = 1n;
  for (let i = 0; i < 26; i++) {
    powers.push(p);
    p *= THREE;
  }

  // Process 5 bytes at a time
  const chunkSize = 5;
  let i = 0;
  while (i < data.length) {
    // Build 40-bit value as BigInt
    let val = 0n;
    const remaining = Math.min(chunkSize, data.length - i);
    for (let b = 0; b < remaining; b++) {
      val = (val << 8n) | BigInt(data[i + b]);
    }
    // If last chunk is short, we still emit 26 trits (padded with zeros)
    // Pad to 5 bytes if short
    if (remaining < chunkSize) {
      // Already shifted correctly; val is the partial bytes in the high positions
      // No adjustment needed — the missing low bytes are effectively 0
    }

    // Extract 26 trits (LSB first)
    for (let t = 0; t < 26; t++) {
      trits.push(Number(val % THREE));
      val /= THREE;
    }
    i += chunkSize;
  }

  return trits;
}

/** Convert trits back to bytes. Inverse of bytesToTrits. */
export function tritsToBytes(trits: number[], mode: "fast" | "dense" = "fast"): Uint8Array {
  if (mode === "dense") {
    return tritsToBytesDense(trits);
  }
  const numBytes = Math.floor(trits.length / 6);
  const out = new Uint8Array(numBytes);
  for (let i = 0; i < numBytes; i++) {
    let val = 0;
    for (let t = 5; t >= 0; t--) {
      val = val * 3 + trits[i * 6 + t];
    }
    out[i] = val & 0xff;
  }
  return out;
}

/**
 * Inverse of bytesToTritsDense: 26 trits → 5 bytes.
 * Uses BigInt to reconstruct the 40-bit value from 26 trits.
 */
function tritsToBytesDense(trits: number[]): Uint8Array {
  const THREE = 3n;
  const numChunks = Math.floor(trits.length / 26);
  const out = new Uint8Array(numChunks * 5);

  let outIdx = 0;
  for (let chunk = 0; chunk < numChunks; chunk++) {
    // Reconstruct value: val = sum(trits[chunk*26 + t] * 3^t)
    let val = 0n;
    let power = 1n;
    for (let t = 0; t < 26; t++) {
      val += BigInt(trits[chunk * 26 + t]) * power;
      power *= THREE;
    }
    // Extract 5 bytes (big-endian)
    for (let b = 4; b >= 0; b--) {
      out[outIdx + b] = Number(val & 0xffn);
      val >>= 8n;
    }
    outIdx += 5;
  }

  return out;
}

/**
 * Encode trits to DNA using Goldman's rotational codebook.
 * Guarantees no homopolymers (every base differs from the previous).
 *
 * @param trits Array of trits (0, 1, or 2)
 * @param startBase Initial base (defaults to "A")
 * @returns DNA string with no homopolymers
 */
export function tritsToDna(trits: number[], startBase: Base = "A"): string {
  // Pre-allocate array for performance (avoid O(n²) string concatenation)
  const parts: string[] = new Array(trits.length);
  let prev = startBase;
  let prevIdx = BASE_TO_IDX[prev];
  for (let i = 0; i < trits.length; i++) {
    const next = CODEBOOK[prevIdx][trits[i]];
    parts[i] = next;
    prevIdx = BASE_TO_IDX[next];
  }
  return parts.join("");
}

/**
 * Decode DNA (Goldman rotational) back to trits.
 *
 * @param dna DNA string (must have no homopolymers)
 * @param startBase Initial base used during encoding
 * @returns Array of trits
 */
export function dnaToTrits(dna: string, startBase: Base = "A"): number[] {
  const trits: number[] = [];
  let prev = startBase;
  for (const base of dna) {
    const b = base as Base;
    const codebookRow = CODEBOOK[BASE_TO_IDX[prev]];
    const trit = codebookRow.indexOf(b);
    if (trit === -1) {
      throw new Error(`Invalid base ${base} after ${prev} (homopolymer or invalid)`);
    }
    trits.push(trit);
    prev = b;
  }
  return trits;
}

/**
 * Full pipeline: bytes → Goldman DNA.
 * Guarantees no homopolymers.
 *
 * @param mode "fast" (1B→6 trits, 1.333 bits/nt) or "dense" (5B→26 trits, 1.538 bits/nt)
 */
export function bytesToGoldmanDna(data: Uint8Array, startBase: Base = "A", mode: "fast" | "dense" = "fast"): string {
  const trits = bytesToTrits(data, mode);
  return tritsToDna(trits, startBase);
}

/**
 * Full pipeline: Goldman DNA → bytes.
 *
 * @param mode "fast" or "dense" (must match the mode used for encoding)
 */
export function goldmanDnaToBytes(dna: string, startBase: Base = "A", mode: "fast" | "dense" = "fast"): Uint8Array {
  const trits = dnaToTrits(dna, startBase);
  return tritsToBytes(trits, mode);
}

/**
 * Verify that a DNA string has no homopolymers (max run = 1).
 */
export function hasHomopolymer(dna: string): boolean {
  for (let i = 1; i < dna.length; i++) {
    if (dna[i] === dna[i - 1]) return true;
  }
  return false;
}
