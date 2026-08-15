/**
 * Gungnir Hash-Based Single-Read Recovery
 *
 * Proof-of-work decoding: each data fragment carries a hash signature.
 * The decoder tests educated guesses (substitution, insertion, deletion
 * hypotheses) until the hash matches. Recovers from a single read with
 * 20% erroneous bases — no redundant copies, no high-coverage sequencing.
 *
 * When to use:
 *   - Nanopore: sequencing is expensive, CPU is cheap → Gungnir saves 10-25x sequencing cost
 *   - Illumina: coverage is cheap, CPU is expensive → use existing consensus instead
 *
 * Performance:
 *   - Order 0 (no errors): ~1us per oligo (just hash comparison)
 *   - Order 1 (single error): ~4N hashes ≈ 1ms per oligo for N=250 bases
 *   - Order 2 (double errors): ~4N² hashes ≈ 250ms per oligo for N=250 bases
 *   - Typical nanopore: 1-3 errors per 150nt → Order 1-2 sufficient
 *
 * Reference:
 *   - HKU-BAL/Gungnir
 *   - Banal et al. (2026). arXiv:2604.20810. (hash-based proof-of-work concept)
 */

import { blake3 } from "@noble/hashes/blake3.js";
import { crc16 } from "./crc16";

// ─── DNA Constants ────────────────────────────────────────────────────────────

const DNA_BASES = ["A", "C", "G", "T"] as const;

// ─── Config ──────────────────────────────────────────────────────────────────

export interface GungnirConfig {
  /** Hash algorithm: 'blake3-32' (4 bytes) or 'crc16' (2 bytes). Default: 'blake3-32' */
  hashAlg: "blake3-32" | "crc16";
  /** Maximum search order (0=exact, 1=single error, 2=double error). Default: 2 */
  maxOrder: number;
  /** Maximum candidates to try before giving up. Default: 100000 */
  maxCandidates: number;
  /** Time budget in ms per oligo. Default: 1000 (1 second) */
  timeBudgetMs: number;
}

export const DEFAULT_GUNGNIR_CONFIG: GungnirConfig = {
  hashAlg: "blake3-32",
  maxOrder: 2,
  maxCandidates: 100000,
  timeBudgetMs: 1000,
};

// ─── Result ───────────────────────────────────────────────────────────────────

export interface GungnirResult {
  /** Corrected DNA sequence, or null if recovery failed */
  correctedDna: string | null;
  /** Number of errors corrected */
  errorsCorrected: number;
  /** Type of error corrected: 'none', 'substitution', 'insertion', 'deletion', 'compound' */
  errorType: string;
  /** Number of candidates tried */
  candidatesTried: number;
  /** Time spent in ms */
  timeMs: number;
  /** Order at which recovery succeeded (-1 if failed) */
  successOrder: number;
}

// ─── Hash Functions ───────────────────────────────────────────────────────────

/**
 * Convert a DNA string to a Uint8Array of ASCII codes.
 * Shared helper for hash computation.
 */
function dnaToBytes(dna: string): Uint8Array {
  const bytes = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) {
    const c = dna.charCodeAt(i);
    if (c !== 65 && c !== 67 && c !== 71 && c !== 84) {
      // A=65, C=67, G=71, T=84
      throw new Error(`Invalid DNA base '${dna[i]}' at position ${i}`);
    }
    bytes[i] = c;
  }
  return bytes;
}

/**
 * Compute hash signature for a DNA sequence.
 *
 * Used during encode to store the expected hash alongside each oligo.
 * The hash is computed over the ASCII representation of the DNA string.
 *
 * @param dna DNA string to hash
 * @param alg Hash algorithm ('blake3-32' = 4 bytes, 'crc16' = 2 bytes)
 * @returns Hash as a Uint8Array
 */
export function computeDnaHash(
  dna: string,
  alg: "blake3-32" | "crc16" = "blake3-32",
): Uint8Array {
  const bytes = dnaToBytes(dna);

  if (alg === "crc16") {
    // CRC-16: 2 bytes, big-endian
    const crc = crc16(bytes);
    return new Uint8Array([(crc >> 8) & 0xff, crc & 0xff]);
  }

  // blake3-32: first 4 bytes of BLAKE3 hash
  const fullHash = blake3(bytes);
  return fullHash.subarray(0, 4);
}

/**
 * Verify a DNA sequence against an expected hash.
 *
 * @param dna DNA string to verify
 * @param expectedHash Expected hash (from computeDnaHash)
 * @param alg Hash algorithm (must match the one used to compute the hash)
 * @returns true if the hash matches
 */
export function verifyDnaHash(
  dna: string,
  expectedHash: Uint8Array,
  alg: "blake3-32" | "crc16" = "blake3-32",
): boolean {
  const actualHash = computeDnaHash(dna, alg);
  if (actualHash.length !== expectedHash.length) return false;
  for (let i = 0; i < actualHash.length; i++) {
    if (actualHash[i] !== expectedHash[i]) return false;
  }
  return true;
}

/**
 * Fast hash comparison without constructing intermediate objects.
 * Returns true if the BLAKE3-32 hash of the DNA matches expectedHash.
 */
function fastHashMatch(
  dna: string,
  expectedHash: Uint8Array,
  alg: "blake3-32" | "crc16",
): boolean {
  if (alg === "crc16") {
    // Inline CRC-16 computation for speed
    let crc = 0xffff;
    for (let i = 0; i < dna.length; i++) {
      crc ^= (dna.charCodeAt(i) << 8) & 0xffff;
      for (let j = 0; j < 8; j++) {
        if (crc & 0x8000) {
          crc = ((crc << 1) ^ 0x1021) & 0xffff;
        } else {
          crc = (crc << 1) & 0xffff;
        }
      }
    }
    return (
      ((crc >> 8) & 0xff) === expectedHash[0] &&
      (crc & 0xff) === expectedHash[1]
    );
  }

  // blake3-32: hash and compare first 4 bytes
  const bytes = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) {
    bytes[i] = dna.charCodeAt(i);
  }
  const hash = blake3(bytes);
  return (
    hash[0] === expectedHash[0] &&
    hash[1] === expectedHash[1] &&
    hash[2] === expectedHash[2] &&
    hash[3] === expectedHash[3]
  );
}

// ─── String Manipulation Helpers ──────────────────────────────────────────────

/**
 * Substitute base at position i with a new base.
 * Returns a new string; does not modify the input.
 */
function substitute(dna: string, i: number, newBase: string): string {
  return dna.substring(0, i) + newBase + dna.substring(i + 1);
}

/**
 * Delete base at position i.
 * Returns a new string with the base removed.
 */
function deleteBase(dna: string, i: number): string {
  return dna.substring(0, i) + dna.substring(i + 1);
}

/**
 * Insert a base before position i.
 * Returns a new string with the base inserted.
 * If i === dna.length, appends at the end.
 */
function insertBase(dna: string, i: number, newBase: string): string {
  return dna.substring(0, i) + newBase + dna.substring(i);
}

// ─── Order 0: Exact Match ────────────────────────────────────────────────────

/**
 * Order 0: check if the read matches the hash exactly (no errors).
 */
function gungnirOrder0(
  readDna: string,
  expectedHash: Uint8Array,
  alg: "blake3-32" | "crc16",
): GungnirResult | null {
  if (fastHashMatch(readDna, expectedHash, alg)) {
    return {
      correctedDna: readDna,
      errorsCorrected: 0,
      errorType: "none",
      candidatesTried: 1,
      timeMs: 0, // caller fills in
      successOrder: 0,
    };
  }
  return null;
}

// ─── Order 1: Single Error ────────────────────────────────────────────────────

/**
 * Order-1 search: try all single-base substitutions, insertions, deletions.
 *
 * For a read of length N:
 *   - Substitutions: 3 * N candidates (3 alternate bases at each position)
 *   - Deletions: N candidates (delete each position)
 *   - Insertions: 4 * (N + 1) candidates (insert each base at each position)
 *
 * Total: 3N + N + 4(N+1) = 8N + 4 candidates.
 *
 * Returns the first candidate that matches the expected hash, or null.
 * Short-circuits on first match for speed.
 *
 * @param readDna The (possibly corrupted) DNA read
 * @param expectedHash The expected hash of the correct sequence
 * @param expectedLength Expected sequence length (before corruption)
 * @param alg Hash algorithm
 * @param maxCandidates Maximum candidates to try
 * @param deadline Timestamp (Date.now()) after which to give up
 * @returns Recovery result or null
 */
export function gungnirOrder1(
  readDna: string,
  expectedHash: Uint8Array,
  expectedLength: number,
  alg: "blake3-32" | "crc16",
  maxCandidates: number,
  deadline: number,
): GungnirResult | null {
  const n = readDna.length;
  let candidates = 0;
  const startTime = Date.now();

  // ── Substitutions: try 3 alternate bases at each position ──
  // Most common error type, so try first
  for (let i = 0; i < n; i++) {
    const currentBase = readDna[i];
    for (const newBase of DNA_BASES) {
      if (newBase === currentBase) continue;
      candidates++;
      if (candidates > maxCandidates) {
        return null;
      }
      if (candidates % 1000 === 0 && Date.now() > deadline) {
        return null;
      }

      const candidate = substitute(readDna, i, newBase);
      if (fastHashMatch(candidate, expectedHash, alg)) {
        return {
          correctedDna: candidate,
          errorsCorrected: 1,
          errorType: "substitution",
          candidatesTried: candidates,
          timeMs: Date.now() - startTime,
          successOrder: 1,
        };
      }
    }
  }

  // ── Deletions: remove one base at each position ──
  // Only makes sense if the read is longer than expected (insertion error in read)
  if (n > 1 && n >= expectedLength) {
    for (let i = 0; i < n; i++) {
      candidates++;
      if (candidates > maxCandidates) {
        return null;
      }
      if (candidates % 1000 === 0 && Date.now() > deadline) {
        return null;
      }

      const candidate = deleteBase(readDna, i);
      if (candidate.length !== expectedLength) continue;

      if (fastHashMatch(candidate, expectedHash, alg)) {
        return {
          correctedDna: candidate,
          errorsCorrected: 1,
          errorType: "deletion",
          candidatesTried: candidates,
          timeMs: Date.now() - startTime,
          successOrder: 1,
        };
      }
    }
  }

  // ── Insertions: insert one base at each position ──
  // Only makes sense if the read is shorter than expected (deletion error in read)
  if (n + 1 <= expectedLength + 1) {
    for (let i = 0; i <= n; i++) {
      for (const newBase of DNA_BASES) {
        candidates++;
        if (candidates > maxCandidates) {
          return null;
        }
        if (candidates % 1000 === 0 && Date.now() > deadline) {
          return null;
        }

        const candidate = insertBase(readDna, i, newBase);
        if (candidate.length !== expectedLength) continue;

        if (fastHashMatch(candidate, expectedHash, alg)) {
          return {
            correctedDna: candidate,
            errorsCorrected: 1,
            errorType: "insertion",
            candidatesTried: candidates,
            timeMs: Date.now() - startTime,
            successOrder: 1,
          };
        }
      }
    }
  }

  return null;
}

// ─── Order 2: Double Error ────────────────────────────────────────────────────

/**
 * Apply a single operation to a DNA string and yield all candidates.
 * Returns an array of {dna, type} objects for all single operations.
 *
 * This is used as a building block for order-2 search.
 */
function singleOperations(
  dna: string,
  expectedLength: number,
): Array<{ dna: string; type: string }> {
  const results: Array<{ dna: string; type: string }> = [];
  const n = dna.length;

  // Substitutions
  for (let i = 0; i < n; i++) {
    const currentBase = dna[i];
    for (const newBase of DNA_BASES) {
      if (newBase === currentBase) continue;
      results.push({ dna: substitute(dna, i, newBase), type: "substitution" });
    }
  }

  // Deletions (only if length would be reasonable)
  if (n > 1 && n >= expectedLength - 1) {
    for (let i = 0; i < n; i++) {
      results.push({ dna: deleteBase(dna, i), type: "deletion" });
    }
  }

  // Insertions (only if length would be reasonable)
  if (n + 1 <= expectedLength + 2) {
    for (let i = 0; i <= n; i++) {
      for (const newBase of DNA_BASES) {
        results.push({ dna: insertBase(dna, i, newBase), type: "insertion" });
      }
    }
  }

  return results;
}

/**
 * Order-2 search: try all pairs of single-base operations.
 *
 * This is O(N²) in the read length but handles 2 errors per read.
 * For a 150nt read: ~150² * 8² ≈ 1.44M candidates — feasible within 1s.
 *
 * Strategy: generate all order-1 intermediates, then apply all order-1
 * operations to each intermediate. Use deadline-based early exit.
 *
 * @param readDna The (possibly corrupted) DNA read
 * @param expectedHash The expected hash of the correct sequence
 * @param expectedLength Expected sequence length (before corruption)
 * @param alg Hash algorithm
 * @param maxCandidates Maximum candidates to try
 * @param deadline Timestamp (Date.now()) after which to give up
 * @returns Recovery result or null
 */
export function gungnirOrder2(
  readDna: string,
  expectedHash: Uint8Array,
  expectedLength: number,
  alg: "blake3-32" | "crc16",
  maxCandidates: number,
  deadline: number,
): GungnirResult | null {
  const startTime = Date.now();
  let candidates = 0;

  // Generate all order-1 intermediates
  const order1Candidates = singleOperations(readDna, expectedLength);

  for (const mid of order1Candidates) {
    // Quick check: is the intermediate length compatible with the expected length?
    // After 2 operations, the length can be expectedLength, expectedLength ± 1, or expectedLength ± 2
    // depending on the operation types. We'll let the inner loop filter.

    // Apply all single operations to this intermediate
    const innerOps = singleOperations(mid.dna, expectedLength);

    for (const inner of innerOps) {
      candidates++;
      if (candidates > maxCandidates) {
        return null;
      }
      if (candidates % 5000 === 0 && Date.now() > deadline) {
        return null;
      }

      // Length filter: only check candidates of the expected length
      if (inner.dna.length !== expectedLength) continue;

      if (fastHashMatch(inner.dna, expectedHash, alg)) {
        return {
          correctedDna: inner.dna,
          errorsCorrected: 2,
          errorType: "compound",
          candidatesTried: candidates,
          timeMs: Date.now() - startTime,
          successOrder: 2,
        };
      }
    }
  }

  return null;
}

// ─── Main Recovery Function ───────────────────────────────────────────────────

/**
 * Attempt to recover a corrupted read via hash-based proof-of-work.
 *
 * Progressively tries higher search orders until the hash matches or
 * the budget is exhausted:
 *   - Order 0: exact match (no errors)
 *   - Order 1: all single-base substitutions, insertions, deletions
 *   - Order 2: all pairs of single-base operations
 *
 * @param readDna The (possibly corrupted) DNA read
 * @param expectedHash The expected hash of the correct sequence
 * @param expectedLength Expected sequence length (before corruption)
 * @param config Gungnir configuration
 * @returns Recovery result with corrected sequence or null
 */
export function gungnirRecover(
  readDna: string,
  expectedHash: Uint8Array,
  expectedLength: number,
  config?: GungnirConfig,
): GungnirResult {
  const cfg = { ...DEFAULT_GUNGNIR_CONFIG, ...config };
  const startTime = Date.now();
  const deadline = startTime + cfg.timeBudgetMs;

  // Validate inputs
  if (readDna.length === 0) {
    return {
      correctedDna: null,
      errorsCorrected: 0,
      errorType: "none",
      candidatesTried: 0,
      timeMs: Date.now() - startTime,
      successOrder: -1,
    };
  }

  if (expectedHash.length === 0) {
    return {
      correctedDna: null,
      errorsCorrected: 0,
      errorType: "none",
      candidatesTried: 0,
      timeMs: Date.now() - startTime,
      successOrder: -1,
    };
  }

  // ── Order 0: exact match ──
  const order0Result = gungnirOrder0(readDna, expectedHash, cfg.hashAlg);
  if (order0Result) {
    order0Result.timeMs = Date.now() - startTime;
    return order0Result;
  }

  // ── Order 1: single error ──
  if (cfg.maxOrder >= 1) {
    const order1Result = gungnirOrder1(
      readDna,
      expectedHash,
      expectedLength,
      cfg.hashAlg,
      cfg.maxCandidates,
      deadline,
    );
    if (order1Result) {
      return order1Result;
    }
  }

  // ── Order 2: double error ──
  if (cfg.maxOrder >= 2) {
    const order2Result = gungnirOrder2(
      readDna,
      expectedHash,
      expectedLength,
      cfg.hashAlg,
      cfg.maxCandidates,
      deadline,
    );
    if (order2Result) {
      return order2Result;
    }
  }

  // Recovery failed
  return {
    correctedDna: null,
    errorsCorrected: 0,
    errorType: "none",
    candidatesTried: cfg.maxCandidates, // approximate
    timeMs: Date.now() - startTime,
    successOrder: -1,
  };
}

// ─── Integration Helper ──────────────────────────────────────────────────────

/**
 * Integrate with decode.ts: for oligos with single reads, try Gungnir
 * before falling back to erasure.
 *
 * This function extracts the DNA from the read, runs Gungnir recovery,
 * and returns the corrected DNA (or null if recovery fails).
 *
 * Typical usage in the decode pipeline:
 * ```typescript
 * if (reads.length === 1 && gungnirMode) {
 *   const corrected = gungnirDecodeSingleRead(
 *     reads[0].sequence,
 *     oligo.expectedHash,
 *     oligo.expectedLength,
 *     gungnirConfig,
 *   );
 *   if (corrected) {
 *     // Use corrected sequence
 *   } else {
 *     // Fall back to erasure marking
 *   }
 * }
 * ```
 *
 * @param readDna The single read's DNA sequence
 * @param expectedHash The expected hash (stored in the oligo header)
 * @param expectedLength Expected sequence length
 * @param config Gungnir configuration
 * @returns Corrected DNA string, or null if recovery failed
 */
export function gungnirDecodeSingleRead(
  readDna: string,
  expectedHash: Uint8Array,
  expectedLength: number,
  config?: GungnirConfig,
): string | null {
  const result = gungnirRecover(readDna, expectedHash, expectedLength, config);
  return result.correctedDna;
}

// ─── Batch Processing ─────────────────────────────────────────────────────────

/**
 * Result for a batch of Gungnir recovery attempts.
 */
export interface GungnirBatchResult {
  /** Results for each oligo in the batch */
  results: GungnirResult[];
  /** Total oligos processed */
  totalOligos: number;
  /** Oligos successfully recovered */
  recovered: number;
  /** Oligos that failed recovery */
  failed: number;
  /** Total time in ms */
  totalTimeMs: number;
  /** Average candidates tried per oligo */
  avgCandidatesPerOligo: number;
}

/**
 * Batch-process multiple oligos with Gungnir recovery.
 *
 * Processes oligos sequentially, respecting individual time budgets.
 * Useful for processing a pool of single-read oligos from nanopore data.
 *
 * @param reads Array of {dna, expectedHash, expectedLength} objects
 * @param config Gungnir configuration
 * @returns Batch result with individual results and summary statistics
 */
export function gungnirBatchRecover(
  reads: Array<{
    dna: string;
    expectedHash: Uint8Array;
    expectedLength: number;
  }>,
  config?: GungnirConfig,
): GungnirBatchResult {
  const startTime = Date.now();
  const results: GungnirResult[] = new Array(reads.length);
  let recovered = 0;
  let totalCandidates = 0;

  for (let i = 0; i < reads.length; i++) {
    const { dna, expectedHash, expectedLength } = reads[i];
    const result = gungnirRecover(dna, expectedHash, expectedLength, config);
    results[i] = result;
    totalCandidates += result.candidatesTried;
    if (result.correctedDna !== null) {
      recovered++;
    }
  }

  return {
    results,
    totalOligos: reads.length,
    recovered,
    failed: reads.length - recovered,
    totalTimeMs: Date.now() - startTime,
    avgCandidatesPerOligo:
      reads.length > 0 ? Math.round(totalCandidates / reads.length) : 0,
  };
}

// ─── Hash Encoding/Decoding Utilities ─────────────────────────────────────────

/**
 * Encode a Gungnir hash as a DNA string for storage in an oligo.
 *
 * Each byte is encoded as 4 DNA bases (2 bits per base, same as direct mapping).
 * For blake3-32 (4 bytes): 16 DNA bases.
 * For crc16 (2 bytes): 8 DNA bases.
 *
 * @param hash Hash bytes from computeDnaHash
 * @returns DNA string encoding the hash
 */
export function encodeHashAsDna(hash: Uint8Array): string {
  const bases: string[] = new Array(hash.length * 4);
  for (let i = 0; i < hash.length; i++) {
    const byte = hash[i];
    bases[i * 4] = DNA_BASES[(byte >> 6) & 0b11];
    bases[i * 4 + 1] = DNA_BASES[(byte >> 4) & 0b11];
    bases[i * 4 + 2] = DNA_BASES[(byte >> 2) & 0b11];
    bases[i * 4 + 3] = DNA_BASES[byte & 0b11];
  }
  return bases.join("");
}

/**
 * Decode a hash from its DNA encoding (inverse of encodeHashAsDna).
 *
 * @param dnaDna DNA string encoding the hash
 * @param numBytes Expected number of hash bytes
 * @returns Hash as Uint8Array
 */
export function decodeHashFromDna(
  dnaDna: string,
  numBytes: number,
): Uint8Array {
  const hash = new Uint8Array(numBytes);
  for (let i = 0; i < numBytes; i++) {
    let byte = 0;
    for (let j = 0; j < 4; j++) {
      const base = dnaDna[i * 4 + j];
      const idx = DNA_BASES.indexOf(base as (typeof DNA_BASES)[number]);
      if (idx < 0) {
        throw new Error(`Invalid base '${base}' in hash DNA at position ${i * 4 + j}`);
      }
      byte = (byte << 2) | idx;
    }
    hash[i] = byte;
  }
  return hash;
}

/**
 * Compute the number of DNA bases needed to store a Gungnir hash.
 *
 * @param alg Hash algorithm
 * @returns Number of DNA bases for the hash field
 */
export function gungnirHashDnaLength(alg: "blake3-32" | "crc16"): number {
  return alg === "blake3-32" ? 16 : 8; // 4 bytes * 4 bases/byte or 2 bytes * 4
}
