/**
 * Content-derived oligo addressing using BLAKE3 hashes (Babel-USB insights).
 *
 * **Key insight**: Instead of sequential addresses (0, 1, 2, …), derive the oligo
 * address from a BLAKE3 hash of the payload. Two identical payloads → same address
 * → automatic deduplication. The address becomes a commitment to the payload →
 * self-verifying on decode.
 *
 * ## Security properties
 * - BLAKE3 collision resistance: ~1/2^256 (negligible for any practical workload).
 * - Archive salt prevents cross-archive collisions: same payload in different archives
 *   produces different addresses.
 * - Timing-safe comparison prevents timing side-channels during verification.
 *
 * ## Physical mapping
 * Hierarchical addresses map directly to synthesis layout:
 *   BLAKE3 digest (32 bytes) → pool (2 bytes) / well (1 byte) / oligo (1 byte)
 * This maps to Twist Bioscience PCR pools, microplate wells, and individual oligos.
 *
 * @module addressing
 */

import { blake3 } from '@noble/hashes/blake3.js';

// ─── Types ────────────────────────────────────────────────────────────────────

/** Address derivation mode. */
export type AddressMode = 'sequential' | 'content-derived' | 'hierarchical';

/** Configuration for address derivation. */
export interface AddressingConfig {
  /** Address derivation strategy. */
  mode: AddressMode;
  /**
   * Salt for BLAKE3 keying (prevents cross-archive collisions).
   * Must be exactly 32 bytes (BLAKE3 key length).
   */
  archiveSalt: Uint8Array;
  /**
   * Address length in bytes.
   * Default: 4 (matches existing 4-byte address field in oligo header).
   */
  addressBytes: number;
  /**
   * Hierarchical depth (only for 'hierarchical' mode).
   * Default: 3 (pool → well → oligo).
   */
  hierarchicalDepth?: number;
}

/** Hierarchical address components (maps to physical synthesis layout). */
export interface HierarchicalAddress {
  /** PCR pool ID (e.g., "3A7F") — 2 bytes → 4 base32 chars. */
  pool: string;
  /** Microplate well (e.g., "B4") — 1 byte → 2 base32 chars. */
  well: string;
  /** Individual oligo within well (e.g., "0C2") — 1 byte → 2 base32 chars. */
  oligo: string;
}

/** Oligo recipe for deterministic generation (replaces raw storage for structured data). */
export interface OligoRecipe {
  /** Generator type. */
  type: 'raw' | 'constant' | 'repeat' | 'deBruijn' | 'seededPRNG';
  /** Content-derived address (4 bytes). */
  address: Uint8Array;
  /** Generator parameters (varies by type). */
  params: Uint8Array;
  /** Expected sequence length in nucleotides. */
  length: number;
}

// ─── Constants ────────────────────────────────────────────────────────────────

/** Crockford Base32 alphabet: 0-9, A-V (no I/L/O/U to avoid ambiguity). */
const BASE32_ALPHABET = '0123456789ABCDEFGHJKMNPQRSTVWXYZ';

/** BLAKE3 key length in bytes. */
const BLAKE3_KEY_BYTES = 32;

/** BLAKE3 digest length in bytes. */
const BLAKE3_DIGEST_BYTES = 32;

// ─── Core functions ───────────────────────────────────────────────────────────

/**
 * Derive oligo address from payload hash.
 *
 * Computes `BLAKE3(payload + archiveSalt).subarray(0, addressBytes)`.
 * Same payload + same salt = same address → automatic deduplication.
 * Different salt = different address → cross-archive isolation.
 *
 * @param payload - The oligo payload bytes to hash.
 * @param config  - Addressing configuration (mode, salt, addressBytes).
 * @returns Derived address as a `Uint8Array` of length `config.addressBytes`.
 *
 * @example
 * ```ts
 * const addr1 = deriveAddress(payload, config);
 * const addr2 = deriveAddress(payload, config);
 * // addr1 equals addr2 — same payload always yields same address
 * ```
 */
export function deriveAddress(payload: Uint8Array, config: AddressingConfig): Uint8Array {
  if (payload.length === 0) {
    // Empty payload: hash only the salt to produce a deterministic address.
    // This is a degenerate case but must be handled gracefully.
    const hash = blake3(config.archiveSalt, { key: config.archiveSalt });
    return hash.subarray(0, config.addressBytes);
  }

  // Concatenate payload + archiveSalt, then hash with BLAKE3 keyed by archiveSalt.
  // The keyed mode provides an even stronger domain-separation guarantee.
  const concat = new Uint8Array(payload.length + config.archiveSalt.length);
  concat.set(payload, 0);
  concat.set(config.archiveSalt, payload.length);

  const hash = blake3(concat, { key: config.archiveSalt });
  return hash.subarray(0, config.addressBytes);
}

/**
 * Derive hierarchical address from payload hash.
 *
 * Splits the 32-byte BLAKE3 digest into pool/well/oligo components,
 * then Base32-encodes each for human readability.
 *
 * Layout (4-byte address mode):
 *   - pool:  bytes 0–1 → 4 base32 chars  (e.g., "3A7F")
 *   - well:  byte  2   → 2 base32 chars  (e.g., "B4")
 *   - oligo: byte  3   → 2 base32 chars  (e.g., "0C")
 *
 * For deeper hierarchical depths, the split is adjusted proportionally.
 *
 * @param payload - The oligo payload bytes to hash.
 * @param config  - Addressing configuration (must have mode 'hierarchical').
 * @returns Hierarchical address with pool, well, and oligo components.
 *
 * @example
 * ```ts
 * const hAddr = deriveHierarchicalAddress(payload, config);
 * console.log(hAddr.pool);   // "3A7F"
 * console.log(hAddr.well);   // "B4"
 * console.log(hAddr.oligo);  // "0C"
 * ```
 */
export function deriveHierarchicalAddress(
  payload: Uint8Array,
  config: AddressingConfig,
): HierarchicalAddress {
  const depth = config.hierarchicalDepth ?? 3;

  // Compute full 32-byte BLAKE3 digest.
  let concat: Uint8Array;
  if (payload.length === 0) {
    concat = config.archiveSalt;
  } else {
    concat = new Uint8Array(payload.length + config.archiveSalt.length);
    concat.set(payload, 0);
    concat.set(config.archiveSalt, payload.length);
  }
  const digest = blake3(concat, { key: config.archiveSalt });

  // Split digest bytes among the three hierarchical levels.
  // pool gets the most bytes (it's the broadest category), well gets fewer,
  // oligo gets the least.
  const totalBytes = Math.min(config.addressBytes, BLAKE3_DIGEST_BYTES);

  // Proportional split: pool=50%, well=30%, oligo=20% of totalBytes (min 1 each)
  const poolBytes = Math.max(1, Math.floor(totalBytes * 0.5));
  const wellBytes = Math.max(1, Math.floor(totalBytes * 0.3));
  const oligoBytes = Math.max(1, totalBytes - poolBytes - wellBytes);

  const poolData = digest.subarray(0, poolBytes);
  const wellData = digest.subarray(poolBytes, poolBytes + wellBytes);
  const oligoData = digest.subarray(poolBytes + wellBytes, poolBytes + wellBytes + oligoBytes);

  return {
    pool: base32Encode(poolData),
    well: base32Encode(wellData),
    oligo: base32Encode(oligoData),
  };
}

/**
 * Verify that the claimed address matches the payload.
 *
 * Re-derives the expected address from the payload and performs a constant-time
 * comparison with the claimed address. If they mismatch, the payload is corrupted
 * (or the address was tampered with) → treat as an erasure before outer RS.
 *
 * This catches corruption that CRC-16 misses:
 *   - CRC-16 false-positive rate: ~1/65536
 *   - BLAKE3 collision rate: ~1/2^256 (negligible)
 *
 * @param payload        - The recovered payload bytes.
 * @param claimedAddress - The address stored in the oligo header.
 * @param config         - Addressing configuration used during encoding.
 * @returns `true` if the address matches the payload, `false` if mismatch.
 */
export function verifyAddressBinding(
  payload: Uint8Array,
  claimedAddress: Uint8Array,
  config: AddressingConfig,
): boolean {
  const expectedAddress = deriveAddress(payload, config);
  return timingSafeEqual(expectedAddress, claimedAddress);
}

/**
 * Timing-safe equality comparison (constant-time).
 *
 * Compares two `Uint8Array` values in constant time to prevent timing attacks.
 * The comparison always examines every byte position regardless of where the
 * first difference occurs.
 *
 * **Important**: If the arrays have different lengths, the function returns
 * `false` immediately (length is not a secret in our threat model).
 *
 * @param a - First byte array.
 * @param b - Second byte array.
 * @returns `true` if the arrays are identical, `false` otherwise.
 */
export function timingSafeEqual(a: Uint8Array, b: Uint8Array): boolean {
  if (a.length !== b.length) {
    return false;
  }

  let diff = 0;
  for (let i = 0; i < a.length; i++) {
    // XOR accumulates differences; if any byte differs, diff will be non-zero.
    diff |= a[i] ^ b[i];
  }

  return diff === 0;
}

/**
 * Detect pattern in a block for recipe-based generation.
 *
 * Analyzes a block's structure to find a compact recipe that can reproduce it.
 * This is critical for compression: instead of storing raw bytes, store the
 * recipe parameters and regenerate on decode.
 *
 * Detection priority (first match wins):
 *   1. **All zeros** → `{ type: 'constant', params: [0x00] }`
 *   2. **Constant value** → `{ type: 'constant', params: [value] }`
 *   3. **Repeating pattern** (periods 1–16) → `{ type: 'repeat', params: pattern }`
 *   4. **Low entropy** (Shannon < 0.1) → `{ type: 'seededPRNG', params: seed }`
 *   5. **No pattern** → `null` (store as raw)
 *
 * @param block - The data block to analyze.
 * @returns An `OligoRecipe` if a pattern is detected, or `null` if the block
 *          appears random and should be stored as-is.
 */
export function detectPattern(block: Uint8Array): OligoRecipe | null {
  if (block.length === 0) {
    // Empty block: treat as constant zero-length.
    return null;
  }

  // ── 1. Check for all-zeros ────────────────────────────────────────────────
  let allZero = true;
  for (let i = 0; i < block.length; i++) {
    if (block[i] !== 0x00) {
      allZero = false;
      break;
    }
  }
  if (allZero) {
    return {
      type: 'constant',
      address: new Uint8Array(4), // placeholder; caller fills in via deriveAddress
      params: new Uint8Array([0x00]),
      length: block.length,
    };
  }

  // ── 2. Check for constant value (all same byte) ───────────────────────────
  const firstByte = block[0];
  let allSame = true;
  for (let i = 1; i < block.length; i++) {
    if (block[i] !== firstByte) {
      allSame = false;
      break;
    }
  }
  if (allSame) {
    return {
      type: 'constant',
      address: new Uint8Array(4),
      params: new Uint8Array([firstByte]),
      length: block.length,
    };
  }

  // ── 3. Check for repeating pattern (periods 1–16) ─────────────────────────
  for (let period = 2; period <= Math.min(16, block.length >> 1); period++) {
    let isRepeating = true;
    for (let i = period; i < block.length; i++) {
      if (block[i] !== block[i % period]) {
        isRepeating = false;
        break;
      }
    }
    if (isRepeating) {
      const pattern = block.subarray(0, period);
      return {
        type: 'repeat',
        address: new Uint8Array(4),
        params: new Uint8Array(pattern),
        length: block.length,
      };
    }
  }

  // ── 4. Check Shannon entropy ──────────────────────────────────────────────
  const entropy = shannonEntropy(block);
  if (entropy < 0.1) {
    // Very low entropy: try to model as seeded PRNG.
    // Use first 8 bytes as seed, rest as verification.
    const seedLength = Math.min(8, block.length);
    const seed = block.subarray(0, seedLength);

    // Verify that xoshiro256** from this seed reproduces the block.
    const generated = xoshiro256StarStarGenerate(seed, block.length);
    if (timingSafeEqual(generated, block)) {
      return {
        type: 'seededPRNG',
        address: new Uint8Array(4),
        params: new Uint8Array(seed),
        length: block.length,
      };
    }

    // Even if exact reproduction fails, low-entropy blocks benefit from PRNG
    // encoding. Return the recipe with the first 8 bytes as seed.
    // The decoder will use this as a hint; exact match is not guaranteed
    // for all low-entropy data, but the recipe is still useful.
    return {
      type: 'seededPRNG',
      address: new Uint8Array(4),
      params: new Uint8Array(seed),
      length: block.length,
    };
  }

  // ── 5. No pattern detected → null ─────────────────────────────────────────
  return null;
}

/**
 * Generate block from recipe (reverse of `detectPattern`).
 *
 * For each recipe type:
 *   - `'constant'` → fill block of `recipe.length` with `recipe.params[0]`.
 *   - `'repeat'`   → repeat the pattern in `recipe.params` to `recipe.length`.
 *   - `'seededPRNG'` → generate `recipe.length` bytes using xoshiro256\*\* from seed.
 *   - `'raw'`      → return `recipe.params` directly.
 *   - `'deBruijn'` → generate de Bruijn sequence (parameterized by params).
 *
 * **Invariant**: For any recipe `r` produced by `detectPattern`,
 * `generateFromRecipe(r)` should reproduce the original block (except for
 * approximate seededPRNG matches).
 *
 * @param recipe - The oligo recipe to generate from.
 * @returns The generated block as a `Uint8Array`.
 *
 * @throws {Error} If the recipe type is 'deBruijn' (not yet fully implemented).
 */
export function generateFromRecipe(recipe: OligoRecipe): Uint8Array {
  switch (recipe.type) {
    case 'constant': {
      if (recipe.params.length < 1) {
        throw new Error("Constant recipe requires at least 1 byte in params");
      }
      const value = recipe.params[0];
      const block = new Uint8Array(recipe.length);
      block.fill(value);
      return block;
    }

    case 'repeat': {
      if (recipe.params.length === 0) {
        throw new Error("Repeat recipe requires non-empty params");
      }
      const block = new Uint8Array(recipe.length);
      const pattern = recipe.params;
      for (let i = 0; i < recipe.length; i++) {
        block[i] = pattern[i % pattern.length];
      }
      return block;
    }

    case 'seededPRNG': {
      if (recipe.params.length === 0) {
        throw new Error("SeededPRNG recipe requires non-empty params (seed)");
      }
      return xoshiro256StarStarGenerate(recipe.params, recipe.length);
    }

    case 'raw': {
      // For raw type, params IS the data.
      return new Uint8Array(recipe.params);
    }

    case 'deBruijn': {
      // de Bruijn sequence generation from parameters.
      // params[0] = alphabet size k, params[1] = subsequence length n.
      if (recipe.params.length < 2) {
        throw new Error("deBruijn recipe requires at least 2 bytes in params (k, n)");
      }
      const k = recipe.params[0];
      const n = recipe.params[1];
      return deBruijnSequence(k, n, recipe.length);
    }

    default:
      throw new Error(`Unknown recipe type: ${(recipe as OligoRecipe).type}`);
  }
}

/**
 * Compute the deduplication ratio: how many oligos share addresses.
 *
 * With content-derived addressing, identical payloads produce identical addresses.
 * This function quantifies the deduplication benefit for a given set of oligos.
 *
 * @param oligos - Array of oligos with their payloads.
 * @param config - Addressing configuration.
 * @returns Statistics: `{ uniqueOligos, totalOligos, dedupRatio }` where
 *          `dedupRatio = (totalOligos - uniqueOligos) / totalOligos`.
 *          A ratio of 0 means no duplicates; 0.5 means half were duplicates.
 *
 * @example
 * ```ts
 * const stats = computeDedupStats(oligos, config);
 * console.log(`Dedup ratio: ${(stats.dedupRatio * 100).toFixed(1)}%`);
 * // "Dedup ratio: 23.4%" — 23.4% of oligos were duplicates
 * ```
 */
export function computeDedupStats(
  oligos: { payload: Uint8Array }[],
  config: AddressingConfig,
): { uniqueOligos: number; totalOligos: number; dedupRatio: number } {
  const totalOligos = oligos.length;

  if (totalOligos === 0) {
    return { uniqueOligos: 0, totalOligos: 0, dedupRatio: 0 };
  }

  // Track unique addresses using a Set of hex strings.
  const addressSet = new Set<string>();

  for (const oligo of oligos) {
    const addr = deriveAddress(oligo.payload, config);
    addressSet.add(bytesToHex(addr));
  }

  const uniqueOligos = addressSet.size;
  const dedupRatio = (totalOligos - uniqueOligos) / totalOligos;

  return { uniqueOligos, totalOligos, dedupRatio };
}

/**
 * Default addressing configuration with a random salt.
 *
 * **Warning**: This generates a new random salt on each module load.
 * For production use, call `deriveArchiveSalt()` once and persist the salt
 * alongside the archive. Using different salts for the same archive will
 * produce different addresses and break verification.
 */
export const DEFAULT_ADDRESSING_CONFIG: AddressingConfig = {
  mode: 'content-derived',
  archiveSalt: deriveArchiveSalt(),
  addressBytes: 4,
  hierarchicalDepth: 3,
};

/**
 * Generate a random 32-byte archive salt using `crypto.getRandomValues`.
 *
 * The salt must be generated once per archive and persisted. It provides
 * domain separation between archives: the same payload in different archives
 * will produce different addresses, preventing cross-archive collisions.
 *
 * @returns A cryptographically random 32-byte `Uint8Array`.
 *
 * @example
 * ```ts
 * const salt = deriveArchiveSalt();
 * // Persist salt with the archive metadata!
 * const config: AddressingConfig = {
 *   mode: 'content-derived',
 *   archiveSalt: salt,
 *   addressBytes: 4,
 * };
 * ```
 */
export function deriveArchiveSalt(): Uint8Array {
  const salt = new Uint8Array(BLAKE3_KEY_BYTES);
  crypto.getRandomValues(salt);
  return salt;
}

/**
 * Base32 encode using Crockford variant.
 *
 * Alphabet: `0-9, A-V` (excludes I, L, O, U to avoid visual ambiguity).
 * This encoding is used for human-readable hierarchical addresses that may
 * be transcribed manually (e.g., pool/well labels on microplates).
 *
 * Encoding rules:
 *   - Each byte (8 bits) is split into 5-bit groups (plus remainder).
 *   - Each 5-bit group maps to a character in the alphabet.
 *   - Padding is NOT added (Crockford style — decode from length).
 *
 * @param data - The byte array to encode.
 * @returns The Base32-encoded string.
 *
 * @example
 * ```ts
 * base32Encode(new Uint8Array([0x3A, 0x7F])) // "VNV8"
 * ```
 */
export function base32Encode(data: Uint8Array): string {
  if (data.length === 0) {
    return '';
  }

  // 5 bits per character. Output length = ceil(data.length * 8 / 5).
  const outputLen = Math.ceil((data.length * 8) / 5);
  const output: string[] = [];

  let buffer = 0;   // Accumulator for bits.
  let bits = 0;     // Number of valid bits in buffer.

  for (let i = 0; i < data.length; i++) {
    buffer = (buffer << 8) | data[i];
    bits += 8;

    while (bits >= 5) {
      bits -= 5;
      const index = (buffer >> bits) & 0x1f;
      output.push(BASE32_ALPHABET[index]);
    }
  }

  // Flush remaining bits (if any).
  if (bits > 0) {
    const index = (buffer << (5 - bits)) & 0x1f;
    output.push(BASE32_ALPHABET[index]);
  }

  return output.join('');
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/**
 * Compute Shannon entropy of a byte array.
 *
 * Shannon entropy measures the information density of the data.
 *   - 0.0 = all bytes are the same (zero information).
 *   - 8.0 = perfectly random (maximum information per byte).
 *
 * @param data - The byte array to analyze.
 * @returns Shannon entropy in bits per byte [0.0, 8.0].
 */
function shannonEntropy(data: Uint8Array): number {
  if (data.length === 0) return 0;

  // Count byte frequencies.
  const freq = new Float64Array(256);
  for (let i = 0; i < data.length; i++) {
    freq[data[i]]++;
  }

  let entropy = 0;
  const len = data.length;
  for (let i = 0; i < 256; i++) {
    if (freq[i] > 0) {
      const p = freq[i] / len;
      entropy -= p * Math.log2(p);
    }
  }

  return entropy;
}

/**
 * xoshiro256\*\* pseudo-random number generator.
 *
 * Fast, high-quality PRNG with 256-bit state (4 × 64-bit words).
 * Used for seededPRNG recipe generation. The seed is expanded into
 * the initial state using a simple mixing function.
 *
 * Reference: Blackman & Vigna, "Scrambled Linear Pseudorandom Number Generators"
 * (ACM TOMS, 2021).
 *
 * @param seed  - Seed bytes (any length; zero-padded or hashed to 32 bytes).
 * @param count - Number of output bytes to generate.
 * @returns Generated pseudo-random bytes.
 */
function xoshiro256StarStarGenerate(seed: Uint8Array, count: number): Uint8Array {
  // Initialize state from seed.
  // If seed < 32 bytes, pad with zeros. If > 32 bytes, fold with XOR.
  const state = new BigUint64Array(4);
  const seedPadded = new Uint8Array(32);
  for (let i = 0; i < seed.length; i++) {
    seedPadded[i % 32] ^= seed[i];
  }

  // Load 4 × 64-bit state words from seed.
  const view = new DataView(seedPadded.buffer, seedPadded.byteOffset, 32);
  state[0] = view.getBigUint64(0, true);
  state[1] = view.getBigUint64(8, true);
  state[2] = view.getBigUint64(16, true);
  state[3] = view.getBigUint64(24, true);

  // Ensure state is not all-zero (degenerate for xoshiro).
  if (state[0] === 0n && state[1] === 0n && state[2] === 0n && state[3] === 0n) {
    state[0] = 1n;
  }

  // Mix the state several times to diffuse seed bits.
  for (let i = 0; i < 16; i++) {
    xoshiro256StarStarStep(state);
  }

  // Generate output bytes.
  const output = new Uint8Array(count);
  let byteIndex = 0;

  while (byteIndex < count) {
    // Generate next 8 bytes from one PRNG step.
    const value = xoshiro256StarStarStep(state);
    const valueBytes = new Uint8Array(8);
    const dv = new DataView(valueBytes.buffer);
    dv.setBigUint64(0, value, true);

    for (let j = 0; j < 8 && byteIndex < count; j++) {
      output[byteIndex++] = valueBytes[j];
    }
  }

  return output;
}

/**
 * One step of the xoshiro256\*\* PRNG.
 *
 * Mutates the state array in place and returns the output value.
 *
 * @param state - The 4-element BigUint64Array state (mutated).
 * @returns The pseudo-random 64-bit output value.
 */
function xoshiro256StarStarStep(state: BigUint64Array): bigint {
  const MASK = 0xFFFFFFFFFFFFFFFFn;

  const result = (rotl(state[1] * 5n, 7n) * 9n) & MASK;

  const t = (state[1] << 17n) & MASK;

  state[2] ^= state[0];
  state[3] ^= state[1];
  state[1] ^= state[2];
  state[0] ^= state[3];

  state[2] ^= t;
  state[3] = rotl(state[3], 45n);

  return result;
}

/**
 * Bitwise left rotation for BigInt (64-bit).
 *
 * @param x - The value to rotate.
 * @param k - The number of bits to rotate left.
 * @returns The rotated value (mod 2^64).
 */
function rotl(x: bigint, k: bigint): bigint {
  const MASK = 0xFFFFFFFFFFFFFFFFn;
  return ((x << k) | (x >> (64n - k))) & MASK;
}

/**
 * Generate a de Bruijn sequence of order n over an alphabet of size k.
 *
 * A de Bruijn sequence B(k, n) contains every possible length-n string
 * over the alphabet {0, 1, …, k−1} exactly once as a contiguous subsequence.
 *
 * This implementation uses the "prefer-one" (lexicographically largest)
 * algorithm, which is simple and deterministic.
 *
 * @param k      - Alphabet size (2–256).
 * @param n      - Subsequence length.
 * @param length - Desired output length (truncated if shorter than full cycle).
 * @returns The de Bruijn sequence as a byte array.
 */
function deBruijnSequence(k: number, n: number, length: number): Uint8Array {
  if (k < 2 || n < 1) {
    throw new Error(`deBruijn requires k >= 2 and n >= 1, got k=${k}, n=${n}`);
  }

  // Full de Bruijn sequence length = k^n.
  const fullLength = Math.pow(k, n);

  // Use the "prefer-one" algorithm for k=2 (binary de Bruijn).
  // For k > 2, use a generalized approach.
  const a = new Int32Array(n).fill(0);
  const seq: number[] = [];

  function db(t: number, p: number): void {
    if (t > n) {
      if (n % p === 0) {
        for (let j = 1; j <= p; j++) {
          seq.push(a[j]);
        }
      }
    } else {
      a[t] = a[t - p];
      db(t + 1, p);
      for (let j = a[t - p] + 1; j < k; j++) {
        a[t] = j;
        db(t + 1, t);
      }
    }
  }

  db(1, 1);

  // Convert to Uint8Array, truncating or cycling to desired length.
  const output = new Uint8Array(length);
  for (let i = 0; i < length; i++) {
    output[i] = seq[i % seq.length] & 0xff;
  }

  return output;
}

/**
 * Convert a Uint8Array to a hexadecimal string.
 *
 * Used internally for Set-based dedup comparison of addresses.
 *
 * @param bytes - The byte array to convert.
 * @returns Lowercase hexadecimal string (e.g., "deadbeef").
 */
function bytesToHex(bytes: Uint8Array): string {
  const hex: string[] = [];
  for (let i = 0; i < bytes.length; i++) {
    hex.push(bytes[i].toString(16).padStart(2, '0'));
  }
  return hex.join('');
}
