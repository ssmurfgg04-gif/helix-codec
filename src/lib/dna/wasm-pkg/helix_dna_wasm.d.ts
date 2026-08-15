/* tslint:disable */
/* eslint-disable */

export class WasmLdpcCode {
    free(): void;
    [Symbol.dispose](): void;
    decode(recv: Uint8Array): Uint8Array;
    encode(data: Uint8Array): Uint8Array;
    constructor(n: number, k: number);
}

export function batch_decode_all(all_reads: Uint8Array, read_offsets: Uint8Array, read_lengths: Uint8Array, fwd_primer: Uint8Array, rev_primer: Uint8Array, oligo_count: number, inner_n: number, inner_k: number, total_inner_bytes: number): Uint8Array;

/**
 * Full decode pipeline: ALL reads → recovered file bytes.
 *
 * Single WASM call runs:
 *   1. Primer trimming + Hamming-distance primer match (≤2 mismatches)
 *   2. DNA→bytes (4 nt → 1 byte) with SIMD128-batched XOR in LDPC encode
 *   3. Cluster by oligo address (XOR-unwhitened)
 *   4. Per-oligo LDPC decode + CRC-16 verification + address verification
 *      + seed-based payload unwhitening
 *   5. Outer RS GF(2^16) pure-erasure recovery (Rust implementation)
 *   6. Concatenate payloads
 *   7. DEFLATE inflate via miniz_oxide (pako-compatible)
 *   8. Trim to file_size
 *
 * Returns: recovered file bytes. Empty vec on total failure.
 */
export function full_decode(all_reads: Uint8Array, read_offsets: Uint8Array, read_lengths: Uint8Array, fwd_primer: Uint8Array, rev_primer: Uint8Array, oligo_count: number, inner_n: number, inner_k: number, total_inner_bytes: number, outer_n: number, outer_k: number, payload_bytes: number, file_size: number, use_deflate: boolean): Uint8Array;

export function full_decode_arithmetic(all_reads: Uint8Array, read_offsets: Uint8Array, read_lengths: Uint8Array, fwd_primer: Uint8Array, rev_primer: Uint8Array, oligo_count: number, inner_n: number, inner_k: number, total_inner_bytes: number, outer_n: number, outer_k: number, payload_bytes: number, file_size: number, use_deflate: boolean, max_homopolymer: number, block_size: number): Uint8Array;

/**
 * Full decode pipeline with ARITHMETIC mapping mode.
 *
 * For arithmetic mode, the DNA is arithmetic-encoded (stateful), so we CANNOT
 * extract the address from specific DNA positions (unlike direct mapping).
 * Instead, we decode each read individually:
 *   1. Trim primers
 *   2. Arithmetic decode inner DNA → bytes
 *   3. LDPC decode → data bytes
 *   4. CRC verify
 *   5. Extract address from data bytes 0-3
 *   6. First successful decode per oligo wins
 *   7. Outer RS + DEFLATE + trim
 *
 * Note: arithmetic coding has unbounded error propagation (1 DNA error → many
 * bit errors), so LDPC correction is less effective. At high coverage (10x+)
 * with low error rates, most reads decode successfully. At lower coverage or
 * higher error rates, more oligos are erased (recovered by outer RS).
 */
export function full_decode_arithmetic_crc(all_reads: Uint8Array, read_offsets: Uint8Array, read_lengths: Uint8Array, fwd_primer: Uint8Array, rev_primer: Uint8Array, oligo_count: number, inner_n: number, inner_k: number, total_inner_bytes: number, outer_n: number, outer_k: number, payload_bytes: number, file_size: number, use_deflate: boolean, max_homopolymer: number, block_size: number): Uint8Array;

/**
 * Full decode pipeline with interleaving support.
 *
 * Same as full_decode but deinterleaves consensus blocks across groups of
 * `interleave_depth` oligos before LDPC decode. This spreads burst errors:
 * a burst in one oligo becomes 1 error per codeword across `depth` oligos.
 *
 * The address bytes (first 4) are NOT interleaved — they remain at fixed
 * positions for clustering to work. Only the payload+parity region
 * (bytes 4..total_inner_bytes) is deinterleaved.
 */
export function full_decode_interleaved(all_reads: Uint8Array, read_offsets: Uint8Array, read_lengths: Uint8Array, fwd_primer: Uint8Array, rev_primer: Uint8Array, oligo_count: number, inner_n: number, inner_k: number, total_inner_bytes: number, outer_n: number, outer_k: number, payload_bytes: number, file_size: number, use_deflate: boolean, interleave_depth: number): Uint8Array;

/**
 * Debug: return the pre-deflate bytes from full_decode (without truncation).
 * This lets us compare Rust-side total vs JS-side totalPayload byte-by-byte.
 */
export function full_decode_pre_deflate(all_reads: Uint8Array, read_offsets: Uint8Array, read_lengths: Uint8Array, fwd_primer: Uint8Array, rev_primer: Uint8Array, oligo_count: number, inner_n: number, inner_k: number, total_inner_bytes: number, outer_n: number, outer_k: number, payload_bytes: number): Uint8Array;

/**
 * Indel-tolerant Viterbi decoder for convolutional codes.
 *
 * Implements an augmented-state Viterbi trellis where each state
 * includes both the encoder shift-register state AND a running drift
 * counter (net insertions − deletions).  This allows the Viterbi
 * algorithm to track and correct insertion/deletion errors that would
 * defeat a standard decoder.
 *
 * ## Trellis structure
 *
 * The trellis is processed one **received symbol** at a time (one
 * symbol = *n* bits for a rate 1/*n* code).  At each depth three
 * transition types are considered:
 *
 * 1. **Match / Substitution** – the received symbol corresponds to an
 *    encoder output.  Advance encoder by one step, consume one
 *    received symbol.  Branch metric = Hamming distance between
 *    expected and received symbol.
 *
 * 2. **Insertion** – the received symbol is spurious (not produced by
 *    the encoder).  Encoder state unchanged, consume one received
 *    symbol.  Branch metric = insertion penalty.
 *
 * 3. **Deletion + Match** – one encoder output was deleted (never
 *    received), and the current received symbol matches the *next*
 *    encoder output.  Advance encoder by two steps, consume one
 *    received symbol.  Branch metric = deletion penalty + Hamming
 *    distance for the matched symbol.
 *
 * ## Parameters
 *
 * - `received_bits`    – hard-decision bits from the demodulator
 * - `num_info_bits`    – number of information bits to decode (excl. tail)
 * - `generators`       – generator polynomials (MSB = input tap),
 *                        e.g. `[0o561, 0o753]` for the NASA K=9 code
 * - `constraint_length`– constraint length *K*
 * - `max_drift`        – maximum net drift (insertions − deletions) to track
 *
 * ## Returns
 *
 * Decoded information bits as a `Vec<u8>` (each element is 0 or 1).
 */
export function indel_viterbi_decode(received_bits: Uint8Array, num_info_bits: number, generators: Uint32Array, constraint_length: number, max_drift: number): Uint8Array;

export function outer_rs_inflate(batch_result: Uint8Array, oligo_count: number, inner_k: number, payload_bytes: number, outer_n: number, outer_k: number, file_size: number, use_deflate: boolean): Uint8Array;

/**
 * Test arithmetic decode: takes DNA bytes + max_homopolymer + expected_bytes.
 * Returns decoded bytes. Used to verify Rust decoder matches JS encoder.
 */
export function test_arithmetic_decode(dna: Uint8Array, max_homopolymer: number, expected_bytes: number): Uint8Array;

/**
 * Test block-wise arithmetic decode. Used to verify blocked decoder matches
 * blocked JS encoder.
 */
export function test_arithmetic_decode_blocked(dna: Uint8Array, max_homopolymer: number, expected_bytes: number, block_size: number): Uint8Array;

/**
 * Test CRC-marker arithmetic decode: returns decoded bytes + erasure bitmap
 * (packed as 1 byte per position: 0 = OK, 1 = erased).
 */
export function test_arithmetic_decode_crc(dna: Uint8Array, max_homopolymer: number, expected_bytes: number, block_size: number): Uint8Array;

/**
 * Test inflate function — returns (success, output_len, error_status)
 * Used to debug DEFLATE issues.
 */
export function test_inflate(data: Uint8Array): string;

/**
 * Test Rust RS216 decoder. Input: codeword (n u16 LE) + erased positions (u32 LE).
 * Output: decoded data (k u16 LE) or empty vec on failure.
 */
export function test_rs216_decode(codeword: Uint8Array, erased: Uint8Array, n: number, k: number): Uint8Array;
