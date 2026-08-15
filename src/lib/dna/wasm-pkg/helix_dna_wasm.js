/* @ts-self-types="./helix_dna_wasm.d.ts" */

class WasmLdpcCode {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WasmLdpcCodeFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_wasmldpccode_free(ptr, 0);
    }
    /**
     * @param {Uint8Array} recv
     * @returns {Uint8Array}
     */
    decode(recv) {
        const ptr0 = passArray8ToWasm0(recv, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.wasmldpccode_decode(this.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v2;
    }
    /**
     * @param {Uint8Array} data
     * @returns {Uint8Array}
     */
    encode(data) {
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.wasmldpccode_encode(this.__wbg_ptr, ptr0, len0);
        var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v2;
    }
    /**
     * @param {number} n
     * @param {number} k
     */
    constructor(n, k) {
        const ret = wasm.wasmldpccode_new(n, k);
        this.__wbg_ptr = ret;
        WasmLdpcCodeFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
}
if (Symbol.dispose) WasmLdpcCode.prototype[Symbol.dispose] = WasmLdpcCode.prototype.free;
exports.WasmLdpcCode = WasmLdpcCode;

/**
 * @param {Uint8Array} all_reads
 * @param {Uint8Array} read_offsets
 * @param {Uint8Array} read_lengths
 * @param {Uint8Array} fwd_primer
 * @param {Uint8Array} rev_primer
 * @param {number} oligo_count
 * @param {number} inner_n
 * @param {number} inner_k
 * @param {number} total_inner_bytes
 * @returns {Uint8Array}
 */
function batch_decode_all(all_reads, read_offsets, read_lengths, fwd_primer, rev_primer, oligo_count, inner_n, inner_k, total_inner_bytes) {
    const ptr0 = passArray8ToWasm0(all_reads, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(read_offsets, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArray8ToWasm0(read_lengths, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArray8ToWasm0(fwd_primer, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArray8ToWasm0(rev_primer, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.batch_decode_all(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, oligo_count, inner_n, inner_k, total_inner_bytes);
    var v6 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v6;
}
exports.batch_decode_all = batch_decode_all;

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
 * @param {Uint8Array} all_reads
 * @param {Uint8Array} read_offsets
 * @param {Uint8Array} read_lengths
 * @param {Uint8Array} fwd_primer
 * @param {Uint8Array} rev_primer
 * @param {number} oligo_count
 * @param {number} inner_n
 * @param {number} inner_k
 * @param {number} total_inner_bytes
 * @param {number} outer_n
 * @param {number} outer_k
 * @param {number} payload_bytes
 * @param {number} file_size
 * @param {boolean} use_deflate
 * @returns {Uint8Array}
 */
function full_decode(all_reads, read_offsets, read_lengths, fwd_primer, rev_primer, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes, file_size, use_deflate) {
    const ptr0 = passArray8ToWasm0(all_reads, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(read_offsets, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArray8ToWasm0(read_lengths, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArray8ToWasm0(fwd_primer, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArray8ToWasm0(rev_primer, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.full_decode(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes, file_size, use_deflate);
    var v6 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v6;
}
exports.full_decode = full_decode;

/**
 * @param {Uint8Array} all_reads
 * @param {Uint8Array} read_offsets
 * @param {Uint8Array} read_lengths
 * @param {Uint8Array} fwd_primer
 * @param {Uint8Array} rev_primer
 * @param {number} oligo_count
 * @param {number} inner_n
 * @param {number} inner_k
 * @param {number} total_inner_bytes
 * @param {number} outer_n
 * @param {number} outer_k
 * @param {number} payload_bytes
 * @param {number} file_size
 * @param {boolean} use_deflate
 * @param {number} max_homopolymer
 * @param {number} block_size
 * @returns {Uint8Array}
 */
function full_decode_arithmetic(all_reads, read_offsets, read_lengths, fwd_primer, rev_primer, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes, file_size, use_deflate, max_homopolymer, block_size) {
    const ptr0 = passArray8ToWasm0(all_reads, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(read_offsets, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArray8ToWasm0(read_lengths, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArray8ToWasm0(fwd_primer, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArray8ToWasm0(rev_primer, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.full_decode_arithmetic(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes, file_size, use_deflate, max_homopolymer, block_size);
    var v6 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v6;
}
exports.full_decode_arithmetic = full_decode_arithmetic;

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
 * @param {Uint8Array} all_reads
 * @param {Uint8Array} read_offsets
 * @param {Uint8Array} read_lengths
 * @param {Uint8Array} fwd_primer
 * @param {Uint8Array} rev_primer
 * @param {number} oligo_count
 * @param {number} inner_n
 * @param {number} inner_k
 * @param {number} total_inner_bytes
 * @param {number} outer_n
 * @param {number} outer_k
 * @param {number} payload_bytes
 * @param {number} file_size
 * @param {boolean} use_deflate
 * @param {number} max_homopolymer
 * @param {number} block_size
 * @returns {Uint8Array}
 */
function full_decode_arithmetic_crc(all_reads, read_offsets, read_lengths, fwd_primer, rev_primer, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes, file_size, use_deflate, max_homopolymer, block_size) {
    const ptr0 = passArray8ToWasm0(all_reads, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(read_offsets, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArray8ToWasm0(read_lengths, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArray8ToWasm0(fwd_primer, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArray8ToWasm0(rev_primer, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.full_decode_arithmetic_crc(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes, file_size, use_deflate, max_homopolymer, block_size);
    var v6 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v6;
}
exports.full_decode_arithmetic_crc = full_decode_arithmetic_crc;

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
 * @param {Uint8Array} all_reads
 * @param {Uint8Array} read_offsets
 * @param {Uint8Array} read_lengths
 * @param {Uint8Array} fwd_primer
 * @param {Uint8Array} rev_primer
 * @param {number} oligo_count
 * @param {number} inner_n
 * @param {number} inner_k
 * @param {number} total_inner_bytes
 * @param {number} outer_n
 * @param {number} outer_k
 * @param {number} payload_bytes
 * @param {number} file_size
 * @param {boolean} use_deflate
 * @param {number} interleave_depth
 * @returns {Uint8Array}
 */
function full_decode_interleaved(all_reads, read_offsets, read_lengths, fwd_primer, rev_primer, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes, file_size, use_deflate, interleave_depth) {
    const ptr0 = passArray8ToWasm0(all_reads, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(read_offsets, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArray8ToWasm0(read_lengths, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArray8ToWasm0(fwd_primer, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArray8ToWasm0(rev_primer, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.full_decode_interleaved(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes, file_size, use_deflate, interleave_depth);
    var v6 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v6;
}
exports.full_decode_interleaved = full_decode_interleaved;

/**
 * Debug: return the pre-deflate bytes from full_decode (without truncation).
 * This lets us compare Rust-side total vs JS-side totalPayload byte-by-byte.
 * @param {Uint8Array} all_reads
 * @param {Uint8Array} read_offsets
 * @param {Uint8Array} read_lengths
 * @param {Uint8Array} fwd_primer
 * @param {Uint8Array} rev_primer
 * @param {number} oligo_count
 * @param {number} inner_n
 * @param {number} inner_k
 * @param {number} total_inner_bytes
 * @param {number} outer_n
 * @param {number} outer_k
 * @param {number} payload_bytes
 * @returns {Uint8Array}
 */
function full_decode_pre_deflate(all_reads, read_offsets, read_lengths, fwd_primer, rev_primer, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes) {
    const ptr0 = passArray8ToWasm0(all_reads, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(read_offsets, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArray8ToWasm0(read_lengths, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    const ptr3 = passArray8ToWasm0(fwd_primer, wasm.__wbindgen_malloc);
    const len3 = WASM_VECTOR_LEN;
    const ptr4 = passArray8ToWasm0(rev_primer, wasm.__wbindgen_malloc);
    const len4 = WASM_VECTOR_LEN;
    const ret = wasm.full_decode_pre_deflate(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3, ptr4, len4, oligo_count, inner_n, inner_k, total_inner_bytes, outer_n, outer_k, payload_bytes);
    var v6 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v6;
}
exports.full_decode_pre_deflate = full_decode_pre_deflate;

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
 * @param {Uint8Array} received_bits
 * @param {number} num_info_bits
 * @param {Uint32Array} generators
 * @param {number} constraint_length
 * @param {number} max_drift
 * @returns {Uint8Array}
 */
function indel_viterbi_decode(received_bits, num_info_bits, generators, constraint_length, max_drift) {
    const ptr0 = passArray8ToWasm0(received_bits, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray32ToWasm0(generators, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.indel_viterbi_decode(ptr0, len0, num_info_bits, ptr1, len1, constraint_length, max_drift);
    var v3 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v3;
}
exports.indel_viterbi_decode = indel_viterbi_decode;

/**
 * @param {Uint8Array} batch_result
 * @param {number} oligo_count
 * @param {number} inner_k
 * @param {number} payload_bytes
 * @param {number} outer_n
 * @param {number} outer_k
 * @param {number} file_size
 * @param {boolean} use_deflate
 * @returns {Uint8Array}
 */
function outer_rs_inflate(batch_result, oligo_count, inner_k, payload_bytes, outer_n, outer_k, file_size, use_deflate) {
    const ptr0 = passArray8ToWasm0(batch_result, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.outer_rs_inflate(ptr0, len0, oligo_count, inner_k, payload_bytes, outer_n, outer_k, file_size, use_deflate);
    var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v2;
}
exports.outer_rs_inflate = outer_rs_inflate;

/**
 * Test arithmetic decode: takes DNA bytes + max_homopolymer + expected_bytes.
 * Returns decoded bytes. Used to verify Rust decoder matches JS encoder.
 * @param {Uint8Array} dna
 * @param {number} max_homopolymer
 * @param {number} expected_bytes
 * @returns {Uint8Array}
 */
function test_arithmetic_decode(dna, max_homopolymer, expected_bytes) {
    const ptr0 = passArray8ToWasm0(dna, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.test_arithmetic_decode(ptr0, len0, max_homopolymer, expected_bytes);
    var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v2;
}
exports.test_arithmetic_decode = test_arithmetic_decode;

/**
 * Test block-wise arithmetic decode. Used to verify blocked decoder matches
 * blocked JS encoder.
 * @param {Uint8Array} dna
 * @param {number} max_homopolymer
 * @param {number} expected_bytes
 * @param {number} block_size
 * @returns {Uint8Array}
 */
function test_arithmetic_decode_blocked(dna, max_homopolymer, expected_bytes, block_size) {
    const ptr0 = passArray8ToWasm0(dna, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.test_arithmetic_decode_blocked(ptr0, len0, max_homopolymer, expected_bytes, block_size);
    var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v2;
}
exports.test_arithmetic_decode_blocked = test_arithmetic_decode_blocked;

/**
 * Test CRC-marker arithmetic decode: returns decoded bytes + erasure bitmap
 * (packed as 1 byte per position: 0 = OK, 1 = erased).
 * @param {Uint8Array} dna
 * @param {number} max_homopolymer
 * @param {number} expected_bytes
 * @param {number} block_size
 * @returns {Uint8Array}
 */
function test_arithmetic_decode_crc(dna, max_homopolymer, expected_bytes, block_size) {
    const ptr0 = passArray8ToWasm0(dna, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.test_arithmetic_decode_crc(ptr0, len0, max_homopolymer, expected_bytes, block_size);
    var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v2;
}
exports.test_arithmetic_decode_crc = test_arithmetic_decode_crc;

/**
 * Test inflate function — returns (success, output_len, error_status)
 * Used to debug DEFLATE issues.
 * @param {Uint8Array} data
 * @returns {string}
 */
function test_inflate(data) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.test_inflate(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}
exports.test_inflate = test_inflate;

/**
 * Test Rust RS216 decoder. Input: codeword (n u16 LE) + erased positions (u32 LE).
 * Output: decoded data (k u16 LE) or empty vec on failure.
 * @param {Uint8Array} codeword
 * @param {Uint8Array} erased
 * @param {number} n
 * @param {number} k
 * @returns {Uint8Array}
 */
function test_rs216_decode(codeword, erased, n, k) {
    const ptr0 = passArray8ToWasm0(codeword, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(erased, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.test_rs216_decode(ptr0, len0, ptr1, len1, n, k);
    var v3 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v3;
}
exports.test_rs216_decode = test_rs216_decode;
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_throw_bb96b2010945f0bc: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbindgen_cast_0000000000000001: function(arg0, arg1) {
            // Cast intrinsic for `Ref(String) -> Externref`.
            const ret = getStringFromWasm0(arg0, arg1);
            return ret;
        },
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        "./helix_dna_wasm_bg.js": import0,
    };
}

const WasmLdpcCodeFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_wasmldpccode_free(ptr, 1));

function getArrayU8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
}

function getStringFromWasm0(ptr, len) {
    return decodeText(ptr >>> 0, len);
}

let cachedUint32ArrayMemory0 = null;
function getUint32ArrayMemory0() {
    if (cachedUint32ArrayMemory0 === null || cachedUint32ArrayMemory0.byteLength === 0) {
        cachedUint32ArrayMemory0 = new Uint32Array(wasm.memory.buffer);
    }
    return cachedUint32ArrayMemory0;
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function passArray32ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 4, 4) >>> 0;
    getUint32ArrayMemory0().set(arg, ptr / 4);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passArray8ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 1, 1) >>> 0;
    getUint8ArrayMemory0().set(arg, ptr / 1);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_externrefs.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
function decodeText(ptr, len) {
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

let WASM_VECTOR_LEN = 0;

const wasmPath = `${__dirname}/helix_dna_wasm_bg.wasm`;
const wasmBytes = require('fs').readFileSync(wasmPath);
const wasmModule = new WebAssembly.Module(wasmBytes);
let wasmInstance = new WebAssembly.Instance(wasmModule, __wbg_get_imports());
let wasm = wasmInstance.exports;
wasm.__wbindgen_start();
