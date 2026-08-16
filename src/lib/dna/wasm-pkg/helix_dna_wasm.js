/* @ts-self-types="./helix_dna_wasm.d.ts" */

/**
 * Result of pack operation with optional error.
 * Use getter functions to access fields from JS.
 */
export class PackResult {
    static __wrap(ptr) {
        const obj = Object.create(PackResult.prototype);
        obj.__wbg_ptr = ptr;
        PackResultFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        PackResultFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_packresult_free(ptr, 0);
    }
    /**
     * @returns {Uint8Array}
     */
    data() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.packresult_data(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayU8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export2(r0, r1 * 1, 1);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * @returns {string}
     */
    error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.packresult_error(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            deferred1_0 = r0;
            deferred1_1 = r1;
            return getStringFromWasm0(r0, r1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export2(deferred1_0, deferred1_1, 1);
        }
    }
}
if (Symbol.dispose) PackResult.prototype[Symbol.dispose] = PackResult.prototype.free;

/**
 * Result of RS decode. Use getter functions to access fields from JS.
 */
export class RsDecodeResult {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        RsDecodeResultFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_rsdecoderesult_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    corrected() {
        const ret = wasm.rsdecoderesult_corrected(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {Uint8Array}
     */
    data() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.rsdecoderesult_data(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayU8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export2(r0, r1 * 1, 1);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * @returns {number}
     */
    erased() {
        const ret = wasm.rsdecoderesult_erased(this.__wbg_ptr);
        return ret >>> 0;
    }
}
if (Symbol.dispose) RsDecodeResult.prototype[Symbol.dispose] = RsDecodeResult.prototype.free;

/**
 * Sequencing error profile.
 */
export class SequencingProfile {
    static __wrap(ptr) {
        const obj = Object.create(SequencingProfile.prototype);
        obj.__wbg_ptr = ptr;
        SequencingProfileFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SequencingProfileFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_sequencingprofile_free(ptr, 0);
    }
    /**
     * Deletion rate added by sequencing.
     * @returns {number}
     */
    get del_rate() {
        const ret = wasm.__wbg_get_sequencingprofile_del_rate(this.__wbg_ptr);
        return ret;
    }
    /**
     * Insertion rate added by sequencing.
     * @returns {number}
     */
    get ins_rate() {
        const ret = wasm.__wbg_get_sequencingprofile_ins_rate(this.__wbg_ptr);
        return ret;
    }
    /**
     * Substitution rate added by sequencing.
     * @returns {number}
     */
    get sub_rate() {
        const ret = wasm.__wbg_get_sequencingprofile_sub_rate(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {SequencingProfile}
     */
    static illumina() {
        const ret = wasm.sequencingprofile_illumina();
        return SequencingProfile.__wrap(ret);
    }
    /**
     * @returns {SequencingProfile}
     */
    static nanopore() {
        const ret = wasm.sequencingprofile_nanopore();
        return SequencingProfile.__wrap(ret);
    }
    constructor() {
        const ret = wasm.sequencingprofile_new();
        this.__wbg_ptr = ret;
        SequencingProfileFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * @returns {SequencingProfile}
     */
    static pacbio_hifi() {
        const ret = wasm.sequencingprofile_pacbio_hifi();
        return SequencingProfile.__wrap(ret);
    }
    /**
     * Deletion rate added by sequencing.
     * @param {number} arg0
     */
    set del_rate(arg0) {
        wasm.__wbg_set_sequencingprofile_del_rate(this.__wbg_ptr, arg0);
    }
    /**
     * Insertion rate added by sequencing.
     * @param {number} arg0
     */
    set ins_rate(arg0) {
        wasm.__wbg_set_sequencingprofile_ins_rate(this.__wbg_ptr, arg0);
    }
    /**
     * Substitution rate added by sequencing.
     * @param {number} arg0
     */
    set sub_rate(arg0) {
        wasm.__wbg_set_sequencingprofile_sub_rate(this.__wbg_ptr, arg0);
    }
}
if (Symbol.dispose) SequencingProfile.prototype[Symbol.dispose] = SequencingProfile.prototype.free;

/**
 * Synthesis error profile.
 */
export class SynthesisProfile {
    static __wrap(ptr) {
        const obj = Object.create(SynthesisProfile.prototype);
        obj.__wbg_ptr = ptr;
        SynthesisProfileFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SynthesisProfileFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_synthesisprofile_free(ptr, 0);
    }
    /**
     * Deletion rate per position.
     * @returns {number}
     */
    get del_rate() {
        const ret = wasm.__wbg_get_synthesisprofile_del_rate(this.__wbg_ptr);
        return ret;
    }
    /**
     * 5' end degradation factor (errors increase toward 5' end).
     * @returns {number}
     */
    get five_prime_bias() {
        const ret = wasm.__wbg_get_synthesisprofile_five_prime_bias(this.__wbg_ptr);
        return ret;
    }
    /**
     * Insertion rate per position.
     * @returns {number}
     */
    get ins_rate() {
        const ret = wasm.__wbg_get_synthesisprofile_ins_rate(this.__wbg_ptr);
        return ret;
    }
    /**
     * Base substitution rate (uniform across positions).
     * @returns {number}
     */
    get sub_rate() {
        const ret = wasm.__wbg_get_synthesisprofile_sub_rate(this.__wbg_ptr);
        return ret;
    }
    /**
     * 3' end degradation factor.
     * @returns {number}
     */
    get three_prime_bias() {
        const ret = wasm.__wbg_get_synthesisprofile_three_prime_bias(this.__wbg_ptr);
        return ret;
    }
    /**
     * Deletion rate per position.
     * @param {number} arg0
     */
    set del_rate(arg0) {
        wasm.__wbg_set_synthesisprofile_del_rate(this.__wbg_ptr, arg0);
    }
    /**
     * 5' end degradation factor (errors increase toward 5' end).
     * @param {number} arg0
     */
    set five_prime_bias(arg0) {
        wasm.__wbg_set_synthesisprofile_five_prime_bias(this.__wbg_ptr, arg0);
    }
    /**
     * Insertion rate per position.
     * @param {number} arg0
     */
    set ins_rate(arg0) {
        wasm.__wbg_set_synthesisprofile_ins_rate(this.__wbg_ptr, arg0);
    }
    /**
     * Base substitution rate (uniform across positions).
     * @param {number} arg0
     */
    set sub_rate(arg0) {
        wasm.__wbg_set_synthesisprofile_sub_rate(this.__wbg_ptr, arg0);
    }
    /**
     * 3' end degradation factor.
     * @param {number} arg0
     */
    set three_prime_bias(arg0) {
        wasm.__wbg_set_synthesisprofile_three_prime_bias(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {SynthesisProfile}
     */
    static illuminina() {
        const ret = wasm.synthesisprofile_illuminina();
        return SynthesisProfile.__wrap(ret);
    }
    /**
     * @returns {SynthesisProfile}
     */
    static nanopore() {
        const ret = wasm.synthesisprofile_nanopore();
        return SynthesisProfile.__wrap(ret);
    }
    constructor() {
        const ret = wasm.synthesisprofile_new();
        this.__wbg_ptr = ret;
        SynthesisProfileFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
}
if (Symbol.dispose) SynthesisProfile.prototype[Symbol.dispose] = SynthesisProfile.prototype.free;

/**
 * Compress DNA bytes using order-0 adaptive arithmetic coding.
 * Input should be raw bytes (not 2-bit packed). For DNA sequences,
 * the caller should first convert to a byte representation.
 * @param {Uint8Array} data
 * @returns {Uint8Array}
 */
export function arith_compress(data) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.arith_compress(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Decompress data compressed by arith_compress.
 * @param {Uint8Array} data
 * @returns {Uint8Array}
 */
export function arith_decompress(data) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.arith_decompress(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * BHE k>1 decode using FSM with u128 arithmetic.
 * @param {Uint8Array} dna
 * @param {number} num_bytes
 * @param {number} max_run
 * @returns {Uint8Array}
 */
export function bhe_decode_fsm(dna, num_bytes, max_run) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(dna, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.bhe_decode_fsm(retptr, ptr0, len0, num_bytes, max_run);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * BHE k=1 decode: convert DNA back to bytes.
 * @param {Uint8Array} dna
 * @param {number} num_bytes
 * @returns {Uint8Array}
 */
export function bhe_decode_k1(dna, num_bytes) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(dna, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.bhe_decode_k1(retptr, ptr0, len0, num_bytes);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * BHE k>1 encode using FSM with u128 arithmetic.
 *
 * For inputs ≤ 16 bytes: direct single-value encoding.
 * For inputs > 16 bytes: chunked encoding with a header storing chunk DNA lengths.
 * @param {Uint8Array} bytes
 * @param {number} max_run
 * @returns {Uint8Array}
 */
export function bhe_encode_fsm(bytes, max_run) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(bytes, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.bhe_encode_fsm(retptr, ptr0, len0, max_run);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * BHE k=1 encode: convert bytes to DNA with no homopolymers.
 *
 * Uses u128 arithmetic for ~50× speedup over JS BigInt.
 * For inputs ≤ 16 bytes: direct single-value encoding.
 * For inputs > 16 bytes: chunked encoding with junction homopolymer avoidance.
 * @param {Uint8Array} bytes
 * @returns {Uint8Array}
 */
export function bhe_encode_k1(bytes) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(bytes, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.bhe_encode_k1(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Validate BHE-encoded DNA satisfies homopolymer constraint.
 * @param {Uint8Array} dna
 * @param {number} max_run
 * @returns {boolean}
 */
export function bhe_validate(dna, max_run) {
    const ptr0 = passArray8ToWasm0(dna, wasm.__wbindgen_export);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.bhe_validate(ptr0, len0, max_run);
    return ret !== 0;
}

/**
 * Compute Hamming distance between two packed DNA arrays.
 * Uses popcount(XOR) / 2 since each mismatching base contributes 2 set bits.
 * @param {Uint8Array} a
 * @param {Uint8Array} b
 * @returns {number}
 */
export function bit_parallel_hamming(a, b) {
    const ptr0 = passArray8ToWasm0(a, wasm.__wbindgen_export);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(b, wasm.__wbindgen_export);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.bit_parallel_hamming(ptr0, len0, ptr1, len1);
    return ret >>> 0;
}

/**
 * Find all exact occurrences of pattern in text using shift-and algorithm.
 * Pattern length limit: 32 bases (fits in u32).
 * @param {Uint8Array} pattern
 * @param {Uint8Array} text
 * @returns {Uint32Array}
 */
export function bit_parallel_match(pattern, text) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(pattern, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArray8ToWasm0(text, wasm.__wbindgen_export);
        const len1 = WASM_VECTOR_LEN;
        wasm.bit_parallel_match(retptr, ptr0, len0, ptr1, len1);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v3 = getArrayU32FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 4, 4);
        return v3;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Compute complement of packed DNA (A↔T, C↔G) via XOR 0xFF per byte.
 * @param {Uint8Array} bits
 * @returns {Uint8Array}
 */
export function complement(bits) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(bits, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.complement(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Compress 2-bit packed DNA using order-1 context modeling.
 * Each base (A=0, C=1, G=2, T=3) is encoded with a context-dependent
 * model where the context is the previous base.
 * @param {Uint8Array} packed
 * @param {number} num_bases
 * @returns {Uint8Array}
 */
export function dna_compress_order1(packed, num_bases) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(packed, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.dna_compress_order1(retptr, ptr0, len0, num_bases);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Decompress DNA compressed by dna_compress_order1.
 * @param {Uint8Array} data
 * @returns {Uint8Array}
 */
export function dna_decompress_order1(data) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.dna_decompress_order1(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * @param {number} n_bits
 * @param {number} k_bits
 * @param {number} dv
 * @returns {number}
 */
export function ldpc_create(n_bits, k_bits, dv) {
    const ret = wasm.ldpc_create(n_bits, k_bits, dv);
    return ret >>> 0;
}

/**
 * @param {number} handle
 * @param {Float64Array} llr
 * @param {number} max_iter
 * @returns {Uint8Array}
 */
export function ldpc_decode(handle, llr, max_iter) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArrayF64ToWasm0(llr, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.ldpc_decode(retptr, handle, ptr0, len0, max_iter);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * @param {number} handle
 * @param {Uint8Array} info
 * @returns {Uint8Array}
 */
export function ldpc_encode(handle, info) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(info, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.ldpc_encode(retptr, handle, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * @param {number} handle
 */
export function ldpc_free(handle) {
    wasm.ldpc_free(handle);
}

/**
 * Pack a DNA string into 2-bit bytes (4 bases per byte, MSB-first).
 *
 * Uses SIMD v128 for the inner loop when available (WASM).
 * For WASM, we process 16 chars at a time.
 * @param {string} dna
 * @returns {Uint8Array}
 */
export function pack_dna_to_bits(dna) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passStringToWasm0(dna, wasm.__wbindgen_export, wasm.__wbindgen_export3);
        const len0 = WASM_VECTOR_LEN;
        wasm.pack_dna_to_bits(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Pack DNA with validation. Returns (packed_bytes, error_message).
 * If error_message is non-empty, the input had invalid bases.
 * @param {string} dna
 * @returns {PackResult}
 */
export function pack_dna_validated(dna) {
    const ptr0 = passStringToWasm0(dna, wasm.__wbindgen_export, wasm.__wbindgen_export3);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.pack_dna_validated(ptr0, len0);
    return PackResult.__wrap(ret);
}

/**
 * Compute Rabin-Karp rolling hashes over packed DNA.
 * Returns array of (numBases - windowSize + 1) 32-bit hashes.
 * @param {Uint8Array} bits
 * @param {number} window_size
 * @returns {Uint32Array}
 */
export function rolling_hash(bits, window_size) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(bits, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.rolling_hash(retptr, ptr0, len0, window_size);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU32FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 4, 4);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * @param {number} n
 * @param {number} k
 * @returns {number}
 */
export function rs_create(n, k) {
    const ret = wasm.rs_create(n, k);
    return ret >>> 0;
}

/**
 * @param {number} handle
 * @param {Uint8Array} recv
 * @returns {Uint8Array}
 */
export function rs_decode(handle, recv) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(recv, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.rs_decode(retptr, handle, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * @param {number} handle
 * @param {Uint8Array} recv
 * @param {Uint32Array} erase_pos
 * @returns {Uint8Array}
 */
export function rs_decode_erasures(handle, recv, erase_pos) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(recv, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArray32ToWasm0(erase_pos, wasm.__wbindgen_export);
        const len1 = WASM_VECTOR_LEN;
        wasm.rs_decode_erasures(retptr, handle, ptr0, len0, ptr1, len1);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v3 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v3;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * @param {number} handle
 * @param {Uint8Array} data
 * @returns {Uint8Array}
 */
export function rs_encode(handle, data) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.rs_encode(retptr, handle, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * @param {number} handle
 */
export function rs_free(handle) {
    wasm.rs_free(handle);
}

/**
 * Simulate multiple oligos through the wetlab pipeline.
 *
 * Each oligo is processed independently with a unique RNG seed
 * derived from the base seed + oligo index.
 *
 * Returns a flat array: for each oligo, [length_u32_le, base1, base2, ...].
 * @param {Uint8Array} oligos
 * @param {Uint32Array} oligo_offsets
 * @param {number} synth_sub
 * @param {number} synth_ins
 * @param {number} synth_del
 * @param {number} seq_sub
 * @param {number} seq_ins
 * @param {number} seq_del
 * @param {number} base_seed
 * @returns {Uint8Array}
 */
export function simulate_batch(oligos, oligo_offsets, synth_sub, synth_ins, synth_del, seq_sub, seq_ins, seq_del, base_seed) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(oligos, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArray32ToWasm0(oligo_offsets, wasm.__wbindgen_export);
        const len1 = WASM_VECTOR_LEN;
        wasm.simulate_batch(retptr, ptr0, len0, ptr1, len1, synth_sub, synth_ins, synth_del, seq_sub, seq_ins, seq_del, base_seed);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v3 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v3;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Simulate a single oligo. Returns the read as ASCII bytes.
 * @param {Uint8Array} oligo
 * @param {number} synth_sub
 * @param {number} synth_ins
 * @param {number} synth_del
 * @param {number} seq_sub
 * @param {number} seq_ins
 * @param {number} seq_del
 * @param {number} seed
 * @returns {Uint8Array}
 */
export function simulate_single(oligo, synth_sub, synth_ins, synth_del, seq_sub, seq_ins, seq_del, seed) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(oligo, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.simulate_single(retptr, ptr0, len0, synth_sub, synth_ins, synth_del, seq_sub, seq_ins, seq_del, seed);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Unpack multiple packed DNA arrays in a single WASM call.
 *
 * This eliminates per-call JS↔WASM boundary overhead, which is the dominant
 * cost for small arrays (<500K bases). Instead of N separate calls each with
 * its own malloc/copy/invoke/free cycle, we:
 *   1. Receive all inputs as a single contiguous buffer (`packed_data`)
 *   2. Use `offsets` to locate each individual array within the buffer
 *   3. Use `num_bases` to know how many bases to output per array
 *   4. Write all results into a single contiguous output buffer
 *
 * The JS wrapper splits the flat output back into individual arrays.
 *
 * # Arguments
 * * `packed_data` — All packed arrays concatenated into one flat buffer
 * * `offsets` — Start offset of each array within `packed_data` (in packed bytes)
 * * `num_bases` — Number of nucleotides to output for each array
 * * `total_bases` — Sum of all `num_bases` (pre-computed for output allocation)
 *
 * # Returns
 * Flat Vec<u8> of all ASCII results concatenated. JS splits by num_bases.
 * @param {Uint8Array} packed_data
 * @param {Uint32Array} offsets
 * @param {Uint32Array} num_bases
 * @param {number} total_bases
 * @returns {Uint8Array}
 */
export function unpack_batch(packed_data, offsets, num_bases, total_bases) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(packed_data, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArray32ToWasm0(offsets, wasm.__wbindgen_export);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passArray32ToWasm0(num_bases, wasm.__wbindgen_export);
        const len2 = WASM_VECTOR_LEN;
        wasm.unpack_batch(retptr, ptr0, len0, ptr1, len1, ptr2, len2, total_bases);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v4 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v4;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Unpack 2-bit bytes into ASCII bytes (A=65, C=67, G=71, T=84).
 *
 * This is the hot path — processes 16 packed bytes → 64 ASCII bytes
 * per SIMD iteration using v128.load/store + i8x16 operations.
 *
 * Returns a Uint8Array of ASCII bytes. The JS wrapper converts to string.
 * @param {Uint8Array} packed
 * @param {number} num_bases
 * @returns {Uint8Array}
 */
export function unpack_bits_to_ascii(packed, num_bases) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(packed, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.unpack_bits_to_ascii(retptr, ptr0, len0, num_bases);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var v2 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export2(r0, r1 * 1, 1);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Unpack 2-bit bytes directly into a DNA string.
 * @param {Uint8Array} packed
 * @param {number} num_bases
 * @returns {string}
 */
export function unpack_bits_to_dna(packed, num_bases) {
    let deferred2_0;
    let deferred2_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(packed, wasm.__wbindgen_export);
        const len0 = WASM_VECTOR_LEN;
        wasm.unpack_bits_to_dna(retptr, ptr0, len0, num_bases);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        deferred2_0 = r0;
        deferred2_1 = r1;
        return getStringFromWasm0(r0, r1);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export2(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Package version
 * @returns {string}
 */
export function version() {
    let deferred1_0;
    let deferred1_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        wasm.version(retptr);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        deferred1_0 = r0;
        deferred1_1 = r1;
        return getStringFromWasm0(r0, r1);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export2(deferred1_0, deferred1_1, 1);
    }
}
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_throw_bb96b2010945f0bc: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
    };
    return {
        __proto__: null,
        "./helix_dna_wasm_bg.js": import0,
    };
}

const PackResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_packresult_free(ptr, 1));
const RsDecodeResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_rsdecoderesult_free(ptr, 1));
const SequencingProfileFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_sequencingprofile_free(ptr, 1));
const SynthesisProfileFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_synthesisprofile_free(ptr, 1));

function getArrayU32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayU8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
}

let cachedDataViewMemory0 = null;
function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null || cachedDataViewMemory0.buffer.detached === true || (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

let cachedFloat64ArrayMemory0 = null;
function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
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

function passArrayF64ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 8, 8) >>> 0;
    getFloat64ArrayMemory0().set(arg, ptr / 8);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passStringToWasm0(arg, malloc, realloc) {
    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }
    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = cachedTextEncoder.encodeInto(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

const cachedTextEncoder = new TextEncoder();

if (!('encodeInto' in cachedTextEncoder)) {
    cachedTextEncoder.encodeInto = function (arg, view) {
        const buf = cachedTextEncoder.encode(arg);
        view.set(buf);
        return {
            read: arg.length,
            written: buf.length
        };
    };
}

let WASM_VECTOR_LEN = 0;

let wasmModule, wasmInstance, wasm;
function __wbg_finalize_init(instance, module) {
    wasmInstance = instance;
    wasm = instance.exports;
    wasmModule = module;
    cachedDataViewMemory0 = null;
    cachedFloat64ArrayMemory0 = null;
    cachedUint32ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (!module.ok) {
            throw new Error(`failed to fetch Wasm: ${module.status} ${module.statusText} fetching '${module.url}'`);
        }

        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('helix_dna_wasm_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
