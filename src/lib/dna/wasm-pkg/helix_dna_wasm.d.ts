/* tslint:disable */
/* eslint-disable */

/**
 * Result of pack operation with optional error.
 * Use getter functions to access fields from JS.
 */
export class PackResult {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    data(): Uint8Array;
    error(): string;
}

/**
 * Result of RS decode. Use getter functions to access fields from JS.
 */
export class RsDecodeResult {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    corrected(): number;
    data(): Uint8Array;
    erased(): number;
}

/**
 * Sequencing error profile.
 */
export class SequencingProfile {
    free(): void;
    [Symbol.dispose](): void;
    static illumina(): SequencingProfile;
    static nanopore(): SequencingProfile;
    constructor();
    static pacbio_hifi(): SequencingProfile;
    /**
     * Deletion rate added by sequencing.
     */
    del_rate: number;
    /**
     * Insertion rate added by sequencing.
     */
    ins_rate: number;
    /**
     * Substitution rate added by sequencing.
     */
    sub_rate: number;
}

/**
 * Synthesis error profile.
 */
export class SynthesisProfile {
    free(): void;
    [Symbol.dispose](): void;
    static illuminina(): SynthesisProfile;
    static nanopore(): SynthesisProfile;
    constructor();
    /**
     * Deletion rate per position.
     */
    del_rate: number;
    /**
     * 5' end degradation factor (errors increase toward 5' end).
     */
    five_prime_bias: number;
    /**
     * Insertion rate per position.
     */
    ins_rate: number;
    /**
     * Base substitution rate (uniform across positions).
     */
    sub_rate: number;
    /**
     * 3' end degradation factor.
     */
    three_prime_bias: number;
}

/**
 * Compress DNA bytes using order-0 adaptive arithmetic coding.
 * Input should be raw bytes (not 2-bit packed). For DNA sequences,
 * the caller should first convert to a byte representation.
 */
export function arith_compress(data: Uint8Array): Uint8Array;

/**
 * Decompress data compressed by arith_compress.
 */
export function arith_decompress(data: Uint8Array): Uint8Array;

/**
 * BHE k>1 decode using FSM with u128 arithmetic.
 */
export function bhe_decode_fsm(dna: Uint8Array, num_bytes: number, max_run: number): Uint8Array;

/**
 * BHE k=1 decode: convert DNA back to bytes.
 */
export function bhe_decode_k1(dna: Uint8Array, num_bytes: number): Uint8Array;

/**
 * BHE k>1 encode using FSM with u128 arithmetic.
 *
 * For inputs ≤ 16 bytes: direct single-value encoding.
 * For inputs > 16 bytes: chunked encoding with a header storing chunk DNA lengths.
 */
export function bhe_encode_fsm(bytes: Uint8Array, max_run: number): Uint8Array;

/**
 * BHE k=1 encode: convert bytes to DNA with no homopolymers.
 *
 * Uses u128 arithmetic for ~50× speedup over JS BigInt.
 * For inputs ≤ 16 bytes: direct single-value encoding.
 * For inputs > 16 bytes: chunked encoding with junction homopolymer avoidance.
 */
export function bhe_encode_k1(bytes: Uint8Array): Uint8Array;

/**
 * Validate BHE-encoded DNA satisfies homopolymer constraint.
 */
export function bhe_validate(dna: Uint8Array, max_run: number): boolean;

/**
 * Compute Hamming distance between two packed DNA arrays.
 * Uses popcount(XOR) / 2 since each mismatching base contributes 2 set bits.
 */
export function bit_parallel_hamming(a: Uint8Array, b: Uint8Array): number;

/**
 * Find all exact occurrences of pattern in text using shift-and algorithm.
 * Pattern length limit: 32 bases (fits in u32).
 */
export function bit_parallel_match(pattern: Uint8Array, text: Uint8Array): Uint32Array;

/**
 * Compute complement of packed DNA (A↔T, C↔G) via XOR 0xFF per byte.
 */
export function complement(bits: Uint8Array): Uint8Array;

/**
 * Convolutional encode with K=7 Voyager code.
 */
export function conv_k7_encode(info_bytes: Uint8Array): Uint8Array;

/**
 * Convolutional encode with K=9 NASA code.
 */
export function conv_k9_encode(info_bytes: Uint8Array): Uint8Array;

/**
 * Compress 2-bit packed DNA using order-1 context modeling.
 * Each base (A=0, C=1, G=2, T=3) is encoded with a context-dependent
 * model where the context is the previous base.
 */
export function dna_compress_order1(packed: Uint8Array, num_bases: number): Uint8Array;

/**
 * Decompress DNA compressed by dna_compress_order1.
 */
export function dna_decompress_order1(data: Uint8Array): Uint8Array;

export function ldpc_create(n_bits: number, k_bits: number, dv: number): number;

export function ldpc_decode(handle: number, llr: Float64Array, max_iter: number): Uint8Array;

export function ldpc_encode(handle: number, info: Uint8Array): Uint8Array;

export function ldpc_free(handle: number): void;

/**
 * Pack a DNA string into 2-bit bytes (4 bases per byte, MSB-first).
 *
 * Uses SIMD v128 for the inner loop when available (WASM).
 * For WASM, we process 16 chars at a time.
 */
export function pack_dna_to_bits(dna: string): Uint8Array;

/**
 * Pack DNA with validation. Returns (packed_bytes, error_message).
 * If error_message is non-empty, the input had invalid bases.
 */
export function pack_dna_validated(dna: string): PackResult;

/**
 * Compute Rabin-Karp rolling hashes over packed DNA.
 * Returns array of (numBases - windowSize + 1) 32-bit hashes.
 */
export function rolling_hash(bits: Uint8Array, window_size: number): Uint32Array;

export function rs_create(n: number, k: number): number;

export function rs_decode(handle: number, recv: Uint8Array): Uint8Array;

export function rs_decode_erasures(handle: number, recv: Uint8Array, erase_pos: Uint32Array): Uint8Array;

export function rs_encode(handle: number, data: Uint8Array): Uint8Array;

export function rs_free(handle: number): void;

/**
 * Simulate multiple oligos through the wetlab pipeline.
 *
 * Each oligo is processed independently with a unique RNG seed
 * derived from the base seed + oligo index.
 *
 * Returns a flat array: for each oligo, [length_u32_le, base1, base2, ...].
 */
export function simulate_batch(oligos: Uint8Array, oligo_offsets: Uint32Array, synth_sub: number, synth_ins: number, synth_del: number, seq_sub: number, seq_ins: number, seq_del: number, base_seed: number): Uint8Array;

/**
 * Simulate a single oligo. Returns the read as ASCII bytes.
 */
export function simulate_single(oligo: Uint8Array, synth_sub: number, synth_ins: number, synth_del: number, seq_sub: number, seq_ins: number, seq_del: number, seed: number): Uint8Array;

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
 */
export function unpack_batch(packed_data: Uint8Array, offsets: Uint32Array, num_bases: Uint32Array, total_bases: number): Uint8Array;

/**
 * Unpack 2-bit bytes into ASCII bytes (A=65, C=67, G=71, T=84).
 *
 * This is the hot path — processes 16 packed bytes → 64 ASCII bytes
 * per SIMD iteration using v128.load/store + i8x16 operations.
 *
 * Returns a Uint8Array of ASCII bytes. The JS wrapper converts to string.
 */
export function unpack_bits_to_ascii(packed: Uint8Array, num_bases: number): Uint8Array;

/**
 * Unpack 2-bit bytes directly into a DNA string.
 */
export function unpack_bits_to_dna(packed: Uint8Array, num_bases: number): string;

/**
 * Package version
 */
export function version(): string;

/**
 * K=7 (Voyager) Indel-Tolerant Viterbi decode (byte-oriented).
 */
export function viterbi_k7_decode(received_bytes: Uint8Array, llr_f64: Uint8Array, num_info_bits: number, max_drift: number, insertion_penalty_x10: number, deletion_penalty_x10: number): Uint8Array;

/**
 * K=9 Indel-Tolerant Viterbi decode (byte-oriented).
 *
 * # Arguments
 * * `received_bytes` — received bytes (hard decisions)
 * * `llr_bytes` — packed LLR values as IEEE 754 f64 little-endian bytes (8 bytes per LLR, one per bit)
 * * `num_info_bits` — number of information bits
 * * `max_drift` — maximum drift to track (default 15)
 * * `insertion_penalty_x10` — insertion penalty × 10 (e.g. 15 for 1.5)
 * * `deletion_penalty_x10` — deletion penalty × 10 (e.g. 10 for 1.0)
 *
 * # Returns
 * Decoded bytes
 */
export function viterbi_k9_decode(received_bytes: Uint8Array, llr_f64: Uint8Array, num_info_bits: number, max_drift: number, insertion_penalty_x10: number, deletion_penalty_x10: number): Uint8Array;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_get_sequencingprofile_del_rate: (a: number) => number;
    readonly __wbg_get_sequencingprofile_ins_rate: (a: number) => number;
    readonly __wbg_get_sequencingprofile_sub_rate: (a: number) => number;
    readonly __wbg_get_synthesisprofile_five_prime_bias: (a: number) => number;
    readonly __wbg_get_synthesisprofile_three_prime_bias: (a: number) => number;
    readonly __wbg_packresult_free: (a: number, b: number) => void;
    readonly __wbg_rsdecoderesult_free: (a: number, b: number) => void;
    readonly __wbg_sequencingprofile_free: (a: number, b: number) => void;
    readonly __wbg_set_sequencingprofile_del_rate: (a: number, b: number) => void;
    readonly __wbg_set_sequencingprofile_ins_rate: (a: number, b: number) => void;
    readonly __wbg_set_sequencingprofile_sub_rate: (a: number, b: number) => void;
    readonly __wbg_set_synthesisprofile_five_prime_bias: (a: number, b: number) => void;
    readonly __wbg_set_synthesisprofile_three_prime_bias: (a: number, b: number) => void;
    readonly __wbg_synthesisprofile_free: (a: number, b: number) => void;
    readonly arith_compress: (a: number, b: number, c: number) => void;
    readonly arith_decompress: (a: number, b: number, c: number) => void;
    readonly bhe_decode_fsm: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly bhe_decode_k1: (a: number, b: number, c: number, d: number) => void;
    readonly bhe_encode_fsm: (a: number, b: number, c: number, d: number) => void;
    readonly bhe_encode_k1: (a: number, b: number, c: number) => void;
    readonly bhe_validate: (a: number, b: number, c: number) => number;
    readonly bit_parallel_hamming: (a: number, b: number, c: number, d: number) => number;
    readonly bit_parallel_match: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly complement: (a: number, b: number, c: number) => void;
    readonly conv_k7_encode: (a: number, b: number, c: number) => void;
    readonly conv_k9_encode: (a: number, b: number, c: number) => void;
    readonly dna_compress_order1: (a: number, b: number, c: number, d: number) => void;
    readonly dna_decompress_order1: (a: number, b: number, c: number) => void;
    readonly ldpc_create: (a: number, b: number, c: number) => number;
    readonly ldpc_decode: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly ldpc_encode: (a: number, b: number, c: number, d: number) => void;
    readonly pack_dna_to_bits: (a: number, b: number, c: number) => void;
    readonly pack_dna_validated: (a: number, b: number) => number;
    readonly packresult_data: (a: number, b: number) => void;
    readonly packresult_error: (a: number, b: number) => void;
    readonly rolling_hash: (a: number, b: number, c: number, d: number) => void;
    readonly rs_create: (a: number, b: number) => number;
    readonly rs_decode: (a: number, b: number, c: number, d: number) => void;
    readonly rs_decode_erasures: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
    readonly rs_encode: (a: number, b: number, c: number, d: number) => void;
    readonly rs_free: (a: number) => void;
    readonly rsdecoderesult_corrected: (a: number) => number;
    readonly rsdecoderesult_data: (a: number, b: number) => void;
    readonly rsdecoderesult_erased: (a: number) => number;
    readonly sequencingprofile_illumina: () => number;
    readonly sequencingprofile_nanopore: () => number;
    readonly sequencingprofile_pacbio_hifi: () => number;
    readonly simulate_batch: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number) => void;
    readonly simulate_single: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number) => void;
    readonly synthesisprofile_illuminina: () => number;
    readonly synthesisprofile_nanopore: () => number;
    readonly synthesisprofile_new: () => number;
    readonly unpack_batch: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => void;
    readonly unpack_bits_to_ascii: (a: number, b: number, c: number, d: number) => void;
    readonly unpack_bits_to_dna: (a: number, b: number, c: number, d: number) => void;
    readonly version: (a: number) => void;
    readonly viterbi_k7_decode: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => void;
    readonly viterbi_k9_decode: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number) => void;
    readonly __wbg_get_synthesisprofile_del_rate: (a: number) => number;
    readonly __wbg_get_synthesisprofile_ins_rate: (a: number) => number;
    readonly __wbg_get_synthesisprofile_sub_rate: (a: number) => number;
    readonly __wbg_set_synthesisprofile_del_rate: (a: number, b: number) => void;
    readonly __wbg_set_synthesisprofile_ins_rate: (a: number, b: number) => void;
    readonly __wbg_set_synthesisprofile_sub_rate: (a: number, b: number) => void;
    readonly ldpc_free: (a: number) => void;
    readonly sequencingprofile_new: () => number;
    readonly __wbindgen_add_to_stack_pointer: (a: number) => number;
    readonly __wbindgen_export: (a: number, b: number) => number;
    readonly __wbindgen_export2: (a: number, b: number, c: number) => void;
    readonly __wbindgen_export3: (a: number, b: number, c: number, d: number) => number;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
