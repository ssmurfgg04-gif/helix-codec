/**
 * Unified API for DNA compressors compiled to WASM via Emscripten.
 *
 * Exposes C functions for:
 *   - GeCo2 (algorithm 0): Order-3 context model + arithmetic coding
 *   - DNA-Diff (algorithm 1): LZ77 delta compression + arithmetic coding
 *   - DNA-QLS (algorithm 2): Order-2 context model + quasi-lexicographic sorting
 *
 * Usage from JS:
 *   const len = _dna_compress(algo, inputPtr, inputLen, outputPtr, outputCap);
 *   const len = _dna_decompress(algo, inputPtr, inputLen, outputPtr, outputCap);
 *   const name = UTF8ToString(_dna_compressor_name(algo));
 *   const count = _dna_compressor_count();
 */

#include "geco.cpp"
#include "dna_diff.cpp"
#include "dna_qls.cpp"

extern "C" {

/**
 * Compress DNA using specified algorithm.
 * @param algo   Algorithm: 0=GeCo, 1=DNA-Diff, 2=DNA-QLS
 * @param input  Pointer to input DNA ASCII bytes
 * @param input_len Length of input
 * @param output Pointer to output buffer
 * @param output_cap Capacity of output buffer
 * @return Number of bytes written, or -1 on error
 */
int dna_compress(int algo, const uint8_t* input, int input_len,
                 uint8_t* output, int output_cap) {
    switch (algo) {
        case 0: return geco::compress(input, input_len, output, output_cap);
        case 1: return dnadiff::compress(input, input_len, output, output_cap);
        case 2: return dnaqls::compress(input, input_len, output, output_cap);
        default: return -1;
    }
}

/**
 * Decompress DNA using specified algorithm.
 * @param algo   Algorithm: 0=GeCo, 1=DNA-Diff, 2=DNA-QLS
 * @param input  Pointer to compressed data
 * @param input_len Length of compressed data
 * @param output Pointer to output buffer
 * @param output_cap Capacity of output buffer
 * @return Number of bytes written, or -1 on error
 */
int dna_decompress(int algo, const uint8_t* input, int input_len,
                   uint8_t* output, int output_cap) {
    switch (algo) {
        case 0: return geco::decompress(input, input_len, output, output_cap);
        case 1: return dnadiff::decompress(input, input_len, output, output_cap);
        case 2: return dnaqls::decompress(input, input_len, output, output_cap);
        default: return -1;
    }
}

/**
 * Get compressor name.
 * @param algo Algorithm index
 * @return Static string with algorithm name
 */
const char* dna_compressor_name(int algo) {
    switch (algo) {
        case 0: return "GeCo2";
        case 1: return "DNA-Diff";
        case 2: return "DNA-QLS";
        default: return "unknown";
    }
}

/**
 * Get number of available compressors.
 * @return 3
 */
int dna_compressor_count() {
    return 3;
}

}  // extern "C"
