/**
 * SIMD-accelerated DNA 2-bit unpacking module — V4 DOUBLE-PUMP.
 *
 * Key optimizations over V1:
 *   1. wasm_i8x16_swizzle for 2-bit→ASCII (replaces 19 ops with 1)
 *   2. Loop-invariant constant loads hoisted
 *   3. wasm_v8x16_shuffle for in-WASM interleaving
 *   4. Double-pump: 2 v128 loads per iteration (32 packed → 128 output bytes)
 *      for better instruction-level parallelism
 *   5. Use wasm_u8x16_shr (unsigned shift) to avoid sign-extension issues
 */

#include <stdint.h>
#include <string.h>
#include <emscripten.h>
#include <wasm_simd128.h>

static uint8_t LUT[4];

static const uint8_t MASK_03[16] __attribute__((aligned(16))) = {
    0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03,
    0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03
};

static const uint8_t LUT_VEC[16] __attribute__((aligned(16))) = {
    0x41, 0x43, 0x47, 0x54, 0x41, 0x43, 0x47, 0x54,
    0x41, 0x43, 0x47, 0x54, 0x41, 0x43, 0x47, 0x54
};

EMSCRIPTEN_KEEPALIVE
void init_lut(void) {
    LUT[0] = 0x41; LUT[1] = 0x43; LUT[2] = 0x47; LUT[3] = 0x54;
}

/**
 * Helper: interleave 4 position vectors into 4 output chunks.
 * Each chunk is 16 bytes: [p0[i],p1[i],p2[i],p3[i]] for i in a group of 4.
 */
#if defined(__wasm_simd128__)
static inline void interleave_and_store(
    v128_t pos0, v128_t pos1, v128_t pos2, v128_t pos3,
    uint8_t* out_ptr, int base) {

    /* Chunk 0: lanes 0-3 from each position */
    v128_t p01_0 = wasm_v8x16_shuffle(pos0, pos1,
        0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23);
    v128_t p23_0 = wasm_v8x16_shuffle(pos2, pos3,
        0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23);
    v128_t chunk0 = wasm_v8x16_shuffle(p01_0, p23_0,
        0, 1, 16, 17, 2, 3, 18, 19, 4, 5, 20, 21, 6, 7, 22, 23);
    wasm_v128_store((v128_t*)(out_ptr + base), chunk0);

    /* Chunk 1: lanes 4-7 */
    v128_t p01_1 = wasm_v8x16_shuffle(pos0, pos1,
        4, 20, 5, 21, 6, 22, 7, 23, 8, 24, 9, 25, 10, 26, 11, 27);
    v128_t p23_1 = wasm_v8x16_shuffle(pos2, pos3,
        4, 20, 5, 21, 6, 22, 7, 23, 8, 24, 9, 25, 10, 26, 11, 27);
    v128_t chunk1 = wasm_v8x16_shuffle(p01_1, p23_1,
        0, 1, 16, 17, 2, 3, 18, 19, 4, 5, 20, 21, 6, 7, 22, 23);
    wasm_v128_store((v128_t*)(out_ptr + base + 16), chunk1);

    /* Chunk 2: lanes 8-11 */
    v128_t p01_2 = wasm_v8x16_shuffle(pos0, pos1,
        8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31);
    v128_t p23_2 = wasm_v8x16_shuffle(pos2, pos3,
        8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31);
    v128_t chunk2 = wasm_v8x16_shuffle(p01_2, p23_2,
        0, 1, 16, 17, 2, 3, 18, 19, 4, 5, 20, 21, 6, 7, 22, 23);
    wasm_v128_store((v128_t*)(out_ptr + base + 32), chunk2);

    /* Chunk 3: lanes 12-15 */
    v128_t p01_3 = wasm_v8x16_shuffle(pos0, pos1,
        12, 28, 13, 29, 14, 30, 15, 31, 12, 28, 13, 29, 14, 30, 15, 31);
    v128_t p23_3 = wasm_v8x16_shuffle(pos2, pos3,
        12, 28, 13, 29, 14, 30, 15, 31, 12, 28, 13, 29, 14, 30, 15, 31);
    v128_t chunk3 = wasm_v8x16_shuffle(p01_3, p23_3,
        0, 1, 16, 17, 2, 3, 18, 19, 4, 5, 20, 21, 6, 7, 22, 23);
    wasm_v128_store((v128_t*)(out_ptr + base + 48), chunk3);
}
#endif

EMSCRIPTEN_KEEPALIVE
void unpack_simd(const uint8_t* ptr, uint8_t* out_ptr, int num_bytes) {
#if defined(__wasm_simd128__)
    v128_t mask03 = wasm_v128_load(MASK_03);
    v128_t lut    = wasm_v128_load(LUT_VEC);

    int i = 0;
    for (; i + 16 <= num_bytes; i += 16) {
        v128_t packed = wasm_v128_load((const v128_t*)(ptr + i));

        v128_t pos0 = wasm_i8x16_swizzle(lut, wasm_u8x16_shr(packed, 6));
        v128_t pos1 = wasm_i8x16_swizzle(lut, wasm_v128_and(wasm_u8x16_shr(packed, 4), mask03));
        v128_t pos2 = wasm_i8x16_swizzle(lut, wasm_v128_and(wasm_u8x16_shr(packed, 2), mask03));
        v128_t pos3 = wasm_i8x16_swizzle(lut, wasm_v128_and(packed, mask03));

        int base = i * 4;
        wasm_v128_store((v128_t*)(out_ptr + base),      pos0);
        wasm_v128_store((v128_t*)(out_ptr + base + 16), pos1);
        wasm_v128_store((v128_t*)(out_ptr + base + 32), pos2);
        wasm_v128_store((v128_t*)(out_ptr + base + 48), pos3);
    }
#endif

    for (int j =
#if defined(__wasm_simd128__)
        i
#else
        0
#endif
        ; j < num_bytes; j++) {
        uint8_t byte = ptr[j];
        out_ptr[j * 4 + 0] = LUT[(byte >> 6) & 0x03];
        out_ptr[j * 4 + 1] = LUT[(byte >> 4) & 0x03];
        out_ptr[j * 4 + 2] = LUT[(byte >> 2) & 0x03];
        out_ptr[j * 4 + 3] = LUT[byte & 0x03];
    }
}

EMSCRIPTEN_KEEPALIVE
void unpack_scalar(const uint8_t* ptr, uint8_t* out_ptr, int num_bytes) {
    for (int i = 0; i < num_bytes; i++) {
        uint8_t byte = ptr[i];
        out_ptr[i * 4 + 0] = LUT[(byte >> 6) & 0x03];
        out_ptr[i * 4 + 1] = LUT[(byte >> 4) & 0x03];
        out_ptr[i * 4 + 2] = LUT[(byte >> 2) & 0x03];
        out_ptr[i * 4 + 3] = LUT[byte & 0x03];
    }
}

EMSCRIPTEN_KEEPALIVE
void unpack_simd_interleaved(const uint8_t* ptr, uint8_t* out_ptr, int num_bytes) {
#if defined(__wasm_simd128__)
    v128_t mask03 = wasm_v128_load(MASK_03);
    v128_t lut    = wasm_v128_load(LUT_VEC);

    int i = 0;

    /* Double-pump: process 32 packed bytes → 128 output bytes per iteration.
     * This improves instruction-level parallelism by allowing the CPU
     * to execute operations on block B while block A's shuffles complete. */
    for (; i + 32 <= num_bytes; i += 32) {
        /* Block A: bytes i..i+15 */
        v128_t packed_a = wasm_v128_load((const v128_t*)(ptr + i));
        v128_t pos0_a = wasm_i8x16_swizzle(lut, wasm_u8x16_shr(packed_a, 6));
        v128_t pos1_a = wasm_i8x16_swizzle(lut, wasm_v128_and(wasm_u8x16_shr(packed_a, 4), mask03));
        v128_t pos2_a = wasm_i8x16_swizzle(lut, wasm_v128_and(wasm_u8x16_shr(packed_a, 2), mask03));
        v128_t pos3_a = wasm_i8x16_swizzle(lut, wasm_v128_and(packed_a, mask03));

        /* Block B: bytes i+16..i+31 */
        v128_t packed_b = wasm_v128_load((const v128_t*)(ptr + i + 16));
        v128_t pos0_b = wasm_i8x16_swizzle(lut, wasm_u8x16_shr(packed_b, 6));
        v128_t pos1_b = wasm_i8x16_swizzle(lut, wasm_v128_and(wasm_u8x16_shr(packed_b, 4), mask03));
        v128_t pos2_b = wasm_i8x16_swizzle(lut, wasm_v128_and(wasm_u8x16_shr(packed_b, 2), mask03));
        v128_t pos3_b = wasm_i8x16_swizzle(lut, wasm_v128_and(packed_b, mask03));

        /* Interleave and store block A (64 bytes) */
        interleave_and_store(pos0_a, pos1_a, pos2_a, pos3_a, out_ptr, i * 4);

        /* Interleave and store block B (64 bytes) */
        interleave_and_store(pos0_b, pos1_b, pos2_b, pos3_b, out_ptr, (i + 16) * 4);
    }

    /* Single-pump for remaining full 16-byte blocks */
    for (; i + 16 <= num_bytes; i += 16) {
        v128_t packed = wasm_v128_load((const v128_t*)(ptr + i));
        v128_t pos0 = wasm_i8x16_swizzle(lut, wasm_u8x16_shr(packed, 6));
        v128_t pos1 = wasm_i8x16_swizzle(lut, wasm_v128_and(wasm_u8x16_shr(packed, 4), mask03));
        v128_t pos2 = wasm_i8x16_swizzle(lut, wasm_v128_and(wasm_u8x16_shr(packed, 2), mask03));
        v128_t pos3 = wasm_i8x16_swizzle(lut, wasm_v128_and(packed, mask03));
        interleave_and_store(pos0, pos1, pos2, pos3, out_ptr, i * 4);
    }
#endif

    for (int j =
#if defined(__wasm_simd128__)
        i
#else
        0
#endif
        ; j < num_bytes; j++) {
        uint8_t byte = ptr[j];
        out_ptr[j * 4 + 0] = LUT[(byte >> 6) & 0x03];
        out_ptr[j * 4 + 1] = LUT[(byte >> 4) & 0x03];
        out_ptr[j * 4 + 2] = LUT[(byte >> 2) & 0x03];
        out_ptr[j * 4 + 3] = LUT[byte & 0x03];
    }
}
