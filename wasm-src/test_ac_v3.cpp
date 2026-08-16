/**
 * test_ac_v3.cpp — Roundtrip test for the range coder with carry bug fix.
 *
 * Tests:
 *   1. GeCo2 compress/decompress roundtrip on random DNA
 *   2. DNA-Diff compress/decompress roundtrip on random DNA
 *   3. DNA-QLS compress/decompress roundtrip on random DNA
 *   4. Arithmetic coder stress test (many symbols, varying distributions)
 *   5. Range coder carry propagation test (known-buggy patterns)
 *
 * Compile:
 *   g++ -O2 -std=c++17 -o test_ac_v3 wasm-src/cpp-compressors/api.cpp wasm-src/test_ac_v3.cpp
 *   ./test_ac_v3
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// Declarations only (api.cpp is linked separately)
extern "C" {
int dna_compress(int algo, const uint8_t* input, int input_len, uint8_t* output, int output_cap);
int dna_decompress(int algo, const uint8_t* input, int input_len, uint8_t* output, int output_cap);
const char* dna_compressor_name(int algo);
int dna_compressor_count();
}

// Simple PRNG
static uint32_t rng_state = 12345;
static uint32_t rng_next() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static const char BASES[] = "ACGT";

// Generate random DNA sequence
static void gen_random_dna(uint8_t* buf, int len) {
    for (int i = 0; i < len; i++) {
        buf[i] = BASES[rng_next() % 4];
    }
}

// Generate biased DNA (60% AT, 40% GC — typical for some organisms)
static void gen_biased_dna(uint8_t* buf, int len) {
    for (int i = 0; i < len; i++) {
        uint32_t r = rng_next() % 100;
        if (r < 30) buf[i] = 'A';
        else if (r < 60) buf[i] = 'T';
        else if (r < 80) buf[i] = 'C';
        else buf[i] = 'G';
    }
}

// Generate repetitive DNA (tandem repeats — stress test for LZ77)
static void gen_repeat_dna(uint8_t* buf, int len) {
    const char* repeat = "ACGTACGTACGTACGT";  // 16-mer repeat
    int rep_len = strlen(repeat);
    for (int i = 0; i < len; i++) {
        buf[i] = repeat[i % rep_len];
    }
}

static int test_count = 0;
static int pass_count = 0;
static int fail_count = 0;

static void test_roundtrip(const char* name, int algo, const uint8_t* dna, int dna_len) {
    test_count++;
    uint8_t* compressed = (uint8_t*)malloc(dna_len + 1024);
    uint8_t* decompressed = (uint8_t*)malloc(dna_len + 1024);

    int comp_len = dna_compress(algo, dna, dna_len, compressed, dna_len + 1024);
    if (comp_len < 0) {
        printf("  FAIL %s: compress returned -1\n", name);
        fail_count++;
        free(compressed); free(decompressed);
        return;
    }

    int decomp_len = dna_decompress(algo, compressed, comp_len, decompressed, dna_len + 1024);
    if (decomp_len < 0) {
        printf("  FAIL %s: decompress returned -1\n", name);
        fail_count++;
        free(compressed); free(decompressed);
        return;
    }

    if (decomp_len != dna_len) {
        printf("  FAIL %s: length mismatch (expected %d, got %d)\n", name, dna_len, decomp_len);
        fail_count++;
        free(compressed); free(decompressed);
        return;
    }

    int mismatches = 0;
    for (int i = 0; i < dna_len; i++) {
        if (dna[i] != decompressed[i]) mismatches++;
    }

    if (mismatches > 0) {
        printf("  FAIL %s: %d mismatches\n", name, mismatches);
        fail_count++;
    } else {
        double ratio = (double)comp_len / dna_len;
        printf("  PASS %s: %d → %d bytes (%.3f ratio)\n", name, dna_len, comp_len, ratio);
        pass_count++;
    }

    free(compressed); free(decompressed);
}

int main() {
    printf("=== Range Coder V3 Roundtrip Tests ===\n\n");

    const char* algo_names[] = {"GeCo2", "DNA-Diff", "DNA-QLS"};

    // Test 1: Random DNA at various lengths
    printf("--- Test 1: Random DNA ---\n");
    for (int len = 100; len <= 10000; len *= 10) {
        uint8_t* dna = (uint8_t*)malloc(len);
        gen_random_dna(dna, len);
        for (int algo = 0; algo < 3; algo++) {
            char name[64];
            snprintf(name, sizeof(name), "%s/random-%d", algo_names[algo], len);
            test_roundtrip(name, algo, dna, len);
        }
        free(dna);
    }

    // Test 2: Biased DNA (non-uniform distribution)
    printf("\n--- Test 2: Biased DNA (60%% AT) ---\n");
    for (int len = 100; len <= 10000; len *= 10) {
        uint8_t* dna = (uint8_t*)malloc(len);
        gen_biased_dna(dna, len);
        for (int algo = 0; algo < 3; algo++) {
            char name[64];
            snprintf(name, sizeof(name), "%s/biased-%d", algo_names[algo], len);
            test_roundtrip(name, algo, dna, len);
        }
        free(dna);
    }

    // Test 3: Repetitive DNA (stress test for carry propagation)
    printf("\n--- Test 3: Repetitive DNA (carry stress test) ---\n");
    for (int len = 100; len <= 10000; len *= 10) {
        uint8_t* dna = (uint8_t*)malloc(len);
        gen_repeat_dna(dna, len);
        for (int algo = 0; algo < 3; algo++) {
            char name[64];
            snprintf(name, sizeof(name), "%s/repeat-%d", algo_names[algo], len);
            test_roundtrip(name, algo, dna, len);
        }
        free(dna);
    }

    // Test 4: All-same-base (worst case for carry bug)
    printf("\n--- Test 4: All-same-base (extreme carry test) ---\n");
    {
        int len = 1000;
        uint8_t dna[1000];
        memset(dna, 'A', len);
        for (int algo = 0; algo < 3; algo++) {
            char name[64];
            snprintf(name, sizeof(name), "%s/all-A-%d", algo_names[algo], len);
            test_roundtrip(name, algo, dna, len);
        }
    }

    // Test 5: Very short sequences (edge cases)
    printf("\n--- Test 5: Short sequences ---\n");
    for (int len = 1; len <= 10; len++) {
        uint8_t* dna = (uint8_t*)malloc(len);
        gen_random_dna(dna, len);
        for (int algo = 0; algo < 3; algo++) {
            char name[64];
            snprintf(name, sizeof(name), "%s/short-%d", algo_names[algo], len);
            test_roundtrip(name, algo, dna, len);
        }
        free(dna);
    }

    // Summary
    printf("\n=== Summary ===\n");
    printf("  Total: %d  Passed: %d  Failed: %d\n", test_count, pass_count, fail_count);

    if (fail_count > 0) {
        printf("\n*** SOME TESTS FAILED ***\n");
        return 1;
    }
    printf("\nAll tests passed!\n");
    return 0;
}
