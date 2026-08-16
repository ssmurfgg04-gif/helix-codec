/**
 * DNA-Diff — LZ77 delta compression + context prediction.
 */
#include "arithmetic_coder.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

namespace dnadiff {

static uint8_t char_to_2bit(uint8_t c) {
    switch (c) { case 'A': case 'a': return 0; case 'C': case 'c': return 1; case 'G': case 'g': return 2; case 'T': case 't': return 3; default: return 255; }
}
static uint8_t BIT2CHAR[4] = {'A', 'C', 'G', 'T'};
static const int MIN_MATCH = 3;
static const int WINDOW_SIZE = 16384;
static const int MAX_MATCH = 64;

static int find_match(const uint8_t* data, int pos, int len, int* outOff) {
    int best = 0, bestOff = 0;
    int start = pos > WINDOW_SIZE ? pos - WINDOW_SIZE : 0;
    for (int j = start; j < pos; j++) {
        int l = 0;
        while (l < MAX_MATCH && pos + l < len && data[j + l] == data[pos + l]) l++;
        if (l > best) { best = l; bestOff = pos - j; }
    }
    *outOff = bestOff;
    return best;
}

int compress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len == 0) return 0;
    if (output_cap < 8) return -1;
    output[0]='D'; output[1]='F'; output[2]='D'; output[3]=0x02;
    output[4]=(input_len>>0)&0xFF; output[5]=(input_len>>8)&0xFF;
    output[6]=(input_len>>16)&0xFF; output[7]=(input_len>>24)&0xFF;

    arith::Enc enc; enc.init(output + 8, output_cap - 8);

    uint8_t* packed = (uint8_t*)malloc(input_len);
    if (!packed) return -1;
    for (int i = 0; i < input_len; i++) { packed[i] = char_to_2bit(input[i]); if (packed[i]==255) packed[i]=0; }

    int i = 0;
    while (i < input_len) {
        int off;
        int ml = find_match(packed, i, input_len, &off);
        if (ml >= MIN_MATCH) {
            enc.put_bit(1);  // match flag
            enc.put_bits(off - 1, 14);  // offset in 14 bits (up to 16K)
            enc.put_bits(ml - MIN_MATCH, 6);  // length in 6 bits (up to 64+3)
            i += ml;
        } else {
            enc.put_bit(0);  // literal flag
            enc.put_bits(packed[i], 2);  // 2-bit literal
            i++;
        }
    }

    free(packed);
    return 8 + enc.fin();
}

int decompress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len < 8) return -1;
    if (input[0]!='D'||input[1]!='F'||input[2]!='D'||input[3]!=0x02) return -1;
    uint32_t orig_len = (uint32_t)input[4]|((uint32_t)input[5]<<8)|((uint32_t)input[6]<<16)|((uint32_t)input[7]<<24);
    if ((int)orig_len > output_cap) return -1;

    arith::Dec dec; dec.init(input + 8, input_len - 8);
    uint8_t* packed = (uint8_t*)malloc(orig_len);
    if (!packed) return -1;

    int out = 0;
    while (out < (int)orig_len) {
        int flag = dec.read_bit();
        if (flag == 0) {
            packed[out++] = (uint8_t)dec.read_bits(2);
        } else {
            uint32_t off = dec.read_bits(14) + 1;
            uint32_t len = dec.read_bits(6) + MIN_MATCH;
            if (off > (uint32_t)out) { free(packed); return -1; }
            for (uint32_t j = 0; j < len && out < (int)orig_len; j++) {
                packed[out] = packed[out - off]; out++;
            }
        }
    }

    for (int i = 0; i < out; i++) output[i] = BIT2CHAR[packed[i]];
    free(packed);
    return out;
}

}  // namespace dnadiff
