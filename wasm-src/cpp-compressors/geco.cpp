/**
 * GeCo2 — DNA compression using order-3 context prediction + VLC.
 */
#include "arithmetic_coder.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

namespace geco {

static uint8_t char_to_2bit(uint8_t c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 255;
    }
}
static uint8_t BIT2CHAR[4] = {'A', 'C', 'G', 'T'};
static const int NUM_CTX = 64;  // 4^3

// Context model: predict next symbol based on previous 3 symbols
struct CtxModel {
    uint32_t count[NUM_CTX][4];  // count[ctx][sym] = frequency

    void init() { memset(count, 0, sizeof(count)); }

    // Predict most likely symbol given context
    uint32_t predict(int ctx) {
        uint32_t best = 0, best_count = count[ctx][0];
        for (uint32_t i = 1; i < 4; i++) {
            if (count[ctx][i] > best_count) { best = i; best_count = count[ctx][i]; }
        }
        return best;
    }

    void update(int ctx, uint32_t sym) { count[ctx][sym]++; }
};

int compress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len == 0) return 0;
    if (output_cap < 8) return -1;

    output[0] = 'G'; output[1] = 'E'; output[2] = 'C'; output[3] = 0x02;
    output[4] = (input_len >> 0) & 0xFF;
    output[5] = (input_len >> 8) & 0xFF;
    output[6] = (input_len >> 16) & 0xFF;
    output[7] = (input_len >> 24) & 0xFF;

    arith::Enc enc; enc.init(output + 8, output_cap - 8);
    CtxModel model; model.init();
    int ctx = 0;

    for (int i = 0; i < input_len; i++) {
        uint8_t sym = char_to_2bit(input[i]);
        if (sym == 255) continue;

        uint32_t pred = model.predict(ctx);
        if (sym == pred) {
            enc.put_bit(0);  // prediction correct: 1 bit
        } else {
            enc.put_bit(1);  // prediction wrong: 1 + 2 bits
            enc.put_bits(sym, 2);
        }
        model.update(ctx, sym);
        ctx = ((ctx << 2) | sym) & (NUM_CTX - 1);
    }

    int enc_len = enc.fin();
    return 8 + enc_len;
}

int decompress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len < 8) return -1;
    if (input[0] != 'G' || input[1] != 'E' || input[2] != 'C' || input[3] != 0x02) return -1;

    uint32_t orig_len = (uint32_t)input[4] | ((uint32_t)input[5] << 8) |
                         ((uint32_t)input[6] << 16) | ((uint32_t)input[7] << 24);
    if ((int)orig_len > output_cap) return -1;

    arith::Dec dec; dec.init(input + 8, input_len - 8);
    CtxModel model; model.init();
    int ctx = 0, out_pos = 0;

    for (uint32_t i = 0; i < orig_len; i++) {
        uint32_t pred = model.predict(ctx);
        int flag = dec.read_bit();
        uint32_t sym;
        if (flag == 0) {
            sym = pred;
        } else {
            sym = dec.read_bits(2);
        }
        if (sym > 3) return -1;
        output[out_pos++] = BIT2CHAR[sym];
        model.update(ctx, sym);
        ctx = ((ctx << 2) | sym) & (NUM_CTX - 1);
    }

    return out_pos;
}

}  // namespace geco
