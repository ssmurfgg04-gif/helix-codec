/**
 * DNA-QLS — Order-2 context model + quasi-lexicographic sorting + arithmetic coding.
 *
 * REWRITE: Replaced VLC with adaptive arithmetic coding, matching the
 * proper QLS algorithm. The quasi-lexicographic sorting reorders reads
 * by similarity before compression, improving context model accuracy.
 *
 * For now, the QLS sorting is simplified (identity permutation) but
 * the arithmetic coding is real and functional.
 */
#include "arithmetic_coder.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

namespace dnaqls {

static uint8_t char_to_2bit(uint8_t c) {
    switch (c) { case 'A': case 'a': return 0; case 'C': case 'c': return 1; case 'G': case 'g': return 2; case 'T': case 't': return 3; default: return 255; }
}
static uint8_t BIT2CHAR[4] = {'A', 'C', 'G', 'T'};
static const int NUM_CTX = 16;  // 4^2

/**
 * Order-2 context model with per-context adaptive frequency table.
 */
struct Ctx2Model {
    uint32_t count[NUM_CTX][4];
    uint32_t cum[NUM_CTX][5];
    uint32_t total[NUM_CTX];

    void init() {
        memset(count, 0, sizeof(count));
        memset(total, 0, sizeof(total));
        for (int ctx = 0; ctx < NUM_CTX; ctx++) {
            for (int s = 0; s < 4; s++) count[ctx][s] = 1;
            total[ctx] = 4;
            recompute(ctx);
        }
    }

    void recompute(int ctx) {
        cum[ctx][0] = 0;
        for (int s = 0; s < 4; s++) cum[ctx][s + 1] = cum[ctx][s] + count[ctx][s];
        total[ctx] = cum[ctx][4];
    }

    void update(int ctx, uint32_t sym) {
        count[ctx][sym]++;
        total[ctx]++;
        if (total[ctx] > arith::MAX_TOTAL) {
            total[ctx] = 0;
            for (int s = 0; s < 4; s++) {
                count[ctx][s] = (count[ctx][s] + 1) >> 1;
                if (count[ctx][s] == 0) count[ctx][s] = 1;
                total[ctx] += count[ctx][s];
            }
        }
        recompute(ctx);
    }
};

int compress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len == 0) return 0;
    if (output_cap < 12) return -1;

    output[0]='Q'; output[1]='L'; output[2]='S'; output[3]=0x03; // v3 = arithmetic
    output[4]=(input_len>>0)&0xFF; output[5]=(input_len>>8)&0xFF;
    output[6]=(input_len>>16)&0xFF; output[7]=(input_len>>24)&0xFF;
    int packed_len = 0;
    for (int i = 0; i < input_len; i++) if (char_to_2bit(input[i]) != 255) packed_len++;
    output[8]=(packed_len>>0)&0xFF; output[9]=(packed_len>>8)&0xFF;
    output[10]=(packed_len>>16)&0xFF; output[11]=(packed_len>>24)&0xFF;

    // Use range coder
    arith::BitOut bout;
    bout.init(output + 12, output_cap - 12);

    arith::RangeEncoder enc;
    enc.init(&bout);

    Ctx2Model model;
    model.init();
    int ctx = 0;

    for (int i = 0; i < input_len; i++) {
        uint8_t sym = char_to_2bit(input[i]);
        if (sym == 255) continue;

        // Arithmetic encode using the context's frequency distribution
        enc.encode(model.cum[ctx][sym], model.count[ctx][sym], model.total[ctx]);
        model.update(ctx, sym);
        ctx = ((ctx << 2) | sym) & (NUM_CTX - 1);
    }

    enc.finish();
    return 12 + bout.fin();
}

int decompress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len < 12) return -1;
    if (input[0]!='Q'||input[1]!='L'||input[2]!='S') return -1;
    int version = input[3];

    uint32_t orig_len = (uint32_t)input[4]|((uint32_t)input[5]<<8)|((uint32_t)input[6]<<16)|((uint32_t)input[7]<<24);
    uint32_t packed_len = (uint32_t)input[8]|((uint32_t)input[9]<<8)|((uint32_t)input[10]<<16)|((uint32_t)input[11]<<24);
    if ((int)packed_len > output_cap) return -1;

    if (version == 0x03) {
        // v3: arithmetic coding
        arith::BitIn bin;
        bin.init(input + 12, input_len - 12);

        arith::RangeDecoder dec;
        dec.init(&bin);

        Ctx2Model model;
        model.init();
        int ctx = 0;

        for (uint32_t i = 0; i < packed_len; i++) {
            uint32_t target = dec.decode_target(model.total[ctx]);
            uint32_t sym = 0;
            for (uint32_t s = 0; s < 4; s++) {
                if (model.cum[ctx][s + 1] > target) { sym = s; break; }
            }
            dec.decode(model.cum[ctx][sym], model.count[ctx][sym], model.total[ctx]);
            output[i] = BIT2CHAR[sym];
            model.update(ctx, sym);
            ctx = ((ctx << 2) | sym) & (NUM_CTX - 1);
        }

        return (int)packed_len;

    } else if (version == 0x02) {
        // v2: backward-compatible VLC decoding
        arith::Dec dec;
        dec.init(input + 12, input_len - 12);

        uint32_t pred_count[NUM_CTX][4];
        memset(pred_count, 0, sizeof(pred_count));
        int ctx = 0;

        for (uint32_t i = 0; i < packed_len; i++) {
            uint32_t best = 0, best_c = pred_count[ctx][0];
            for (uint32_t s = 1; s < 4; s++) {
                if (pred_count[ctx][s] > best_c) { best = s; best_c = pred_count[ctx][s]; }
            }
            int flag = dec.read_bit();
            uint32_t sym = (flag == 0) ? best : dec.read_bits(2);
            if (sym > 3) return -1;
            output[i] = BIT2CHAR[sym];
            pred_count[ctx][sym]++;
            ctx = ((ctx << 2) | sym) & (NUM_CTX - 1);
        }

        return (int)packed_len;
    }

    return -1;
}

}  // namespace dnaqls
