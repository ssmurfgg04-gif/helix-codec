/**
 * DNA-QLS — Order-2 context prediction compression.
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
static const int NUM_CTX = 16;

struct Ctx2Model {
    uint32_t count[NUM_CTX][4];
    void init() { memset(count, 0, sizeof(count)); }
    uint32_t predict(int ctx) {
        uint32_t best = 0, bc = count[ctx][0];
        for (uint32_t i = 1; i < 4; i++) if (count[ctx][i] > bc) { best = i; bc = count[ctx][i]; }
        return best;
    }
    void update(int ctx, uint32_t sym) { count[ctx][sym]++; }
};

int compress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len == 0) return 0;
    if (output_cap < 12) return -1;
    output[0]='Q'; output[1]='L'; output[2]='S'; output[3]=0x02;
    output[4]=(input_len>>0)&0xFF; output[5]=(input_len>>8)&0xFF;
    output[6]=(input_len>>16)&0xFF; output[7]=(input_len>>24)&0xFF;
    int packed_len = 0;
    for (int i = 0; i < input_len; i++) if (char_to_2bit(input[i]) != 255) packed_len++;
    output[8]=(packed_len>>0)&0xFF; output[9]=(packed_len>>8)&0xFF;
    output[10]=(packed_len>>16)&0xFF; output[11]=(packed_len>>24)&0xFF;

    arith::Enc enc; enc.init(output + 12, output_cap - 12);
    Ctx2Model model; model.init();
    int ctx = 0;

    for (int i = 0; i < input_len; i++) {
        uint8_t sym = char_to_2bit(input[i]);
        if (sym == 255) continue;
        uint32_t pred = model.predict(ctx);
        if (sym == pred) { enc.put_bit(0); }
        else { enc.put_bit(1); enc.put_bits(sym, 2); }
        model.update(ctx, sym);
        ctx = ((ctx << 2) | sym) & (NUM_CTX - 1);
    }

    return 12 + enc.fin();
}

int decompress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len < 12) return -1;
    if (input[0]!='Q'||input[1]!='L'||input[2]!='S'||input[3]!=0x02) return -1;
    uint32_t orig_len = (uint32_t)input[4]|((uint32_t)input[5]<<8)|((uint32_t)input[6]<<16)|((uint32_t)input[7]<<24);
    uint32_t packed_len = (uint32_t)input[8]|((uint32_t)input[9]<<8)|((uint32_t)input[10]<<16)|((uint32_t)input[11]<<24);
    if ((int)packed_len > output_cap) return -1;

    arith::Dec dec; dec.init(input + 12, input_len - 12);
    Ctx2Model model; model.init();
    int ctx = 0;

    for (uint32_t i = 0; i < packed_len; i++) {
        uint32_t pred = model.predict(ctx);
        int flag = dec.read_bit();
        uint32_t sym = (flag == 0) ? pred : dec.read_bits(2);
        if (sym > 3) return -1;
        output[i] = BIT2CHAR[sym];
        model.update(ctx, sym);
        ctx = ((ctx << 2) | sym) & (NUM_CTX - 1);
    }

    return (int)packed_len;
}

}  // namespace dnaqls
