/**
 * DNA-Diff — LZ77 delta compression + context prediction + arithmetic coding.
 *
 * REWRITE: Replaced VLC with adaptive arithmetic coding for the literal
 * and match tokens. The LZ77 tokens are now entropy-coded using the
 * range coder, achieving better compression for the offset/length fields
 * which have non-uniform distributions (short offsets are more common).
 *
 * Token format (arithmetic-coded):
 *   - Flag token: 0 = literal, 1 = match (2-symbol alphabet)
 *   - Literal token: 4-symbol alphabet {A,C,G,T}
 *   - Match offset: 16-symbol alphabet for the 14-bit offset (grouped)
 *   - Match length: 8-symbol alphabet for the 6-bit length (grouped)
 *
 * The LZ77 engine remains the same — only the entropy coding layer changed.
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

    // Header: [magic 4B] [orig_len 4B LE]
    output[0]='D'; output[1]='F'; output[2]='D'; output[3]=0x03; // v3 = arithmetic
    output[4]=(input_len>>0)&0xFF; output[5]=(input_len>>8)&0xFF;
    output[6]=(input_len>>16)&0xFF; output[7]=(input_len>>24)&0xFF;

    // Use range coder for all tokens
    arith::BitOut bout;
    bout.init(output + 8, output_cap - 8);

    arith::RangeEncoder enc;
    enc.init(&bout);

    // Models for different token types
    arith::AdaptiveModel flagModel;   // 2 symbols: literal(0), match(1)
    arith::AdaptiveModel litModel;    // 4 symbols: A,C,G,T
    arith::AdaptiveModel offModel;    // 16 symbols for offset (grouped 4-bit chunks)
    arith::AdaptiveModel lenModel;    // 8 symbols for length (grouped 3-bit chunks)
    flagModel.init(2);
    litModel.init(4);
    offModel.init(16);
    lenModel.init(8);

    uint8_t* packed = (uint8_t*)malloc(input_len);
    if (!packed) return -1;
    for (int i = 0; i < input_len; i++) { packed[i] = char_to_2bit(input[i]); if (packed[i]==255) packed[i]=0; }

    int i = 0;
    while (i < input_len) {
        int off;
        int ml = find_match(packed, i, input_len, &off);
        if (ml >= MIN_MATCH) {
            // Encode match flag
            enc.encode(flagModel.getCumFreq(1), flagModel.getFreq(1), flagModel.getTotal());
            flagModel.update(1);

            // Encode offset in 4-bit groups (14 bits = 4+4+4+2)
            // Group 1: bits 10-13 (most significant)
            uint32_t g1 = (off - 1) >> 10;
            enc.encode(offModel.getCumFreq(g1), offModel.getFreq(g1), offModel.getTotal());
            offModel.update(g1);
            // Group 2: bits 6-9
            uint32_t g2 = ((off - 1) >> 6) & 0xF;
            enc.encode(offModel.getCumFreq(g2), offModel.getFreq(g2), offModel.getTotal());
            offModel.update(g2);
            // Group 3: bits 2-5
            uint32_t g3 = ((off - 1) >> 2) & 0xF;
            enc.encode(offModel.getCumFreq(g3), offModel.getFreq(g3), offModel.getTotal());
            offModel.update(g3);
            // Group 4: bits 0-1 (encoded as 4-symbol for simplicity)
            uint32_t g4 = (off - 1) & 0x3;
            enc.encode(litModel.getCumFreq(g4), litModel.getFreq(g4), litModel.getTotal());
            litModel.update(g4);

            // Encode length in 3-bit groups (6 bits = 3+3)
            uint32_t l1 = (ml - MIN_MATCH) >> 3;
            enc.encode(lenModel.getCumFreq(l1), lenModel.getFreq(l1), lenModel.getTotal());
            lenModel.update(l1);
            uint32_t l2 = (ml - MIN_MATCH) & 0x7;
            enc.encode(lenModel.getCumFreq(l2), lenModel.getFreq(l2), lenModel.getTotal());
            lenModel.update(l2);

            i += ml;
        } else {
            // Encode literal flag
            enc.encode(flagModel.getCumFreq(0), flagModel.getFreq(0), flagModel.getTotal());
            flagModel.update(0);

            // Encode literal symbol
            enc.encode(litModel.getCumFreq(packed[i]), litModel.getFreq(packed[i]), litModel.getTotal());
            litModel.update(packed[i]);
            i++;
        }
    }

    free(packed);
    enc.finish();
    return 8 + bout.fin();
}

int decompress(const uint8_t* input, int input_len, uint8_t* output, int output_cap) {
    if (input_len < 8) return -1;
    if (input[0]!='D'||input[1]!='F'||input[2]!='D') return -1;
    int version = input[3];

    uint32_t orig_len = (uint32_t)input[4]|((uint32_t)input[5]<<8)|((uint32_t)input[6]<<16)|((uint32_t)input[7]<<24);
    if ((int)orig_len > output_cap) return -1;

    if (version == 0x03) {
        // v3: arithmetic coding
        arith::BitIn bin;
        bin.init(input + 8, input_len - 8);

        arith::RangeDecoder dec;
        dec.init(&bin);

        arith::AdaptiveModel flagModel, litModel, offModel, lenModel;
        flagModel.init(2);
        litModel.init(4);
        offModel.init(16);
        lenModel.init(8);

        uint8_t* packed = (uint8_t*)malloc(orig_len);
        if (!packed) return -1;

        int out = 0;
        while (out < (int)orig_len) {
            // Decode flag
            uint32_t target = dec.decode_target(flagModel.getTotal());
            uint32_t flag = flagModel.findSymbol(target);
            dec.decode(flagModel.getCumFreq(flag), flagModel.getFreq(flag), flagModel.getTotal());
            flagModel.update(flag);

            if (flag == 0) {
                // Literal
                target = dec.decode_target(litModel.getTotal());
                uint32_t sym = litModel.findSymbol(target);
                dec.decode(litModel.getCumFreq(sym), litModel.getFreq(sym), litModel.getTotal());
                litModel.update(sym);
                packed[out++] = (uint8_t)sym;
            } else {
                // Match: decode offset in 4-bit groups
                target = dec.decode_target(offModel.getTotal());
                uint32_t g1 = offModel.findSymbol(target);
                dec.decode(offModel.getCumFreq(g1), offModel.getFreq(g1), offModel.getTotal());
                offModel.update(g1);

                target = dec.decode_target(offModel.getTotal());
                uint32_t g2 = offModel.findSymbol(target);
                dec.decode(offModel.getCumFreq(g2), offModel.getFreq(g2), offModel.getTotal());
                offModel.update(g2);

                target = dec.decode_target(offModel.getTotal());
                uint32_t g3 = offModel.findSymbol(target);
                dec.decode(offModel.getCumFreq(g3), offModel.getFreq(g3), offModel.getTotal());
                offModel.update(g3);

                target = dec.decode_target(litModel.getTotal());
                uint32_t g4 = litModel.findSymbol(target);
                dec.decode(litModel.getCumFreq(g4), litModel.getFreq(g4), litModel.getTotal());
                litModel.update(g4);

                uint32_t off = (g1 << 10) | (g2 << 6) | (g3 << 2) | g4;
                off += 1;

                // Decode length in 3-bit groups
                target = dec.decode_target(lenModel.getTotal());
                uint32_t l1 = lenModel.findSymbol(target);
                dec.decode(lenModel.getCumFreq(l1), lenModel.getFreq(l1), lenModel.getTotal());
                lenModel.update(l1);

                target = dec.decode_target(lenModel.getTotal());
                uint32_t l2 = lenModel.findSymbol(target);
                dec.decode(lenModel.getCumFreq(l2), lenModel.getFreq(l2), lenModel.getTotal());
                lenModel.update(l2);

                uint32_t len = (l1 << 3) | l2;
                len += MIN_MATCH;

                if (off > (uint32_t)out) { free(packed); return -1; }
                for (uint32_t j = 0; j < len && out < (int)orig_len; j++) {
                    packed[out] = packed[out - off]; out++;
                }
            }
        }

        for (int i = 0; i < out; i++) output[i] = BIT2CHAR[packed[i]];
        free(packed);
        return out;

    } else if (version == 0x02) {
        // v2: backward-compatible VLC decoding
        arith::Dec dec;
        dec.init(input + 8, input_len - 8);

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

    return -1;
}

}  // namespace dnadiff
