#include "cpp-compressors/arithmetic_coder.h"
#include <stdio.h>
#include <string.h>

using namespace arith;

// Manual encode/decode without model to isolate the issue
int main() {
    printf("=== Minimal Range Coder Debug ===\n\n");

    // Encode symbols [0, 1, 2, 3] with uniform distribution (total=4)
    uint8_t compressed[256];
    
    // Manual encoding
    RangeEncoder enc;
    enc.init(compressed, sizeof(compressed));
    
    // Symbol 0: cf=0, f=1, t=4
    printf("Before encode(0): Low=%08x Range=%08x\n", enc.Low, enc.Range);
    enc.Encode(0, 1, 4);
    printf("After  encode(0): Low=%08x Range=%08x\n", enc.Low, enc.Range);
    
    // Symbol 1: cf=1, f=1, t=4
    enc.Encode(1, 1, 4);
    printf("After  encode(1): Low=%08x Range=%08x\n", enc.Low, enc.Range);
    
    // Symbol 2
    enc.Encode(2, 1, 4);
    printf("After  encode(2): Low=%08x Range=%08x\n", enc.Low, enc.Range);
    
    // Symbol 3
    enc.Encode(3, 1, 4);
    printf("After  encode(3): Low=%08x Range=%08x\n", enc.Low, enc.Range);
    
    int comp_len = enc.Finish();
    printf("Compressed %d bytes:", comp_len);
    for (int i = 0; i < comp_len; i++) printf(" %02x", compressed[i]);
    printf("\n\n");
    
    // Now decode
    RangeDecoder dec;
    dec.init(compressed, comp_len);
    printf("Decoder: Code=%08x Range=%08x\n", dec.Code, dec.Range);
    
    uint32_t cf = dec.GetCumFreq(4);
    printf("GetCumFreq(4) = %u (Range after = %08x)\n", cf, dec.Range);
    
    // cf should be 0 (symbol 0)
    dec.DecodeSymbol(0, 1);
    printf("After DecodeSymbol(0,1): Code=%08x Range=%08x\n", dec.Code, dec.Range);
    
    cf = dec.GetCumFreq(4);
    printf("GetCumFreq(4) = %u\n", cf);
    dec.DecodeSymbol(1, 1);
    printf("After DecodeSymbol(1,1): Code=%08x Range=%08x\n", dec.Code, dec.Range);
    
    cf = dec.GetCumFreq(4);
    printf("GetCumFreq(4) = %u\n", cf);
    dec.DecodeSymbol(2, 1);
    printf("After DecodeSymbol(2,1): Code=%08x Range=%08x\n", dec.Code, dec.Range);
    
    cf = dec.GetCumFreq(4);
    printf("GetCumFreq(4) = %u\n", cf);
    dec.DecodeSymbol(3, 1);
    printf("After DecodeSymbol(3,1): Code=%08x Range=%08x\n", dec.Code, dec.Range);
    
    return 0;
}
