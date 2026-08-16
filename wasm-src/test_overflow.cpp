#include "cpp-compressors/arithmetic_coder.h"
#include <stdio.h>
#include <stdlib.h>
using namespace arith;
int main() {
    // Check if cache + top_byte ever overflows
    const int N = 10000;
    uint32_t* syms = new uint32_t[N];
    uint32_t st = 42;
    for (int i = 0; i < N; i++) { st ^= st << 13; st ^= st >> 17; st ^= st << 5; syms[i] = (st >> 16) % 4; }
    
    // Manually encode and check overflow
    uint8_t comp[16384];
    Enc enc; enc.init(comp, sizeof(comp));
    uint32_t fr[5]; AModel m; m.init(4, fr);
    
    for (int i = 0; i < N; i++) {
        m.esym(enc, syms[i]);
        // Check if cache + top_byte could overflow
        // cache is 0..0xFE, top_byte is 0..0xFE
        // max = 0xFE + 0xFE = 0x1FC
    }
    int cl = enc.fin();
    printf("Encoded %d -> %d bytes\n", N, cl);
    
    // Decode
    Dec dec; dec.init(comp, cl);
    uint32_t fr2[5]; AModel m2; m2.init(4, fr2);
    int ok = 1, ff = -1;
    for (int i = 0; i < N; i++) {
        uint32_t s = m2.dsym(dec);
        if (s != syms[i]) { if (ff == -1) ff = i; ok = 0; }
    }
    printf("Result: %s", ok ? "PASS" : "FAIL");
    if (!ok) printf(" (fail@%d)", ff);
    printf("\n");
    
    // If failing, try to understand what the first few encoded bytes are
    printf("First 16 bytes:");
    for (int i = 0; i < 16 && i < cl; i++) printf(" %02x", comp[i]);
    printf("\n");
    
    delete[] syms;
    return 0;
}
