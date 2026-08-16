#include "cpp-compressors/arithmetic_coder.h"
#include <stdio.h>
#include <stdlib.h>
using namespace arith;
int main() {
    // Generate same sequence
    const int N = 20;
    uint32_t syms[N];
    uint32_t st = 42;
    for (int i = 0; i < N; i++) { st ^= st << 13; st ^= st >> 17; st ^= st << 5; syms[i] = (st >> 16) % 4; }

    printf("Symbols: ");
    for (int i = 0; i < N; i++) printf("%u ", syms[i]);
    printf("\n");

    // Encode
    uint8_t comp[256];
    Enc enc; enc.init(comp, sizeof(comp));
    uint32_t fr[5]; AModel m; m.init(4, fr);
    for (int i = 0; i < N; i++) {
        uint64_t old_lo = enc.lo;
        m.esym(enc, syms[i]);
        printf("After sym %u (pos %d): lo=%010llx rn=%08x (lo_bits40=%05llx)\n",
               syms[i], i, (unsigned long long)enc.lo, enc.rn, (unsigned long long)(enc.lo >> 24));
    }
    int cl = enc.fin();
    printf("Compressed %d bytes:", cl);
    for (int i = 0; i < cl; i++) printf(" %02x", comp[i]);
    printf("\n\n");

    // Decode
    Dec dec; dec.init(comp, cl);
    uint32_t fr2[5]; AModel m2; m2.init(4, fr2);
    for (int i = 0; i < N; i++) {
        uint32_t c = dec.gcf(m2.tt);
        printf("Decode pos %d: cd=%08x rn=%08x gcf(%u)=%u", i, dec.cd, dec.rn, m2.tt, c);
        
        // Find symbol
        uint32_t lo2=0, hi2=m2.n;
        while(lo2<hi2){uint32_t m2v=(lo2+hi2)>>1; if(fr2[m2v+1]<=c) lo2=m2v+1; else hi2=m2v;}
        uint32_t s = lo2;
        printf(" -> sym %u (expected %u) %s\n", s, syms[i], s==syms[i]?"OK":"FAIL");
        
        dec.dec(fr2[s], fr2[s+1]-fr2[s]);
        m2.up(s);
    }
    return 0;
}
