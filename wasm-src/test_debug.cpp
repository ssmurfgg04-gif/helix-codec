#include "cpp-compressors/arithmetic_coder.h"
#include <stdio.h>
#include <string>
#include <stdlib.h>

using namespace arith;

int main() {
    // Encode [0, 1, 2, 3] with uniform model
    uint8_t comp[1024];
    Enc enc; enc.init(comp, sizeof(comp));
    uint32_t fr[5]; AModel m; m.init(4, fr);
    
    printf("Encode symbol 0\n"); m.esym(enc, 0);
    printf("Encode symbol 1\n"); m.esym(enc, 1);
    printf("Encode symbol 2\n"); m.esym(enc, 2);
    printf("Encode symbol 3\n"); m.esym(enc, 3);
    int cl = enc.fin();
    
    printf("Compressed %d bytes:", cl);
    for (int i=0;i<cl;i++) printf(" %02x", comp[i]);
    printf("\n\n");
    
    // Decode
    Dec dec; dec.init(comp, cl);
    uint32_t fr2[5]; AModel m2; m2.init(4, fr2);
    
    printf("Decoder: code=%08x range=%08x\n", dec.cd, dec.rn);
    
    for (int i=0;i<4;i++) {
        uint32_t c = dec.gcf(m2.tt);
        printf("  gcf(%u) = %u (range_after=%08x)\n", m2.tt, c, dec.rn);
        
        // Find symbol
        uint32_t lo=0,hi=m2.n;
        while(lo<hi){uint32_t mid=(lo+hi)>>1; if(fr2[mid+1]<=c) lo=mid+1; else hi=mid;}
        uint32_t sym=lo;
        printf("  -> symbol %u\n", sym);
        
        dec.dec(fr2[sym], fr2[sym+1]-fr2[sym]);
        m2.up(sym);
        printf("  after dec: code=%08x range=%08x\n", dec.cd, dec.rn);
    }
    
    return 0;
}
