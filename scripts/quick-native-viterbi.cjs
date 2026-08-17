/**
 * Minimal Nanopore validation — native Viterbi only (no MSA)
 * Tests K=7 Viterbi + LDPC + RS pipeline at 9% IDS
 */
const { getCachedLDPCInner } = require('../src/lib/dna/ldpc-codec');
const { ReedSolomon } = require('../src/lib/dna/reedsolomon');
const { crc16Bytes } = require('../src/lib/dna/crc16');

// Load native addon
const mod = { exports: {} };
process.dlopen(mod, require('path').resolve(__dirname, '../rust/helix-dna-napi/target/release/libhelix_dna_napi.so'));
const addon = mod.exports;
console.log(addon.napiVersion());

const BASES = 'ACGT';
class Rng { constructor(s) { this.s = (s>>>0)||1; } next() { this.s^=this.s<<13;this.s^=this.s>>>17;this.s^=this.s<<5;this.s=this.s>>>0;return this.s/0x100000000; } nextInt(m) { return Math.floor(this.next()*m); } }

function noisyChannel(dna, sub, ins, del, rng) {
  const r=[]; let s=0,ni=0,d=0;
  for(let p=0;p<dna.length;p++){
    if(rng.next()<del){d++;continue;}
    let b=dna[p]; if(rng.next()<sub){let nb;do{nb=BASES[rng.nextInt(4)];}while(nb===b);b=nb;s++;}
    r.push(b); if(rng.next()<ins){r.push(BASES[rng.nextInt(4)]);ni++;}
  }
  return {noisy:r.join(''),s,i:ni,d};
}

function dnaToBytes(dna) {
  const bits=[]; for(const c of dna){const code='ACGT'.indexOf(c);bits.push((code>>1)&1);bits.push(code&1);}
  const bytes=new Uint8Array(Math.floor(bits.length/8));
  for(let b=0;b<bytes.length*8&&b<bits.length;b++) bytes[b>>3]|=bits[b]<<(7-(b&7));
  return bytes;
}
function bytesToDna(data) {
  const dna=[]; for(const byte of data) for(let bit=7;bit>=1;bit-=2){const code=((byte>>bit)&1)<<1|((byte>>(bit-1))&1);dna.push(BASES[code]);}
  return dna.join('');
}
function pluralityConsensus(reads, len) {
  const r=[]; for(let p=0;p<len;p++){const v=[0,0,0,0];for(const rd of reads){if(p<rd.length){const i='ACGT'.indexOf(rd[p]);if(i>=0)v[i]++;}}let b=0;for(let i=1;i<4;i++)if(v[i]>v[b])b=i;r.push(BASES[b]);}
  return r.join('');
}

const rng = new Rng(42);
const numOligos = 15;
const payloadBytes = 30;
const idsRate = 0.09;
const delR = idsRate*0.45, insR = idsRate*0.30, subR = idsRate*0.25;

let ldpcCode; try { ldpcCode = getCachedLDPCInner(payloadBytes+8, payloadBytes); } catch(e) { console.log('LDPC err:', e.message); }
let rsCode; try { rsCode = new ReedSolomon({n:255,k:223}); } catch {}

const configs = [
  { cov: 10, par: 8, k7: true },
  { cov: 20, par: 8, k7: true },
  { cov: 20, par: 10, k7: true },
  { cov: 30, par: 8, k7: true },
  { cov: 20, par: 8, k7: false },
  { cov: 30, par: 10, k7: false },
];

for (const cfg of configs) {
  let ldpc; try { ldpc = getCachedLDPCInner(payloadBytes+cfg.par, payloadBytes); } catch { continue; }
  let recovered=0, viterbiOk=0, crcOk=0, ldpcOk=0;
  const t0 = Date.now();

  for(let oligo=0; oligo<numOligos; oligo++){
    const payload = new Uint8Array(payloadBytes);
    for(let j=0;j<payloadBytes;j++) payload[j]=rng.nextInt(256);
    let ldpcCW = payload;
    try { ldpcCW=ldpc.encode(payload); } catch {}
    const withCrc = new Uint8Array(ldpcCW.length+2);
    withCrc.set(ldpcCW,0);
    const crc=crc16Bytes(ldpcCW); withCrc[ldpcCW.length]=crc[0]; withCrc[ldpcCW.length+1]=crc[1];
    const convOut = cfg.k7 ? addon.convK7Encode(withCrc) : addon.convK9Encode(withCrc);
    const dna = bytesToDna(convOut);

    const reads = [];
    for(let r=0;r<cfg.cov;r++) reads.push(noisyChannel(dna,subR,insR,delR,rng).noisy);
    const cons = pluralityConsensus(reads, dna.length);
    const consBytes = dnaToBytes(cons);

    let afterConv;
    try {
      afterConv = new Uint8Array(cfg.k7 ? addon.viterbiK7Decode(consBytes, {maxDrift:10,insertionPenalty:1.5,deletionPenalty:1.0}) : addon.viterbiK9Decode(consBytes, {maxDrift:10,insertionPenalty:1.5,deletionPenalty:1.0}));
      viterbiOk++;
    } catch { afterConv=consBytes; }

    let crcPass = false;
    if(afterConv.length>=2){const d=afterConv.slice(0,afterConv.length-2);const c=crc16Bytes(d);crcPass=afterConv[afterConv.length-2]===c[0]&&afterConv[afterConv.length-1]===c[1];if(crcPass)crcOk++;}

    let decoded=null;
    if(ldpc&&afterConv.length>=payloadBytes+cfg.par){try{const{data}=ldpc.decode(afterConv.slice(0,payloadBytes+cfg.par));if(data.length===payloadBytes){decoded=data;ldpcOk++;}}catch{}}
    if(!decoded&&crcPass&&afterConv.length>=payloadBytes) decoded=afterConv.slice(0,payloadBytes);

    if(decoded){let m=true;for(let b=0;b<payloadBytes;b++)if(decoded[b]!==payload[b]){m=false;break;}if(m)recovered++;}
  }

  const ms = Date.now()-t0;
  console.log(
    `${cfg.k7?'K=7':'K=9'} ${cfg.cov}× ${cfg.par}B: ` +
    `${(recovered/numOligos*100).toFixed(1).padStart(6)}% ` +
    `(V:${viterbiOk} C:${crcOk} L:${ldpcOk}) ` +
    `[${ms}ms, ${(ms/numOligos).toFixed(1)}ms/oligo]`
  );
}
