const path = require('path');
const mod = { exports: {} };
process.dlopen(mod, path.resolve(__dirname, '../rust/helix-dna-napi/target/release/libhelix_dna_napi.so'));
const addon = mod.exports;
console.log(addon.napiVersion());

const payloadBytes = 30;
const BASES = 'ACGT';

class Rng { constructor(s) { this.s = (s>>>0)||1; } next() { this.s^=this.s<<13;this.s^=this.s>>>17;this.s^=this.s<<5;this.s=this.s>>>0;return this.s/0x100000000; } nextInt(m) { return Math.floor(this.next()*m); } }

function bytesToDna(data) { const d=[]; for(const b of data)for(let bit=7;bit>=1;bit-=2)d.push(BASES[((b>>bit)&1)<<1|((b>>(bit-1))&1)]);return d.join(''); }
function dnaToBytes(dna) { const bits=[]; for(const c of dna){const code='ACGT'.indexOf(c);bits.push((code>>1)&1);bits.push(code&1);}const bytes=new Uint8Array(Math.floor(bits.length/8));for(let b=0;b<bytes.length*8&&b<bits.length;b++)bytes[b>>3]|=bits[b]<<(7-(b&7));return bytes; }
function plurality(reads,len) { const r=[];for(let p=0;p<len;p++){const v=[0,0,0,0];for(const rd of reads){if(p<rd.length){const i='ACGT'.indexOf(rd[p]);if(i>=0)v[i]++;}}let b=0;for(let i=1;i<4;i++)if(v[i]>v[b])b=i;r.push(BASES[b]);}return r.join(''); }

// Clean roundtrip
const payload = Buffer.alloc(payloadBytes);
for (let i = 0; i < payloadBytes; i++) payload[i] = (i * 37 + 17) & 0xFF;
const enc = addon.convK7Encode(payload);
const dec = addon.viterbiK7Decode(enc, {maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.0, expectedLength: 30});
let match = true;
for (let i = 0; i < payloadBytes; i++) if (dec[i] !== payload[i]) match = false;
console.log('Clean roundtrip:', match ? 'PASS ✓' : 'FAIL');

// Noisy channel test
const rng = new Rng(42);
let recovered = 0, viterbiOk = 0;
const t0 = Date.now();

for (let oligo = 0; oligo < 20; oligo++) {
  const pl = Buffer.alloc(payloadBytes);
  for (let i = 0; i < payloadBytes; i++) pl[i] = rng.nextInt(256);
  const enc = addon.convK7Encode(pl);
  const dna = bytesToDna(enc);

  // Generate 20 noisy reads at 9% IDS
  const reads = [];
  for (let r = 0; r < 20; r++) {
    let noisy = '';
    for (let p = 0; p < dna.length; p++) {
      if (rng.next() < 0.04) continue; // deletion 4%
      let base = dna[p];
      if (rng.next() < 0.0225) { let nb; do { nb = BASES[rng.nextInt(4)]; } while (nb === base); base = nb; } // sub 2.25%
      noisy += base;
      if (rng.next() < 0.0275) noisy += BASES[rng.nextInt(4)]; // insertion 2.75%
    }
    reads.push(noisy);
  }

  const cons = plurality(reads, dna.length);
  const consBytes = dnaToBytes(cons);

  // Viterbi decode with expectedLength
  try {
    const dec = addon.viterbiK7Decode(consBytes, {maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.0, expectedLength: 30});
    viterbiOk++;
    let m = true;
    for (let i = 0; i < payloadBytes; i++) if (dec[i] !== pl[i]) m = false;
    if (m) recovered++;
  } catch(e) {
    // Viterbi failed
  }
}

console.log('\n=== Nanopore 9% IDS — K=7 Viterbi + Plurality Consensus ===');
console.log('Recovered:', recovered + '/20 (' + (recovered/20*100).toFixed(1) + '%)');
console.log('Viterbi ok:', viterbiOk + '/20');
console.log('Time:', Date.now()-t0, 'ms (' + ((Date.now()-t0)/20).toFixed(1) + 'ms/oligo)');

// Also test with MSA (if available)
try {
  const { msaConsensus } = require(path.resolve(__dirname, '../src/lib/dna/msa-consensus'));
  console.log('\nMSA module available, testing MSA consensus...');

  let msaRecovered = 0;
  const rng2 = new Rng(123);
  for (let oligo = 0; oligo < 5; oligo++) {
    const pl = Buffer.alloc(payloadBytes);
    for (let i = 0; i < payloadBytes; i++) pl[i] = rng2.nextInt(256);
    const enc = addon.convK7Encode(pl);
    const dna = bytesToDna(enc);
    const reads = [];
    const qualities = [];
    for (let r = 0; r < 10; r++) {
      let noisy = '';
      const q = [];
      for (let p = 0; p < dna.length; p++) {
        if (rng2.next() < 0.04) continue;
        let base = dna[p];
        if (rng2.next() < 0.0225) { let nb; do { nb = BASES[rng2.nextInt(4)]; } while (nb === base); base = nb; q.push(15); }
        else { q.push(30); }
        noisy += base;
        if (rng2.next() < 0.0275) { noisy += BASES[rng2.nextInt(4)]; q.push(10); }
      }
      reads.push(noisy);
      qualities.push(new Uint8Array(q));
    }
    const msaResult = msaConsensus(reads, qualities, { iterations: 1 });
    const consBytes = dnaToBytes(msaResult.consensus);
    try {
      const dec = addon.viterbiK7Decode(consBytes, {maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.0, expectedLength: 30});
      let m = true;
      for (let i = 0; i < payloadBytes; i++) if (dec[i] !== pl[i]) m = false;
      if (m) msaRecovered++;
    } catch {}
  }
  console.log('MSA + K=7 recovered:', msaRecovered + '/5');
} catch(e) {
  console.log('MSA not directly testable from CJS (TS module):', e.message);
}
