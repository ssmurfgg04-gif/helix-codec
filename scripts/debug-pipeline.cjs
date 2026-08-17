const path = require('path');
const { getCachedLDPCInner } = require(path.resolve(__dirname, '../src/lib/dna/ldpc-codec'));
const { crc16Bytes } = require(path.resolve(__dirname, '../src/lib/dna/crc16'));

const mod = { exports: {} };
process.dlopen(mod, path.resolve(__dirname, '../rust/helix-dna-napi/target/release/libhelix_dna_napi.so'));
const addon = mod.exports;
console.log(addon.napiVersion());

const payloadBytes = 30, ldpcParity = 8;
const ldpcCode = getCachedLDPCInner(payloadBytes + ldpcParity, payloadBytes);

const payload = new Uint8Array(payloadBytes);
for (let i = 0; i < payloadBytes; i++) payload[i] = (i * 37 + 17) & 0xFF;

const ldpcCW = ldpcCode.encode(payload);
console.log('Payload:', payloadBytes, 'B → LDPC CW:', ldpcCW.length, 'B');

const withCrc = new Uint8Array(ldpcCW.length + 2);
withCrc.set(ldpcCW, 0);
const crc = crc16Bytes(ldpcCW);
withCrc[ldpcCW.length] = crc[0]; withCrc[ldpcCW.length + 1] = crc[1];
console.log('With CRC:', withCrc.length, 'B (last 2:', withCrc[withCrc.length-2].toString(16), withCrc[withCrc.length-1].toString(16), ')');

const convOut = addon.convK7Encode(withCrc);
console.log('Conv encoded:', convOut.length, 'B');

// Standard decode (no noise, no indels)
const stdDec = addon.viterbiK9DecodeStandard(addon.convK9Encode(withCrc));
console.log('K=9 standard roundtrip:', stdDec.length, 'B');
let stdMatch = true;
for (let i = 0; i < withCrc.length; i++) if (stdDec[i] !== withCrc[i]) stdMatch = false;
console.log('K=9 standard match:', stdMatch ? 'PASS' : 'FAIL');

// Indel decode (no noise — should recover perfectly)
const indelDec7 = addon.viterbiK7Decode(convOut, {maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.0});
console.log('\nK=7 indel decode (clean):', indelDec7.length, 'B');
let match7 = true, diff7 = 0;
for (let i = 0; i < Math.min(withCrc.length, indelDec7.length); i++) if (withCrc[i] !== indelDec7[i]) { match7 = false; diff7++; }
console.log('K=7 match:', match7 ? 'PASS' : `FAIL (${diff7} byte diffs, expected ${withCrc.length}, got ${indelDec7.length})`);

const indelDec9 = addon.viterbiK9Decode(addon.convK9Encode(withCrc), {maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.0});
console.log('K=9 indel decode (clean):', indelDec9.length, 'B');
let match9 = true, diff9 = 0;
for (let i = 0; i < Math.min(withCrc.length, indelDec9.length); i++) if (withCrc[i] !== indelDec9[i]) { match9 = false; diff9++; }
console.log('K=9 match:', match9 ? 'PASS' : `FAIL (${diff9} byte diffs)`);

// Now test with DNA conversion
const BASES = 'ACGT';
function bytesToDna(data) {
  const dna=[]; for(const byte of data) for(let bit=7;bit>=1;bit-=2){dna.push(BASES[((byte>>bit)&1)<<1|((byte>>(bit-1))&1)]);}
  return dna.join('');
}
function dnaToBytes(dna) {
  const bits=[]; for(const c of dna){const code='ACGT'.indexOf(c);bits.push((code>>1)&1);bits.push(code&1);}
  const bytes=new Uint8Array(Math.floor(bits.length/8));
  for(let b=0;b<bytes.length*8&&b<bits.length;b++) bytes[b>>3]|=bits[b]<<(7-(b&7));
  return bytes;
}

const dna = bytesToDna(convOut);
const backBytes = dnaToBytes(dna);
console.log('\nDNA:', dna.length, 'nt → bytes:', backBytes.length, 'B');
let dnaMatch = true;
for (let i = 0; i < convOut.length; i++) if (convOut[i] !== backBytes[i]) dnaMatch = false;
console.log('DNA roundtrip:', dnaMatch ? 'PASS' : 'FAIL');

const vitFromDna = addon.viterbiK7Decode(backBytes, {maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.0});
console.log('K=7 from DNA bytes:', vitFromDna.length, 'B');

// CRC check
if (vitFromDna.length >= 2) {
  const dataPart = vitFromDna.slice(0, vitFromDna.length - 2);
  const computedCrc = crc16Bytes(dataPart);
  const crcPass = vitFromDna[vitFromDna.length-2] === computedCrc[0] && vitFromDna[vitFromDna.length-1] === computedCrc[1];
  console.log('CRC check:', crcPass ? 'PASS' : 'FAIL');
  console.log('Expected CRC:', withCrc[withCrc.length-2].toString(16).padStart(2,'0'), withCrc[withCrc.length-1].toString(16).padStart(2,'0'));
  console.log('Got CRC:', computedCrc[0].toString(16).padStart(2,'0'), computedCrc[1].toString(16).padStart(2,'0'));
}

// LDPC decode
if (vitFromDna.length >= payloadBytes + ldpcParity) {
  try {
    const { data } = ldpcCode.decode(vitFromDna.slice(0, payloadBytes + ldpcParity));
    let ldpcMatch = true;
    for (let i = 0; i < payloadBytes; i++) if (data[i] !== payload[i]) ldpcMatch = false;
    console.log('LDPC decode:', ldpcMatch ? 'PASS' : 'FAIL');
  } catch(e) { console.log('LDPC decode error:', e.message); }
}
