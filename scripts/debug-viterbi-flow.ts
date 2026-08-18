import * as path from "node:path";
import { getCachedLDPCInner } from "../src/lib/dna/ldpc-codec";
import { crc16Bytes } from "../src/lib/dna/crc16";

const addonPath = path.resolve(process.cwd(), "rust/helix-dna-napi/target/release/libhelix_dna_napi.so");
const mod = { exports: {} } as any;
(process as any).dlopen(mod, addonPath);
const addon = mod.exports;
console.log(`Native Viterbi: ${addon.napiVersion()}`);

const data = new Uint8Array(30);
for (let i = 0; i < 30; i++) data[i] = i * 7 + 13;

const ldpc = getCachedLDPCInner(30 + 8, 30);
const ldpcCW = ldpc.encode(data);
console.log('LDPC codeword:', ldpcCW.length, 'bytes');

const withCrc = new Uint8Array(ldpcCW.length + 2);
withCrc.set(ldpcCW, 0);
const crc = crc16Bytes(ldpcCW);
withCrc[ldpcCW.length] = crc[0]; withCrc[ldpcCW.length + 1] = crc[1];
console.log('With CRC:', withCrc.length, 'bytes');

const convOut = new Uint8Array(addon.convK9Encode(withCrc));
console.log('Conv encoded:', convOut.length, 'bytes', convOut.length * 8, 'bits');

// Generate 20 noisy reads at 9% IDS (bit-level)
const delR = 0.09 * 0.45, insR = 0.09 * 0.30, subR = 0.09 * 0.25;
const reads: number[][] = [];
for (let r = 0; r < 20; r++) {
  const bits: number[] = [];
  for (let i = 0; i < convOut.length * 8; i++) bits.push((convOut[Math.floor(i/8)] >> (7 - (i%8))) & 1);
  const result: number[] = [];
  let pos = 0;
  while (pos < bits.length) {
    if (Math.random() < delR) { pos++; continue; }
    if (Math.random() < subR) result.push(bits[pos] ^ 1); else result.push(bits[pos]);
    pos++;
    if (Math.random() < insR) result.push(Math.random() < 0.5 ? 0 : 1);
  }
  reads.push(result);
}

// Consensus: majority vote per bit position
const avgLen = Math.round(reads.reduce((s, r) => s + r.length, 0) / reads.length);
console.log('Average read length:', avgLen, '(original:', convOut.length * 8, ')');

const consBits: number[] = [];
for (let pos = 0; pos < avgLen; pos++) {
  let ones = 0, zeros = 0;
  for (const r of reads) { if (pos < r.length) { if (r[pos] === 1) ones++; else zeros++; } }
  consBits.push(ones >= zeros ? 1 : 0);
}
const consBytes = new Uint8Array(Math.ceil(consBits.length / 8));
for (let b = 0; b < consBits.length; b++) consBytes[Math.floor(b/8)] |= consBits[b] << (7 - (b%8));

console.log('Consensus:', consBytes.length, 'bytes', consBits.length, 'bits');

// Count errors in consensus vs original
let consensusErrors = 0;
for (let i = 0; i < Math.min(consBits.length, convOut.length * 8); i++) {
  const origBit = (convOut[Math.floor(i/8)] >> (7 - (i%8))) & 1;
  if (consBits[i] !== origBit) consensusErrors++;
}
console.log('Consensus bit errors:', consensusErrors, '/', convOut.length * 8, '(' + (consensusErrors / (convOut.length * 8) * 100).toFixed(2) + '%)');

// Viterbi decode the consensus
const vitDec = new Uint8Array(addon.viterbiK9Decode(Buffer.from(consBytes), { maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.5, numInfoBits: withCrc.length * 8 }));
console.log('Viterbi decoded:', vitDec.length, 'bytes');
console.log('Expected inner:', withCrc.length, 'bytes');

// Compare with original
let vitErrors = 0;
for (let i = 0; i < Math.min(vitDec.length, withCrc.length); i++) {
  if (vitDec[i] !== withCrc[i]) vitErrors++;
}
console.log('Viterbi output byte errors:', vitErrors, '/', withCrc.length);

// CRC check
if (vitDec.length >= withCrc.length) {
  const dp = vitDec.slice(0, vitDec.length - 2);
  const c = crc16Bytes(dp);
  const crcOk = vitDec[vitDec.length - 2] === c[0] && vitDec[vitDec.length - 1] === c[1];
  console.log('CRC check:', crcOk);
}

// LDPC decode
if (vitDec.length >= 38) {
  try {
    const { data: ldpcDec } = ldpc.decode(vitDec.slice(0, 38));
    console.log('LDPC decoded:', ldpcDec.length, 'bytes');
    console.log('LDPC match:', Buffer.from(ldpcDec).equals(Buffer.from(data)));
  } catch(e: any) { console.log('LDPC failed:', e.message?.slice(0, 100)); }
}

// Also try: direct Viterbi on each read individually (no consensus)
console.log('\n--- Per-read decode (no consensus) ---');
let perReadRecovered = 0;
for (let r = 0; r < reads.length; r++) {
  const readBytes = new Uint8Array(Math.ceil(reads[r].length / 8));
  for (let b = 0; b < reads[r].length; b++) readBytes[Math.floor(b/8)] |= reads[r][b] << (7 - (b%8));
  try {
    const dec = new Uint8Array(addon.viterbiK9Decode(Buffer.from(readBytes), { maxDrift: 10, insertionPenalty: 1.5, deletionPenalty: 1.5, numInfoBits: withCrc.length * 8 }));
    if (dec.length >= 38) {
      try {
        const { data: ldpcDec } = ldpc.decode(dec.slice(0, 38));
        if (ldpcDec.length === 30 && Buffer.from(ldpcDec).equals(Buffer.from(data))) perReadRecovered++;
      } catch {}
    }
  } catch {}
}
console.log('Per-read recovered:', perReadRecovered, '/', reads.length);
