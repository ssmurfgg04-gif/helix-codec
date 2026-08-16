import {
  compressWithNAF,
  decompressWithNAF,
} from '../src/lib/dna/dna-compress-real';

// Let me trace the NAF pipeline step by step
const NUCLEOTIDE_2BIT: Record<number, number> = {
  0x41: 0, 0x43: 1, 0x47: 2, 0x54: 3,
  0x61: 0, 0x63: 1, 0x67: 2, 0x74: 3,
};
const BIT2_NUCLEOTIDE = [0x41, 0x43, 0x47, 0x54];

const dna = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const dnaBytes = new Uint8Array(dna.length);
for (let i = 0; i < dna.length; i++) dnaBytes[i] = dna.charCodeAt(i);

// Step 1: extract 2-bit
const dna2bit = new Uint8Array(dna.length);
for (let i = 0; i < dna.length; i++) dna2bit[i] = NUCLEOTIDE_2BIT[dnaBytes[i]];
console.log("2-bit values (first 20):", Array.from(dna2bit.slice(0, 20)));

// Step 2: pack 2-bit
const count = dna2bit.length;
const packedLen = Math.ceil(count / 4);
const packed = new Uint8Array(4 + packedLen);
new DataView(packed.buffer).setUint32(0, count, true);
for (let i = 0; i < count; i++) {
  packed[4 + (i >> 2)] |= (dna2bit[i] & 0b11) << (6 - (i % 4) * 2);
}
console.log("Packed length:", packed.length, "bytes (4 header +", packedLen, "data)");
console.log("Packed (first 10):", Array.from(packed.slice(0, 10)));

// Step 3: RLE
const rle: number[] = [];
let ri = 0;
while (ri < packed.length) {
  const val = packed[ri];
  let cnt = 1;
  while (ri + cnt < packed.length && packed[ri + cnt] === val && cnt < 255) cnt++;
  rle.push(val, cnt);
  ri += cnt;
}
const rleBytes = new Uint8Array(rle);
console.log("RLE length:", rleBytes.length, "bytes");
console.log("RLE (all):", Array.from(rleBytes));

// Now compress and then decompress
const compressed = compressWithNAF(dnaBytes, 6);
console.log("\nCompressed:", compressed.length, "bytes");

// Parse header
const v = new DataView(compressed.buffer, compressed.byteOffset);
const seqLen = v.getUint32(4, true);
const rleLen = v.getUint32(8, true);
const modelLen = v.getUint32(12, true);
console.log("seqLen:", seqLen, "rleLen:", rleLen, "modelLen:", modelLen);

// Manually decode the arithmetic stream
const model = compressed.slice(16, 16 + modelLen);
const compData = compressed.slice(16 + modelLen);
console.log("Model length:", model.length, "Compressed data length:", compData.length);

// Try decode
const { ArithmeticDecoder, AdaptiveFrequencyModel } = await import('../src/lib/dna/arithmetic-coder');
const dec = new ArithmeticDecoder(compData);
const arithModel = AdaptiveFrequencyModel.deserialize(model, 0).model;
const decodedRle = new Uint8Array(rleLen);
for (let i = 0; i < rleLen; i++) {
  const ft = arithModel.getFrequencyTable();
  decodedRle[i] = dec.decode(ft);
  arithModel.update(decodedRle[i]);
}
console.log("Decoded RLE:", Array.from(decodedRle));
console.log("RLE match:", JSON.stringify(Array.from(decodedRle)) === JSON.stringify(Array.from(rleBytes)));

// RLE decode
const rleDecoded: number[] = [];
for (let i = 0; i < decodedRle.length; i += 2) {
  for (let j = 0; j < decodedRle[i + 1]; j++) rleDecoded.push(decodedRle[i]);
}
const rleDecodedBytes = new Uint8Array(rleDecoded);
console.log("RLE decoded length:", rleDecodedBytes.length);
console.log("RLE decoded (first 10):", Array.from(rleDecodedBytes.slice(0, 10)));

// 2-bit unpack
const numNuc = new DataView(rleDecodedBytes.buffer, rleDecodedBytes.byteOffset, rleDecodedBytes.byteLength).getUint32(0, true);
console.log("Num nucleotides:", numNuc);
const values = new Uint8Array(numNuc);
for (let i = 0; i < numNuc; i++) {
  values[i] = (rleDecodedBytes[4 + (i >> 2)] >> (6 - (i % 4) * 2)) & 0b11;
}
console.log("2-bit values (first 20):", Array.from(values.slice(0, 20)));

// Convert to DNA bytes
const out = new Uint8Array(numNuc);
for (let i = 0; i < numNuc; i++) out[i] = BIT2_NUCLEOTIDE[values[i]];
console.log("Output DNA (first 20):", String.fromCharCode(...out.slice(0, 20)));
console.log("Expected DNA (first 20):", dna.substring(0, 20));
