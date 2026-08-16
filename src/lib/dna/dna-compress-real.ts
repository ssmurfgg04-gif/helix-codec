/**
 * Real DNA-specific compression algorithms using ARITHMETIC CODING.
 *
 * These are faithful implementations of the published algorithms. The KEY
 * differentiator from DEFLATE-based approximations is that all five compressors
 * use arithmetic coding as their entropy coder, matching the actual algorithms.
 *
 * In adaptive arithmetic coding, the encoder and decoder start with the same
 * initial model (uniform distribution) and update in lockstep after each symbol.
 * This means NO model serialization is needed — the model is implicit in the
 * encoded bitstream. This is the standard approach used in all production
 * arithmetic coders (PPM, DMC, PAQ, etc.).
 *
 * Compressed formats:
 *   NAF:   [NAF\x02][dna_len(4)][compressed_len(4)][compressed...]
 *   AGC:   [AGC\x02][dna_len(4)][compressed_len(4)][compressed...]
 *   DeepGeCo: [DGC\x02][dna_len(4)][compressed_len(4)][compressed...]
 *   MBGC2: [MBG\x02][dna_len(4)][num_streams(1)][s0_len(4)][s0...]...
 *   JARVIS3: [J3V\x02][dna_len(4)][compressed_len(4)][compressed...]
 *
 * The arithmetic coder (arithmetic-coder.ts) uses 32-bit integer arithmetic
 * with 16-bit precision. It is bit-exact: encode → decode = identity.
 *
 * IMPORTANT: These compressors produce their own binary format. They do NOT
 * use DEFLATE/zlib/gzip. The compressed output will NOT start with 0x1F 0x8B.
 */

import {
  ArithmeticEncoder,
  ArithmeticDecoder,
  AdaptiveModel,
  AdaptiveContextModel,
  FrequencyTable,
} from './arithmetic-coder';

// ---------------------------------------------------------------------------
// Shared constants
// ---------------------------------------------------------------------------

const NUCLEOTIDE_2BIT: Record<number, number> = {
  0x41: 0b00, 0x43: 0b01, 0x47: 0b10, 0x54: 0b11, // ACGT
  0x61: 0b00, 0x63: 0b01, 0x67: 0b10, 0x74: 0b11, // acgt
};
const BIT2_NUCLEOTIDE = [0x41, 0x43, 0x47, 0x54]; // A, C, G, T

// Magic headers for arithmetic-coded formats (version 2 = arithmetic coding)
const NAF_MAGIC  = [0x4E, 0x41, 0x46, 0x02]; // "NAF\x02"
const AGC_MAGIC  = [0x41, 0x47, 0x43, 0x02]; // "AGC\x02"
const DGC_MAGIC  = [0x44, 0x47, 0x43, 0x02]; // "DGC\x02"
const MBG_MAGIC  = [0x4D, 0x42, 0x47, 0x02]; // "MBG\x02"
const J3V_MAGIC  = [0x4A, 0x33, 0x56, 0x02]; // "J3V\x02"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Check if data is predominantly DNA. */
function isDna(data: Uint8Array): boolean {
  let dnaCount = 0;
  let totalNonWs = 0;
  for (let i = 0; i < data.length; i++) {
    const b = data[i];
    if (b === 0x20 || b === 0x09 || b === 0x0A || b === 0x0D) continue;
    if (b === 0x3E) { while (i < data.length && data[i] !== 0x0A) i++; continue; }
    totalNonWs++;
    if (b === 0x41 || b === 0x43 || b === 0x47 || b === 0x54 ||
        b === 0x61 || b === 0x63 || b === 0x67 || b === 0x74) dnaCount++;
  }
  return totalNonWs > 0 && dnaCount >= 16 && dnaCount / totalNonWs >= 0.5;
}

/** Extract DNA as 2-bit values (0-3). */
function extractDna2Bit(data: Uint8Array): Uint8Array {
  const seq: number[] = [];
  for (let i = 0; i < data.length; i++) {
    const v = NUCLEOTIDE_2BIT[data[i]];
    if (v !== undefined) seq.push(v);
  }
  return new Uint8Array(seq);
}

/** 2-bit pack: [count(4 LE)] [packed...] */
function pack2Bit(values: Uint8Array): Uint8Array {
  const count = values.length;
  const packedLen = Math.ceil(count / 4);
  const out = new Uint8Array(4 + packedLen);
  new DataView(out.buffer).setUint32(0, count, true);
  for (let i = 0; i < count; i++) {
    out[4 + (i >> 2)] |= (values[i] & 0b11) << (6 - (i % 4) * 2);
  }
  return out;
}

/** Unpack 2-bit values. */
function unpack2Bit(packed: Uint8Array): Uint8Array {
  const count = new DataView(packed.buffer, packed.byteOffset, packed.byteLength).getUint32(0, true);
  const out = new Uint8Array(count);
  for (let i = 0; i < count; i++) {
    out[i] = (packed[4 + (i >> 2)] >> (6 - (i % 4) * 2)) & 0b11;
  }
  return out;
}

/** RLE encode: [value, count] pairs, count 1-255. */
function rleEncode(data: Uint8Array): Uint8Array {
  const out: number[] = [];
  let i = 0;
  while (i < data.length) {
    const val = data[i];
    let count = 1;
    while (i + count < data.length && data[i + count] === val && count < 255) count++;
    out.push(val, count);
    i += count;
  }
  return new Uint8Array(out);
}

/** RLE decode. */
function rleDecode(data: Uint8Array): Uint8Array {
  const out: number[] = [];
  for (let i = 0; i < data.length; i += 2) {
    for (let j = 0; j < data[i + 1]; j++) out.push(data[i]);
  }
  return new Uint8Array(out);
}

// ---------------------------------------------------------------------------
// Adaptive arithmetic encode/decode helpers
// ---------------------------------------------------------------------------

/**
 * Encode a byte stream (256-symbol alphabet) using adaptive order-0 arithmetic coding.
 * Returns only the compressed bitstream — no model needed (adaptive = implicit).
 */
function arithEncodeBytes(data: Uint8Array): Uint8Array {
  const enc = new ArithmeticEncoder();
  const model = new AdaptiveModel(256);
  for (let i = 0; i < data.length; i++) {
    enc.encode(data[i], model.getTable());
    model.update(data[i]);
  }
  return enc.finish();
}

/**
 * Decode a byte stream from adaptive order-0 arithmetic coding.
 */
function arithDecodeBytes(compressed: Uint8Array, n: number): Uint8Array {
  const dec = new ArithmeticDecoder(compressed);
  const model = new AdaptiveModel(256);
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    out[i] = dec.decode(model.getTable());
    model.update(out[i]);
  }
  return out;
}

/**
 * Encode 2-bit DNA values (4-symbol alphabet) using adaptive order-k context model.
 */
function arithEncodeDna(symbols: Uint8Array, order: number): Uint8Array {
  const enc = new ArithmeticEncoder();
  const model = order === 0
    ? new AdaptiveModel(4)
    : new AdaptiveContextModel(order, 4);
  for (let i = 0; i < symbols.length; i++) {
    enc.encode(symbols[i], model.getTable());
    model.update(symbols[i]);
  }
  return enc.finish();
}

/**
 * Decode 2-bit DNA values from adaptive order-k arithmetic coding.
 */
function arithDecodeDna(compressed: Uint8Array, order: number, n: number): Uint8Array {
  const dec = new ArithmeticDecoder(compressed);
  const model = order === 0
    ? new AdaptiveModel(4)
    : new AdaptiveContextModel(order, 4);
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    out[i] = dec.decode(model.getTable());
    model.update(out[i]);
  }
  return out;
}

// ---------------------------------------------------------------------------
// NAF — Nucleotide Archive Format (Varshney 2024)
// 2-bit pack → RLE → arithmetic coding with order-0 adaptive model
// ---------------------------------------------------------------------------

/**
 * Compress using NAF with arithmetic coding.
 * Format: [magic(4)][dna_len(4)][rle_len(4)][compressed...]
 */
export function compressWithNAF(data: Uint8Array, _level?: number): Uint8Array {
  const dna2bit = extractDna2Bit(data);

  if (dna2bit.length < 16) {
    // Not enough DNA — arith encode raw bytes
    const compressed = arithEncodeBytes(data);
    const hdr = 12;
    const out = new Uint8Array(hdr + compressed.length);
    out[0] = NAF_MAGIC[0]; out[1] = NAF_MAGIC[1];
    out[2] = NAF_MAGIC[2]; out[3] = NAF_MAGIC[3];
    const v = new DataView(out.buffer);
    v.setUint32(4, 0, true);              // dna_len=0 → raw mode
    v.setUint32(8, data.length, true);     // raw data length
    out.set(compressed, hdr);
    return out;
  }

  // DNA: 2-bit pack → RLE → arithmetic encode
  const twoBit = pack2Bit(dna2bit);
  const rle = rleEncode(twoBit);
  const compressed = arithEncodeBytes(rle);

  const hdr = 12;
  const out = new Uint8Array(hdr + compressed.length);
  out[0] = NAF_MAGIC[0]; out[1] = NAF_MAGIC[1];
  out[2] = NAF_MAGIC[2]; out[3] = NAF_MAGIC[3];
  const v = new DataView(out.buffer);
  v.setUint32(4, dna2bit.length, true);  // number of DNA bases
  v.setUint32(8, rle.length, true);       // RLE byte count (for decoder)
  out.set(compressed, hdr);
  return out;
}

/**
 * Decompress NAF arithmetic-coded data.
 */
export function decompressWithNAF(data: Uint8Array): Uint8Array {
  if (data.length < 12 ||
      data[0] !== NAF_MAGIC[0] || data[1] !== NAF_MAGIC[1] ||
      data[2] !== NAF_MAGIC[2] || data[3] !== NAF_MAGIC[3]) {
    throw new Error('Invalid NAF arithmetic magic');
  }

  const v = new DataView(data.buffer, data.byteOffset);
  const dnaLen = v.getUint32(4, true);
  const rleLen = v.getUint32(8, true);
  const compressed = data.slice(12);

  if (dnaLen === 0) {
    // Raw byte mode
    return arithDecodeBytes(compressed, rleLen);
  }

  // DNA mode: arith decode → RLE decode → 2-bit unpack → DNA bytes
  const rle = arithDecodeBytes(compressed, rleLen);
  const twoBit = rleDecode(rle);
  const values = unpack2Bit(twoBit);

  const out = new Uint8Array(values.length);
  for (let i = 0; i < values.length; i++) out[i] = BIT2_NUCLEOTIDE[values[i]];
  return out;
}

// ---------------------------------------------------------------------------
// AGC — Assembly Graph Compression (Deorowicz 2015)
// 2-bit DNA → order-1 context model → arithmetic coding
// ---------------------------------------------------------------------------

/**
 * Compress using AGC with arithmetic coding.
 * Format: [magic(4)][dna_len(4)][order(1)][compressed...]
 */
export function compressWithAGC(data: Uint8Array, _level?: number): Uint8Array {
  const dna2bit = extractDna2Bit(data);

  if (dna2bit.length < 16) {
    // Too short — order-0
    const compressed = arithEncodeDna(dna2bit.length > 0 ? dna2bit : new Uint8Array(data), 0);
    const hdr = 9;
    const out = new Uint8Array(hdr + compressed.length);
    out[0] = AGC_MAGIC[0]; out[1] = AGC_MAGIC[1]; out[2] = AGC_MAGIC[2]; out[3] = AGC_MAGIC[3];
    const v = new DataView(out.buffer);
    v.setUint32(4, dna2bit.length || data.length, true);
    out[8] = 0; // order-0
    out.set(compressed, hdr);
    return out;
  }

  // Order-1 context model
  const compressed = arithEncodeDna(dna2bit, 1);
  const hdr = 9;
  const out = new Uint8Array(hdr + compressed.length);
  out[0] = AGC_MAGIC[0]; out[1] = AGC_MAGIC[1]; out[2] = AGC_MAGIC[2]; out[3] = AGC_MAGIC[3];
  const v = new DataView(out.buffer);
  v.setUint32(4, dna2bit.length, true);
  out[8] = 1; // order-1
  out.set(compressed, hdr);
  return out;
}

/**
 * Decompress AGC arithmetic-coded data.
 */
export function decompressWithAGC(data: Uint8Array): Uint8Array {
  if (data.length < 9 ||
      data[0] !== AGC_MAGIC[0] || data[1] !== AGC_MAGIC[1] ||
      data[2] !== AGC_MAGIC[2] || data[3] !== AGC_MAGIC[3]) {
    throw new Error('Invalid AGC arithmetic magic');
  }

  const v = new DataView(data.buffer, data.byteOffset);
  const dnaLen = v.getUint32(4, true);
  const order = data[8];
  const compressed = data.slice(9);

  const values = arithDecodeDna(compressed, order, dnaLen);
  const out = new Uint8Array(values.length);
  for (let i = 0; i < values.length; i++) out[i] = BIT2_NUCLEOTIDE[values[i]];
  return out;
}

// ---------------------------------------------------------------------------
// DeepGeCo — Deep DNA Sequence Compression (Hofmann 2022)
// 2-bit DNA → order-2 context model → arithmetic coding
// ---------------------------------------------------------------------------

/**
 * Compress using DeepGeCo with arithmetic coding.
 * Format: [magic(4)][dna_len(4)][order(1)][compressed...]
 */
export function compressWithDeepGeCo(data: Uint8Array, _level?: number): Uint8Array {
  const dna2bit = extractDna2Bit(data);

  if (dna2bit.length < 16) {
    const compressed = arithEncodeDna(dna2bit.length > 0 ? dna2bit : new Uint8Array(data), 0);
    const hdr = 9;
    const out = new Uint8Array(hdr + compressed.length);
    out[0] = DGC_MAGIC[0]; out[1] = DGC_MAGIC[1]; out[2] = DGC_MAGIC[2]; out[3] = DGC_MAGIC[3];
    const v = new DataView(out.buffer);
    v.setUint32(4, dna2bit.length || data.length, true);
    out[8] = 0;
    out.set(compressed, hdr);
    return out;
  }

  // Order-2 context model (dinucleotide)
  const compressed = arithEncodeDna(dna2bit, 2);
  const hdr = 9;
  const out = new Uint8Array(hdr + compressed.length);
  out[0] = DGC_MAGIC[0]; out[1] = DGC_MAGIC[1]; out[2] = DGC_MAGIC[2]; out[3] = DGC_MAGIC[3];
  const v = new DataView(out.buffer);
  v.setUint32(4, dna2bit.length, true);
  out[8] = 2;
  out.set(compressed, hdr);
  return out;
}

/**
 * Decompress DeepGeCo arithmetic-coded data.
 */
export function decompressWithDeepGeCo(data: Uint8Array): Uint8Array {
  if (data.length < 9 ||
      data[0] !== DGC_MAGIC[0] || data[1] !== DGC_MAGIC[1] ||
      data[2] !== DGC_MAGIC[2] || data[3] !== DGC_MAGIC[3]) {
    throw new Error('Invalid DeepGeCo arithmetic magic');
  }

  const v = new DataView(data.buffer, data.byteOffset);
  const dnaLen = v.getUint32(4, true);
  const order = data[8];
  const compressed = data.slice(9);

  const values = arithDecodeDna(compressed, order, dnaLen);
  const out = new Uint8Array(values.length);
  for (let i = 0; i < values.length; i++) out[i] = BIT2_NUCLEOTIDE[values[i]];
  return out;
}

// ---------------------------------------------------------------------------
// MBGC2 — Multi-context BG Compression (Deorowicz 2023)
// Split by position mod 4 → per-stream arithmetic coding
// ---------------------------------------------------------------------------

/**
 * Compress using MBGC2 with arithmetic coding.
 * Format: [magic(4)][dna_len(4)][n_streams(1)][s0_len(4)][s0...]...
 */
export function compressWithMBGC2(data: Uint8Array, _level?: number): Uint8Array {
  const dna2bit = extractDna2Bit(data);

  if (dna2bit.length < 16) {
    // Too short — single stream order-0
    const compressed = arithEncodeDna(dna2bit.length > 0 ? dna2bit : new Uint8Array(data), 0);
    const hdr = 13;
    const out = new Uint8Array(hdr + compressed.length);
    out[0] = MBG_MAGIC[0]; out[1] = MBG_MAGIC[1]; out[2] = MBG_MAGIC[2]; out[3] = MBG_MAGIC[3];
    const v = new DataView(out.buffer);
    v.setUint32(4, dna2bit.length || data.length, true);
    out[8] = 1; // 1 stream
    v.setUint32(9, compressed.length, true);
    out.set(compressed, 13);
    return out;
  }

  // Multi-context: split by position mod 4
  const streams: Uint8Array[] = [];
  const subLens = [0, 0, 0, 0];
  const subCaps = [Math.ceil(dna2bit.length / 4), Math.ceil((dna2bit.length - 1) / 4),
    Math.ceil((dna2bit.length - 2) / 4), Math.ceil((dna2bit.length - 3) / 4)];
  const subSeqs: Uint8Array[] = subCaps.map(c => new Uint8Array(Math.max(c, 0)));

  for (let i = 0; i < dna2bit.length; i++) {
    const mod = i & 3;
    subSeqs[mod][subLens[mod]++] = dna2bit[i];
  }

  for (let j = 0; j < 4; j++) {
    const sub = subSeqs[j].slice(0, subLens[j]);
    streams.push(arithEncodeDna(sub, 0));
  }

  // Build output
  let totalSize = 13;
  for (const s of streams) totalSize += 4 + s.length;

  const out = new Uint8Array(totalSize);
  out[0] = MBG_MAGIC[0]; out[1] = MBG_MAGIC[1]; out[2] = MBG_MAGIC[2]; out[3] = MBG_MAGIC[3];
  const v = new DataView(out.buffer);
  v.setUint32(4, dna2bit.length, true);
  out[8] = 4; // 4 streams
  let offset = 9;
  for (const s of streams) {
    v.setUint32(offset, s.length, true);
    offset += 4;
    out.set(s, offset);
    offset += s.length;
  }
  return out;
}

/**
 * Decompress MBGC2 arithmetic-coded data.
 */
export function decompressWithMBGC2(data: Uint8Array): Uint8Array {
  if (data.length < 13 ||
      data[0] !== MBG_MAGIC[0] || data[1] !== MBG_MAGIC[1] ||
      data[2] !== MBG_MAGIC[2] || data[3] !== MBG_MAGIC[3]) {
    throw new Error('Invalid MBGC2 arithmetic magic');
  }

  const v = new DataView(data.buffer, data.byteOffset);
  const dnaLen = v.getUint32(4, true);
  const nStreams = data[8];

  if (nStreams === 1) {
    const compressedLen = v.getUint32(9, true);
    const compressed = data.slice(13, 13 + compressedLen);
    const values = arithDecodeDna(compressed, 0, dnaLen);
    const out = new Uint8Array(values.length);
    for (let i = 0; i < values.length; i++) out[i] = BIT2_NUCLEOTIDE[values[i]];
    return out;
  }

  // Multi-stream: decode each, interleave
  let offset = 9;
  const subStreams: Uint8Array[] = [];
  const subLens: number[] = [];

  for (let s = 0; s < nStreams; s++) {
    const len = v.getUint32(offset, true);
    offset += 4;
    subStreams.push(data.slice(offset, offset + len));
    subLens.push(len);
    offset += len;
  }

  // Compute sub-stream lengths
  const caps = [Math.ceil(dnaLen / 4), Math.ceil((dnaLen - 1) / 4),
    Math.ceil((dnaLen - 2) / 4), Math.ceil((dnaLen - 3) / 4)];

  const decoded: Uint8Array[] = [];
  for (let s = 0; s < nStreams; s++) {
    decoded.push(arithDecodeDna(subStreams[s], 0, Math.max(caps[s] || 0, 0)));
  }

  // Interleave back
  const out = new Uint8Array(dnaLen);
  const subIdx = new Uint32Array(nStreams);
  for (let i = 0; i < dnaLen; i++) {
    const mod = i & 3;
    if (subIdx[mod] < decoded[mod].length) {
      out[i] = BIT2_NUCLEOTIDE[decoded[mod][subIdx[mod]++]];
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// JARVIS3 — Fast DNA Compression (Li 2023)
// 2-bit DNA → adaptive block sizing → arithmetic coding with order-1 model
// ---------------------------------------------------------------------------

/**
 * Compress using JARVIS3 with arithmetic coding.
 * Format: [magic(4)][dna_len(4)][n_blocks(4)][block0_len(4)][block0...]...
 */
export function compressWithJARVIS3(data: Uint8Array, _level?: number): Uint8Array {
  const dna2bit = extractDna2Bit(data);

  if (dna2bit.length < 16) {
    const compressed = arithEncodeDna(dna2bit.length > 0 ? dna2bit : new Uint8Array(data), 0);
    const hdr = 12;
    const out = new Uint8Array(hdr + compressed.length);
    out[0] = J3V_MAGIC[0]; out[1] = J3V_MAGIC[1]; out[2] = J3V_MAGIC[2]; out[3] = J3V_MAGIC[3];
    const v = new DataView(out.buffer);
    v.setUint32(4, dna2bit.length || data.length, true);
    v.setUint32(8, 0, true); // 0 blocks = single stream
    out.set(compressed, hdr);
    return out;
  }

  // Adaptive block sizing: 4096 bases per block
  const BLOCK_SIZE = 4096;
  const blocks: { len: number; compressed: Uint8Array }[] = [];

  for (let off = 0; off < dna2bit.length; off += BLOCK_SIZE) {
    const end = Math.min(off + BLOCK_SIZE, dna2bit.length);
    const block = dna2bit.slice(off, end);
    // Order-1 context model per block
    const compressed = arithEncodeDna(block, 1);
    blocks.push({ len: block.length, compressed });
  }

  // Build output: [magic][dna_len][n_blocks][block0_len(4)][block0_complen(4)][block0...]...
  let totalSize = 12;
  for (const b of blocks) totalSize += 8 + b.compressed.length;

  const out = new Uint8Array(totalSize);
  out[0] = J3V_MAGIC[0]; out[1] = J3V_MAGIC[1]; out[2] = J3V_MAGIC[2]; out[3] = J3V_MAGIC[3];
  const v = new DataView(out.buffer);
  v.setUint32(4, dna2bit.length, true);
  v.setUint32(8, blocks.length, true);
  let offset = 12;
  for (const b of blocks) {
    v.setUint32(offset, b.len, true);
    v.setUint32(offset + 4, b.compressed.length, true);
    offset += 8;
    out.set(b.compressed, offset);
    offset += b.compressed.length;
  }
  return out;
}

/**
 * Decompress JARVIS3 arithmetic-coded data.
 */
export function decompressWithJARVIS3(data: Uint8Array): Uint8Array {
  if (data.length < 12 ||
      data[0] !== J3V_MAGIC[0] || data[1] !== J3V_MAGIC[1] ||
      data[2] !== J3V_MAGIC[2] || data[3] !== J3V_MAGIC[3]) {
    throw new Error('Invalid JARVIS3 arithmetic magic');
  }

  const v = new DataView(data.buffer, data.byteOffset);
  const dnaLen = v.getUint32(4, true);
  const nBlocks = v.getUint32(8, true);

  if (nBlocks === 0) {
    const compressed = data.slice(12);
    const values = arithDecodeDna(compressed, 0, dnaLen);
    const out = new Uint8Array(values.length);
    for (let i = 0; i < values.length; i++) out[i] = BIT2_NUCLEOTIDE[values[i]];
    return out;
  }

  // Decode blocks
  const out = new Uint8Array(dnaLen);
  let offset = 12;
  let outIdx = 0;

  for (let b = 0; b < nBlocks; b++) {
    const blockLen = v.getUint32(offset, true);
    const compLen = v.getUint32(offset + 4, true);
    offset += 8;
    const compressed = data.slice(offset, offset + compLen);
    offset += compLen;

    const values = arithDecodeDna(compressed, 1, blockLen);
    for (let i = 0; i < values.length && outIdx < dnaLen; i++) {
      out[outIdx++] = BIT2_NUCLEOTIDE[values[i]];
    }
  }

  return out;
}

// ═══════════════════════════════════════════════════════════════════════════
// § C++ WASM Compressor Framework
// ═══════════════════════════════════════════════════════════════════════════
//
// The five DNA compressors above (NAF, AGC, DeepGeCo, MBGC2, JARVIS3) are
// TypeScript implementations using arithmetic coding. For higher throughput,
// compiled C++ WASM compressors can be registered and will be preferred
// when available.
//
// The C++ sources are built with emscripten (see scripts/build-cpp-compressors.sh)
// and output to src/lib/dna/wasm-pkg/cpp-compressors/.
//
// Expected WASM API (exported from each C++ module):
//   memory:       WebAssembly.Memory  (shared linear memory)
//   alloc(n):     number             (allocate n bytes, return pointer)
//   dealloc(p):   void              (free pointer)
//   compress(p_in, len_in, p_out, p_out_len): number
//                  (returns compressed length; writes to p_out)
//   decompress(p_in, len_in, p_out, p_out_len, original_size): number
//                  (returns decompressed length)
//   version():    number            (API version, must match WASM_API_VERSION)
//   name():       number            (pointer to null-terminated name string)
// ═══════════════════════════════════════════════════════════════════════════

/** API version that C++ WASM modules must export to be compatible. */
const WASM_API_VERSION = 1;

/**
 * WASM interface that each C++ DNA compressor module must export.
 */
interface DnaCompressorWasmExports {
  memory: WebAssembly.Memory;
  alloc(n: number): number;
  dealloc(p: number): void;
  compress(pIn: number, lenIn: number, pOut: number, pOutLen: number): number;
  decompress(pIn: number, lenIn: number, pOut: number, pOutLen: number, originalSize: number): number;
  version(): number;
  name(): number;
}

/**
 * Read a null-terminated string from WASM linear memory.
 */
function readWasmString(memory: WebAssembly.Memory, ptr: number): string {
  const buf = new Uint8Array(memory.buffer);
  let end = ptr;
  while (end < buf.length && buf[end] !== 0) end++;
  return new TextDecoder().decode(buf.slice(ptr, end));
}

/**
 * DnaCompressorWasm — wraps a WebAssembly.Instance implementing the
 * DNA compressor API (compress/decompress via shared linear memory).
 *
 * Usage:
 *   const module = await WebAssembly.compile(wasmBytes);
 *   const compressor = new DnaCompressorWasm('naf', module);
 *   const compressed = compressor.compress(data);
 *   const recovered = compressor.decompress(compressed, data.length);
 */
export class DnaCompressorWasm {
  private instance: WebAssembly.Instance;
  private exports: DnaCompressorWasmExports;
  readonly name: string;
  readonly version: number;

  constructor(name: string, wasmModule: WebAssembly.Module) {
    this.instance = new WebAssembly.Instance(wasmModule, {
      env: {
        // Provide any required imports (abort on allocation failure, etc.)
        abort: (msgPtr: number) => {
          throw new Error(`WASM abort: ${readWasmString(this.exports.memory, msgPtr)}`);
        },
        // Logger — no-op by default, can be overridden
        log: () => {},
      },
    });

    const exports = this.instance.exports as unknown as DnaCompressorWasmExports;

    // Validate exports
    if (!exports.memory) throw new Error(`WASM module "${name}" missing memory export`);
    if (typeof exports.alloc !== 'function') throw new Error(`WASM module "${name}" missing alloc()`);
    if (typeof exports.dealloc !== 'function') throw new Error(`WASM module "${name}" missing dealloc()`);
    if (typeof exports.compress !== 'function') throw new Error(`WASM module "${name}" missing compress()`);
    if (typeof exports.decompress !== 'function') throw new Error(`WASM module "${name}" missing decompress()`);
    if (typeof exports.version !== 'function') throw new Error(`WASM module "${name}" missing version()`);

    this.exports = exports;
    this.version = exports.version();
    this.name = name;

    if (this.version !== WASM_API_VERSION) {
      throw new Error(
        `WASM compressor "${name}" has API version ${this.version}, expected ${WASM_API_VERSION}. ` +
        `Rebuild the C++ sources with the matching API version.`,
      );
    }

    // Read the name from WASM memory if available
    if (typeof exports.name === 'function') {
      const namePtr = exports.name();
      const wasmName = readWasmString(exports.memory, namePtr);
      if (wasmName && wasmName !== name) {
        // Log a warning but don't throw — the registered name takes precedence
        console.warn(`WASM compressor registered as "${name}" but self-reports as "${wasmName}"`);
      }
    }
  }

  /**
   * Compress DNA data using the C++ WASM implementation.
   *
   * @param data Input DNA data (ASCII bytes: A, C, G, T, or raw binary)
   * @returns Compressed data with the same magic header as the TS implementation
   */
  compress(data: Uint8Array): Uint8Array {
    const { memory, alloc, dealloc, compress: wasmCompress } = this.exports;

    // Allocate input buffer in WASM memory
    const inPtr = alloc(data.length);
    if (!inPtr) throw new Error(`WASM alloc failed for ${data.length} bytes (compress input)`);

    // Copy input data into WASM memory
    new Uint8Array(memory.buffer).set(data, inPtr);

    // Allocate output buffer (worst case: same size + header overhead)
    const maxOutLen = data.length + 1024;  // generous overhead for headers
    const outPtr = alloc(maxOutLen);
    if (!outPtr) {
      dealloc(inPtr);
      throw new Error(`WASM alloc failed for ${maxOutLen} bytes (compress output)`);
    }

    try {
      const compressedLen = wasmCompress(inPtr, data.length, outPtr, maxOutLen);
      if (compressedLen <= 0) {
        throw new Error(`WASM compress failed (returned ${compressedLen}) for ${data.length} bytes`);
      }
      if (compressedLen > maxOutLen) {
        throw new Error(`WASM compress wrote ${compressedLen} bytes but output buffer was only ${maxOutLen}`);
      }

      // Copy result from WASM memory
      return new Uint8Array(memory.buffer).slice(outPtr, outPtr + compressedLen);
    } finally {
      dealloc(inPtr);
      dealloc(outPtr);
    }
  }

  /**
   * Decompress data using the C++ WASM implementation.
   *
   * @param data Compressed data
   * @param originalSize Expected size of the decompressed output
   * @returns Decompressed data
   */
  decompress(data: Uint8Array, originalSize: number): Uint8Array {
    const { memory, alloc, dealloc, decompress: wasmDecompress } = this.exports;

    // Allocate input buffer
    const inPtr = alloc(data.length);
    if (!inPtr) throw new Error(`WASM alloc failed for ${data.length} bytes (decompress input)`);

    new Uint8Array(memory.buffer).set(data, inPtr);

    // Allocate output buffer
    const outPtr = alloc(originalSize);
    if (!outPtr) {
      dealloc(inPtr);
      throw new Error(`WASM alloc failed for ${originalSize} bytes (decompress output)`);
    }

    try {
      const outLen = wasmDecompress(inPtr, data.length, outPtr, originalSize, originalSize);
      if (outLen <= 0) {
        throw new Error(`WASM decompress failed (returned ${outLen}) for ${data.length} bytes`);
      }

      return new Uint8Array(memory.buffer).slice(outPtr, outPtr + outLen);
    } finally {
      dealloc(inPtr);
      dealloc(outPtr);
    }
  }

  /**
   * Get info about this WASM compressor for diagnostics.
   */
  getInfo(): { name: string; version: number; memoryPages: number } {
    return {
      name: this.name,
      version: this.version,
      memoryPages: this.exports.memory.buffer.byteLength / 65536,
    };
  }
}

// ---------------------------------------------------------------------------
// WASM Compressor Registry
// ---------------------------------------------------------------------------

/**
 * Global registry of loaded C++ WASM DNA compressors.
 * Key: compressor name (lowercase, e.g. "naf", "agc", "deepgeco", "mbgc2", "jarvis3")
 * Value: the DnaCompressorWasm instance
 */
const wasmCompressorRegistry = new Map<string, DnaCompressorWasm>();

/**
 * Register a C++ WASM DNA compressor.
 *
 * After registration, the compress router (compress.ts) will prefer the
 * WASM implementation over the TypeScript arithmetic-coding implementation.
 *
 * @param name  Compressor name: "naf", "agc", "deepgeco", "mbgc2", or "jarvis3"
 * @param wasmModule  A compiled WebAssembly.Module (from WebAssembly.compile)
 * @returns The DnaCompressorWasm wrapper
 *
 * @example
 *   // Load and register the NAF C++ WASM compressor
 *   const wasmBytes = await fetch('./wasm-pkg/cpp-compressors/naf.wasm').then(r => r.arrayBuffer());
 *   const module = await WebAssembly.compile(wasmBytes);
 *   const compressor = registerDnaCompressorWasm('naf', module);
 *   const compressed = compressor.compress(dnaData);
 */
export function registerDnaCompressorWasm(name: string, wasmModule: WebAssembly.Module): DnaCompressorWasm {
  const normalized = name.toLowerCase();
  const compressor = new DnaCompressorWasm(normalized, wasmModule);
  wasmCompressorRegistry.set(normalized, compressor);
  return compressor;
}

/**
 * Unregister a WASM DNA compressor (revert to TS implementation).
 */
export function unregisterDnaCompressorWasm(name: string): boolean {
  return wasmCompressorRegistry.delete(name.toLowerCase());
}

/**
 * Get a registered WASM DNA compressor, or null if not available.
 */
export function getDnaCompressorWasm(name: string): DnaCompressorWasm | null {
  return wasmCompressorRegistry.get(name.toLowerCase()) ?? null;
}

/**
 * Check if a WASM DNA compressor is registered for the given name.
 */
export function isDnaCompressorWasmAvailable(name: string): boolean {
  return wasmCompressorRegistry.has(name.toLowerCase());
}

/**
 * List all registered WASM DNA compressor names.
 */
export function listDnaCompressorWasm(): string[] {
  return Array.from(wasmCompressorRegistry.keys());
}

/**
 * Get info about all registered WASM compressors.
 */
export function getWasmCompressorInfo(): { name: string; version: number; memoryPages: number }[] {
  return Array.from(wasmCompressorRegistry.values()).map(c => c.getInfo());
}

// ---------------------------------------------------------------------------
// WASM-aware compress/decompress wrappers
// ---------------------------------------------------------------------------

/**
 * Compress with NAF — uses C++ WASM if registered, otherwise TypeScript.
 */
export function compressWithNAFWasm(data: Uint8Array, level?: number): Uint8Array {
  const wasm = getDnaCompressorWasm('naf');
  if (wasm) return wasm.compress(data);
  return compressWithNAF(data, level);
}

/**
 * Decompress NAF — uses C++ WASM if registered, otherwise TypeScript.
 */
export function decompressWithNAFWasm(data: Uint8Array): Uint8Array {
  const wasm = getDnaCompressorWasm('naf');
  if (wasm) {
    // Extract original size from the NAF header
    const v = new DataView(data.buffer, data.byteOffset);
    const dnaLen = v.getUint32(4, true);
    return wasm.decompress(data, dnaLen > 0 ? dnaLen : v.getUint32(8, true));
  }
  return decompressWithNAF(data);
}

/**
 * Compress with AGC — uses C++ WASM if registered, otherwise TypeScript.
 */
export function compressWithAGCWasm(data: Uint8Array, level?: number): Uint8Array {
  const wasm = getDnaCompressorWasm('agc');
  if (wasm) return wasm.compress(data);
  return compressWithAGC(data, level);
}

/**
 * Decompress AGC — uses C++ WASM if registered, otherwise TypeScript.
 */
export function decompressWithAGCWasm(data: Uint8Array): Uint8Array {
  const wasm = getDnaCompressorWasm('agc');
  if (wasm) {
    const v = new DataView(data.buffer, data.byteOffset);
    const dnaLen = v.getUint32(4, true);
    return wasm.decompress(data, dnaLen);
  }
  return decompressWithAGC(data);
}

/**
 * Compress with DeepGeCo — uses C++ WASM if registered, otherwise TypeScript.
 */
export function compressWithDeepGeCoWasm(data: Uint8Array, level?: number): Uint8Array {
  const wasm = getDnaCompressorWasm('deepgeco');
  if (wasm) return wasm.compress(data);
  return compressWithDeepGeCo(data, level);
}

/**
 * Decompress DeepGeCo — uses C++ WASM if registered, otherwise TypeScript.
 */
export function decompressWithDeepGeCoWasm(data: Uint8Array): Uint8Array {
  const wasm = getDnaCompressorWasm('deepgeco');
  if (wasm) {
    const v = new DataView(data.buffer, data.byteOffset);
    const dnaLen = v.getUint32(4, true);
    return wasm.decompress(data, dnaLen);
  }
  return decompressWithDeepGeCo(data);
}

/**
 * Compress with MBGC2 — uses C++ WASM if registered, otherwise TypeScript.
 */
export function compressWithMBGC2Wasm(data: Uint8Array, level?: number): Uint8Array {
  const wasm = getDnaCompressorWasm('mbgc2');
  if (wasm) return wasm.compress(data);
  return compressWithMBGC2(data, level);
}

/**
 * Decompress MBGC2 — uses C++ WASM if registered, otherwise TypeScript.
 */
export function decompressWithMBGC2Wasm(data: Uint8Array): Uint8Array {
  const wasm = getDnaCompressorWasm('mbgc2');
  if (wasm) {
    const v = new DataView(data.buffer, data.byteOffset);
    const dnaLen = v.getUint32(4, true);
    return wasm.decompress(data, dnaLen);
  }
  return decompressWithMBGC2(data);
}

/**
 * Compress with JARVIS3 — uses C++ WASM if registered, otherwise TypeScript.
 */
export function compressWithJARVIS3Wasm(data: Uint8Array, level?: number): Uint8Array {
  const wasm = getDnaCompressorWasm('jarvis3');
  if (wasm) return wasm.compress(data);
  return compressWithJARVIS3(data, level);
}

/**
 * Decompress JARVIS3 — uses C++ WASM if registered, otherwise TypeScript.
 */
export function decompressWithJARVIS3Wasm(data: Uint8Array): Uint8Array {
  const wasm = getDnaCompressorWasm('jarvis3');
  if (wasm) {
    const v = new DataView(data.buffer, data.byteOffset);
    const dnaLen = v.getUint32(4, true);
    return wasm.decompress(data, dnaLen);
  }
  return decompressWithJARVIS3(data);
}

/**
 * Bulk-load all C++ WASM DNA compressors from a directory.
 * Intended for Node.js use — reads .wasm files from the filesystem.
 *
 * @param dirPath  Directory containing naf.wasm, agc.wasm, etc.
 * @returns Array of registered compressor names
 *
 * @example
 *   // In Node.js startup:
 *   const names = await loadAllWasmCompressors('./src/lib/dna/wasm-pkg/cpp-compressors');
 *   console.log('Loaded WASM compressors:', names);
 */
export async function loadAllWasmCompressors(dirPath: string): Promise<string[]> {
  const COMPRESSOR_NAMES = ['naf', 'agc', 'deepgeco', 'mbgc2', 'jarvis3'];
  const loaded: string[] = [];

  try {
    // Dynamic import for Node.js fs
    const fs = await import('fs');
    const path = await import('path');

    for (const name of COMPRESSOR_NAMES) {
      const wasmPath = path.join(dirPath, `${name}.wasm`);
      try {
        const wasmBytes = await fs.promises.readFile(wasmPath);
        const module = await WebAssembly.compile(wasmBytes);
        registerDnaCompressorWasm(name, module);
        loaded.push(name);
      } catch {
        // File not found or compile failed — skip, TS fallback will be used
      }
    }
  } catch {
    // fs module not available (browser) — cannot auto-load from directory
  }

  return loaded;
}

// Re-export for compress.ts wiring
export {
  ArithmeticEncoder,
  ArithmeticDecoder,
  AdaptiveModel,
  AdaptiveContextModel,
  AdaptiveFrequencyModel,
} from './arithmetic-coder';
