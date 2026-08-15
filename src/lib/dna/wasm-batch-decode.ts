/**
 * WASM Decode — Full-pipeline single-call decode
 *
 * Two entry points:
 *   - batchDecodeAll (legacy): runs only the per-read LDPC/CRC loop in WASM,
 *     returns per-oligo payloads. JS still does outer RS + decompress.
 *   - fullDecode (v5+): runs the ENTIRE pipeline in WASM — primer trimming,
 *     LDPC decode, CRC, outer RS GF(2^16) erasure recovery, DEFLATE inflate,
 *     file-size trim. Returns the final recovered file bytes.
 *
 * fullDecode eliminates ALL of:
 *   - 50K+ JS↔WASM boundary crossings (per-read calls)
 *   - JS outer RS overhead (~400ms for 256KB at 10x)
 *   - pako DEFLATE overhead (~300ms for 256KB compressed)
 */

// @ts-ignore
const wasm = require("./wasm-pkg/helix_dna_wasm.js");

let wasmReady = false;
function ensureWasm() {
  if (wasmReady) return;
  wasmReady = true;
}

export interface BatchDecodeResult {
  /** Flat array: for each oligo, 1 byte status + k bytes payload. */
  data: Uint8Array;
  /** Number of oligos decoded successfully. */
  decoded: number;
  /** Number of oligos erased. */
  erased: number;
}

/**
 * Batch decode ALL reads in a single WASM call (legacy mode — only LDPC/CRC).
 * Use fullDecode() for the complete pipeline including outer RS + DEFLATE.
 */
export function batchDecodeAll(
  reads: any[],
  fwdPrimer: string,
  revPrimer: string,
  oligoCount: number,
  innerN: number,
  innerK: number,
  totalInnerBytes: number,
): BatchDecodeResult {
  ensureWasm();

  let totalLen = 0;
  for (const read of reads) totalLen += read.sequence.length;

  const allReads = new Uint8Array(totalLen);
  const readOffsets = new Uint8Array(reads.length * 4);
  const readLengths = new Uint8Array(reads.length * 4);

  let offset = 0;
  for (let i = 0; i < reads.length; i++) {
    const seq = reads[i].sequence;
    const view = new DataView(readOffsets.buffer, i * 4, 4);
    view.setUint32(0, offset, true);
    const view2 = new DataView(readLengths.buffer, i * 4, 4);
    view2.setUint32(0, seq.length, true);
    for (let j = 0; j < seq.length; j++) {
      allReads[offset + j] = seq.charCodeAt(j);
    }
    offset += seq.length;
  }

  const fwdPrimerBytes = new Uint8Array(fwdPrimer.length);
  for (let i = 0; i < fwdPrimer.length; i++) fwdPrimerBytes[i] = fwdPrimer.charCodeAt(i);
  const revPrimerBytes = new Uint8Array(revPrimer.length);
  for (let i = 0; i < revPrimer.length; i++) revPrimerBytes[i] = revPrimer.charCodeAt(i);

  const result = wasm.batch_decode_all(
    allReads, readOffsets, readLengths,
    fwdPrimerBytes, revPrimerBytes,
    oligoCount, innerN, innerK, totalInnerBytes,
  );

  let decoded = 0;
  let erased = 0;
  // Stride is 1+innerK (1 status byte + innerK bytes per oligo, with only
  // payload_bytes used — the remaining 4 bytes are zero padding for alignment
  // with the JS consumer's stride computation).
  const stride = 1 + innerK;
  for (let i = 0; i < oligoCount; i++) {
    const status = result[i * stride];
    if (status === 1) decoded++;
    else erased++;
  }
  return { data: result, decoded, erased };
}

/**
 * Full pipeline in a single WASM call.
 *
 * Runs (entirely in WASM):
 *   1. Primer trimming + Hamming match (≤2 mismatches)
 *   2. DNA→bytes conversion (4 nt → 1 byte)
 *   3. Cluster reads by oligo index (XOR-unwhitened address)
 *   4. Per-oligo LDPC decode + CRC-16 + address verification + payload unwhiten
 *   5. Outer RS GF(2^16) pure-erasure recovery (Rust impl, matches JS ReedSolomon216)
 *   6. Payload concatenation
 *   7. DEFLATE inflate (miniz_oxide, pako-compatible zlib stream)
 *   8. Trim to file_size
 *
 * Returns: recovered file bytes (Uint8Array).
 */
export function fullDecode(
  reads: any[],
  fwdPrimer: string,
  revPrimer: string,
  oligoCount: number,
  innerN: number,
  innerK: number,
  totalInnerBytes: number,
  outerN: number,
  outerK: number,
  payloadBytes: number,
  fileSize: number,
  useDeflate: boolean,
): Uint8Array {
  ensureWasm();

  // Flatten all reads — v58: use Buffer.from(seq, 'latin1') for native-speed
  // ASCII string → Uint8Array copy (10× faster than charCodeAt loop).
  let totalLen = 0;
  for (const read of reads) totalLen += read.sequence.length;

  const allReads = new Uint8Array(totalLen);
  const readOffsets = new Uint8Array(reads.length * 4);
  const readLengths = new Uint8Array(reads.length * 4);

  // Use a single DataView over the entire offsets/lengths buffers for fast writes
  const offView = new DataView(readOffsets.buffer);
  const lenView = new DataView(readLengths.buffer);

  let offset = 0;
  const isNode = typeof Buffer !== 'undefined';
  for (let i = 0; i < reads.length; i++) {
    const seq = reads[i].sequence;
    offView.setUint32(i * 4, offset, true);
    lenView.setUint32(i * 4, seq.length, true);
    if (isNode) {
      // Node/Bun: Buffer.from(seq, 'latin1') copies ASCII bytes natively
      allReads.set(Buffer.from(seq, 'latin1'), offset);
    } else {
      // Browser fallback: charCodeAt loop (slow path)
      for (let j = 0; j < seq.length; j++) allReads[offset + j] = seq.charCodeAt(j);
    }
    offset += seq.length;
  }

  const fwdPrimerBytes = isNode
    ? Buffer.from(fwdPrimer, 'latin1')
    : Uint8Array.from(fwdPrimer, c => c.charCodeAt(0));
  const revPrimerBytes = isNode
    ? Buffer.from(revPrimer, 'latin1')
    : Uint8Array.from(revPrimer, c => c.charCodeAt(0));

  // Single WASM call — entire pipeline runs in Rust
  return wasm.full_decode(
    allReads, readOffsets, readLengths,
    fwdPrimerBytes, revPrimerBytes,
    oligoCount, innerN, innerK, totalInnerBytes,
    outerN, outerK, payloadBytes, fileSize, useDeflate,
  );
}

/**
 * Full pipeline with interleaving support (single WASM call).
 *
 * Same as fullDecode but deinterleaves consensus blocks across groups of
 * `interleaveDepth` oligos before LDPC decode. Spreads burst errors.
 */
export function fullDecodeInterleaved(
  reads: any[],
  fwdPrimer: string,
  revPrimer: string,
  oligoCount: number,
  innerN: number,
  innerK: number,
  totalInnerBytes: number,
  outerN: number,
  outerK: number,
  payloadBytes: number,
  fileSize: number,
  useDeflate: boolean,
  interleaveDepth: number,
): Uint8Array {
  ensureWasm();

  let totalLen = 0;
  for (const read of reads) totalLen += read.sequence.length;

  const allReads = new Uint8Array(totalLen);
  const readOffsets = new Uint8Array(reads.length * 4);
  const readLengths = new Uint8Array(reads.length * 4);

  let offset = 0;
  const isNodeFlatten = typeof Buffer !== 'undefined';
  for (let i = 0; i < reads.length; i++) {
    const seq = reads[i].sequence;
    new DataView(readOffsets.buffer, i * 4, 4).setUint32(0, offset, true);
    new DataView(readLengths.buffer, i * 4, 4).setUint32(0, seq.length, true);
    if (isNodeFlatten) {
      allReads.set(Buffer.from(seq, 'latin1'), offset);
    } else {
      for (let j = 0; j < seq.length; j++) allReads[offset + j] = seq.charCodeAt(j);
    }
    offset += seq.length;
  }

  const fwdPrimerBytes = isNodeFlatten ? Buffer.from(fwdPrimer, 'latin1') : Uint8Array.from(fwdPrimer, c => c.charCodeAt(0));
  const revPrimerBytes = isNodeFlatten ? Buffer.from(revPrimer, 'latin1') : Uint8Array.from(revPrimer, c => c.charCodeAt(0));

  return wasm.full_decode_interleaved(
    allReads, readOffsets, readLengths,
    fwdPrimerBytes, revPrimerBytes,
    oligoCount, innerN, innerK, totalInnerBytes,
    outerN, outerK, payloadBytes, fileSize, useDeflate, interleaveDepth,
  );
}

/**
 * Full pipeline with ARITHMETIC mapping mode + per-block CRC sync markers.
 * Uses full_decode_arithmetic_crc — the DNA-Aeon approach with per-block
 * CRC-8 verification. Blocks with failed CRC are skipped, and the read is
 * treated as erased (other reads of the same oligo provide the data).
 */
export function fullDecodeArithmetic(
  reads: any[],
  fwdPrimer: string,
  revPrimer: string,
  oligoCount: number,
  innerN: number,
  innerK: number,
  totalInnerBytes: number,
  outerN: number,
  outerK: number,
  payloadBytes: number,
  fileSize: number,
  useDeflate: boolean,
  maxHomopolymer: number,
  blockSize: number,
): Uint8Array {
  ensureWasm();

  let totalLen = 0;
  for (const read of reads) totalLen += read.sequence.length;

  const allReads = new Uint8Array(totalLen);
  const readOffsets = new Uint8Array(reads.length * 4);
  const readLengths = new Uint8Array(reads.length * 4);

  let offset = 0;
  const isNodeArith = typeof Buffer !== 'undefined';
  for (let i = 0; i < reads.length; i++) {
    const seq = reads[i].sequence;
    new DataView(readOffsets.buffer, i * 4, 4).setUint32(0, offset, true);
    new DataView(readLengths.buffer, i * 4, 4).setUint32(0, seq.length, true);
    if (isNodeArith) {
      allReads.set(Buffer.from(seq, 'latin1'), offset);
    } else {
      for (let j = 0; j < seq.length; j++) allReads[offset + j] = seq.charCodeAt(j);
    }
    offset += seq.length;
  }

  const fwdPrimerBytes = isNodeArith ? Buffer.from(fwdPrimer, 'latin1') : Uint8Array.from(fwdPrimer, c => c.charCodeAt(0));
  const revPrimerBytes = isNodeArith ? Buffer.from(revPrimer, 'latin1') : Uint8Array.from(revPrimer, c => c.charCodeAt(0));

  return wasm.full_decode_arithmetic_crc(
    allReads, readOffsets, readLengths,
    fwdPrimerBytes, revPrimerBytes,
    oligoCount, innerN, innerK, totalInnerBytes,
    outerN, outerK, payloadBytes, fileSize, useDeflate, maxHomopolymer, blockSize,
  );
}

/**
 * Phase 3-6 only: takes batch_decode_all output + does outer RS + DEFLATE + trim.
 * Used by the parallel decode pipeline (workers run batch_decode_all, main thread
 * merges results and calls this to finish).
 */
export function outerRsInflate(
  batchResult: Uint8Array,
  oligoCount: number,
  innerK: number,
  payloadBytes: number,
  outerN: number,
  outerK: number,
  fileSize: number,
  useDeflate: boolean,
): Uint8Array {
  ensureWasm();
  return wasm.outer_rs_inflate(
    batchResult, oligoCount, innerK, payloadBytes,
    outerN, outerK, fileSize, useDeflate,
  );
}

/**
 * Flatten reads into Uint8Array buffers for sending to workers.
 * Returns the flat arrays + offsets/lengths.
 */
export function flattenReads(reads: any[]): {
  allReads: Uint8Array;
  readOffsets: Uint8Array;
  readLengths: Uint8Array;
} {
  let totalLen = 0;
  for (const read of reads) totalLen += read.sequence.length;

  const allReads = new Uint8Array(totalLen);
  const readOffsets = new Uint8Array(reads.length * 4);
  const readLengths = new Uint8Array(reads.length * 4);

  let offset = 0;
  const isNodeFlat = typeof Buffer !== 'undefined';
  for (let i = 0; i < reads.length; i++) {
    const seq = reads[i].sequence;
    new DataView(readOffsets.buffer, i * 4, 4).setUint32(0, offset, true);
    new DataView(readLengths.buffer, i * 4, 4).setUint32(0, seq.length, true);
    if (isNodeFlat) {
      allReads.set(Buffer.from(seq, 'latin1'), offset);
    } else {
      for (let j = 0; j < seq.length; j++) allReads[offset + j] = seq.charCodeAt(j);
    }
    offset += seq.length;
  }

  return { allReads, readOffsets, readLengths };
}

/**
 * Run batch_decode_all on a subset of reads (for parallel decode).
 * This is the per-worker function.
 */
export function batchDecodeChunk(
  allReads: Uint8Array,
  readOffsets: Uint8Array,
  readLengths: Uint8Array,
  fwdPrimerBytes: Uint8Array,
  revPrimerBytes: Uint8Array,
  oligoCount: number,
  innerN: number,
  innerK: number,
  totalInnerBytes: number,
): Uint8Array {
  ensureWasm();
  return wasm.batch_decode_all(
    allReads, readOffsets, readLengths,
    fwdPrimerBytes, revPrimerBytes,
    oligoCount, innerN, innerK, totalInnerBytes,
  );
}

export { ensureWasm };
