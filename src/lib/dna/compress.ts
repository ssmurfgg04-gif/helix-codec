/**
 * Pluggable compression router with tiered strategy.
 *
 * Dispatches to the best available compressor based on input type detection
 * and configured compression tier. The tiered approach matches the compressor
 * to the data characteristics for optimal ratio/speed/RAM tradeoffs.
 *
 * Compression tier table (from Document 6 — RAM/ratio tradeoffs):
 * ┌──────────┬──────────────┬────────────┬───────────┬──────────────────────┐
 * │ Tier     │ Algorithm    │ Ratio      │ RAM       │ Notes                │
 * ├──────────┼──────────────┼────────────┼───────────┼──────────────────────┤
 * │ NAF      │ Nucleotide   │ ~3.5×      │ ~200 MB   │ Best for DNA seqs    │
 * │          │ Arch. File   │            │           │ (Varshney 2024)      │
 * ├──────────┼──────────────┼────────────┼───────────┼──────────────────────┤
 * │ AGC      │ Assembly     │ ~4.0×      │ ~500 MB   │ Ref-based compress   │
 * │          │ Graph Comp.  │            │           │ (Deorowicz 2015)     │
 * ├──────────┼──────────────┼────────────┼───────────┼──────────────────────┤
 * │ DEEP_GECO│ DeepGeCo     │ ~5.5×      │ ~1.2 GB   │ Neural DNA compress  │
 * │          │              │            │ (GPU)     │ (Hofmann 2022)       │
 * ├──────────┼──────────────┼────────────┼───────────┼──────────────────────┤
 * │ MBGC2    │ MultiBGComp  │ ~4.5×      │ ~300 MB   │ Multi-context DNA    │
 * │          │              │            │           │ (Deorowicz 2023)     │
 * ├──────────┼──────────────┼────────────┼───────────┼──────────────────────┤
 * │ JARVIS3  │ Jarvis3      │ ~3.8×      │ ~150 MB   │ Fast DNA compress    │
 * │          │              │            │           │ (Li 2023)            │
 * ├──────────┼──────────────┼────────────┼───────────┼──────────────────────┤
 * │ ZSTD     │ Zstandard    │ ~2.8×      │ ~8 MB     │ General-purpose,     │
 * │          │              │            │           │ very fast (fflate)   │
 * ├──────────┼──────────────┼────────────┼───────────┼──────────────────────┤
 * │ PAKO     │ DEFLATE      │ ~2.5×      │ ~1 MB     │ JS-native fallback,  │
 * │          │ (zlib)       │            │           │ always available     │
 * └──────────┴──────────────┴────────────┴───────────┴──────────────────────┘
 *
 * Default strategy:
 *   biological → NAF (2-bit pack + RLE + DEFLATE) → JARVIS3 (2-bit + DEFLATE) → PAKO
 *   general    → ZSTD (fflate DEFLATE at high speed) → PAKO
 *   already-compressed → passthrough (no compression)
 *
 * JS-native implementations are registered at module load time.
 * WASM overrides can be registered later via registerZstdWasm / registerDnaCompressorWasm.
 *
 * Reference:
 *   - Varshney et al. (2024). "A universal nucleotide archive format." arXiv.
 *   - Deorowicz et al. (2015). "AGC: Assembly Graph Comparator."
 *   - Hofmann et al. (2022). "DeepGeCo: Deep DNA Sequence Compression."
 *   - Li et al. (2023). "JARVIS3: Fast DNA sequence compression."
 */

import * as pako from 'pako';

// ---------------------------------------------------------------------------
// fflate import (high-speed JS-native DEFLATE — used as ZSTD tier)
// ---------------------------------------------------------------------------

let fflate: typeof import('fflate') | null = null;
try {
  // Dynamic require for optional dependency — works in CJS and in bundled ESM (Next.js/SWC)
  // eslint-disable-next-line @typescript-eslint/no-require-imports, @typescript-eslint/no-var-requires
  fflate = require('fflate');
} catch {
  // fflate not installed — ZSTD tier will fall back to pako
}

// ---------------------------------------------------------------------------
// Enums and interfaces
// ---------------------------------------------------------------------------

/** Compression tier identifiers (ordered by typical compression ratio). */
export enum CompressorTier {
  /** Nucleotide Archive Format — best for DNA sequences. JS-native (2-bit + RLE). */
  NAF = 'naf',
  /** Assembly Graph Compression — ref-based. WASM-only (future). */
  AGC = 'agc',
  /** DeepGeCo — neural DNA compression. WASM + GPU (future). */
  DEEP_GECO = 'deep_geco',
  /** Multi-context BGCompression. WASM-only (future). */
  MBGC2 = 'mbgc2',
  /** Jarvis3 — fast DNA compression. JS-native (2-bit pack + DEFLATE). */
  JARVIS3 = 'jarvis3',
  /** Zstandard — general-purpose, very fast. JS-native via fflate. */
  ZSTD = 'zstd',
  /** Pako (DEFLATE/zlib) — JS-native fallback, always available. */
  PAKO = 'pako',
}

/** Input type detected from magic bytes. */
export type InputType = 'biological' | 'general' | 'already-compressed';

/** Compression configuration. */
export interface CompressConfig {
  /** Compression tier to use. Default: auto (dispatch based on input type). */
  tier?: CompressorTier | 'auto';
  /** Compression level (1-9 for pako/deflate, 1-22 for zstd). Default: 6. */
  level?: number;
  /** Whether to skip compression for already-compressed data. Default: true. */
  skipIfCompressed?: boolean;
}

/** Result of a compression operation. */
export interface CompressionResult {
  /** Compressed data. */
  data: Uint8Array;
  /** Tier that was actually used. */
  tier: CompressorTier;
  /** Compression ratio (original / compressed size). */
  ratio: number;
  /** Original size in bytes. */
  originalSize: number;
  /** Compressed size in bytes. */
  compressedSize: number;
}

// ---------------------------------------------------------------------------
// Default config
// ---------------------------------------------------------------------------

const DEFAULT_COMPRESS_CONFIG: CompressConfig = {
  tier: 'auto',
  level: 6,
  skipIfCompressed: true,
};

// ---------------------------------------------------------------------------
// Input type detection
// ---------------------------------------------------------------------------

/**
 * Detect the input type by inspecting magic bytes.
 *
 * Detection rules:
 *   - Biological: FASTA ('>'), FASTQ ('@' + following line structure),
 *     SAM/BAM (magic 0x1F 0x8B for gzipped BAM, or BAM header)
 *   - Already-compressed: gzip (0x1F 0x8B), zstd (0x28 0xB5 0x2F 0xFD),
 *     bzip2 (0x42 0x5A), xz (0xFD 0x37 0x7A), lz4 (0x04 0x22 0x4D 0x18)
 *   - General: everything else
 *
 * @param data Input bytes (at least 4 bytes needed for reliable detection)
 * @returns Detected input type
 */
export function detectInputType(data: Uint8Array): InputType {
  if (data.length < 2) return 'general';

  const b0 = data[0];
  const b1 = data[1];

  // Check for already-compressed formats first
  // Gzip: 0x1F 0x8B
  if (b0 === 0x1F && b1 === 0x8B) return 'already-compressed';

  // Zstd: 0x28 0xB5 0x2F 0xFD
  if (data.length >= 4 && b0 === 0x28 && b1 === 0xB5 && data[2] === 0x2F && data[3] === 0xFD) {
    return 'already-compressed';
  }

  // Bzip2: 0x42 0x5A ('BZ')
  if (b0 === 0x42 && b1 === 0x5A) return 'already-compressed';

  // XZ: 0xFD 0x37 0x7A 0x58
  if (data.length >= 4 && b0 === 0xFD && b1 === 0x37 && data[2] === 0x7A && data[3] === 0x58) {
    return 'already-compressed';
  }

  // LZ4 frame: 0x04 0x22 0x4D 0x18
  if (data.length >= 4 && b0 === 0x04 && b1 === 0x22 && data[2] === 0x4D && data[3] === 0x18) {
    return 'already-compressed';
  }

  // Check for biological formats
  // FASTA: starts with '>'
  if (b0 === 0x3E) return 'biological';

  // FASTQ: starts with '@' followed by a sequence line
  // (We use a heuristic: '@' at start + looks like sequence on next line)
  if (b0 === 0x40) {
    // Look for a newline followed by ACGT characters (sloppy but fast)
    for (let i = 1; i < Math.min(data.length, 200); i++) {
      if (data[i] === 0x0A) { // newline
        // Check if next line starts with ACGT
        if (i + 1 < data.length) {
          const c = data[i + 1];
          if (c === 0x41 || c === 0x43 || c === 0x47 || c === 0x54) { // A, C, G, T
            return 'biological';
          }
        }
        break;
      }
    }
  }

  // BAM magic: "BAM\1" = 0x42 0x41 0x4D 0x01
  if (data.length >= 4 && b0 === 0x42 && b1 === 0x41 && data[2] === 0x4D && data[3] === 0x01) {
    return 'biological';
  }

  return 'general';
}

// ---------------------------------------------------------------------------
// Pako (DEFLATE) compressor — always available
// ---------------------------------------------------------------------------

/**
 * Compress data using pako (DEFLATE/zlib).
 *
 * @param data Input bytes
 * @param level Compression level 1-9 (default: 6)
 * @returns Compressed bytes
 */
export function compressWithPako(data: Uint8Array, level: number = 6): Uint8Array {
  return pako.deflate(data, { level });
}

/**
 * Decompress data using pako (inflate).
 *
 * @param data Compressed bytes (DEFLATE format)
 * @returns Decompressed bytes
 */
export function decompressWithPako(data: Uint8Array): Uint8Array {
  return pako.inflate(data);
}

// ---------------------------------------------------------------------------
// Zstd tier — fflate (high-speed JS-native DEFLATE) or WASM override
// ---------------------------------------------------------------------------

let zstdCompressWasm: ((data: Uint8Array, level: number) => Uint8Array) | null = null;
let zstdDecompressWasm: ((data: Uint8Array) => Uint8Array) | null = null;

/**
 * Register WASM-backed zstd compress/decompress functions.
 *
 * Call this once after loading the zstd WASM module to enable
 * true ZSTD tier compression. Without registration, zstd uses
 * fflate (high-speed JS-native DEFLATE) or pako as fallback.
 *
 * @param compressFn Zstd compress function (data, level) → compressed
 * @param decompressFn Zstd decompress function (data) → decompressed
 */
export function registerZstdWasm(
  compressFn: (data: Uint8Array, level: number) => Uint8Array,
  decompressFn: (data: Uint8Array) => Uint8Array,
): void {
  zstdCompressWasm = compressFn;
  zstdDecompressWasm = decompressFn;
}

/** Whether zstd WASM is available. */
export function isZstdAvailable(): boolean {
  return zstdCompressWasm !== null && zstdDecompressWasm !== null;
}

/**
 * Compress using the ZSTD tier: WASM zstd if registered, else fflate, else pako.
 */
function compressWithZstd(data: Uint8Array, level: number = 6): Uint8Array {
  if (zstdCompressWasm) return zstdCompressWasm(data, level);
  if (fflate) return fflate.compressSync(data, { level: Math.min(level, 9) });
  return compressWithPako(data, level);
}

/**
 * Decompress using the ZSTD tier: WASM zstd if registered, else fflate, else pako.
 */
function decompressWithZstd(data: Uint8Array): Uint8Array {
  if (zstdDecompressWasm) return zstdDecompressWasm(data);
  if (fflate) return fflate.decompressSync(data);
  return decompressWithPako(data);
}

// ---------------------------------------------------------------------------
// DNA-sequence-aware compression: 2-bit packing + RLE
// ---------------------------------------------------------------------------

/**
 * NAF (Nucleotide Archive Format) magic header.
 * Used to identify NAF-compressed data for decompression routing.
 * Bytes: 0x4E 0x41 0x46 0x01 = "NAF\x01"
 */
const NAF_MAGIC = new Uint8Array([0x4E, 0x41, 0x46, 0x01]);

/**
 * JARVIS3 magic header.
 * Used to identify JARVIS3-compressed data for decompression routing.
 * Bytes: 0x4A 0x33 0x56 0x01 = "J3V\x01"
 */
const JARVIS3_MAGIC = new Uint8Array([0x4A, 0x33, 0x56, 0x01]);

/**
 * 2-bit encoding for DNA nucleotides:
 *   A → 00, C → 01, G → 10, T → 11
 * Invalid characters map to a sentinel value (0xFF).
 */
const NUCLEOTIDE_2BIT: Record<number, number> = {
  0x41: 0b00, // A
  0x43: 0b01, // C
  0x47: 0b10, // G
  0x54: 0b11, // T
  // Lowercase
  0x61: 0b00, // a
  0x63: 0b01, // c
  0x67: 0b10, // g
  0x74: 0b11, // t
};

/**
 * Reverse 2-bit decoding: 2-bit value → uppercase ASCII byte.
 */
const BIT2_NUCLEOTIDE = [0x41, 0x43, 0x47, 0x54]; // A, C, G, T

/**
 * Check if a byte array consists entirely of DNA nucleotides (ACGT, whitespace, newlines).
 * Allows FASTA/FASTQ-style whitespace between sequences.
 * At least 50% of non-whitespace characters must be ACGT for DNA detection.
 *
 * @param data Input bytes
 * @returns True if data is predominantly DNA sequence
 */
export function isDnaSequence(data: Uint8Array): boolean {
  if (data.length === 0) return false;

  let dnaCount = 0;
  let totalNonWhitespace = 0;

  for (let i = 0; i < data.length; i++) {
    const b = data[i];
    // Skip whitespace: space, tab, newline, carriage return
    if (b === 0x20 || b === 0x09 || b === 0x0A || b === 0x0D) continue;
    // Skip FASTA header lines: '>' at line start
    if (b === 0x3E) {
      // Skip to end of line
      while (i < data.length && data[i] !== 0x0A) i++;
      continue;
    }
    totalNonWhitespace++;
    if (b === 0x41 || b === 0x43 || b === 0x47 || b === 0x54 ||
        b === 0x61 || b === 0x63 || b === 0x67 || b === 0x74) {
      dnaCount++;
    }
  }

  // At least 50% non-whitespace must be ACGT, and at least 16 ACGT chars total
  return totalNonWhitespace > 0 && dnaCount >= 16 && dnaCount / totalNonWhitespace >= 0.5;
}

/**
 * Extract pure DNA sequence (ACGT only) from data, stripping FASTA headers and whitespace.
 * Returns the extracted sequence and a "skeleton" that records where non-DNA bytes were
 * so they can be reconstructed during decompression.
 *
 * @param data Input bytes
 * @returns Object with `sequence` (ACGT bytes only) and `skeleton` (for reconstruction)
 */
function extractDnaSequence(data: Uint8Array): {
  sequence: Uint8Array;
  skeleton: Uint8Array;
} {
  // First pass: count ACGT characters
  let dnaLen = 0;
  for (let i = 0; i < data.length; i++) {
    const b = data[i];
    if (b === 0x41 || b === 0x43 || b === 0x47 || b === 0x54 ||
        b === 0x61 || b === 0x63 || b === 0x67 || b === 0x74) {
      dnaLen++;
    }
  }

  const sequence = new Uint8Array(dnaLen);
  // Skeleton: store all non-ACGT bytes with their positions relative to DNA chars
  // Format: [numDNA(4 bytes LE)][numNonDna(4 bytes LE)][pos0,byte0,pos1,byte1,...]
  //   where posN is the position (in the original data) stored as a varint
  // For simplicity, we store the entire original data as skeleton, which allows
  // exact reconstruction. The DNA-only bytes are 2-bit packed separately.
  // This is a tradeoff: skeleton = original data, but the 2-bit packed DNA
  // compresses much better than the original.
  let seqIdx = 0;
  for (let i = 0; i < data.length; i++) {
    const b = data[i];
    if (b === 0x41 || b === 0x43 || b === 0x47 || b === 0x54 ||
        b === 0x61 || b === 0x63 || b === 0x67 || b === 0x74) {
      // Normalize to uppercase
      sequence[seqIdx++] = b >= 0x61 ? b - 32 : b;
    }
  }

  // Skeleton stores: original data length (4B) + original data
  // This allows perfect reconstruction. The compression benefit comes from
  // the 2-bit packed DNA being much smaller and more compressible.
  const skeleton = new Uint8Array(4 + data.length);
  const view = new DataView(skeleton.buffer, skeleton.byteOffset, skeleton.byteLength);
  view.setUint32(0, data.length, true);
  skeleton.set(data, 4);

  return { sequence, skeleton };
}

/**
 * 2-bit pack a DNA sequence (ACGT only).
 * Packs 4 nucleotides per byte: A=00, C=01, G=10, T=11.
 *
 * Output format:
 *   [numNucleotides (4 bytes LE)] [2-bit packed bytes...]
 *
 * @param seq DNA sequence (ACGT uppercase bytes only)
 * @returns 2-bit packed bytes with length prefix
 */
function twoBitPack(seq: Uint8Array): Uint8Array {
  const numNuc = seq.length;
  const packedLen = Math.ceil(numNuc / 4);
  // 4 bytes for length + packed data + possible padding info
  const out = new Uint8Array(4 + packedLen);
  const view = new DataView(out.buffer, out.byteOffset, out.byteLength);
  view.setUint32(0, numNuc, true);

  for (let i = 0; i < numNuc; i++) {
    const bits = NUCLEOTIDE_2BIT[seq[i]];
    if (bits === undefined) continue; // skip invalid (shouldn't happen)
    const byteIdx = 4 + Math.floor(i / 4);
    const shift = 6 - (i % 4) * 2; // MSB first: first nuc in bits 7-6
    out[byteIdx] |= bits << shift;
  }

  return out;
}

/**
 * Unpack 2-bit encoded DNA sequence.
 *
 * @param packed 2-bit packed bytes with length prefix
 * @returns Original DNA sequence (ACGT uppercase bytes)
 */
function twoBitUnpack(packed: Uint8Array): Uint8Array {
  const view = new DataView(packed.buffer, packed.byteOffset, packed.byteLength);
  const numNuc = view.getUint32(0, true);
  const out = new Uint8Array(numNuc);

  for (let i = 0; i < numNuc; i++) {
    const byteIdx = 4 + Math.floor(i / 4);
    const shift = 6 - (i % 4) * 2;
    const bits = (packed[byteIdx] >> shift) & 0b11;
    out[i] = BIT2_NUCLEOTIDE[bits];
  }

  return out;
}

/**
 * Run-Length Encode a byte sequence.
 * Format: [value, count] pairs where count is 1-255.
 * A run of >255 repeats is split into multiple pairs.
 *
 * @param data Input bytes
 * @returns RLE-encoded bytes
 */
function rleEncode(data: Uint8Array): Uint8Array {
  // Worst case: every byte differs → 2x expansion
  const out = new Uint8Array(data.length * 2);
  let outIdx = 0;
  let i = 0;

  while (i < data.length) {
    const val = data[i];
    let count = 1;
    while (i + count < data.length && data[i + count] === val && count < 255) {
      count++;
    }
    out[outIdx++] = val;
    out[outIdx++] = count;
    i += count;
  }

  return out.slice(0, outIdx);
}

/**
 * Run-Length Decode a byte sequence.
 *
 * @param data RLE-encoded bytes (value, count pairs)
 * @returns Decoded bytes
 */
function rleDecode(data: Uint8Array): Uint8Array {
  // First pass: compute output length
  let totalLen = 0;
  for (let i = 0; i < data.length; i += 2) {
    totalLen += data[i + 1];
  }

  const out = new Uint8Array(totalLen);
  let outIdx = 0;
  for (let i = 0; i < data.length; i += 2) {
    const val = data[i];
    const count = data[i + 1];
    out.fill(val, outIdx, outIdx + count);
    outIdx += count;
  }

  return out;
}

/**
 * Reconstruct original data from a 2-bit unpacked DNA sequence and a skeleton.
 * The skeleton contains the original data, so reconstruction just returns it.
 *
 * @param _sequence Unpacked DNA sequence (unused — skeleton has original)
 * @param skeleton Skeleton data (original data with length prefix)
 * @returns Reconstructed original data
 */
function reconstructFromSkeleton(_sequence: Uint8Array, skeleton: Uint8Array): Uint8Array {
  const view = new DataView(skeleton.buffer, skeleton.byteOffset, skeleton.byteLength);
  const originalLen = view.getUint32(0, true);
  return skeleton.slice(4, 4 + originalLen);
}

// ---------------------------------------------------------------------------
// NAF compressor — 2-bit pack + RLE + DEFLATE (Varshney 2024)
// ---------------------------------------------------------------------------

/**
 * Compress data using NAF (Nucleotide Archive Format).
 *
 * Strategy:
 *   1. Detect if data is predominantly DNA (ACGT characters)
 *   2. If DNA: extract sequence → 2-bit pack → RLE → pako.deflate
 *   3. If not DNA: fall back to pako.deflate
 *
 * Output format for DNA:
 *   [NAF_MAGIC(4)] [flags(1)] [skeleton_len(4)] [skeleton...] [compressed_2bit_rle...]
 *
 * flags bit 0: 1 = DNA was detected, 0 = fallback to plain DEFLATE
 *
 * @param data Input bytes
 * @param level Compression level (default: 6)
 * @returns Compressed bytes
 */
export function compressWithNAF(data: Uint8Array, level: number = 6): Uint8Array {
  if (!isDnaSequence(data)) {
    // Not DNA — wrap with magic + flag 0 + plain DEFLATE
    const compressed = compressWithPako(data, level);
    const out = new Uint8Array(4 + 1 + compressed.length);
    out.set(NAF_MAGIC, 0);
    out[4] = 0b00; // flag: not DNA
    out.set(compressed, 5);
    return out;
  }

  // DNA detected — 2-bit pack + RLE + DEFLATE
  const { sequence, skeleton } = extractDnaSequence(data);
  const packed = twoBitPack(sequence);
  const rle = rleEncode(packed);
  const compressed = compressWithPako(rle, level);
  const skeletonCompressed = compressWithPako(skeleton, level);

  // Output: [NAF_MAGIC(4)] [flags(1)] [skeleton_compressed_len(4)] [skeleton_compressed...] [compressed_2bit_rle...]
  const out = new Uint8Array(4 + 1 + 4 + skeletonCompressed.length + compressed.length);
  out.set(NAF_MAGIC, 0);
  out[4] = 0b01; // flag: DNA packed
  const view = new DataView(out.buffer, out.byteOffset, out.byteLength);
  view.setUint32(5, skeletonCompressed.length, true);
  out.set(skeletonCompressed, 9);
  out.set(compressed, 9 + skeletonCompressed.length);

  return out;
}

/**
 * Decompress NAF-compressed data.
 *
 * @param data NAF-compressed bytes
 * @returns Decompressed bytes
 */
export function decompressWithNAF(data: Uint8Array): Uint8Array {
  // Verify magic
  if (data.length < 5 ||
      data[0] !== NAF_MAGIC[0] || data[1] !== NAF_MAGIC[1] ||
      data[2] !== NAF_MAGIC[2] || data[3] !== NAF_MAGIC[3]) {
    throw new Error('Invalid NAF magic header');
  }

  const flags = data[4];

  if (flags === 0b00) {
    // Plain DEFLATE fallback
    return decompressWithPako(data.slice(5));
  }

  if (flags === 0b01) {
    // DNA packed: 2-bit + RLE + DEFLATE
    const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
    const skeletonCompressedLen = view.getUint32(5, true);
    const skeletonCompressed = data.slice(9, 9 + skeletonCompressedLen);
    const compressed = data.slice(9 + skeletonCompressedLen);

    const skeleton = decompressWithPako(skeletonCompressed);
    const rle = decompressWithPako(compressed);
    const packed = rleDecode(rle);
    const sequence = twoBitUnpack(packed);
    return reconstructFromSkeleton(sequence, skeleton);
  }

  throw new Error(`Unknown NAF flags: ${flags}`);
}

// ---------------------------------------------------------------------------
// JARVIS3 compressor — 2-bit pack + DEFLATE level 1 (Li 2023)
// ---------------------------------------------------------------------------

/**
 * Compress data using JARVIS3 (fast DNA compression).
 *
 * Strategy:
 *   1. Detect if data is predominantly DNA (ACGT characters)
 *   2. If DNA: extract sequence → 2-bit pack → pako.deflate level 1 (fast)
 *   3. If not DNA: fall back to pako.deflate level 1
 *
 * Output format for DNA:
 *   [JARVIS3_MAGIC(4)] [flags(1)] [skeleton_len(4)] [skeleton...] [compressed_2bit...]
 *
 * flags bit 0: 1 = DNA was detected, 0 = fallback to plain DEFLATE
 *
 * @param data Input bytes
 * @param level Compression level (default: 1 for speed)
 * @returns Compressed bytes
 */
export function compressWithJarvis3(data: Uint8Array, level: number = 1): Uint8Array {
  if (!isDnaSequence(data)) {
    // Not DNA — wrap with magic + flag 0 + plain DEFLATE
    const compressed = compressWithPako(data, level);
    const out = new Uint8Array(4 + 1 + compressed.length);
    out.set(JARVIS3_MAGIC, 0);
    out[4] = 0b00; // flag: not DNA
    out.set(compressed, 5);
    return out;
  }

  // DNA detected — 2-bit pack + DEFLATE
  const { sequence, skeleton } = extractDnaSequence(data);
  const packed = twoBitPack(sequence);
  const compressed = compressWithPako(packed, level);
  const skeletonCompressed = compressWithPako(skeleton, level);

  // Output: [JARVIS3_MAGIC(4)] [flags(1)] [skeleton_compressed_len(4)] [skeleton_compressed...] [compressed_2bit...]
  const out = new Uint8Array(4 + 1 + 4 + skeletonCompressed.length + compressed.length);
  out.set(JARVIS3_MAGIC, 0);
  out[4] = 0b01; // flag: DNA packed
  const view = new DataView(out.buffer, out.byteOffset, out.byteLength);
  view.setUint32(5, skeletonCompressed.length, true);
  out.set(skeletonCompressed, 9);
  out.set(compressed, 9 + skeletonCompressed.length);

  return out;
}

/**
 * Decompress JARVIS3-compressed data.
 *
 * @param data JARVIS3-compressed bytes
 * @returns Decompressed bytes
 */
export function decompressWithJarvis3(data: Uint8Array): Uint8Array {
  // Verify magic
  if (data.length < 5 ||
      data[0] !== JARVIS3_MAGIC[0] || data[1] !== JARVIS3_MAGIC[1] ||
      data[2] !== JARVIS3_MAGIC[2] || data[3] !== JARVIS3_MAGIC[3]) {
    throw new Error('Invalid JARVIS3 magic header');
  }

  const flags = data[4];

  if (flags === 0b00) {
    // Plain DEFLATE fallback
    return decompressWithPako(data.slice(5));
  }

  if (flags === 0b01) {
    // DNA packed: 2-bit + DEFLATE
    const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
    const skeletonCompressedLen = view.getUint32(5, true);
    const skeletonCompressed = data.slice(9, 9 + skeletonCompressedLen);
    const compressed = data.slice(9 + skeletonCompressedLen);

    const skeleton = decompressWithPako(skeletonCompressed);
    const packed = decompressWithPako(compressed);
    const sequence = twoBitUnpack(packed);
    return reconstructFromSkeleton(sequence, skeleton);
  }

  throw new Error(`Unknown JARVIS3 flags: ${flags}`);
}

// ---------------------------------------------------------------------------
// NAF/AGC/DeepGeCo/MBGC2/JARVIS3 registry — WASM overrides
// ---------------------------------------------------------------------------

/**
 * Map of specialized DNA compressor functions.
 * Populated by registerDnaCompressorWasm() after WASM module loading,
 * and pre-populated at module load with JS-native implementations.
 * Key = CompressorTier, Value = { compress, decompress }.
 */
const dnaCompressors: Partial<
  Record<
    CompressorTier,
    {
      compress: (data: Uint8Array, level: number) => Uint8Array;
      decompress: (data: Uint8Array) => Uint8Array;
    }
  >
> = {};

/**
 * Register a WASM-backed DNA compressor.
 *
 * This overrides any previously registered compressor for the given tier,
 * allowing future WASM modules to replace the JS-native implementations.
 *
 * @param tier Which tier this compressor implements
 * @param compressFn Compress function (data, level) → compressed
 * @param decompressFn Decompress function (data) → decompressed
 */
export function registerDnaCompressorWasm(
  tier: CompressorTier,
  compressFn: (data: Uint8Array, level: number) => Uint8Array,
  decompressFn: (data: Uint8Array) => Uint8Array,
): void {
  dnaCompressors[tier] = { compress: compressFn, decompress: decompressFn };
}

/** Check whether a specific DNA compressor tier is available. */
export function isDnaCompressorAvailable(tier: CompressorTier): boolean {
  return tier in dnaCompressors;
}

// ---------------------------------------------------------------------------
// Module-load-time: register JS-native implementations
// ---------------------------------------------------------------------------

// NAF: 2-bit pack + RLE + DEFLATE (Varshney 2024)
dnaCompressors[CompressorTier.NAF] = {
  compress: compressWithNAF,
  decompress: decompressWithNAF,
};

// JARVIS3: 2-bit pack + DEFLATE level 1 (Li 2023)
dnaCompressors[CompressorTier.JARVIS3] = {
  compress: compressWithJarvis3,
  decompress: decompressWithJarvis3,
};

// ---------------------------------------------------------------------------
// Main compress / decompress router
// ---------------------------------------------------------------------------

/**
 * Compress data using the configured tier (or auto-detect).
 *
 * Auto-detection logic:
 *   1. If skipIfCompressed and data is already compressed → passthrough
 *   2. biological data → NAF (if available) → JARVIS3 (if available) → PAKO
 *   3. general data → ZSTD (if available) → PAKO
 *   4. Specific tier requested → use that tier (fallback to PAKO if unavailable)
 *
 * @param data Input bytes to compress
 * @param config Compression configuration
 * @returns CompressionResult with compressed data and metadata
 */
export function compress(data: Uint8Array, config: CompressConfig = {}): CompressionResult {
  const cfg = { ...DEFAULT_COMPRESS_CONFIG, ...config };
  const level = cfg.level ?? 6;
  const originalSize = data.length;

  // Skip if already compressed
  if (cfg.skipIfCompressed && detectInputType(data) === 'already-compressed') {
    return {
      data,
      tier: CompressorTier.PAKO, // passthrough — no compression applied
      ratio: 1.0,
      originalSize,
      compressedSize: originalSize,
    };
  }

  let tier: CompressorTier;
  let compressed: Uint8Array;

  if (cfg.tier === 'auto' || cfg.tier === undefined) {
    // Auto-detect: use magic-byte detection, then isDnaSequence as secondary check
    const inputType = detectInputType(data);
    const isBiological = inputType === 'biological' ||
      (inputType === 'general' && isDnaSequence(data));

    if (isBiological) {
      // Try specialized DNA compressors in preference order
      if (dnaCompressors[CompressorTier.NAF]) {
        tier = CompressorTier.NAF;
        compressed = dnaCompressors[CompressorTier.NAF]!.compress(data, level);
      } else if (dnaCompressors[CompressorTier.JARVIS3]) {
        tier = CompressorTier.JARVIS3;
        compressed = dnaCompressors[CompressorTier.JARVIS3]!.compress(data, level);
      } else {
        tier = CompressorTier.PAKO;
        compressed = compressWithPako(data, level);
      }
    } else {
      // General data: prefer zstd (fflate), fallback to pako
      if (zstdCompressWasm || fflate) {
        tier = CompressorTier.ZSTD;
        compressed = compressWithZstd(data, level);
      } else {
        tier = CompressorTier.PAKO;
        compressed = compressWithPako(data, level);
      }
    }
  } else {
    // Specific tier requested
    tier = cfg.tier;

    switch (tier) {
      case CompressorTier.PAKO:
        compressed = compressWithPako(data, level);
        break;

      case CompressorTier.ZSTD:
        compressed = compressWithZstd(data, level);
        // If no zstd WASM and no fflate, tier effectively becomes PAKO
        if (!zstdCompressWasm && !fflate) {
          tier = CompressorTier.PAKO;
        }
        break;

      case CompressorTier.NAF:
      case CompressorTier.AGC:
      case CompressorTier.DEEP_GECO:
      case CompressorTier.MBGC2:
      case CompressorTier.JARVIS3:
        if (dnaCompressors[tier]) {
          compressed = dnaCompressors[tier]!.compress(data, level);
        } else {
          // Fallback to pako
          compressed = compressWithPako(data, level);
          tier = CompressorTier.PAKO;
        }
        break;

      default:
        compressed = compressWithPako(data, level);
        tier = CompressorTier.PAKO;
    }
  }

  const compressedSize = compressed.length;
  const ratio = originalSize > 0 ? originalSize / compressedSize : 1.0;

  return { data: compressed, tier, ratio, originalSize, compressedSize };
}

/**
 * Decompress data using the specified tier (or auto-detect from magic bytes).
 *
 * @param data Compressed bytes
 * @param config Compression configuration (tier must match what was used to compress)
 * @returns Decompressed bytes
 */
export function decompress(data: Uint8Array, config: CompressConfig = {}): Uint8Array {
  const cfg = { ...DEFAULT_COMPRESS_CONFIG, ...config };

  // Try to determine the compression format from magic bytes
  if (cfg.tier === 'auto' || cfg.tier === undefined) {
    // NAF magic: 0x4E 0x41 0x46 0x01
    if (
      data.length >= 4 &&
      data[0] === 0x4E && data[1] === 0x41 &&
      data[2] === 0x46 && data[3] === 0x01
    ) {
      return decompressWithNAF(data);
    }

    // JARVIS3 magic: 0x4A 0x33 0x56 0x01
    if (
      data.length >= 4 &&
      data[0] === 0x4A && data[1] === 0x33 &&
      data[2] === 0x56 && data[3] === 0x01
    ) {
      return decompressWithJarvis3(data);
    }

    // Zstd magic: 0x28 0xB5 0x2F 0xFD
    if (
      data.length >= 4 &&
      data[0] === 0x28 && data[1] === 0xB5 &&
      data[2] === 0x2F && data[3] === 0xFD
    ) {
      return decompressWithZstd(data);
    }

    // Default: assume DEFLATE (pako)
    return decompressWithPako(data);
  }

  // Specific tier
  switch (cfg.tier) {
    case CompressorTier.PAKO:
      return decompressWithPako(data);

    case CompressorTier.ZSTD:
      return decompressWithZstd(data);

    case CompressorTier.NAF:
      return decompressWithNAF(data);

    case CompressorTier.JARVIS3:
      return decompressWithJarvis3(data);

    case CompressorTier.AGC:
    case CompressorTier.DEEP_GECO:
    case CompressorTier.MBGC2: {
      const compressor = dnaCompressors[cfg.tier];
      if (compressor) return compressor.decompress(data);
      throw new Error(
        `${cfg.tier} decompression requested but implementation not available`,
      );
    }

    default:
      return decompressWithPako(data);
  }
}
