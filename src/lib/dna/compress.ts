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
 * │          │              │            │           │ very fast (WASM)     │
 * ├──────────┼──────────────┼────────────┼───────────┼──────────────────────┤
 * │ PAKO     │ DEFLATE      │ ~2.5×      │ ~1 MB     │ JS-native fallback,  │
 * │          │ (zlib)       │            │           │ always available     │
 * └──────────┴──────────────┴────────────┴───────────┴──────────────────────┘
 *
 * Default strategy:
 *   biological → NAF (if WASM loaded) else PAKO
 *   general    → ZSTD (if WASM loaded) else PAKO
 *   already-compressed → passthrough (no compression)
 *
 * WASM loading is deferred: zstd and NIF/AGC stubs check for a loaded module
 * at first use and fall back to pako if unavailable.
 *
 * Reference:
 *   - Varshney et al. (2024). "A universal nucleotide archive format." arXiv.
 *   - Deorowicz et al. (2015). "AGC: Assembly Graph Comparator."
 *   - Hofmann et al. (2022). "DeepGeCo: Deep DNA Sequence Compression."
 */

import * as pako from 'pako';

// ---------------------------------------------------------------------------
// Enums and interfaces
// ---------------------------------------------------------------------------

/** Compression tier identifiers (ordered by typical compression ratio). */
export enum CompressorTier {
  /** Nucleotide Archive Format — best for DNA sequences. WASM-only. */
  NAF = 'naf',
  /** Assembly Graph Compression — ref-based. WASM-only. */
  AGC = 'agc',
  /** DeepGeCo — neural DNA compression. WASM + GPU. */
  DEEP_GECO = 'deep_geco',
  /** Multi-context BGCompression. WASM-only. */
  MBGC2 = 'mbgc2',
  /** Jarvis3 — fast DNA compression. WASM-only. */
  JARVIS3 = 'jarvis3',
  /** Zstandard — general-purpose, very fast. WASM-only. */
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
// Zstd stubs — future WASM loading
// ---------------------------------------------------------------------------

let zstdCompressWasm: ((data: Uint8Array, level: number) => Uint8Array) | null = null;
let zstdDecompressWasm: ((data: Uint8Array) => Uint8Array) | null = null;

/**
 * Register WASM-backed zstd compress/decompress functions.
 *
 * Call this once after loading the zstd WASM module to enable
 * ZSTD tier compression. Without registration, zstd falls back to pako.
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

// ---------------------------------------------------------------------------
// NAF/AGC/DeepGeCo/MBGC2/JARVIS3 stubs — future WASM loading
// ---------------------------------------------------------------------------

/**
 * Map of specialized DNA compressor functions.
 * Populated by registerDnaCompressorWasm() after WASM module loading.
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
    // Auto-detect
    const inputType = detectInputType(data);

    if (inputType === 'biological') {
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
      // General data: prefer zstd, fallback to pako
      if (zstdCompressWasm) {
        tier = CompressorTier.ZSTD;
        compressed = zstdCompressWasm(data, level);
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
        if (zstdCompressWasm) {
          compressed = zstdCompressWasm(data, level);
        } else {
          // Fallback to pako with a warning (silent in production)
          compressed = compressWithPako(data, level);
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
    // Zstd magic: 0x28 0xB5 0x2F 0xFD
    if (
      data.length >= 4 &&
      data[0] === 0x28 && data[1] === 0xB5 &&
      data[2] === 0x2F && data[3] === 0xFD
    ) {
      if (zstdDecompressWasm) {
        return zstdDecompressWasm(data);
      }
      // Can't decompress zstd without WASM — fall through to pako error
      throw new Error('Zstd-compressed data but zstd WASM not loaded');
    }

    // Default: assume DEFLATE (pako)
    return decompressWithPako(data);
  }

  // Specific tier
  switch (cfg.tier) {
    case CompressorTier.PAKO:
      return decompressWithPako(data);

    case CompressorTier.ZSTD:
      if (zstdDecompressWasm) return zstdDecompressWasm(data);
      throw new Error('Zstd decompression requested but zstd WASM not loaded');

    case CompressorTier.NAF:
    case CompressorTier.AGC:
    case CompressorTier.DEEP_GECO:
    case CompressorTier.MBGC2:
    case CompressorTier.JARVIS3: {
      const compressor = dnaCompressors[cfg.tier];
      if (compressor) return compressor.decompress(data);
      throw new Error(
        `${cfg.tier} decompression requested but WASM module not loaded`,
      );
    }

    default:
      return decompressWithPako(data);
  }
}
