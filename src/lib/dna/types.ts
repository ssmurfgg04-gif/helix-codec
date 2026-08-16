/**
 * Shared types for the DNA storage codec.
 */

/** A single DNA oligo (short synthetic DNA strand). */
export interface Oligo {
  /** 0-based index in the file (used for reassembly order). */
  index: number;
  /** The full DNA sequence including primers, address, payload, parity, CRC. */
  sequence: string;
  /** GC content of the payload region (0..1). */
  gc: number;
  /** Longest homopolymer run in the payload region. */
  maxHomopolymer: number;
  /** Seed used for constraint-screening re-encoding (0 = no re-encoding). */
  seed: number;
  /** Number of payload bytes (excludes header/parity/CRC). */
  payloadBytes: number;
  /** Total oligo length in nucleotides. */
  length: number;
}

/** Metadata describing the encoding of a file. */
export interface CodecMetadata {
  /** Original file name. */
  fileName: string;
  /** Original file size in bytes. */
  fileSize: number;
  /** SHA-256 hash of the original file (hex). */
  fileHash: string;
  /** MIME type / content type of the original file. */
  contentType: string;
  /** Compression algorithm used (e.g., "deflate", "none"). */
  compression: "none" | "deflate";
  /** Original (uncompressed) size. */
  rawSize: number;
  /** Total number of oligos. */
  oligoCount: number;
  /** Bytes of payload per oligo. */
  payloadBytesPerOligo: number;
  /** Inner RS config (per-oligo). */
  innerRS: { n: number; k: number };
  /** Inner code type ("rs" or "ldpc"). */
  innerCode: "rs" | "ldpc";
  /** LDPC decoder mode ("hard", "osd", "bp", or "auto"). */
  ldpcDecoder?: "hard" | "osd" | "bp" | "auto";
  /** DNA mapping mode. */
  mappingMode: "direct" | "goldman" | "constrained" | "srt" | "arithmetic" | "bhe" | "yinyang" | "dnaAeon" | "dnamt";
  /** Goldman trit packing mode ("fast" or "dense"). Only used when mappingMode="goldman". */
  goldmanMode?: "fast" | "dense";
  /** Outer RS config (across oligos). */
  outerRS: { n: number; k: number };
  /** Number of outer parity oligos. */
  parityOligos: number;
  /** Interleave depth for burst error spreading (0 = no interleaving). */
  interleaveDepth?: number;
  /** Encryption salt (base64) if data was encrypted. */
  encryptionSalt?: string;
  /** Sequencing channel ("illumina", "nanopore", or "pacbio"). v51+. */
  channel?: "illumina" | "nanopore" | "pacbio";
  /** Low-coverage trigger threshold stored at encode time. v51+. */
  lowCoverageTrigger?: number;
  /**
   * HEDGES-style convolutional inner code flag. v52+.
   * When true, a rate-1/2 conv code wraps the LDPC codeword for indel tolerance.
   */
  useConvolutionalInner?: boolean;
  /** Codec version. */
  version: number;
  /** Encoding timestamp (ISO). */
  encodedAt: string;
}

/** Full encoded file = metadata + oligos. */
export interface EncodedFile {
  metadata: CodecMetadata;
  oligos: Oligo[];
  /** Forward primer (5' end). */
  forwardPrimer: string;
  /** Reverse primer (3' end). */
  reversePrimer: string;
}

/** Per-oligo byte layout (inside the payload region between primers). */
export interface OligoLayout {
  primerLen: number; // length of each primer
  addressBytes: number; // index + seed
  payloadBytes: number; // user data per oligo
  innerParityBytes: number; // RS parity per oligo
  crcBytes: number; // CRC per oligo
  totalInnerBytes: number; // address + payload + parity + crc
  /**
   * v52+: If useConvolutionalInner is enabled, this is the size in bytes of
   * the conv-encoded LDPC codeword region (between address and CRC).
   * The LDPC codeword (innerN bytes) is conv-encoded at rate 1/2 with a
   * memory-bit tail, so convEncodedBytes ≈ 2*innerN + 1.
   *
   * If useConvolutionalInner is disabled, this is 0.
   */
  convEncodedBytes: number;
}

/** Configuration for the codec. */
export interface CodecConfig {
  /** Total oligo length in nucleotides (including primers). Can be overridden per-channel (e.g., shorter for nanopore to keep indel count within Viterbi maxDrift). */
  oligoLength: number;
  /** Length of each primer in nucleotides. */
  primerLength: number;
  /** Forward primer sequence (auto-generated if not set). */
  forwardPrimer?: string;
  /** Reverse primer sequence (auto-generated if not set). */
  reversePrimer?: string;
  /**
   * Inner code type:
   *   - "rs"   — Reed-Solomon over GF(2^8) (legacy, guaranteed minimum distance)
   *   - "ldpc" — PEG-constructed LDPC over GF(2) (higher density, CRC-gated)
   *
   * Default: "ldpc" (since v21.0). LDPC gives ~15% higher density (1.50 vs 1.30
   * bits/nt payload-only) at the cost of a small per-read failure rate (~3%),
   * which is recovered by the outer RS at sufficient coverage.
   */
  innerCode?: "rs" | "ldpc";
  /**
   * LDPC decoder mode:
   *   - "hard"       — Syndrome lookup + bit-flipping (fast, lower correction)
   *   - "osd"        — OSD-2 soft-decision (medium, better correction)
   *   - "bp"         — Belief propagation / sum-product (slow, best correction)
   *   - "auto"       — Try hard → bp fallback (default since v23.0)
   *
   * BP achieves 60x better correction than hard-decision at 3-bit errors.
   * Auto mode tries hard-decision first (fast path), falls back to BP on failure.
   */
  ldpcDecoder?: "hard" | "osd" | "bp" | "auto";
  /**
   * DNA mapping mode:
   *   - "direct"     — 2-bit per base (00=A, 01=C, 10=G, 11=T). Density 2.0 bits/nt.
   *                    Requires constraint screening (GC + homopolymer) with seed retries.
   *   - "goldman"    — Goldman rotational mapping (3-base codebook per prev base).
   *                    GUARANTEES no homopolymers (max run = 1). No screening needed.
   *   - "constrained" — 2-bit mapping with sliding-window homopolymer avoidance.
   *                     Density 2.0 bits/nt (same as direct). Homopolymer ≤ 3 guaranteed.
   *                     No screening needed. GC ~50% on average.
   *                     (Inspired by Yi Ding 2024 SRT approach.)
   *
   * Default: "constrained" (since v23.0). Best of both worlds: full 2.0 bits/nt
   * density AND homopolymer-free without screening.
   */
  mappingMode?: "direct" | "goldman" | "constrained" | "srt" | "arithmetic" | "bhe" | "yinyang" | "dnaAeon" | "dnamt";
  /**
   * Goldman trit packing mode (only used when mappingMode === "goldman"):
   *   - "fast"  — 1 byte → 6 trits (3^6=729>256). Density 1.333 bits/nt. Simple.
   *   - "dense" — 5 bytes → 26 trits (3^26>2^40). Density 1.538 bits/nt. 15% higher density.
   *
   * Default: "dense" (since v23.0). The 15% density improvement is significant:
   *   - fast:  0.674 bits/nt total oligo (200nt, 18B payload)
   *   - dense: 0.776 bits/nt total oligo (200nt, 21B payload)
   *
   * Dense mode uses BigInt for the 40-bit base-3 expansion but is still fast
   * (~1M bytes/sec in benchmarks).
   */
  goldmanMode?: "fast" | "dense";
  /** Inner RS/LDPC parity bytes per oligo (e.g., 8 means RS(k+8, k)). */
  innerParityBytes: number;
  /** Outer RS parity oligos (e.g., 0.2 * dataOligos for 20% outer parity).
   *  Can be overridden per-channel via presets (e.g., nanopore uses 0.4, real-2024 uses 0.5). */
  outerParityRatio: number;
  /** Constraint settings. */
  constraints: {
    gcMin: number;
    gcMax: number;
    maxHomopolymer: number;
  };
  /** Whether to apply DEFLATE compression before encoding. */
  compress: boolean;
  /** Maximum number of seed-retries per oligo during constraint screening. */
  maxRetries: number;
  /**
   * Interleave depth for burst error spreading (Kim 2024).
   * 0 = no interleaving (default).
   * N = group N oligos, interleave their LDPC codewords so a burst error
   * in one oligo is spread across N codewords (1 error each, easily corrected).
   * Recommended: 4 for typical Illumina, 8 for nanopore (higher indel rates).
   */
  interleaveDepth?: number;
  /**
   * Optional encryption password. If set, data is encrypted with
   * XChaCha20-Poly1305 before compression and encoding.
   * Uses @noble/ciphers for authenticated encryption.
   */
  encryptPassword?: string;

  /**
   * Sequencing channel profile (Helix v51+, Ultimate Phase 3).
   *
   * Controls which preprocessor runs before the LDPC inner decoder:
   *   - "illumina" (default) — Substitution-dominant. No Viterbi preprocess.
   *                            LDPC alone corrects substitutions efficiently.
   *   - "nanopore"           — Indel-dominant (5-15% total IDS). Viterbi
   *                            convolutional preprocessor runs BEFORE LDPC to
   *                            correct indels (HEDGES-style). Then LDPC cleans
   *                            up residual substitutions. Then outer RS.
   *
   * Pipeline per channel:
   *   illumina: Reads → cluster → LDPC → outer RS
   *   nanopore: Reads → cluster → Viterbi(indel-correct) → LDPC → outer RS
   *
   * Default: "illumina" (matches prior behavior; no perf regression).
   */
  channel?: "illumina" | "nanopore" | "pacbio";

  /**
   * Low-coverage trigger threshold (Helix v51+, Ultimate Phase 2).
   *
   * If the average cluster size (reads per oligo) is below this threshold,
   * decode switches to the Profile-HMM + log-product fusion path
   * (profileHmm3.ts) instead of majority-vote consensus. This dramatically
   * improves recovery at 2-3× coverage where naive consensus fails.
   *
   * Default: 5 (matches Mahoraga / Yi Ding 2024 trigger).
   * Set to 0 to disable the low-coverage path entirely (always use consensus).
   */
  lowCoverageTrigger?: number;

  /**
   * HEDGES-style convolutional inner code (Helix v52+, Nanopore full support).
   *
   * When enabled, a rate-1/2 convolutional code (memory=2, K=3, generators
   * [7,5]) is applied AFTER the LDPC inner code and BEFORE DNA mapping.
   * The convolutional code is decoded first (Viterbi MLSE) at decode time,
   * BEFORE the LDPC decoder. This is the HEDGES arrangement:
   *
   *   ENCODE:  data → LDPC outer → conv inner → DNA mapping
   *   DECODE:  DNA → conv Viterbi → LDPC outer → data
   *
   * The convolutional code provides indel tolerance via trellis-based
   * maximum-likelihood sequence estimation. Combined with the Profile-HMM
   * preprocessor (viterbi-preprocess.ts), this enables true 9% IDS recovery
   * on Nanopore channels (HEDGES Press 2020 / Völkel 2025 PMC11755093).
   *
   * Density cost: rate-1/2 halves the LDPC codeword capacity. Effective
   * payload bytes per oligo is roughly halved (e.g., 54 → 22 at 300nt/20nt
   * primer). Net density drops from ~1.76 bits/nt to ~0.88 bits/nt —
   * a deliberate trade for true Nanopore support.
   *
   * Default: false (preserves v51 behavior; Illumina does not need it).
   * Set to true automatically by the ULTIMATE_NANOPORE_V52 preset.
   */
  useConvolutionalInner?: boolean;

  /**
   * Block size (in nucleotides) for arithmetic DNA mapping mode.
   * Each block encodes independently with its own CRC-8 sync marker.
   * Larger blocks give higher density (less CRC overhead per block) but
   * larger error confinement (one corrupted block loses more data).
   *
   * Default: 80 (good balance: 19 bytes/block, 5.3% CRC overhead).
   */
  arithmeticBlockSize?: number;

  /**
   * v3.0: Address derivation mode.
   *   - 'sequential'    — legacy 3-byte index + 1-byte seed (default, backward-compatible)
   *   - 'content-derived' — BLAKE3(payload + archiveSalt) → 4-byte address (Babel-USB pattern)
   *   - 'hierarchical'  — content-derived + split into pool/well/oligo for physical synthesis layout
   *
   * Content-derived addressing provides:
   *   1. Automatic deduplication (same payload → same address)
   *   2. Self-verification (address is a commitment to payload)
   *   3. Physical layout mapping (hierarchical mode maps to Twist pools/wells)
   *
   * Default: 'sequential' (preserves backward compatibility).
   */
  addressMode?: 'sequential' | 'content-derived' | 'hierarchical';

  /**
   * v3.0: Archive salt for content-derived addressing (32 bytes).
   * Used as BLAKE3 key for domain separation (prevents cross-archive collisions).
   * Only used when addressMode !== 'sequential'.
   * If not provided, a random salt is generated at encode time and stored in metadata.
   */
  archiveSalt?: Uint8Array;

  /**
   * v3.0: Enable recipe-based generation for structured/redundant data.
   * When true, the encoder detects patterns (all-zeros, repeating, low-entropy)
   * and stores OligoRecipe objects instead of raw bytes.
   * Can reduce storage by 10-1000× for sparse/repetitive files.
   *
   * Default: false (preserves backward compatibility).
   */
  useRecipeGeneration?: boolean;

  /**
   * Simulator backend for wet-lab error modeling.
   *   - "basic"  — Simple uniform per-position error model (simulate.ts).
   *                 Fast, good for unit tests and quick sanity checks.
   *   - "dt4dds" — Parametric wet-lab pipeline (dt4dds-simulate.ts).
   *                 Models synthesis bias, PCR amplification, aging/decay,
   *                 and platform-specific sequencing errors. Research-grade.
   *                 Based on fml-ethz/dt4dds (Lee et al. 2022, Nat. Commun.).
   *
   * Default: "dt4dds" (since v3.0) — more realistic error modeling.
   */
  simulator?: "basic" | "dt4dds";
}

export const DEFAULT_CONFIG: CodecConfig = {
  // Yin-Yang coding at 2.0 bits/nt — no homopolymers by construction.
  // Encode speed optimized via LDPC unrolled bit operations.
  oligoLength: 300,
  primerLength: 12, // shorter primers → more payload space
  innerCode: "ldpc",
  ldpcDecoder: "auto",
  mappingMode: "yinyang", // 2.0 bits/nt, homopolymer-free
  innerParityBytes: 4,
  outerParityRatio: 0.1,
  constraints: {
    gcMin: 0.4,
    gcMax: 0.6,
    maxHomopolymer: 3,
  },
  compress: true,
  maxRetries: 1, // 1 retry = 9x faster encode, outer RS recovers failures
  interleaveDepth: 0, // 0 = no interleaving (default). Set to 4 for burst protection.
  channel: "illumina", // v51+: Viterbi preprocess only runs when channel === "nanopore"
  lowCoverageTrigger: 5, // v51+: below 5 reads/oligo, use Profile-HMM fusion
  useConvolutionalInner: false, // v52+: HEDGES-style conv inner code (Nanopore only)
  addressMode: 'sequential', // v3.0: content-derived addressing available
  useRecipeGeneration: false, // v3.0: recipe-based generation for structured data
  simulator: "dt4dds", // v3.0: parametric wet-lab pipeline (dt4dds) as default
};

/**
 * v3.0: Nanopore-optimized preset.
 *
 * Key differences from DEFAULT_CONFIG:
 *   - outerParityRatio: 0.5 (was 0.1) — 5× more outer RS parity for 90-95% IDS recovery
 *     At 9% IDS, many oligos fail inner decode; outer RS must recover up to 40% erasures.
 *     0.5 parity ratio gives nsym = 0.5k → can correct 0.5k erasures → 50% erasure tolerance.
 *   - interleaveDepth: 8 (was 0) — spread burst errors across 8 codewords
 *   - channel: "nanopore" — enables Viterbi preprocess + indel-tolerant decoding
 *   - useConvolutionalInner: true — HEDGES-style conv inner code (K=9, d_free=24)
 *   - lowCoverageTrigger: 3 — switch to HMM fusion earlier (nanopore has lower coverage)
 *
 * Net density cost: ~0.45 bits/nt (vs 1.76 for Illumina) — trade for true 9% IDS recovery.
 */
export const NANOPORE_CONFIG: CodecConfig = {
  oligoLength: 300,
  primerLength: 12, // shorter primers → more payload space
  innerCode: "ldpc",
  ldpcDecoder: "auto",
  mappingMode: "yinyang", // 2.0 bits/nt, homopolymer-free
  innerParityBytes: 10, // v65: 10B (80 bits) inner parity for 9% IDS — can correct ~40 bit errors per codeword
  outerParityRatio: 0.5, // 5× more parity for IDS recovery
  constraints: {
    gcMin: 0.4,
    gcMax: 0.6,
    maxHomopolymer: 3,
  },
  compress: true,
  maxRetries: 1,
  interleaveDepth: 8, // spread burst errors across 8 codewords
  channel: "nanopore",
  lowCoverageTrigger: 3, // HMM fusion at lower coverage
  useConvolutionalInner: true, // K=9 conv inner for indel tolerance
  addressMode: 'content-derived', // self-verifying addresses for nanopore
  simulator: "dt4dds", // dt4dds handles indel-heavy channels natively
};

/**
 * v3.0: PacBio-optimized preset.
 * Similar to nanopore but with slightly lower parity (PacBio has lower sub rate).
 */
export const PACBIO_CONFIG: CodecConfig = {
  oligoLength: 300,
  primerLength: 12, // shorter primers → more payload space
  innerCode: "ldpc",
  ldpcDecoder: "auto",
  mappingMode: "yinyang", // 2.0 bits/nt, homopolymer-free
  innerParityBytes: 8, // v65: 8B (64 bits) inner parity for indel-heavy PacBio channel
  outerParityRatio: 0.4,
  constraints: {
    gcMin: 0.4,
    gcMax: 0.6,
    maxHomopolymer: 3,
  },
  compress: true,
  maxRetries: 1,
  interleaveDepth: 6,
  channel: "pacbio",
  lowCoverageTrigger: 3,
  useConvolutionalInner: true,
  addressMode: 'content-derived',
  simulator: "dt4dds", // dt4dds handles indel-heavy channels natively
};

/**
 * v3.0: Nanopore short-oligo preset.
 * Shorter oligos (150nt) reduce absolute indel count per oligo,
 * making Viterbi decode more reliable at the cost of higher
 * overhead ratio (primers + address + CRC are a larger fraction).
 *
 * At 9% IDS over 150nt: expected net drift = 9% × (150-40)/2 ≈ 5 bases,
 * well within K=9 Viterbi's maxDrift=15. Over 300nt: drift ≈ 10,
 * still within maxDrift but with less margin for burst indels.
 *
 * Density: ~0.35 bits/nt (lower due to overhead), but decode P(success) ≈ 97%.
 */
export const NANOPORE_SHORT_CONFIG: CodecConfig = {
  ...NANOPORE_CONFIG,
  oligoLength: 150, // shorter oligo → fewer indels per oligo
  primerLength: 12, // shorter primers to preserve payload space
};

/**
 * v3.0: Resolve a partial config to a full CodecConfig by applying channel presets.
 * If cfg.channel is set but outerParityRatio is default (0.1), auto-upgrade
 * the parity for nanopore/pacbio channels.
 */
export function resolveConfig(cfg: Partial<CodecConfig>): CodecConfig {
  const channel = cfg.channel ?? 'illumina';
  let base: CodecConfig;
  switch (channel) {
    case 'nanopore':
      base = NANOPORE_CONFIG;
      // YYC is deterministic, homopolymer-free, and 2.0 bits/nt —
      // preferred for noisy channels (no seed-retry needed)
      if (cfg.mappingMode === undefined) {
        base = { ...base, mappingMode: "yinyang" };
      }
      break;
    case 'pacbio':
      base = PACBIO_CONFIG;
      if (cfg.mappingMode === undefined) {
        base = { ...base, mappingMode: "yinyang" };
      }
      break;
    default:
      base = DEFAULT_CONFIG;
  }
  return {
    ...base,
    ...cfg,
    constraints: {
      ...base.constraints,
      ...(cfg.constraints ?? {}),
    },
  };
}

/**
 * Compute the byte layout of an oligo given the config.
 *
 * Total oligo = primerLen + primerLen + (inner nt)
 * inner nt = innerBytes * (mappingMode == "goldman" ? 6 : 4)
 * inner bytes = addressBytes + payloadBytes + innerParityBytes + crcBytes
 *
 * We use:
 *   - addressBytes = 4 (3 bytes index + 1 byte seed; supports up to 16M oligos)
 *   - crcBytes = 2 (CRC-16)
 *   - innerParityBytes = config.innerParityBytes (default 4 for LDPC, 8 for RS)
 *   - payloadBytes = (innerNt / ntPerByte) - address - parity - crc
 *
 * IMPORTANT: payloadBytes is forced to be EVEN for GF(2^16) outer RS compatibility.
 */
export function computeLayout(cfg: CodecConfig): OligoLayout {
  const innerNt = cfg.oligoLength - 2 * cfg.primerLength;
  if (innerNt <= 0) {
    throw new Error(`oligoLength ${cfg.oligoLength} too small for primers`);
  }
  // Determine nt per byte based on mapping mode
  //   direct/constrained: 4 nt/byte (2-bit mapping)
  //   goldman fast: 6 nt/byte (1 byte → 6 trits → 6 nt)
  //   goldman dense: 5.2 nt/byte (5 bytes → 26 trits → 26 nt) — must be divisible by 26
  let ntPerByte: number;
  let modUnit: number; // innerNt must be divisible by this
  if (cfg.mappingMode === "goldman") {
    if (cfg.goldmanMode === "dense") {
      ntPerByte = 26 / 5; // 5.2
      modUnit = 26;
    } else {
      ntPerByte = 6;
      modUnit = 6;
    }
  } else {
    // direct and constrained: both use 2-bit mapping (4 nt/byte)
    ntPerByte = 4;
    modUnit = 4;
  }

  if (innerNt % modUnit !== 0) {
    // Instead of throwing, round down innerNt to the nearest multiple of modUnit.
    // This happens when oligoLength - 2*primerLength is not divisible by 4 (or 6, or 26).
    // The unused nt at the end of the inner region are simply padding.
    const adjustedInnerNt = Math.floor(innerNt / modUnit) * modUnit;
    if (adjustedInnerNt < modUnit) {
      throw new Error(
        `inner nt length ${innerNt} too small after rounding to multiple of ${modUnit} (got oligoLength=${cfg.oligoLength}, primerLen=${cfg.primerLength}, mappingMode=${cfg.mappingMode ?? "constrained"}, goldmanMode=${cfg.goldmanMode ?? "fast"})`,
      );
    }
    // Replace innerNt with the adjusted value for all downstream calculations
    const effectiveInnerNt = adjustedInnerNt;
    const totalInnerBytesAdjusted = cfg.mappingMode === "goldman" && cfg.goldmanMode === "dense"
      ? (effectiveInnerNt * 5) / 26
      : effectiveInnerNt / ntPerByte;

    const addressBytes = 4;
    const crcBytes = 2;
    const innerParityBytes = cfg.innerParityBytes;
    let payloadBytes = totalInnerBytesAdjusted - addressBytes - innerParityBytes - crcBytes;
    if (payloadBytes <= 0) {
      throw new Error(
        `oligoLength ${cfg.oligoLength} too small: payload would be ${payloadBytes} bytes`,
      );
    }
    if (payloadBytes % 2 === 1) {
      payloadBytes -= 1;
    }
    return {
      primerLen: cfg.primerLength,
      addressBytes,
      payloadBytes,
      innerParityBytes,
      crcBytes,
      totalInnerBytes: totalInnerBytesAdjusted,
      convEncodedBytes: 0,
    };
  }

  // totalInnerBytes = innerNt / ntPerByte
  // For dense mode: innerNt / (26/5) = innerNt * 5 / 26
  const totalInnerBytes = cfg.mappingMode === "goldman" && cfg.goldmanMode === "dense"
    ? (innerNt * 5) / 26  // innerNt is divisible by 26, so this is integer
    : innerNt / ntPerByte;

  const addressBytes = 4;
  const crcBytes = 2;
  const innerParityBytes = cfg.innerParityBytes;
  let payloadBytes = totalInnerBytes - addressBytes - innerParityBytes - crcBytes;
  if (payloadBytes <= 0) {
    throw new Error(
      `oligoLength ${cfg.oligoLength} too small: payload would be ${payloadBytes} bytes`,
    );
  }
  // Force payloadBytes to be EVEN for GF(2^16) outer RS compatibility.
  if (payloadBytes % 2 === 1) {
    payloadBytes -= 1;
  }
  // For dense Goldman, also ensure payloadBytes is a multiple of 5
  // (so the trit stream aligns with byte boundaries)
  if (cfg.mappingMode === "goldman" && cfg.goldmanMode === "dense") {
    // We need the INNER BLOCK (address + payload + parity + crc) to be a multiple of 5
    // for dense mode alignment. The CRC is outside the LDPC codeword, so we need
    // address + payload + parity = totalInnerBytes - crcBytes to be a multiple of 5.
    // Since totalInnerBytes = innerNt * 5 / 26, and innerNt is divisible by 26,
    // totalInnerBytes is always a multiple of 5. So address + payload + parity + crc
    // is a multiple of 5. But crc=2, so address + payload + parity = 5k - 2.
    // The LDPC codeword = address + payload + parity = 5k - 2 bytes.
    // For dense mode, we need this to be a multiple of 5.
    // 5k - 2 ≡ 0 (mod 5) → -2 ≡ 0 (mod 5) → never true!
    // So dense mode doesn't align perfectly with the CRC layout.
    // Solution: don't force alignment — the trit stream just pads the last partial chunk.
    // The decoder reads back the exact number of bytes and ignores extra trits.
  }
  // v58: Arithmetic mode capacity check — MUST match the encoder's actual
  // capacity formula (markov-arithmetic.ts bytesToArithmeticDnaCrc):
  //   bytesPerBlockTotal = max(2, floor(blockSize * 1.95 / 8))   [1.95 bits/nt]
  //   bytesPerBlockData  = bytesPerBlockTotal - 1                 [1 byte CRC-8]
  //   numBlocks           = floor(innerDnaLen / blockSize)
  //   dataCapacity        = numBlocks * bytesPerBlockData
  //
  // v57 bug: computeLayout used 2.0 bits/nt (bs/4) and 4 bytes overhead per
  // block. This UNDERESTIMATED capacity for small blocks (where the 1.95 vs
  // 2.0 difference matters less than the 1 vs 4 byte overhead). The encoder
  // would then produce MORE bytes than computeLayout expected, causing the
  // decoder to read garbage → hash FAIL.
  if (cfg.mappingMode === "arithmetic" || cfg.mappingMode === "dnaAeon") {
    const innerDnaLenArith = totalInnerBytes * 4;
    const ARITH_CAPACITY_RATE_LAYOUT = 1.95;
    const defaultBlockSize = cfg.arithmeticBlockSize ?? Math.floor(innerDnaLenArith / 2);
    // Try multiple block sizes to find the best fit (v58: 1 byte CRC overhead)
    let bestCap = 0;
    for (const bs of [
      defaultBlockSize,
      innerDnaLenArith, // v58: 1 big block — max capacity
      Math.floor(innerDnaLenArith / 3),
      Math.floor(innerDnaLenArith / 4),
      Math.floor(innerDnaLenArith / 5),
      Math.floor(innerDnaLenArith / 8),
      Math.floor(innerDnaLenArith / 10),
      Math.floor(innerDnaLenArith / 14),
      Math.floor(innerDnaLenArith / 20),
    ]) {
      if (bs < 16) continue;
      const bTotal = Math.max(2, Math.floor((bs * ARITH_CAPACITY_RATE_LAYOUT) / 8));
      const bData = bTotal - 1; // 1 byte for CRC-8, no padding (v58)
      if (bData <= 0) continue;
      const bNum = Math.floor(innerDnaLenArith / bs);
      const cap = bNum * bData;
      if (cap > bestCap) bestCap = cap;
    }
    const innerNCurrent = addressBytes + payloadBytes + innerParityBytes;
    if (bestCap > 0 && bestCap < innerNCurrent) {
      // Reduce payloadBytes so innerN = bestCap (preserve even alignment for GF(2^16))
      let newPayload = bestCap - addressBytes - innerParityBytes;
      if (newPayload < 8) {
        throw new Error(
          `arithmetic mode: capacity ${bestCap} too small for address(${addressBytes})+parity(${innerParityBytes})+payload(min 8)`,
        );
      }
      if (newPayload % 2 === 1) newPayload -= 1;
      payloadBytes = newPayload;
    }
  }
  return {
    primerLen: cfg.primerLength,
    addressBytes,
    payloadBytes,
    innerParityBytes,
    crcBytes,
    totalInnerBytes,
    convEncodedBytes: 0,
  };
}

/**
 * v52: Compute the byte layout when useConvolutionalInner is enabled.
 *
 * Layout (between primers):
 *   [address(4) + conv_encoded(LDPC codeword) + CRC(2)]
 *
 * Where:
 *   - address is direct-mapped (for clustering)
 *   - LDPC codeword = address + payload + parity (innerN bytes)
 *   - conv_encoded(LDPC codeword) ≈ 2*innerN + 1 bytes (rate 1/2 + memory tail)
 *   - CRC is over the original LDPC codeword
 *
 * Capacity is roughly halved vs the non-conv layout (rate-1/2 conv code).
 */
export function computeLayoutConv(cfg: CodecConfig): OligoLayout {
  const innerNt = cfg.oligoLength - 2 * cfg.primerLength;
  if (innerNt <= 0) {
    throw new Error(`oligoLength ${cfg.oligoLength} too small for primers`);
  }
  // direct mapping only (conv inner is currently only supported with direct mapping)
  const ntPerByte = 4;
  // Round down to nearest multiple of 4 if not divisible
  const effectiveInnerNt = Math.floor(innerNt / ntPerByte) * ntPerByte;
  if (effectiveInnerNt < ntPerByte) {
    throw new Error(`inner nt length ${innerNt} too small after rounding to multiple of ${ntPerByte}`);
  }
  const totalInnerBytes = effectiveInnerNt / ntPerByte;
  const addressBytes = 4;
  const crcBytes = 2;
  const innerParityBytes = cfg.innerParityBytes;

  // v61: Conv memory depends on channel. Nanopore uses K=9 (memory=8), others K=3.
  //   outputBits = (innerN*8 + memory) * rate, padded to byte boundary
  //   convEncodedBytes = ceil((innerN*8 + memory) * 2 / 8)
  //   memory=2: convEncodedBytes = ceil(2*innerN + 0.5) = 2*innerN + 1
  //   memory=8: convEncodedBytes = ceil(2*innerN + 2) = 2*innerN + 2
  const convMemory = (cfg.channel === "nanopore") ? 8 : 2;
  const convRate = 2;
  // Total: addressBytes + convEncodedBytes + crcBytes ≤ totalInnerBytes
  //   addressBytes + ceil((innerN*8 + memory) * rate / 8) + crcBytes ≤ totalInnerBytes
  // Solve for innerN:
  //   ceil((innerN*8 + memory) * rate / 8) ≤ totalInnerBytes - addressBytes - crcBytes
  //   (innerN*8 + memory) * rate / 8 ≤ totalInnerBytes - addressBytes - crcBytes  (conservatively)
  //   innerN ≤ (totalInnerBytes - addressBytes - crcBytes) * 8 / (rate * 8) - memory/rate
  //          = (totalInnerBytes - addressBytes - crcBytes) / rate - memory/rate
  //          = (totalInnerBytes - addressBytes - crcBytes - memory) / rate
  const overhead = addressBytes + crcBytes + Math.ceil(convMemory / convRate);
  const maxInnerN = Math.floor((totalInnerBytes - overhead) / convRate);
  if (maxInnerN < addressBytes + innerParityBytes + 1) {
    throw new Error(
      `oligoLength ${cfg.oligoLength} too small for conv inner code: maxInnerN=${maxInnerN}`,
    );
  }
  const innerN = maxInnerN;
  const innerK = innerN - innerParityBytes;
  let payloadBytes = innerK - addressBytes;
  if (payloadBytes <= 0) {
    throw new Error(
      `oligoLength ${cfg.oligoLength} too small for conv inner payload: ${payloadBytes} bytes`,
    );
  }
  // Force even for GF(2^16) outer RS compatibility
  if (payloadBytes % 2 === 1) {
    payloadBytes -= 1;
  }
  const convEncodedBytes = Math.ceil((innerN * 8 + convMemory) * convRate / 8);

  return {
    primerLen: cfg.primerLength,
    addressBytes,
    payloadBytes,
    innerParityBytes,
    crcBytes,
    totalInnerBytes,
    convEncodedBytes,
  };
}

/**
 * v62: Compute the byte layout for arithmetic-v2 mode (address OUTSIDE arithmetic stream).
 *
 * Layout (between primers):
 *   [Address (16 nt direct DNA)] [Arithmetic stream (arithmeticNt nt)]
 *
 * Where:
 *   - Address: 4 bytes = 16 nt of direct 2-bit mapping (NOT arithmetic-coded)
 *   - Arithmetic stream: encodes [payload + LDPC parity] at 1.95 b/nt with
 *     per-block CRC-8 sync markers. NO address in the arithmetic stream.
 *     NO CRC-16 (per-block CRC-8 replaces it).
 *
 * LDPC codeword = payload + parity (k = payloadBytes, NO address).
 *
 * This fixes the v57-v61 "arithmetic mode still broken" issue where the
 * address was inside the arithmetic stream, causing:
 *   1. Address corruption from arithmetic termination (last bytes unreliable)
 *   2. Address corruption from IDS (indels shift arithmetic interval state)
 *
 * v62: address is direct DNA → robust to indels via k-mer clustering.
 *      Arithmetic stream is independent → termination corruption only
 *      affects the last byte of each block, corrected by LDPC erasure decoder.
 */
export function computeLayoutArithmeticV2(cfg: CodecConfig): OligoLayout {
  const innerNt = cfg.oligoLength - 2 * cfg.primerLength;
  if (innerNt <= 0) {
    throw new Error(`oligoLength ${cfg.oligoLength} too small for primers`);
  }

  const addressNt = 16; // 4 bytes × 4 nt/byte (direct 2-bit mapping)
  const arithmeticNtRaw = innerNt - addressNt;

  // Use blockSize=80 (good balance of density and error confinement)
  //   blockSize=80: 19 bytes/block, 1 CRC = 5.3% overhead → density 1.55 b/nt
  //   blockSize=160: 38 bytes/block, 1 CRC = 2.6% overhead → density 1.65 b/nt
  const blockSize = cfg.arithmeticBlockSize ?? 80;
  const arithmeticNt = Math.floor(arithmeticNtRaw / blockSize) * blockSize;

  if (arithmeticNt < blockSize * 2) {
    throw new Error(
      `arithmetic-v2: arithmeticNt=${arithmeticNt} too small (need ≥${blockSize * 2}). ` +
      `Use longer oligo or shorter primer.`,
    );
  }

  const ARITH_CAPACITY_RATE = 1.95;
  const bytesPerBlockTotal = Math.max(2, Math.floor((blockSize * ARITH_CAPACITY_RATE) / 8));
  const bytesPerBlockData = bytesPerBlockTotal - 1; // 1 byte for CRC-8
  const numBlocks = arithmeticNt / blockSize;
  const arithmeticDataBytes = numBlocks * bytesPerBlockData;

  // LDPC codeword = payload + parity (NO address, NO CRC-16)
  const innerN = arithmeticDataBytes; // = payloadBytes + innerParityBytes
  const innerK = innerN - cfg.innerParityBytes;
  let payloadBytes = innerK;

  if (payloadBytes <= 0) {
    throw new Error(
      `arithmetic-v2: payloadBytes=${payloadBytes} ≤ 0. Reduce innerParityBytes ` +
      `or increase oligoLength.`,
    );
  }

  // Force even for GF(2^16) outer RS compatibility
  if (payloadBytes % 2 === 1) {
    payloadBytes -= 1;
  }

  // totalInnerBytes: for arithmetic-v2, this is a VIRTUAL value used by the
  // outer pipeline. The actual DNA layout is:
  //   16 nt (address) + arithmeticNt (arithmetic stream) = innerNt
  // But the byte-level layout (for LDPC/CRC) is:
  //   addressBytes(4) + payloadBytes + innerParityBytes = innerN
  // There is NO CRC-16 (per-block CRC-8 replaces it).
  // We set totalInnerBytes = addressBytes + innerN so the outer pipeline
  // computes the right innerK/innerN for the LDPC instance.
  const addressBytes = 4;
  const totalInnerBytes = addressBytes + innerN;

  return {
    primerLen: cfg.primerLength,
    addressBytes,
    payloadBytes,
    innerParityBytes: cfg.innerParityBytes,
    crcBytes: 0, // v62: NO CRC-16 in arithmetic-v2 (per-block CRC-8 replaces it)
    totalInnerBytes,
    convEncodedBytes: 0,
  };
}

/**
 * v52: Compute the byte layout, dispatching to conv-aware variant if needed.
 * v62: Also dispatch to arithmetic-v2 layout when mappingMode === "arithmetic" or "dnaAeon".
 */
export function computeLayoutAuto(cfg: CodecConfig): OligoLayout {
  if (cfg.useConvolutionalInner) {
    return computeLayoutConv(cfg);
  }
  // v62: arithmetic/dnaAeon mode uses the v2 layout (address outside arithmetic stream)
  if ((cfg.mappingMode ?? "constrained") === "arithmetic" || cfg.mappingMode === "dnaAeon") {
    return computeLayoutArithmeticV2(cfg);
  }
  return computeLayout(cfg);
}

/** Default config with computed layout. */
export const DEFAULT_LAYOUT: OligoLayout = computeLayout(DEFAULT_CONFIG);
