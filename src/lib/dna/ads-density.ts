/**
 * Adaptive Density System (ADS Codex pattern)
 *
 * Tunes encoding parameters to achieve higher density (up to 0.99 bits/nt).
 * The key insight is that the outer RS parity ratio and inner code parameters
 * can be optimized for the specific input size, channel, and synthesis platform.
 *
 * Current Helix: fixed ~0.84 bits/nt (DEFAULT_CONFIG with 10% outer parity)
 * ADS target: up to 0.99 bits/nt (optimized parity + payload tuning)
 *
 * How:
 *   1. For small files (<1KB): reduce outer parity (fewer oligos → RS overhead dominates)
 *   2. For large files (>100KB): increase interleaveDepth (burst protection at scale)
 *   3. For constrained mapping: switch from direct to constrained mode (no retries needed)
 *   4. Tune innerParityBytes based on channel error rates
 *   5. For nanopore: use convolutional inner + shorter oligos (more oligos, smaller errors each)
 *
 * Reference:
 *   - lanl/adscodex — achieves 0.99 bits/nt with lookup-table acceleration
 *   - ADS adapts to synthesis and sequencing requirements
 */

import type { CodecConfig } from "./types";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

export interface DensityConfig {
  /** Target density in bits/nt. Default: 0.95 */
  targetDensity: number;
  /** Minimum acceptable density. Default: 0.80 */
  minDensity: number;
  /** Synthesis platform: 'twist', 'idt', 'custom'. Default: 'twist' */
  synthesisPlatform: "twist" | "idt" | "custom";
  /** Maximum oligo length the synthesizer supports. Default: 300 */
  maxOligoLength: number;
}

export const DEFAULT_DENSITY_CONFIG: DensityConfig = {
  targetDensity: 0.95,
  minDensity: 0.80,
  synthesisPlatform: "twist",
  maxOligoLength: 300,
};

export interface DensityResult {
  /** Optimized CodecConfig */
  config: CodecConfig;
  /** Achieved density in bits/nt */
  achievedDensity: number;
  /** Payload efficiency (payload bytes / total oligo bytes) */
  payloadEfficiency: number;
  /** Estimated outer RS overhead */
  outerRSOverhead: number;
  /** Estimated inner code overhead */
  innerCodeOverhead: number;
  /** Estimated constraint overhead */
  constraintOverhead: number;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Clamp a value to [lo, hi]. */
function clamp(value: number, lo: number, hi: number): number {
  return Math.max(lo, Math.min(hi, value));
}

/**
 * Approximate mapping rate (bits/nt) per mapping mode.
 * Used for density estimation without requiring a full encode.
 */
const MAPPING_RATE: Record<string, number> = {
  direct: 2.0,
  constrained: 2.0,
  srt: 2.0,
  goldman: 1.538,
  arithmetic: 1.9,
};

// ---------------------------------------------------------------------------
// optimalOuterParity
// ---------------------------------------------------------------------------

/**
 * Compute the optimal outer parity ratio for a given file size and channel.
 *
 * The key insight: smaller files need proportionally MORE parity because
 * RS overhead dominates when there are few oligos. Larger files can use
 * LESS parity because the law of large numbers gives reliable recovery
 * even with a thinner erasure code.
 *
 * Channel-specific formulas:
 *   - Illumina:  clamp(0.05 + 500/fileSize, 0.05, 0.15)
 *   - Nanopore:  clamp(0.3  + 5000/fileSize, 0.3, 0.5)
 *   - PacBio:    clamp(0.2  + 3000/fileSize, 0.2, 0.4)
 *
 * @param fileSize - Input file size in bytes
 * @param channel  - Sequencing channel: "illumina" | "nanopore" | "pacbio"
 */
export function optimalOuterParity(
  fileSize: number,
  channel: string,
): number {
  const safeSize = Math.max(1, fileSize); // guard against zero

  switch (channel) {
    case "nanopore":
      return clamp(0.3 + 5000 / safeSize, 0.3, 0.5);

    case "pacbio":
      return clamp(0.2 + 3000 / safeSize, 0.2, 0.4);

    case "illumina":
    default:
      return clamp(0.05 + 500 / safeSize, 0.05, 0.15);
  }
}

// ---------------------------------------------------------------------------
// optimalInnerParity
// ---------------------------------------------------------------------------

/**
 * Compute the optimal inner parity bytes for a given channel.
 *
 * Higher error-rate channels need more inner parity to correct per-oligo
 * errors before the outer RS sees them. The oligo length also matters:
 * longer oligos accumulate more errors per read, so they need more parity.
 *
 * Channel defaults:
 *   - Illumina:  4 bytes  (low sub rate ~0.1%)
 *   - Nanopore:  6 bytes  (high IDS ~5-15%)
 *   - PacBio:    5 bytes  (moderate IDS ~2-5%)
 *
 * For very long oligos (>300nt), we add one extra byte per 100nt above 200nt
 * to account for the increased per-read error count.
 *
 * @param channel    - Sequencing channel
 * @param oligoLength - Total oligo length in nucleotides
 */
export function optimalInnerParity(
  channel: string,
  oligoLength: number,
): number {
  let base: number;

  switch (channel) {
    case "nanopore":
      base = 6;
      break;
    case "pacbio":
      base = 5;
      break;
    case "illumina":
    default:
      base = 4;
      break;
  }

  // Add extra parity for long oligos (more errors accumulate per read)
  if (oligoLength > 200) {
    base += Math.floor((oligoLength - 200) / 100);
  }

  return base;
}

// ---------------------------------------------------------------------------
// computeDensity
// ---------------------------------------------------------------------------

/**
 * Compute the actual density of a given config + file size.
 *
 * Density = fileBits / totalNucleotides
 *
 * where totalNucleotides = numOligos × oligoLength
 *
 * This uses the same analytical model as presets.ts:computeDensity but
 * accounts for file-size-dependent outer RS overhead.
 *
 * @param cfg      - Codec configuration
 * @param fileSize - Input file size in bytes
 * @returns Density in bits/nt
 */
export function computeDensity(cfg: CodecConfig, fileSize: number): number {
  const rate = MAPPING_RATE[cfg.mappingMode ?? "direct"] ?? 2.0;
  const innerNt = cfg.oligoLength - 2 * cfg.primerLength;
  const innerBits = innerNt * rate;
  const innerBytes = innerBits / 8;

  // Per-oligo overhead: 4B address + 2B CRC-16 + innerParityBytes
  const overheadBytes = 4 + 2 + cfg.innerParityBytes;
  const payloadBytes = innerBytes - overheadBytes;

  if (payloadBytes <= 0) return 0; // config cannot carry any payload

  const payloadBits = payloadBytes * 8;

  // Outer RS efficiency: only (1 - outerParityRatio) of oligos carry data
  const outerEfficiency = 1 / (1 + cfg.outerParityRatio);
  const infoBits = payloadBits * outerEfficiency;

  return infoBits / cfg.oligoLength;
}

// ---------------------------------------------------------------------------
// optimizeForDensity
// ---------------------------------------------------------------------------

/**
 * Optimize a CodecConfig for maximum density.
 *
 * Adjusts outerParityRatio, innerParityBytes, mappingMode, and other
 * parameters to approach the target density.
 *
 * Strategy:
 *   1. Use channel-specific optimalOuterParity() for the outer RS ratio
 *   2. Use channel-specific optimalInnerParity() for the inner parity bytes
 *   3. Switch to "constrained" mapping (zero retries, full 2.0 bits/nt)
 *   4. Tune interleaveDepth: 4 for files > 10KB, 8 for nanopore
 *   5. For nanopore: enable convolutional inner code
 *   6. Respect the synthesizer's maxOligoLength
 *
 * @param cfg        - Base codec configuration to optimize
 * @param fileSize   - Input file size in bytes
 * @param densityCfg - Density tuning parameters (optional, uses defaults)
 */
export function optimizeForDensity(
  cfg: CodecConfig,
  fileSize: number,
  densityCfg?: DensityConfig,
): DensityResult {
  const dc = densityCfg ?? DEFAULT_DENSITY_CONFIG;
  const channel = cfg.channel ?? "illumina";

  // --- 1. Outer parity ratio (channel + file-size dependent) ---
  const outerParityRatio = optimalOuterParity(fileSize, channel);

  // --- 2. Inner parity bytes (channel dependent) ---
  const oligoLength = Math.min(cfg.oligoLength, dc.maxOligoLength);
  const innerParityBytes = optimalInnerParity(channel, oligoLength);

  // --- 3. Mapping mode: switch to constrained for density ---
  // Constrained mapping gives full 2.0 bits/nt with zero retries,
  // unlike seed-retry direct mapping which wastes 1 byte per oligo on the seed.
  const mappingMode = "constrained" as const;

  // --- 4. Interleave depth ---
  let interleaveDepth = cfg.interleaveDepth ?? 0;
  if (fileSize > 10_000) {
    interleaveDepth = channel === "nanopore" ? 8 : 4;
  }

  // --- 5. Build optimized config ---
  const optimized: CodecConfig = {
    ...cfg,
    oligoLength,
    outerParityRatio,
    innerParityBytes,
    mappingMode,
    maxRetries: 0, // constrained mode never needs retries
    interleaveDepth,
  };

  // Nanopore-specific: enable convolutional inner code for indel tolerance
  if (channel === "nanopore" && !cfg.useConvolutionalInner) {
    optimized.useConvolutionalInner = true;
    // Shorter oligos for nanopore: keep per-block indel count within
    // K=9 Viterbi maxDrift=15 correction radius
    if (optimized.oligoLength > 150) {
      optimized.oligoLength = 150;
    }
  }

  // --- 6. Compute achieved density ---
  const achievedDensity = computeDensity(optimized, fileSize);

  // --- 7. Compute overhead breakdown ---
  const innerNt = optimized.oligoLength - 2 * optimized.primerLength;
  const rate = MAPPING_RATE[mappingMode] ?? 2.0;
  const innerBytes = (innerNt * rate) / 8;
  const payloadBytes = innerBytes - 4 - 2 - innerParityBytes;

  // Outer RS overhead = fraction of oligos that are parity
  const outerRSOverhead = outerParityRatio / (1 + outerParityRatio);

  // Inner code overhead = innerParityBytes / innerBytes
  const innerCodeOverhead =
    innerBytes > 0 ? innerParityBytes / innerBytes : 0;

  // Constraint overhead: constrained mode has 0 bytes per-oligo overhead
  // (unlike seed-retry which wastes 1 byte on the seed)
  const constraintOverhead =
    cfg.mappingMode === "direct" && cfg.maxRetries > 0
      ? 1 / innerBytes // legacy seed-retry cost
      : 0;

  // Payload efficiency = payload / total inner bytes
  const payloadEfficiency =
    innerBytes > 0 ? Math.max(0, payloadBytes) / innerBytes : 0;

  return {
    config: optimized,
    achievedDensity,
    payloadEfficiency,
    outerRSOverhead,
    innerCodeOverhead,
    constraintOverhead,
  };
}
