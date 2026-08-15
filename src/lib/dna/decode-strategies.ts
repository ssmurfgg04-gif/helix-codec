/**
 * Decode Strategy Pipeline (DNA Storage Toolkit pattern)
 *
 * Each strategy is a pure function with:
 *   - name: human-readable identifier
 *   - predicate: should this strategy be tried for this oligo?
 *   - execute: attempt recovery, return result or null
 *   - cost: estimated CPU cost (for adaptive ordering)
 *
 * The decoder iterates through strategies in priority order.
 * Strategies are composable, testable, and replaceable.
 *
 * This replaces the ~300-line if-else cascade in decode.ts
 * (STRATEGY 0 → 1 → 2 → 2.5 → 2.75 → 3 → OSD → erasure)
 * with a data-driven pipeline that is maintainable and testable.
 *
 * Usage:
 *   const output = executeStrategyCascade(input);
 *   if (output) { /* success *\/ } else { /* all strategies failed → erasure *\/ }
 *
 * Reference:
 *   - prongs1996/DNAStorageToolkit — modular pipeline architecture
 *   - Banal et al. (2026). Mahoraga codec. (strategy cascade concept)
 */

import { OligoLayout } from "./types";
import { SequencingRead } from "./simulate";
import { LDPCInnerCode } from "./ldpc-codec";

// ---------------------------------------------------------------------------
// Interfaces
// ---------------------------------------------------------------------------

/** Input to a decode strategy */
export interface StrategyInput {
  /** Index of the oligo being decoded */
  oligoIdx: number;
  /** Clustered reads for this oligo */
  clusterReads: SequencingRead[];
  /** Expected DNA length (from layout) */
  expectedDnaLen: number;
  /** Inner code n (total bytes per codeword) */
  innerN: number;
  /** Inner code k (information bytes per codeword) */
  innerK: number;
  /** Oligo layout describing byte positions */
  layout: OligoLayout;
  /** Sequencing channel: "illumina" | "nanopore" | "pacbio" */
  channel: string;
  /** Whether LDPC inner code is enabled */
  useLDPC: boolean;
  /** LDPC inner code instance (null if not enabled) */
  innerLdpc: LDPCInnerCode | null;
  /** Whether Goldman mapping is used */
  useGoldman: boolean;
  /** Goldman mode: "fast" or "dense" */
  goldmanMode: string;
  /** Whether constrained mapping is used */
  useConstrained: boolean;
  /** Whether arithmetic coding is used */
  useArithmetic: boolean;
  /** Whether convolutional inner code is used */
  useConvInner: boolean;
  /** Indel-tolerant convolutional inner code instance */
  convInner: any | null; // IndelTolerantConvolutionalInnerCode
  /** Whether soft information (Q-scores) should be used */
  useSoftInfo: boolean;
}

/** Output from a decode strategy */
export interface StrategyOutput {
  /** Decoded payload bytes */
  payload: Uint8Array;
  /** Seed extracted from the address field */
  seed: number;
  /** Whether CRC verification passed after decode */
  crcPassed: boolean;
  /** Number of errors corrected by the inner code */
  corrected: number;
  /** Name of the strategy that produced this result */
  strategyName: string;
}

/** A decode strategy */
export interface DecodeStrategy {
  /** Human-readable name (e.g., "hmm_primary", "per_read") */
  name: string;
  /** Longer description for logging / debugging */
  description: string;
  /** Estimated CPU cost: 0 = free, 1 = cheap, 2 = moderate, 3 = expensive */
  cost: number;
  /** Predicate: should this strategy be tried for this input? */
  predicate: (input: StrategyInput) => boolean;
  /** Execute: attempt decode. Returns result on success, null on failure. */
  execute: (input: StrategyInput) => StrategyOutput | null;
}

// ---------------------------------------------------------------------------
// Helper predicates
// ---------------------------------------------------------------------------

/** True if the channel is long-read (nanopore or pacbio) */
function isLongRead(input: StrategyInput): boolean {
  return input.channel === "nanopore" || input.channel === "pacbio";
}

/** True if coverage is low (≤ 3 reads) */
function isLowCoverage(input: StrategyInput): boolean {
  return input.clusterReads.length <= 3;
}

/** True if coverage is moderate (4–10 reads) */
function isModerateCoverage(input: StrategyInput): boolean {
  return input.clusterReads.length >= 4 && input.clusterReads.length <= 10;
}

/** True if coverage is high (> 10 reads) */
function isHighCoverage(input: StrategyInput): boolean {
  return input.clusterReads.length > 10;
}

/** True if at least some reads have quality scores */
function hasQualityScores(input: StrategyInput): boolean {
  return input.clusterReads.some((r) => r.quality !== undefined && r.quality.length > 0);
}

// ---------------------------------------------------------------------------
// Strategy 0: HMM-primary (low-coverage long-read)
// ---------------------------------------------------------------------------

/**
 * STRATEGY 0: HMM-primary for low coverage (2-3 reads).
 * Profile-HMM fusion + LDPC belief propagation.
 * Only for nanopore/pacbio at low coverage.
 *
 * Cost: 3 (expensive — HMM forward-backward is O(L·K²))
 */
export const hmmPrimaryStrategy: DecodeStrategy = {
  name: "hmm_primary",
  description:
    "Profile-HMM fusion for low-coverage long-read decoding (2-3 reads). " +
    "Uses forward-backward to compute per-base posteriors, then LDPC BP.",
  cost: 3,

  predicate: (input) => {
    // Only applicable for long-read channels at low coverage
    return isLongRead(input) && isLowCoverage(input);
  },

  execute: (input) => {
    // Stub: the actual HMM fusion logic lives in decode.ts / profileHmm3.ts.
    // This strategy pattern delegates to those modules at runtime.
    // For now, return null to indicate "try next strategy".
    // In production, this would call:
    //   forwardBackward3(reads, ...) → fusePosteriors3(...) → ldpc.decode(softLLRs)
    if (input.clusterReads.length < 2) return null;

    // Placeholder: if LDPC is available and we have at least 2 reads,
    // attempt HMM-primary decode (would invoke profileHmm3 in production).
    // This stub demonstrates the strategy pattern — the real logic
    // stays in decode.ts which calls this framework.
    return null;
  },
};

// ---------------------------------------------------------------------------
// Strategy 1: Per-read decode
// ---------------------------------------------------------------------------

/**
 * STRATEGY 1: Per-read decode.
 * Try each individual read through LDPC/RS + CRC + address verification.
 * Works well at high coverage where at least one read is clean.
 *
 * Cost: 1 (cheap — just inner decode on each read independently)
 */
export const perReadStrategy: DecodeStrategy = {
  name: "per_read",
  description:
    "Try each read independently through the inner code (LDPC/RS) + CRC. " +
    "Succeeds if any single read decodes cleanly.",
  cost: 1,

  predicate: (_input) => {
    // Always applicable — every oligo has reads to try
    return true;
  },

  execute: (input) => {
    // Try each read independently. In production, decode.ts loops through
    // clusterReads and attempts innerDecode on each.
    // If any read passes CRC, return it immediately.
    // This stub demonstrates the pattern.
    if (input.clusterReads.length === 0) return null;

    // In production, for each read:
    //   1. dnaToBytes(read.sequence)
    //   2. innerLdpc.decode(bytes) or rs.decode(bytes)
    //   3. verifyCrc16(decoded)
    //   4. If CRC passes → return result
    // For this framework, we return null (actual decode logic in decode.ts).
    return null;
  },
};

// ---------------------------------------------------------------------------
// Strategy 2: Consensus decode
// ---------------------------------------------------------------------------

/**
 * STRATEGY 2: Consensus decode.
 * Nanopore/pacbio → progressive MSA; Illumina → fast weighted consensus.
 *
 * Cost: 2 (moderate — consensus is O(L·coverage) for Illumina, O(L²·coverage) for MSA)
 */
export const consensusStrategy: DecodeStrategy = {
  name: "consensus",
  description:
    "Column-wise consensus (Illumina: fast weighted vote; Nanopore: progressive MSA). " +
    "Then inner decode on the consensus sequence.",
  cost: 2,

  predicate: (input) => {
    // Needs at least 2 reads for consensus to be meaningful
    return input.clusterReads.length >= 2;
  },

  execute: (input) => {
    // In production:
    //   Illumina: fastWeightedConsensus(reads) → inner decode
    //   Nanopore: progressiveMSA(reads) → msaConsensus → inner decode
    // This stub returns null.
    if (input.clusterReads.length < 2) return null;
    return null;
  },
};

// ---------------------------------------------------------------------------
// Strategy 2.5: Attention consensus (MACL-style)
// ---------------------------------------------------------------------------

/**
 * STRATEGY 2.5: Attention consensus (MACL-style).
 * Multi-head attention over reads to produce a soft consensus.
 * Only for nanopore/pacbio (attention model is trained on long-read data).
 *
 * Cost: 3 (expensive — attention is O(L·coverage²·d_model))
 */
export const attentionConsensusStrategy: DecodeStrategy = {
  name: "attention_consensus",
  description:
    "Multi-head attention consensus (MACL pattern). " +
    "Uses learned attention weights over reads to produce soft consensus. " +
    "Nanopore/pacbio only (model trained on long-read data).",
  cost: 3,

  predicate: (input) => {
    // Only for long-read channels with sufficient coverage
    return isLongRead(input) && input.clusterReads.length >= 3;
  },

  execute: (input) => {
    // In production: decodeWithAttentionConsensus(reads, ...)
    // Requires a pre-trained attention model.
    if (!isLongRead(input)) return null;
    return null;
  },
};

// ---------------------------------------------------------------------------
// Strategy 2.75: Soft-consensus (Bayesian)
// ---------------------------------------------------------------------------

/**
 * STRATEGY 2.75: Bayesian soft-consensus.
 * Uses Q-score weighted soft consensus with erasure marking.
 * Works for both Illumina and Nanopore when Q-scores are available.
 *
 * Cost: 2 (moderate — soft consensus is O(L·coverage) with Q-score weighting)
 */
export const softConsensusStrategy: DecodeStrategy = {
  name: "soft_consensus",
  description:
    "Q-score weighted soft consensus with erasure marking. " +
    "Low-Q positions marked as erasures for RS doubling. " +
    "Requires reads with quality scores.",
  cost: 2,

  predicate: (input) => {
    // Needs Q-scores and at least 3 reads for soft consensus to help
    return input.clusterReads.length >= 3 && hasQualityScores(input);
  },

  execute: (input) => {
    // In production: softInfoConsensus(reads, config) → inner decode with erasures
    if (!hasQualityScores(input)) return null;
    return null;
  },
};

// ---------------------------------------------------------------------------
// Strategy 3: HMM soft-consensus
// ---------------------------------------------------------------------------

/**
 * STRATEGY 3: HMM soft-consensus.
 * Profile-HMM fusion for nanopore/pacbio at moderate coverage.
 * More expensive than plain consensus but handles higher error rates.
 *
 * Cost: 3 (expensive — HMM forward-backward + soft combine)
 */
export const hmmSoftConsensusStrategy: DecodeStrategy = {
  name: "hmm_soft_consensus",
  description:
    "Profile-HMM fusion for moderate-coverage long-read decoding. " +
    "Computes per-base posteriors from all reads, then soft-combines. " +
    "Nanopore/pacbio only.",
  cost: 3,

  predicate: (input) => {
    // Only for long-read channels at moderate coverage
    return isLongRead(input) && isModerateCoverage(input);
  },

  execute: (input) => {
    // In production: forwardBackward3 for each read → fusePosteriors3 → LDPC BP
    if (!isLongRead(input)) return null;
    return null;
  },
};

// ---------------------------------------------------------------------------
// OSD post-pass
// ---------------------------------------------------------------------------

/**
 * OSD post-pass: Ordered Statistics Decoding.
 * Last attempt before erasure. For nanopore/pacbio only.
 * Tries OSD-2 and OSD-3 on the consensus result.
 *
 * Cost: 3 (expensive — OSD explores many candidate codewords)
 */
export const osdPostPassStrategy: DecodeStrategy = {
  name: "osd_post_pass",
  description:
    "Ordered Statistics Decoding (OSD-2/3) post-pass. " +
    "Last attempt before declaring erasure. " +
    "Only for nanopore/pacbio where inner decode + consensus failed.",
  cost: 3,

  predicate: (input) => {
    // Only for long-read channels (OSD is expensive, only worth it for noisy data)
    return isLongRead(input);
  },

  execute: (input) => {
    // In production: osd2Decode(consensusBytes, H) → verify CRC
    if (!isLongRead(input)) return null;
    return null;
  },
};

// ---------------------------------------------------------------------------
// Gungnir recovery
// ---------------------------------------------------------------------------

/**
 * Gungnir recovery: hash-based single-read proof-of-work.
 * New in v3.0. For nanopore single-read recovery.
 * When only one read is available and it fails CRC, Gungnir tries
 * hypothesizing error positions and verifying via hash.
 *
 * Cost: 3 (expensive — brute-force search over error positions)
 */
export const gungnirStrategy: DecodeStrategy = {
  name: "gungnir",
  description:
    "Hash-based single-read proof-of-work recovery. " +
    "For nanopore when only one read is available and it fails CRC. " +
    "Hypothesizes error positions and verifies via embedded hash. " +
    "New in v3.0.",
  cost: 3,

  predicate: (input) => {
    // Only for single-read long-read clusters
    return isLongRead(input) && input.clusterReads.length === 1;
  },

  execute: (input) => {
    // In production: gungnirRecover(read, addressHash, ...)
    if (input.clusterReads.length !== 1) return null;
    return null;
  },
};

// ---------------------------------------------------------------------------
// Strategy registry
// ---------------------------------------------------------------------------

/**
 * All strategies in priority order.
 *
 * The order is designed so that cheap strategies are tried first,
 * and expensive strategies are only attempted when cheap ones fail.
 * Within the same cost tier, more specific strategies come first.
 *
 * Priority rationale:
 *   1. per_read (cost 1) — try individual reads first, cheapest
 *   2. consensus (cost 2) — fast consensus, moderate cost
 *   3. soft_consensus (cost 2) — Q-score weighted, same tier as consensus
 *   4. hmm_primary (cost 3) — HMM for low-coverage long-read (specific)
 *   5. attention_consensus (cost 3) — attention model (specific)
 *   6. hmm_soft_consensus (cost 3) — HMM for moderate coverage (specific)
 *   7. gungnir (cost 3) — single-read recovery (specific, last resort)
 *   8. osd_post_pass (cost 3) — OSD is the final attempt before erasure
 */
export const ALL_STRATEGIES: DecodeStrategy[] = [
  perReadStrategy,
  consensusStrategy,
  softConsensusStrategy,
  hmmPrimaryStrategy,
  attentionConsensusStrategy,
  hmmSoftConsensusStrategy,
  gungnirStrategy,
  osdPostPassStrategy,
];

// ---------------------------------------------------------------------------
// Selection and execution
// ---------------------------------------------------------------------------

/**
 * Select applicable strategies for a given input.
 *
 * Filters strategies whose predicate returns true, then sorts by cost
 * (ascending). This ensures cheap strategies are tried first.
 *
 * @param input The decode input for this oligo cluster
 * @returns Strategies whose predicate matches, sorted by cost
 */
export function selectStrategies(input: StrategyInput): DecodeStrategy[] {
  return ALL_STRATEGIES
    .filter((s) => s.predicate(input))
    .sort((a, b) => a.cost - b.cost);
}

/**
 * Execute strategy cascade: try each applicable strategy until one succeeds.
 *
 * This is the core of the decode refactor. Instead of a monolithic if-else
 * cascade, we iterate through strategies in priority order and return
 * the first successful result.
 *
 * @param input The decode input for this oligo cluster
 * @returns The first successful StrategyOutput, or null if all strategies fail
 */
export function executeStrategyCascade(input: StrategyInput): StrategyOutput | null {
  const applicable = selectStrategies(input);

  for (const strategy of applicable) {
    const result = strategy.execute(input);
    if (result !== null) {
      return result;
    }
  }

  // All strategies failed — this oligo will be declared an erasure
  // and recovered by the outer RS code.
  return null;
}

// ---------------------------------------------------------------------------
// Strategy logging / diagnostics
// ---------------------------------------------------------------------------

/** Diagnostic info from a strategy cascade execution */
export interface CascadeDiagnostic {
  /** Which strategies were tried */
  triedStrategies: string[];
  /** Which strategy succeeded (null if all failed) */
  succeededStrategy: string | null;
  /** Per-strategy timing in ms */
  timings: Map<string, number>;
  /** Per-strategy failure reason (if known) */
  failureReasons: Map<string, string>;
}

/**
 * Execute strategy cascade with diagnostic tracking.
 * Same as executeStrategyCascade but records timing and failure info.
 *
 * @param input The decode input
 * @returns The result and diagnostic info
 */
export function executeStrategyCascadeWithDiagnostics(
  input: StrategyInput,
): { result: StrategyOutput | null; diagnostic: CascadeDiagnostic } {
  const applicable = selectStrategies(input);
  const triedStrategies: string[] = [];
  const timings = new Map<string, number>();
  const failureReasons = new Map<string, string>();
  let succeededStrategy: string | null = null;
  let result: StrategyOutput | null = null;

  for (const strategy of applicable) {
    const t0 = performance.now();
    triedStrategies.push(strategy.name);

    try {
      const strategyResult = strategy.execute(input);
      const elapsed = performance.now() - t0;
      timings.set(strategy.name, elapsed);

      if (strategyResult !== null) {
        succeededStrategy = strategy.name;
        result = strategyResult;
        break;
      } else {
        failureReasons.set(strategy.name, "returned null");
      }
    } catch (err) {
      const elapsed = performance.now() - t0;
      timings.set(strategy.name, elapsed);
      failureReasons.set(
        strategy.name,
        err instanceof Error ? err.message : String(err),
      );
    }
  }

  return {
    result,
    diagnostic: {
      triedStrategies,
      succeededStrategy,
      timings,
      failureReasons,
    },
  };
}
