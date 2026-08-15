/**
 * Soft-Information Consensus Decoder
 *
 * Uses per-base Phred quality scores (Q-scores) from simulated reads to build
 * a confidence-aware consensus. Low-Q positions are marked as erasures for
 * Reed-Solomon decoding, doubling correction capacity (erasures cost 1 parity
 * symbol vs. 2 for unknown errors).
 *
 * This is the "cheap version" of Banal et al. 2026's profile-HMM approach:
 *   - Banal: forward-backward on profile HMM → per-position posteriors → LLRs → OSD
 *   - This:  column-wise plurality weighted by Q-scores → per-position confidence → erasure hints
 *
 * The full HMM approach handles indels inside the posterior (no pre-alignment
 * needed). This simpler version assumes reads are roughly aligned (no major
 * indels), which is true for Illumina-style sub-dominant channels.
 *
 * Reference:
 *   - Banal et al. (2026). "Mahoraga." arXiv:2604.20810.
 *   - Schwarz et al. (2024). Q<10 erasure threshold. PMC11066528.
 */

import { SequencingRead } from "./simulate";

export interface ConsensusResult {
  /** Consensus DNA sequence. */
  sequence: string;
  /** Per-base confidence (0..1). Low confidence → erasure hint. */
  confidence: Float32Array;
  /** Positions marked as erasures (Q < threshold). */
  erasurePositions: number[];
  /** Per-position base distribution (for visualization). */
  baseCounts: { A: number; C: number; G: number; T: number }[];
  /** Number of reads used for consensus. */
  readCount: number;
}

export interface SoftInfoConfig {
  /** Q-score threshold below which a position is marked as erasure. Default: 10. */
  erasureQThreshold: number;
  /** Minimum confidence (0..1) for a position to NOT be an erasure. Default: 0.5. */
  minConfidence: number;
  /** Whether to use Q-weighted voting (true) or simple plurality (false). */
  useQualityWeights: boolean;
}

export const DEFAULT_SOFT_INFO_CONFIG: SoftInfoConfig = {
  erasureQThreshold: 10,
  minConfidence: 0.5,
  useQualityWeights: true,
};

/**
 * Convert Phred Q-score to log-likelihood ratio (LLR).
 *
 * LLR = log(P(b=0|r) / P(b=1|r)) = Q * ln(10) / 10 ≈ Q * 0.2303
 *
 * Positive LLR → more likely 0, negative → more likely 1.
 * High Q-score → high |LLR| → high confidence.
 *
 * This is used to convert per-base Q-scores from sequencing into soft-decision
 * inputs for the Viterbi decoder, providing 2-3 dB coding gain over
 * hard-decision Hamming distance.
 */
export function qScoreToLLR(qScore: number): number {
  // LLR = Q * ln(10) / 10 ≈ Q * 0.2303
  return qScore * 0.2303;
}

/**
 * Convert an array of Phred Q-scores to LLRs (Float32Array).
 *
 * For DNA bases mapped to bits, each base contributes 2 bits.
 * The Q-score for a base is applied to both bits (they share
 * the same confidence since they come from the same base call).
 *
 * @param qScores  Per-base Q-scores (one per DNA base).
 * @param bitsPerBase  Number of bits per base (default 2 for 2-bit mapping).
 * @returns Per-bit LLRs (length = qScores.length * bitsPerBase).
 */
export function qScoresToBitLLRs(qScores: ArrayLike<number>, bitsPerBase = 2): Float32Array {
  const llrs = new Float32Array(qScores.length * bitsPerBase);
  for (let i = 0; i < qScores.length; i++) {
    const llr = qScoreToLLR(qScores[i]);
    for (let b = 0; b < bitsPerBase; b++) {
      llrs[i * bitsPerBase + b] = llr;
    }
  }
  return llrs;
}

/**
 * Compute Phred Q-score → probability of error.
 * Q = -10 * log10(P_error), so P_error = 10^(-Q/10).
 */
export function qToErrorProb(q: number): number {
  return Math.pow(10, -q / 10);
}

/**
 * Compute confidence from Q-score: 1 - P_error.
 */
export function qToConfidence(q: number): number {
  return 1 - qToErrorProb(q);
}

/**
 * Build a soft-information consensus from a set of reads.
 *
 * Algorithm:
 *   1. For each position, collect (base, Q-score) from all reads covering it.
 *   2. Weight each vote by Q-score confidence (1 - P_error).
 *   3. Pick the base with highest weighted vote.
 *   4. Compute confidence = winning_vote_weight / total_weight.
 *   5. If confidence < minConfidence OR mean Q < erasureQThreshold, mark as erasure.
 *
 * Erasure positions are passed to the RS decoder, which can correct them at
 * half the cost of unknown errors.
 */
export function softInfoConsensus(
  reads: SequencingRead[],
  cfg: SoftInfoConfig = DEFAULT_SOFT_INFO_CONFIG,
): ConsensusResult {
  if (reads.length === 0) {
    return {
      sequence: "",
      confidence: new Float32Array(0),
      erasurePositions: [],
      baseCounts: [],
      readCount: 0,
    };
  }
  if (reads.length === 1) {
    // Single read: use its Q-scores directly
    const read = reads[0];
    const confidence = new Float32Array(read.sequence.length);
    const erasurePositions: number[] = [];
    const baseCounts = read.sequence.split("").map((base, i) => {
      const q = read.quality?.[i] ?? 30;
      const conf = qToConfidence(q);
      confidence[i] = conf;
      if (q < cfg.erasureQThreshold || conf < cfg.minConfidence) {
        erasurePositions.push(i);
      }
      const counts = { A: 0, C: 0, G: 0, T: 0 };
      counts[base as keyof typeof counts] = 1;
      return counts;
    });
    return {
      sequence: read.sequence,
      confidence,
      erasurePositions,
      baseCounts,
      readCount: 1,
    };
  }

  // Find consensus length (median of read lengths)
  const lengths = reads.map((r) => r.sequence.length);
  const sorted = lengths.slice().sort((a, b) => a - b);
  const medianLen = sorted[Math.floor(sorted.length / 2)];

  const sequence: string[] = [];
  const confidence = new Float32Array(medianLen);
  const erasurePositions: number[] = [];
  const baseCounts: { A: number; C: number; G: number; T: number }[] = [];

  for (let col = 0; col < medianLen; col++) {
    const votes: Record<string, number> = { A: 0, C: 0, G: 0, T: 0 };
    let totalWeight = 0;
    let totalQ = 0;
    let qCount = 0;
    let lowQCount = 0; // count of reads with Q < threshold at this position
    let minQAtCol = 40;
    const counts = { A: 0, C: 0, G: 0, T: 0 };

    for (const read of reads) {
      if (col >= read.sequence.length) continue;
      const base = read.sequence[col];
      if (!(base in votes)) continue;
      counts[base as keyof typeof counts]++;

      let weight: number;
      if (cfg.useQualityWeights && read.quality) {
        const q = read.quality[col] ?? 30;
        weight = qToConfidence(q);
        totalQ += q;
        qCount++;
        if (q < cfg.erasureQThreshold) lowQCount++;
        if (q < minQAtCol) minQAtCol = q;
      } else {
        weight = 1;
      }
      votes[base] += weight;
      totalWeight += weight;
    }

    // Pick base with highest weighted vote
    let bestBase = "A";
    let bestWeight = -1;
    for (const b of ["A", "C", "G", "T"]) {
      if (votes[b] > bestWeight) {
        bestWeight = votes[b];
        bestBase = b;
      }
    }

    sequence.push(bestBase);
    baseCounts.push(counts);

    // Confidence = winning_weight / total_weight
    const conf = totalWeight > 0 ? bestWeight / totalWeight : 0;
    confidence[col] = conf;

    // Mean Q-score for this position
    const meanQ = qCount > 0 ? totalQ / qCount : 30;

    // Mark as erasure if:
    //   - Low confidence (votes disagree), OR
    //   - Min Q at this position is low (at least one read is uncertain), OR
    //   - Mean Q is low
    // The "min Q" criterion catches positions where one read has a substitution
    // (low Q) even if the consensus is correct via majority vote.
    if (conf < cfg.minConfidence || minQAtCol < cfg.erasureQThreshold || meanQ < cfg.erasureQThreshold) {
      erasurePositions.push(col);
    }
  }

  return {
    sequence: sequence.join(""),
    confidence,
    erasurePositions,
    baseCounts,
    readCount: reads.length,
  };
}

/**
 * Cap the erasure list to a maximum size, keeping only the lowest-confidence positions.
 * This prevents exceeding the RS decoder's erasure capacity (which is nsym erasures max).
 */
export function capErasures(
  erasurePositions: number[],
  confidence: Float32Array,
  maxErasures: number,
): number[] {
  if (erasurePositions.length <= maxErasures) return erasurePositions;
  // Sort by confidence ascending (lowest first = most likely to be errors)
  const sorted = erasurePositions.slice().sort((a, b) => confidence[a] - confidence[b]);
  return sorted.slice(0, maxErasures).sort((a, b) => a - b);
}

/**
 * Convert erasure positions (in DNA coordinates) to erasure positions (in byte coordinates).
 *
 * Since we use 2-bit DNA mapping (4 bases = 1 byte), an erasure at any of the 4
 * bases in a byte marks the whole byte as an erasure.
 */
export function dnaErasureToByteErasure(
  dnaErasurePositions: number[],
): number[] {
  const byteErasures = new Set<number>();
  for (const pos of dnaErasurePositions) {
    byteErasures.add(Math.floor(pos / 4));
  }
  return Array.from(byteErasures).sort((a, b) => a - b);
}

/**
 * Map erasure positions from DNA-string coordinates to inner-block byte coordinates.
 *
 * The inner block starts after the forward primer, so we subtract the primer offset.
 */
export function mapErasuresToInnerBlock(
  dnaErasurePositions: number[],
  primerLength: number,
): number[] {
  return dnaErasurePositions
    .map((pos) => pos - primerLength)
    .filter((pos) => pos >= 0)
    .map((pos) => Math.floor(pos / 4));
}
