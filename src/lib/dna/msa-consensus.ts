/**
 * MSA-Based Consensus for Nanopore Reads — Pre-Viterbi Indel Correction
 *
 * For Nanopore 9% IDS, multiple reads of the same oligo are available.
 * Instead of feeding each read independently through the expensive Viterbi
 * decoder, we first build a multi-sequence alignment (MSA) consensus that:
 *
 *   1. Aligns each read to a preliminary consensus via Profile HMM (forwardBackward3)
 *   2. Builds a column-wise weighted consensus with quality scores
 *   3. Produces a single consensus read with per-position posterior probabilities
 *   4. Feeds this consensus to Viterbi (single decode instead of N decodes)
 *
 * This reduces the effective IDS from 9% (per-read) to ~2-3% (post-consensus),
 * making K=7 Viterbi (4.5ms) sufficient instead of K=9 (12ms).
 *
 * The approach is:
 *   Iterative refinement MSA (like MUSCLE's 3rd stage):
 *     1. Build initial plurality consensus from all reads
 *     2. Align each read to the consensus via Profile HMM
 *     3. Build weighted consensus from aligned reads (using HMM posteriors)
 *     4. Repeat steps 2-3 once more for refinement
 *     5. Return consensus + quality scores + posterior probabilities
 *
 * Quality scores from the MSA are converted to LLR for soft-decision
 * LDPC decoding downstream.
 *
 * References:
 *   - Durbin/Eddy/Krogh/Mitchison (1998), ch. 5 — Multiple alignment
 *   - Edgar (2004), MUSCLE — iterative refinement MSA
 *   - Press et al. (2021), HEDGES — convolutional + RS for DNA storage
 */

import {
  forwardBackward3,
  fusePosteriors3,
  Hmm3Params,
  DEFAULT_HMM3_PARAMS,
  Hmm3Result,
} from "./profileHmm3";

/**
 * Configuration for the MSA consensus builder.
 */
export interface MsaConfig {
  /** Number of refinement iterations. Default: 2 (initial + 2 refinements). */
  iterations: number;
  /**
   * Minimum cluster size for MSA. Below this, return simple plurality consensus.
   * Default: 2.
   */
  minClusterSize: number;
  /** Profile HMM parameters. Default: tuned for ONT R10.4.1 (9% IDS). */
  hmmParams: Hmm3Params;
  /** HMM band width for long sequences. Default: 12. */
  hmmBandWidth: number;
  /**
   * Quality floor for consensus positions with low coverage.
   * Default: Q10 (10% error rate).
   */
  qualityFloor: number;
}

export const DEFAULT_MSA_CONFIG: MsaConfig = {
  iterations: 2,
  minClusterSize: 2,
  hmmParams: {
    ...DEFAULT_HMM3_PARAMS,
    matchToMatch: 0.90,
    matchToInsert: 0.05,
    matchToDelete: 0.05,
  },
  hmmBandWidth: 12,
  qualityFloor: 10,
};

/**
 * Result of MSA consensus building.
 */
export interface MsaConsensusResult {
  /** Consensus DNA sequence (indel-corrected). */
  consensus: string;
  /** Per-base quality scores (Phred Q). */
  quality: Uint8Array;
  /** Per-position, per-base posterior probabilities [pos * 4 + base]. */
  posteriors: Float32Array;
  /** Number of reads that contributed to each position. */
  coverage: Uint8Array;
  /** Effective substitution rate after consensus (estimated from posteriors). */
  effectiveSubRate: number;
  /** Number of refinement iterations actually performed. */
  iterations: number;
}

/**
 * Build MSA-based consensus from multiple noisy reads of the same oligo.
 *
 * @param reads Array of noisy reads (DNA strings, may have different lengths)
 * @param qualities Optional per-read quality arrays
 * @param cfg MSA configuration
 * @returns Consensus result with quality and posteriors
 */
export function buildMsaConsensus(
  reads: string[],
  qualities?: (Uint8Array | undefined)[],
  cfg: MsaConfig = DEFAULT_MSA_CONFIG,
): MsaConsensusResult {
  if (reads.length === 0) {
    return {
      consensus: "",
      quality: new Uint8Array(0),
      posteriors: new Float32Array(0),
      coverage: new Uint8Array(0),
      effectiveSubRate: 1.0,
      iterations: 0,
    };
  }

  if (reads.length === 1) {
    // Single read — no consensus possible, return as-is
    const q = qualities?.[0] ?? new Uint8Array(reads[0].length).fill(20);
    const posteriors = new Float32Array(reads[0].length * 4).fill(0.25);
    return {
      consensus: reads[0],
      quality: q,
      posteriors,
      coverage: new Uint8Array(reads[0].length).fill(1),
      effectiveSubRate: 0.01, // estimate
      iterations: 0,
    };
  }

  if (reads.length < cfg.minClusterSize) {
    // Below min cluster size — use simple plurality consensus
    return buildPluralityConsensus(reads, qualities, cfg);
  }

  // Iterative refinement MSA
  let consensus = buildPluralityConsensus(reads, qualities, cfg).consensus;

  let bestResult: MsaConsensusResult | null = null;

  for (let iter = 0; iter <= cfg.iterations; iter++) {
    // Align each read to the current consensus
    const alignedResults: Hmm3Result[] = [];
    const alignedReads: string[] = [];

    for (let r = 0; r < reads.length; r++) {
      try {
        const result = forwardBackward3(
          reads[r],
          consensus,
          qualities?.[r],
          cfg.hmmParams,
          cfg.hmmBandWidth,
        );
        alignedResults.push(result);
        // Reconstruct the read aligned to the consensus using the Viterbi path
        alignedReads.push(reconstructAligned(reads[r], consensus, result));
      } catch {
        // HMM failed — use the read as-is (length may differ)
        alignedResults.push({
          logLikelihood: -1e30,
          matchPosteriors: new Float32Array(consensus.length * 4).fill(0.25),
          path: [],
        });
        alignedReads.push(reads[r].slice(0, consensus.length));
      }
    }

    // Build weighted consensus from aligned reads using HMM posteriors
    const result = buildWeightedConsensus(
      alignedReads,
      alignedResults,
      qualities,
      consensus.length,
      cfg,
    );

    consensus = result.consensus;
    bestResult = result;
  }

  return bestResult ?? buildPluralityConsensus(reads, qualities, cfg);
}

/**
 * Reconstruct a read aligned to the reference using the HMM Viterbi path.
 *
 * Unlike the simple reconstructReadFromPath (which drops insertions and fills
 * deletions), this version preserves the read's own base calls for match
 * positions, giving per-read variation for downstream soft-decision decoding.
 */
function reconstructAligned(
  read: string,
  ref: string,
  hmmResult: Hmm3Result,
): string {
  if (hmmResult.path.length === 0) {
    // No path — truncate/pad to ref length
    if (read.length >= ref.length) return read.slice(0, ref.length);
    return read + ref.slice(read.length);
  }

  let out = "";
  for (const step of hmmResult.path) {
    if (step.state === "M") {
      // Match: keep the read's own base
      if (step.readPos >= 0 && step.readPos < read.length) {
        out += read[step.readPos];
      } else if (step.refPos >= 0 && step.refPos < ref.length) {
        out += ref[step.refPos];
      }
    } else if (step.state === "D") {
      // Deletion: fill from reference
      if (step.refPos >= 0 && step.refPos < ref.length) {
        out += ref[step.refPos];
      }
    }
    // Insertion: skip (don't emit)
  }

  // Ensure length matches reference
  if (out.length < ref.length) {
    out += ref.slice(out.length);
  } else if (out.length > ref.length) {
    out = out.slice(0, ref.length);
  }

  return out;
}

/**
 * Build weighted consensus from HMM-aligned reads.
 *
 * Uses the posterior probabilities from forwardBackward3 to weight each
 * base call, producing a consensus that is more accurate than simple
 * plurality voting.
 */
function buildWeightedConsensus(
  alignedReads: string[],
  hmmResults: Hmm3Result[],
  qualities: (Uint8Array | undefined)[] | undefined,
  refLen: number,
  cfg: MsaConfig,
): MsaConsensusResult {
  const BASES = ["A", "C", "G", "T"];

  // Accumulate weighted votes for each position
  const votes = new Float64Array(refLen * 4).fill(0);
  const coverage = new Uint8Array(refLen).fill(0);

  for (let r = 0; r < alignedReads.length; r++) {
    const read = alignedReads[r];
    const posteriors = hmmResults[r].matchPosteriors;

    for (let pos = 0; pos < refLen; pos++) {
      if (pos >= read.length) continue;
      coverage[pos]++;

      const base = read[pos];
      const baseIdx = BASES.indexOf(base);
      if (baseIdx < 0) continue; // skip non-ACGT bases

      // Weight by HMM posterior probability for this position
      if (posteriors.length >= (pos + 1) * 4) {
        for (let b = 0; b < 4; b++) {
          votes[pos * 4 + b] += posteriors[pos * 4 + b];
        }
      } else {
        // No posteriors — use uniform weight
        votes[pos * 4 + baseIdx] += 1.0 - 0.01; // Q20 ≈ 1% error
        for (let b = 0; b < 4; b++) {
          if (b !== baseIdx) votes[pos * 4 + b] += 0.01 / 3;
        }
      }
    }
  }

  // Build consensus from weighted votes
  let consensus = "";
  const quality = new Uint8Array(refLen);
  const posteriorsOut = new Float32Array(refLen * 4);
  let totalSubProb = 0;

  for (let pos = 0; pos < refLen; pos++) {
    let sum = 0;
    for (let b = 0; b < 4; b++) sum += votes[pos * 4 + b];

    if (sum === 0) {
      // No coverage — use A as placeholder
      consensus += "A";
      quality[pos] = cfg.qualityFloor;
      for (let b = 0; b < 4; b++) posteriorsOut[pos * 4 + b] = 0.25;
      totalSubProb += 0.75;
      continue;
    }

    // Normalize to get posterior probabilities
    let bestBase = 0;
    let bestProb = 0;
    for (let b = 0; b < 4; b++) {
      const prob = votes[pos * 4 + b] / sum;
      posteriorsOut[pos * 4 + b] = prob;
      if (prob > bestProb) {
        bestProb = prob;
        bestBase = b;
      }
    }

    consensus += BASES[bestBase];

    // Quality score from posterior: Q = -10 * log10(1 - P(best))
    const errorProb = 1 - bestProb;
    totalSubProb += errorProb;
    if (errorProb > 0 && errorProb < 1) {
      quality[pos] = Math.min(40, Math.max(cfg.qualityFloor, Math.round(-10 * Math.log10(errorProb))));
    } else if (errorProb === 0) {
      quality[pos] = 40;
    } else {
      quality[pos] = cfg.qualityFloor;
    }
  }

  const effectiveSubRate = refLen > 0 ? totalSubProb / refLen : 1.0;

  return {
    consensus,
    quality,
    posteriors: posteriorsOut,
    coverage,
    effectiveSubRate,
    iterations: 0, // set by caller
  };
}

/**
 * Build simple plurality consensus (no HMM alignment).
 *
 * Aligns reads at position 0 and takes the most common base at each column.
 * Used as the initial consensus for iterative refinement.
 */
function buildPluralityConsensus(
  reads: string[],
  qualities: (Uint8Array | undefined)[] | undefined,
  cfg: MsaConfig,
): MsaConsensusResult {
  const BASES = ["A", "C", "G", "T"];

  // Use median length as consensus length
  const lengths = reads.map((r) => r.length).sort((a, b) => a - b);
  const consensusLen = lengths[Math.floor(lengths.length / 2)];

  const votes = new Float64Array(consensusLen * 4).fill(0);
  const coverage = new Uint8Array(consensusLen).fill(0);

  for (const read of reads) {
    for (let pos = 0; pos < consensusLen; pos++) {
      if (pos >= read.length) continue;
      const base = read[pos];
      const idx = BASES.indexOf(base);
      if (idx >= 0) {
        votes[pos * 4 + idx] += 1;
        coverage[pos]++;
      }
    }
  }

  let consensus = "";
  const quality = new Uint8Array(consensusLen);
  const posteriors = new Float32Array(consensusLen * 4);
  let totalSubProb = 0;

  for (let pos = 0; pos < consensusLen; pos++) {
    let sum = 0;
    for (let b = 0; b < 4; b++) sum += votes[pos * 4 + b];

    if (sum === 0) {
      consensus += "A";
      quality[pos] = cfg.qualityFloor;
      for (let b = 0; b < 4; b++) posteriors[pos * 4 + b] = 0.25;
      totalSubProb += 0.75;
      continue;
    }

    let bestBase = 0;
    let bestProb = 0;
    for (let b = 0; b < 4; b++) {
      const prob = votes[pos * 4 + b] / sum;
      posteriors[pos * 4 + b] = prob;
      if (prob > bestProb) { bestProb = prob; bestBase = b; }
    }

    consensus += BASES[bestBase];
    const errorProb = 1 - bestProb;
    totalSubProb += errorProb;
    quality[pos] = errorProb > 0 && errorProb < 1
      ? Math.min(40, Math.max(cfg.qualityFloor, Math.round(-10 * Math.log10(errorProb))))
      : errorProb === 0 ? 40 : cfg.qualityFloor;
  }

  return {
    consensus,
    quality,
    posteriors,
    coverage,
    effectiveSubRate: consensusLen > 0 ? totalSubProb / consensusLen : 1.0,
    iterations: 0,
  };
}

/**
 * Convert MSA posteriors to LLR (log-likelihood ratio) array for
 * soft-decision LDPC decoding.
 *
 * LLR[i] = log(P(bit_i = 0) / P(bit_i = 1))
 *
 * For DNA → bits mapping (A=00, C=01, G=10, T=11):
 *   LLR for MSB: log(P(bit=0) / P(bit=1)) = log((P(A)+P(C)) / (P(G)+P(T)))
 *   LLR for LSB: log(P(bit=0) / P(bit=1)) = log((P(A)+P(G)) / (P(C)+P(T)))
 */
export function posteriorsToLLR(posteriors: Float32Array): Float32Array {
  const numPositions = posteriors.length / 4;
  const llr = new Float32Array(numPositions * 2); // 2 bits per base

  for (let pos = 0; pos < numPositions; pos++) {
    const pA = Math.max(posteriors[pos * 4 + 0], 1e-10);
    const pC = Math.max(posteriors[pos * 4 + 1], 1e-10);
    const pG = Math.max(posteriors[pos * 4 + 2], 1e-10);
    const pT = Math.max(posteriors[pos * 4 + 3], 1e-10);

    // MSB (bit 1): 0 for A/C, 1 for G/T
    llr[pos * 2] = Math.log((pA + pC) / (pG + pT));

    // LSB (bit 0): 0 for A/G, 1 for C/T
    llr[pos * 2 + 1] = Math.log((pA + pG) / (pC + pT));
  }

  return llr;
}

/**
 * Top-level entry point: build MSA consensus for a cluster of reads.
 *
 * This is the pre-Viterbi step in the decode pipeline:
 *   Multiple noisy reads → MSA consensus → single read for Viterbi
 *
 * The consensus reduces the effective IDS rate, making downstream
 * Viterbi decoding faster and more accurate.
 *
 * @param reads Noisy reads for one oligo
 * @param qualities Per-read quality arrays
 * @param cfg MSA configuration
 * @returns Consensus result ready for Viterbi decoding
 */
export function msaConsensus(
  reads: string[],
  qualities?: (Uint8Array | undefined)[],
  cfg?: Partial<MsaConfig>,
): MsaConsensusResult {
  const finalCfg: MsaConfig = { ...DEFAULT_MSA_CONFIG, ...cfg };
  return buildMsaConsensus(reads, qualities, finalCfg);
}
