/**
 * Bayesian Quality-Weighted Consensus
 *
 * Replaces simple plurality voting with a proper Bayesian posterior:
 *   P(base | reads) ∝ Π_read P(read_i | base) · Q_i
 *
 * Each read contributes a likelihood weighted by its Q-score. The consensus
 * is the argmax base per position. This extends the HMM soft-info to the
 * consensus layer — the natural complement to the 3-state HMM.
 *
 * Algorithm:
 *   For each position j:
 *     For each candidate base b ∈ {A,C,G,T}:
 *       posterior[b] = prior[b] * Π_read P(read_i[j] | b, Q_i[j])
 *     consensus[j] = argmax_b posterior[b]
 *     confidence[j] = max_b posterior[b] / sum_b posterior[b]
 *
 * This is the "full Bayesian" version of what Mahoraga does at the read level
 * (log-product fusion) — extended to the consensus layer.
 *
 * Reference:
 *   - Durbin et al. (1998). Biological Sequence Analysis, Ch. 4-5.
 *   - Banal et al. (2026). arXiv:2604.20810. (log-product fusion)
 *   - Li (2011). "Statistical modeling of sequencing errors." NAR.
 */

const BASES = ["A", "C", "G", "T"];

export interface BayesianConsensusResult {
  /** Consensus DNA sequence. */
  sequence: string;
  /** Per-position confidence (0..1). */
  confidence: Float32Array;
  /** Per-position posterior probabilities for each base. */
  posteriors: Float32Array; // length = seqLen * 4
  /** Positions marked as erasures (confidence < threshold). */
  erasurePositions: number[];
  /** Number of reads used. */
  readCount: number;
}

export interface BayesianConsensusConfig {
  /** Confidence threshold below which a position is marked as erasure. */
  erasureThreshold: number;
  /** Prior distribution over bases (default: uniform). */
  prior: [number, number, number, number]; // P(A), P(C), P(G), P(T)
  /** Whether to use Q-score weighting. */
  useQualityWeights: boolean;
}

export const DEFAULT_BAYESIAN_CONFIG: BayesianConsensusConfig = {
  erasureThreshold: 0.7,
  prior: [0.25, 0.25, 0.25, 0.25],
  useQualityWeights: true,
};

/**
 * Compute a Bayesian quality-weighted consensus from multiple reads.
 *
 * @param reads Array of reads (DNA strings)
 * @param qualities Array of Q-score arrays (one per read)
 * @param config Consensus configuration
 * @returns Consensus + posteriors + erasure positions
 */
export function bayesianConsensus(
  reads: string[],
  qualities: (Uint8Array | undefined)[] | undefined,
  config: BayesianConsensusConfig = DEFAULT_BAYESIAN_CONFIG,
): BayesianConsensusResult {
  if (reads.length === 0) {
    return {
      sequence: "",
      confidence: new Float32Array(0),
      posteriors: new Float32Array(0),
      erasurePositions: [],
      readCount: 0,
    };
  }

  // Use median length as consensus length
  const lengths = reads.map((r) => r.length);
  const sorted = lengths.slice().sort((a, b) => a - b);
  const consensusLen = sorted[Math.floor(sorted.length / 2)];

  const posteriors = new Float32Array(consensusLen * 4);
  const confidence = new Float32Array(consensusLen);
  const erasurePositions: number[] = [];
  let consensus = "";

  for (let j = 0; j < consensusLen; j++) {
    // Start with prior (in log-space for numerical stability)
    const logPost = [Math.log(config.prior[0]), Math.log(config.prior[1]),
                     Math.log(config.prior[2]), Math.log(config.prior[3])];

    let readCount = 0;
    for (let r = 0; r < reads.length; r++) {
      if (j >= reads[r].length) continue;
      const obs = reads[r][j];
      const obsIdx = BASES.indexOf(obs);
      if (obsIdx === -1) continue; // skip invalid bases

      readCount++;
      const q = qualities?.[r]?.[j] ?? 30;
      const pError = Math.pow(10, -q / 10);
      const pCorrect = 1 - pError;

      // P(obs | true_base) for each candidate base
      // If true_base == obs: pCorrect
      // If true_base != obs: pError / 3
      for (let b = 0; b < 4; b++) {
        const likelihood = b === obsIdx ? pCorrect : pError / 3;
        logPost[b] += Math.log(likelihood);
      }
    }

    // Normalize: convert log-posteriors to probabilities
    const maxLog = Math.max(...logPost);
    const expPost = logPost.map((lp) => Math.exp(lp - maxLog));
    const sum = expPost.reduce((s, v) => s + v, 0);
    const probs = expPost.map((v) => v / sum);

    // Store posteriors
    for (let b = 0; b < 4; b++) {
      posteriors[j * 4 + b] = probs[b];
    }

    // Pick best base
    let bestBase = 0;
    let bestProb = probs[0];
    for (let b = 1; b < 4; b++) {
      if (probs[b] > bestProb) {
        bestProb = probs[b];
        bestBase = b;
      }
    }
    consensus += BASES[bestBase];
    confidence[j] = bestProb;

    // Mark as erasure if confidence below threshold
    if (bestProb < config.erasureThreshold) {
      erasurePositions.push(j);
    }
  }

  return {
    sequence: consensus,
    confidence,
    posteriors,
    erasurePositions,
    readCount: reads.length,
  };
}

/**
 * Fuse posteriors from multiple reads using log-product.
 * This is the Mahoraga approach: P(base | all reads) ∝ Π P(base | read_i).
 */
export function fusePosteriorsBayesian(
  posteriorsList: Float32Array[],
): Float32Array {
  if (posteriorsList.length === 0) return new Float32Array(0);
  const len = posteriorsList[0].length;
  const fused = new Float32Array(len);

  for (let i = 0; i < len; i++) {
    let logSum = 0;
    for (const post of posteriorsList) {
      logSum += Math.log(Math.max(post[i], 1e-10));
    }
    fused[i] = Math.exp(logSum);
  }

  // Normalize per position (every 4 values)
  for (let j = 0; j < len; j += 4) {
    let sum = 0;
    for (let b = 0; b < 4; b++) sum += fused[j + b];
    if (sum > 0) {
      for (let b = 0; b < 4; b++) fused[j + b] /= sum;
    } else {
      for (let b = 0; b < 4; b++) fused[j + b] = 0.25;
    }
  }

  return fused;
}
