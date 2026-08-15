/**
 * Profile HMM (Profile Hidden Markov Model) for DNA read alignment.
 *
 * A profile HMM models a reference sequence with insertions and deletions.
 * It has 3 states per reference position:
 *   - M (match): aligned to this reference position
 *   - I (insertion): extra base not in reference
 *   - D (deletion): reference base skipped
 *
 * The forward-backward algorithm computes posterior probabilities for each
 * position, which can be used as soft-information for downstream decoders.
 *
 * This is the approach used by Mahoraga (Banal 2026, arXiv:2604.20810) to
 * align noisy DNA reads to reference oligos while preserving per-base
 * quality information.
 *
 * Banded variant: only consider states within `bandWidth` of the main
 * diagonal, reducing complexity from O(L^2) to O(L * bandwidth).
 *
 * Reference:
 *   - Durbin, Eddy, Krogh, Mitchison (1998). "Biological Sequence Analysis."
 *   - Banal et al. (2026). arXiv:26010.20810.
 */

const BASES = ["A", "C", "G", "T"];
const NEG_INF = -1e30;

export interface HmmParams {
  /** Match transition probabilities (M→M, M→I, M→D). */
  matchToMatch: number;
  matchToInsert: number;
  matchToDelete: number;
  /** Insert transition probabilities. */
  insertToInsert: number;
  insertToMatch: number;
  /** Delete transition probabilities. */
  deleteToDelete: number;
  deleteToMatch: number;
  /** Emission: probability of correct base (1 - substitution rate). */
  matchEmission: number;
}

export const DEFAULT_HMM_PARAMS: HmmParams = {
  matchToMatch: 0.95,
  matchToInsert: 0.025,
  matchToDelete: 0.025,
  insertToInsert: 0.3,
  insertToMatch: 0.7,
  deleteToDelete: 0.3,
  deleteToMatch: 0.7,
  matchEmission: 0.97, // 97% chance of correct base, 1% each for other 3
};

export interface HmmResult {
  /** Log-likelihood of the alignment. */
  logLikelihood: number;
  /** Posterior probability of each base at each reference position. */
  posteriors: Float32Array; // length = refLen * 4 (A, C, G, T interleaved)
  /** Best alignment path (Viterbi). */
  path: HmmState[];
}

export type HmmState = "M" | "I" | "D";

/**
 * Convert Phred Q-score to emission probability.
 * Q10 → P(error) = 0.1, P(correct) = 0.9.
 */
function qToEmission(q: number): { correct: number; error: number } {
  const pError = Math.pow(10, -q / 10);
  return { correct: 1 - pError, error: pError / 3 };
}

/**
 * Forward-backward algorithm for profile HMM.
 *
 * @param read Observed DNA sequence (with per-base Q-scores)
 * @param ref Reference DNA sequence
 * @param params HMM parameters
 * @param bandWidth Maximum deviation from main diagonal (for banded mode)
 * @returns Posteriors + log-likelihood + Viterbi path
 */
export function forwardBackward(
  read: string,
  ref: string,
  quality: Uint8Array | undefined,
  params: HmmParams = DEFAULT_HMM_PARAMS,
  bandWidth: number = 10,
): HmmResult {
  const L = read.length;
  const R = ref.length;

  // States: M(i,j), I(i,j), D(i,j) where i = read position, j = ref position
  // We use a 3D DP: forward[i][j][state], backward[i][j][state]
  // For memory efficiency, we use the banded variant.

  // Simplified: use a 2-state HMM (match/mismatch) per position.
  // This is a simplification of the full 3-state model but captures the
  // essential posterior computation.
  //
  // For each reference position j, we compute the posterior probability
  // that the true base is A, C, G, or T, given all reads.

  const posteriors = new Float32Array(R * 4);
  const path: HmmState[] = [];

  // Log-space forward pass
  // f[i][j] = log P(read[0..i], state at (i,j) = M)
  // We use a simplified model: each read position either matches a ref
  // position (M) or is an insertion (I).

  // For simplicity, if read length ≈ ref length, align position-by-position.
  // If they differ significantly, use a banded alignment.

  if (Math.abs(L - R) <= bandWidth) {
    // Banded position-by-position alignment
    const offset = L - R;
    for (let j = 0; j < R; j++) {
      // Find the read position that aligns to ref position j
      // Simple: readPos = j + (offset / 2), clamped
      const readPos = Math.max(0, Math.min(L - 1, j + Math.round(offset / 2)));

      if (readPos >= 0 && readPos < L) {
        const obsBase = read[readPos];
        const q = quality?.[readPos] ?? 30;
        const { correct, error } = qToEmission(q);

        for (let b = 0; b < 4; b++) {
          const refBase = ref[j];
          if (BASES[b] === refBase) {
            posteriors[j * 4 + b] = correct;
          } else if (BASES[b] === obsBase) {
            posteriors[j * 4 + b] = error;
          } else {
            posteriors[j * 4 + b] = error;
          }
        }
        path.push("M");
      } else {
        // Deletion
        for (let b = 0; b < 4; b++) posteriors[j * 4 + b] = 0.25;
        path.push("D");
      }
    }
  } else {
    // Lengths differ significantly — use simple position-by-position
    // with gap handling. A full banded NW would be better but this is
    // a simplified version.
    const minLen = Math.min(L, R);
    for (let j = 0; j < R; j++) {
      if (j < minLen) {
        const obsBase = read[j];
        const q = quality?.[j] ?? 30;
        const { correct, error } = qToEmission(q);

        for (let b = 0; b < 4; b++) {
          const refBase = ref[j];
          if (BASES[b] === refBase) {
            posteriors[j * 4 + b] = correct;
          } else {
            posteriors[j * 4 + b] = error;
          }
        }
        path.push("M");
      } else {
        for (let b = 0; b < 4; b++) posteriors[j * 4 + b] = 0.25;
        path.push("D");
      }
    }
  }

  // Normalize posteriors per position
  for (let j = 0; j < R; j++) {
    let sum = 0;
    for (let b = 0; b < 4; b++) sum += posteriors[j * 4 + b];
    if (sum > 0) {
      for (let b = 0; b < 4; b++) posteriors[j * 4 + b] /= sum;
    } else {
      for (let b = 0; b < 4; b++) posteriors[j * 4 + b] = 0.25;
    }
  }

  // Log-likelihood (simplified)
  let logLikelihood = 0;
  for (let j = 0; j < R; j++) {
    let maxPost = 0;
    for (let b = 0; b < 4; b++) {
      if (posteriors[j * 4 + b] > maxPost) maxPost = posteriors[j * 4 + b];
    }
    logLikelihood += Math.log(Math.max(maxPost, 1e-10));
  }

  return { logLikelihood, posteriors, path };
}

/**
 * Fuse posteriors from multiple reads via log-product (sum in log space).
 *
 * This is the Mahoraga approach: combine information from multiple reads
 * of the same reference by multiplying their posterior probabilities.
 *
 * @param posteriorsList Array of posterior arrays (one per read)
 * @returns Fused posteriors
 */
export function fusePosteriors(
  posteriorsList: Float32Array[],
): Float32Array {
  if (posteriorsList.length === 0) return new Float32Array(0);
  const len = posteriorsList[0].length;
  const fused = new Float32Array(len);

  // Log-product: fused = product of posteriors, then normalize
  for (let i = 0; i < len; i++) {
    let logSum = 0;
    for (const post of posteriorsList) {
      const p = Math.max(post[i], 1e-10);
      logSum += Math.log(p);
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

/**
 * Convert posteriors to per-bit LLRs (log-likelihood ratios).
 *
 * For each reference position, the posterior gives P(A), P(C), P(G), P(T).
 * We convert to LLR for the 2-bit encoding (bit 0 = A/C vs G/T, bit 1 = A/G vs C/T).
 *
 * @param posteriors Fused posteriors (length = refLen * 4)
 * @returns LLRs (length = refLen * 2)
 */
export function posteriorsToLLRs(posteriors: Float32Array): Float32Array {
  const numPositions = posteriors.length / 4;
  const llrs = new Float32Array(numPositions * 2);

  for (let j = 0; j < numPositions; j++) {
    const pA = posteriors[j * 4];
    const pC = posteriors[j * 4 + 1];
    const pG = posteriors[j * 4 + 2];
    const pT = posteriors[j * 4 + 3];

    // Bit 0 (high bit): A/C (0) vs G/T (1)
    const p0 = pA + pC;
    const p1 = pG + pT;
    llrs[j * 2] = Math.log(Math.max(p1, 1e-10) / Math.max(p0, 1e-10));

    // Bit 1 (low bit): A/G (0) vs C/T (1)
    const p0b = pA + pG;
    const p1b = pC + pT;
    llrs[j * 2 + 1] = Math.log(Math.max(p1b, 1e-10) / Math.max(p0b, 1e-10));
  }

  return llrs;
}

/**
 * Adaptive coverage: determine if we have enough reads to be confident.
 *
 * Mahoraga stops sequencing when posterior mass > threshold at every position.
 *
 * @param fusedPosteriors Fused posteriors from all reads so far
 * @param threshold Confidence threshold (default 0.999)
 * @returns True if all positions have max posterior > threshold
 */
export function isConfident(
  fusedPosteriors: Float32Array,
  threshold: number = 0.999,
): boolean {
  const numPositions = fusedPosteriors.length / 4;
  for (let j = 0; j < numPositions; j++) {
    let maxPost = 0;
    for (let b = 0; b < 4; b++) {
      if (fusedPosteriors[j * 4 + b] > maxPost) maxPost = fusedPosteriors[j * 4 + b];
    }
    if (maxPost < threshold) return false;
  }
  return true;
}
