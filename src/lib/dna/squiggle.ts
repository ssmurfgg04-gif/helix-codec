/**
 * Nanopore Squiggle-Native Decoding
 *
 * Skips basecalling entirely — decodes DNA storage data directly from raw
 * electrical signals (squiggles) produced by Oxford Nanopore sequencers.
 *
 * Traditional pipeline:  raw signal → basecaller → bases → decoder
 * Squiggle-native:       raw signal → squiggle matcher → decoder
 *
 * This avoids basecalling errors (which dominate Nanopore error rates) and
 * enables near-real-time decoding.
 *
 * Approach:
 *   1. Pre-compute reference squiggles for each possible k-mer (using a
 *      nanopore signal model — ONT's k-mer template table)
 *   2. Slide a window over the raw signal and match to reference squiggles
 *   3. Use dynamic time warping (DTW) for elastic matching
 *   4. Build a consensus from multiple squiggle reads
 *
 * Reference:
 *   - Loose, Malla, Stout (2016). "Real-time selective sequencing with
 *     ReadUntil." Nat Methods 13:751-754.
 *   - Bao et al. (2022). "SquiggleNet: real-time, direct classification
 *     of nanopore signals." Genome Biol 23.
 *   - ONT signal model: github.com/nanoporetech/kmer_models
 */

export interface SquiggleConfig {
  /** K-mer size for signal model (default 5 for R9.4 chemistry). */
  kmerSize: number;
  /** Sampling rate in Hz (default 4000 for Nanopore). */
  sampleRate: number;
  /** Window size for signal matching. */
  windowSize: number;
  /** DTW constraint (maximum warp distance). */
  dtwConstraint: number;
}

export const DEFAULT_SQUIGGLE_CONFIG: SquiggleConfig = {
  kmerSize: 5,
  sampleRate: 4000,
  windowSize: 200,
  dtwConstraint: 10,
};

// Simulated nanopore signal model: k-mer → mean signal level (pA)
// In production, load from ONT's kmer_models file
const KMER_SIGNAL_TABLE: Map<string, number> = new Map();

function initSignalModel(): void {
  if (KMER_SIGNAL_TABLE.size > 0) return;
  const bases = "ACGT";
  // Generate pseudo-random signal levels for each k-mer
  let seed = 42;
  const rng = () => {
    seed = (seed * 1103515245 + 12345) & 0x7fffffff;
    return seed / 0x7fffffff;
  };
  // Generate all k-mers of size 5
  const generateKmers = (k: number, prefix = ""): string[] => {
    if (prefix.length === k) return [prefix];
    const result: string[] = [];
    for (const b of bases) result.push(...generateKmers(k, prefix + b));
    return result;
  };
  for (const kmer of generateKmers(5)) {
    // Signal level: 40-90 pA (typical Nanopore range)
    KMER_SIGNAL_TABLE.set(kmer, 40 + rng() * 50);
  }
}

/**
 * Generate a synthetic squiggle (raw signal) from a DNA sequence.
 * Used for testing the squiggle-native decoder.
 */
export function generateSquiggle(
  sequence: string,
  config: SquiggleConfig = DEFAULT_SQUIGGLE_CONFIG,
  noiseLevel: number = 1.0,
): Float32Array {
  initSignalModel();
  const signals: number[] = [];

  for (let i = 0; i <= sequence.length - config.kmerSize; i++) {
    const kmer = sequence.slice(i, i + config.kmerSize);
    const meanLevel = KMER_SIGNAL_TABLE.get(kmer) ?? 60;

    // Each k-mer produces ~sampleRate/450 samples (450 bases/sec typical)
    const samplesPerKmer = Math.floor(config.sampleRate / 450);
    for (let s = 0; s < samplesPerKmer; s++) {
      // Add Gaussian noise
      const noise = (Math.random() - 0.5) * 2 * noiseLevel;
      signals.push(meanLevel + noise);
    }
  }

  return new Float32Array(signals);
}

/**
 * Dynamic Time Warping (DTW) distance between two signal segments.
 * Used for elastic matching of squiggles to reference signals.
 */
export function dtwDistance(
  signal1: Float32Array,
  signal2: Float32Array,
  constraint: number = 10,
): number {
  const n = signal1.length;
  const m = signal2.length;
  const dp: Float32Array[] = Array.from({ length: n + 1 }, () =>
    new Float32Array(m + 1).fill(Infinity),
  );
  dp[0][0] = 0;

  for (let i = 1; i <= n; i++) {
    const jStart = Math.max(1, i - constraint);
    const jEnd = Math.min(m, i + constraint);
    for (let j = jStart; j <= jEnd; j++) {
      const cost = Math.abs(signal1[i - 1] - signal2[j - 1]);
      dp[i][j] = cost + Math.min(
        dp[i - 1][j],     // insertion
        dp[i][j - 1],     // deletion
        dp[i - 1][j - 1], // match
      );
    }
  }

  return dp[n][m];
}

export interface SquiggleMatch {
  /** Reference sequence position. */
  refPos: number;
  /** Signal segment start. */
  signalStart: number;
  /** Signal segment end. */
  signalEnd: number;
  /** DTW distance (lower = better match). */
  distance: number;
  /** Confidence (0-1, derived from distance). */
  confidence: number;
}

/**
 * Match a raw squiggle against a reference DNA sequence.
 * Returns the best-matching positions with confidence scores.
 *
 * This replaces basecalling + alignment with direct signal-to-reference matching.
 */
export function matchSquiggleToReference(
  squiggle: Float32Array,
  reference: string,
  config: SquiggleConfig = DEFAULT_SQUIGGLE_CONFIG,
): SquiggleMatch[] {
  initSignalModel();
  const matches: SquiggleMatch[] = [];
  const windowSize = config.windowSize;

  for (let sigStart = 0; sigStart < squiggle.length; sigStart += windowSize / 2) {
    const sigEnd = Math.min(sigStart + windowSize, squiggle.length);
    const sigSegment = squiggle.slice(sigStart, sigEnd);

    let bestPos = 0;
    let bestDist = Infinity;

    // Slide over reference and find best match
    for (let refPos = 0; refPos <= reference.length - config.kmerSize; refPos++) {
      const refKmer = reference.slice(refPos, refPos + config.kmerSize);
      const refSignal = generateSquiggle(refKmer, config, 0);
      const dist = dtwDistance(sigSegment, refSignal, config.dtwConstraint);
      if (dist < bestDist) {
        bestDist = dist;
        bestPos = refPos;
      }
    }

    const confidence = Math.exp(-bestDist / 100);
    matches.push({
      refPos: bestPos,
      signalStart: sigStart,
      signalEnd: sigEnd,
      distance: bestDist,
      confidence,
    });
  }

  return matches;
}

/**
 * Build a consensus from multiple squiggle matches.
 * Each match contributes a vote weighted by its confidence.
 */
export function squiggleConsensus(matches: SquiggleMatch[]): {
  positions: number[];
  confidence: Float32Array;
} {
  if (matches.length === 0) {
    return { positions: [], confidence: new Float32Array(0) };
  }

  // Group matches by reference position
  const positionVotes = new Map<number, number>();
  for (const match of matches) {
    const current = positionVotes.get(match.refPos) ?? 0;
    positionVotes.set(match.refPos, current + match.confidence);
  }

  // Sort by position
  const sortedPositions = Array.from(positionVotes.keys()).sort((a, b) => a - b);
  const confidence = new Float32Array(sortedPositions.map((p) => positionVotes.get(p)!));

  // Normalize confidence
  const maxConf = Math.max(...confidence);
  if (maxConf > 0) {
    for (let i = 0; i < confidence.length; i++) confidence[i] /= maxConf;
  }

  return { positions: sortedPositions, confidence };
}
