/**
 * Progressive Multiple Sequence Alignment (Kalign-style)
 *
 * Implements the progressive alignment approach used by Kalign, Clustal, and
 * MAFFT:
 *   1. Compute pairwise distances between all sequences (k-mer based, fast)
 *   2. Build a guide tree (UPGMA clustering)
 *   3. Progressively align sequences following the guide tree (profile-profile
 *      alignment via Needleman-Wunsch)
 *
 * This produces much more accurate alignments than the center-star method
 * (which we had before), especially for sequences with many indels.
 *
 * Complexity: O(N² L) for distance matrix, O(N log N) for tree, O(N L²) for
 * progressive alignment.
 *
 * Reference:
 *   - Lassmann & Sonnhammer (2005). "Kalign — an accurate and fast multiple
 *     sequence alignment algorithm." BMC Bioinformatics 6:298.
 *   - Feng & Doolittle (1987). "Progressive sequence alignment."
 *   - Thompson, Higgins, Gibson (1994). "CLUSTAL W." NAR 22:22.
 */

const GAP = "-";

export interface ProgressiveMSAConfig {
  /** K-mer size for distance computation. */
  kmerSize: number;
  /** Gap open penalty. */
  gapOpen: number;
  /** Gap extend penalty. */
  gapExtend: number;
  /** Match score. */
  matchScore: number;
  /** Mismatch penalty. */
  mismatchPenalty: number;
}

export const DEFAULT_MSA_CONFIG: ProgressiveMSAConfig = {
  kmerSize: 6,
  gapOpen: -10,
  gapExtend: -1,
  matchScore: 2,
  mismatchPenalty: -1,
};

interface GuideTreeNode {
  left: GuideTreeNode | null;
  right: GuideTreeNode | null;
  sequences: number[]; // indices of sequences in this subtree
  height: number;
}

/**
 * Compute k-mer distance between two sequences.
 * Distance = 1 - (shared k-mers / total unique k-mers)
 */
function kmerDistance(seq1: string, seq2: string, k: number): number {
  const kmers1 = new Set<string>();
  const kmers2 = new Set<string>();
  for (let i = 0; i <= seq1.length - k; i++) kmers1.add(seq1.slice(i, i + k));
  for (let i = 0; i <= seq2.length - k; i++) kmers2.add(seq2.slice(i, i + k));

  let shared = 0;
  for (const kmer of kmers1) if (kmers2.has(kmer)) shared++;

  const total = kmers1.size + kmers2.size - shared;
  return total === 0 ? 1 : 1 - shared / total;
}

/**
 * Build a distance matrix using k-mer distances.
 */
function buildDistanceMatrix(sequences: string[], k: number): number[][] {
  const n = sequences.length;
  const matrix: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const dist = kmerDistance(sequences[i], sequences[j], k);
      matrix[i][j] = dist;
      matrix[j][i] = dist;
    }
  }
  return matrix;
}

/**
 * Build a guide tree using UPGMA (Unweighted Pair Group Method with Arithmetic mean).
 */
function buildGuideTree(distances: number[][]): GuideTreeNode {
  const n = distances.length;
  const nodes: GuideTreeNode[] = [];
  const activeIndices: number[] = [];

  // Initialize: each sequence is a leaf
  for (let i = 0; i < n; i++) {
    nodes.push({ left: null, right: null, sequences: [i], height: 0 });
    activeIndices.push(i);
  }

  // Iteratively merge closest pairs
  while (activeIndices.length > 1) {
    let minDist = Infinity;
    let minI = -1, minJ = -1;
    for (let i = 0; i < activeIndices.length; i++) {
      for (let j = i + 1; j < activeIndices.length; j++) {
        const idx1 = activeIndices[i];
        const idx2 = activeIndices[j];
        // Average distance between clusters
        let totalDist = 0;
        let count = 0;
        for (const s1 of nodes[idx1].sequences) {
          for (const s2 of nodes[idx2].sequences) {
            totalDist += distances[s1][s2];
            count++;
          }
        }
        const avgDist = totalDist / count;
        if (avgDist < minDist) {
          minDist = avgDist;
          minI = idx1;
          minJ = idx2;
        }
      }
    }

    // Merge minI and minJ
    const merged: GuideTreeNode = {
      left: nodes[minI],
      right: nodes[minJ],
      sequences: [...nodes[minI].sequences, ...nodes[minJ].sequences],
      height: minDist / 2,
    };
    nodes.push(merged);
    activeIndices.splice(activeIndices.indexOf(minI), 1);
    activeIndices.splice(activeIndices.indexOf(minJ), 1);
    activeIndices.push(nodes.length - 1);
  }

  return nodes[activeIndices[0]];
}

/**
 * Needleman-Wunsch alignment of two profiles (multiple sequences).
 * Each profile is an array of aligned sequences (same length, with gaps).
 */
function alignProfiles(
  profile1: string[],
  profile2: string[],
  config: ProgressiveMSAConfig,
): { aligned1: string[]; aligned2: string[] } {
  const len1 = profile1[0]?.length ?? 0;
  const len2 = profile2[0]?.length ?? 0;

  if (len1 === 0) {
    const gapSeq = GAP.repeat(len2);
    return { aligned1: profile1.map(() => gapSeq), aligned2: profile2 };
  }
  if (len2 === 0) {
    const gapSeq = GAP.repeat(len1);
    return { aligned1: profile1, aligned2: profile2.map(() => gapSeq) };
  }

  // Scoring: average score over all pairs
  const score = (i: number, j: number): number => {
    let total = 0;
    for (const s1 of profile1) {
      for (const s2 of profile2) {
        const c1 = s1[i];
        const c2 = s2[j];
        if (c1 === GAP || c2 === GAP) total += config.gapExtend;
        else if (c1 === c2) total += config.matchScore;
        else total += config.mismatchPenalty;
      }
    }
    return total / (profile1.length * profile2.length);
  };

  // DP matrix
  const dp: number[][] = Array.from({ length: len1 + 1 }, () => new Array(len2 + 1).fill(0));
  const trace: number[][] = Array.from({ length: len1 + 1 }, () => new Array(len2 + 1).fill(0));

  for (let i = 0; i <= len1; i++) { dp[i][0] = i * config.gapOpen; trace[i][0] = 1; }
  for (let j = 0; j <= len2; j++) { dp[0][j] = j * config.gapOpen; trace[0][j] = 2; }
  trace[0][0] = 0;

  for (let i = 1; i <= len1; i++) {
    for (let j = 1; j <= len2; j++) {
      const diag = dp[i - 1][j - 1] + score(i - 1, j - 1);
      const up = dp[i - 1][j] + config.gapOpen;
      const left = dp[i][j - 1] + config.gapOpen;
      if (diag >= up && diag >= left) { dp[i][j] = diag; trace[i][j] = 0; }
      else if (up >= left) { dp[i][j] = up; trace[i][j] = 1; }
      else { dp[i][j] = left; trace[i][j] = 2; }
    }
  }

  // Traceback
  const aligned1: string[] = profile1.map(() => "");
  const aligned2: string[] = profile2.map(() => "");
  let i = len1, j = len2;
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && trace[i][j] === 0) {
      for (let k = 0; k < profile1.length; k++) aligned1[k] = profile1[k][i - 1] + aligned1[k];
      for (let k = 0; k < profile2.length; k++) aligned2[k] = profile2[k][j - 1] + aligned2[k];
      i--; j--;
    } else if (i > 0 && (j === 0 || trace[i][j] === 1)) {
      for (let k = 0; k < profile1.length; k++) aligned1[k] = profile1[k][i - 1] + aligned1[k];
      for (let k = 0; k < profile2.length; k++) aligned2[k] = GAP + aligned2[k];
      i--;
    } else {
      for (let k = 0; k < profile1.length; k++) aligned1[k] = GAP + aligned1[k];
      for (let k = 0; k < profile2.length; k++) aligned2[k] = profile2[k][j - 1] + aligned2[k];
      j--;
    }
  }

  return { aligned1, aligned2 };
}

/**
 * Progressively align sequences following the guide tree.
 */
function progressiveAlign(
  node: GuideTreeNode,
  sequences: string[],
  config: ProgressiveMSAConfig,
): string[] {
  if (node.left === null && node.right === null) {
    // Leaf — return single sequence
    return [sequences[node.sequences[0]]];
  }

  const leftAligned = node.left ? progressiveAlign(node.left, sequences, config) : [];
  const rightAligned = node.right ? progressiveAlign(node.right, sequences, config) : [];

  if (leftAligned.length === 0) return rightAligned;
  if (rightAligned.length === 0) return leftAligned;

  const { aligned1, aligned2 } = alignProfiles(leftAligned, rightAligned, config);
  return [...aligned1, ...aligned2];
}

/**
 * Perform progressive multiple sequence alignment.
 *
 * @param sequences Array of DNA strings to align
 * @param config MSA configuration
 * @returns Aligned sequences (all same length, with gaps)
 */
export function progressiveMSA(
  sequences: string[],
  config: ProgressiveMSAConfig = DEFAULT_MSA_CONFIG,
): string[] {
  if (sequences.length <= 1) return sequences.slice();

  const distances = buildDistanceMatrix(sequences, config.kmerSize);
  const tree = buildGuideTree(distances);
  return progressiveAlign(tree, sequences, config);
}

/**
 * Build a consensus from a progressive MSA.
 * Column-wise plurality vote, skipping gaps.
 */
export function msaConsensus(aligned: string[]): {
  sequence: string;
  confidence: Float32Array;
  coverage: Float32Array;
} {
  if (aligned.length === 0) {
    return { sequence: "", confidence: new Float32Array(0), coverage: new Float32Array(0) };
  }

  const msaLen = aligned[0].length;
  let consensus = "";
  const confidence: number[] = [];
  const coverage: number[] = [];

  for (let col = 0; col < msaLen; col++) {
    const counts: Record<string, number> = { A: 0, C: 0, G: 0, T: 0 };
    let gapCount = 0;
    let total = 0;

    for (const seq of aligned) {
      const base = seq[col];
      if (base === GAP) { gapCount++; continue; }
      if (base in counts) { counts[base]++; total++; }
    }

    if (gapCount > aligned.length / 2) continue; // skip gap-heavy columns

    let bestBase = "A";
    let bestCount = 0;
    for (const b of ["A", "C", "G", "T"]) {
      if (counts[b] > bestCount) { bestCount = counts[b]; bestBase = b; }
    }

    consensus += bestBase;
    confidence.push(total > 0 ? bestCount / total : 0);
    coverage.push(total / aligned.length);
  }

  return {
    sequence: consensus,
    confidence: new Float32Array(confidence),
    coverage: new Float32Array(coverage),
  };
}
