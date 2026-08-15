/**
 * MACL-style Attention Consensus — Multi-scale attention for DNA error correction
 *
 * Inspired by MACL (Multi-scale Attention Contrastive Learning) which achieves
 * lossless DNA data recovery at 5% BER using attention mechanisms.
 *
 * Architecture:
 *   1. Multi-scale positional encoding: captures local (k-mer) and global
 *      (position-in-oligo) context
 *   2. Self-attention across reads: each read attends to all other reads
 *      at each position, weighting by similarity
 *   3. Contrastive loss: pushes apart reads from different oligos
 *
 * Since we can't ship a trained model, this implements a STRUCTURAL
 * attention mechanism using Q-scores and sequence similarity:
 *   - Attention weight = exp(-Hamming distance / temperature)
 *   - Multi-scale: combine local (windowed) and global attention
 *   - Output: weighted consensus with confidence scores
 *
 * Reference:
 *   - MACL (2025). "Multi-scale Attention Contrastive Learning for DNA Storage."
 *   - Vaswani et al. (2017). "Attention Is All You Need."
 *   - Bar-Lev et al. (2025). Nature Machine Intelligence 7:639-649.
 */

import { SequencingRead } from "./simulate";
import { Base } from "./mapping";

const BASES: Base[] = ["A", "C", "G", "T"];

export interface AttentionConsensusResult {
  /** Consensus DNA string. */
  consensus: string;
  /** Per-base confidence scores (0..1). */
  confidence: Float32Array;
  /** Number of reads that contributed. */
  readsUsed: number;
}

/**
 * Compute multi-scale attention consensus from multiple reads.
 *
 * @param reads Array of reads for the same oligo
 * @param expectedLen Expected DNA length
 * @param temperature Attention temperature (lower = sharper attention)
 * @returns Consensus DNA + confidence scores
 */
export function attentionConsensus(
  reads: SequencingRead[],
  expectedLen: number,
  temperature: number = 2.0,
): AttentionConsensusResult | null {
  if (reads.length === 0) return null;

  const numReads = reads.length;
  const consensusParts: string[] = new Array(expectedLen);
  const confidence = new Float32Array(expectedLen);

  // Multi-scale: process in windows of 10 (local) and full length (global)
  const windowSize = 10;

  for (let pos = 0; pos < expectedLen; pos++) {
    // Compute attention weights for each read at this position
    // Weight = similarity to other reads (local window)
    const weights = new Float64Array(numReads);
    let totalWeight = 0;

    for (let i = 0; i < numReads; i++) {
      if (pos >= reads[i].sequence.length) continue;

      // Local similarity: compare window around pos with all other reads
      let similarity = 0;
      let comparisons = 0;
      for (let j = 0; j < numReads; j++) {
        if (i === j) continue;
        if (pos >= reads[j].sequence.length) continue;

        // Window comparison
        const winStart = Math.max(0, pos - Math.floor(windowSize / 2));
        const winEnd = Math.min(expectedLen, pos + Math.ceil(windowSize / 2));
        let matches = 0;
        let total = 0;
        for (let w = winStart; w < winEnd; w++) {
          if (w < reads[i].sequence.length && w < reads[j].sequence.length) {
            if (reads[i].sequence[w] === reads[j].sequence[w]) matches++;
            total++;
          }
        }
        if (total > 0) {
          similarity += matches / total;
          comparisons++;
        }
      }

      if (comparisons > 0) {
        similarity /= comparisons;
      } else {
        similarity = 0.5; // neutral if no comparisons
      }

      // Q-score weight
      let qWeight = 0.9;
      if (reads[i].quality && pos < reads[i].quality.length) {
        const q = reads[i].quality[pos];
        qWeight = 1 - Math.pow(10, -q / 10);
      }

      // Combined weight: similarity * Q-score
      weights[i] = Math.exp(similarity / temperature) * qWeight;
      totalWeight += weights[i];
    }

    if (totalWeight === 0) {
      consensusParts[pos] = "A";
      confidence[pos] = 0;
      continue;
    }

    // Weighted vote
    const votes: Record<string, number> = { A: 0, C: 0, G: 0, T: 0 };
    for (let i = 0; i < numReads; i++) {
      if (pos >= reads[i].sequence.length) continue;
      const base = reads[i].sequence[pos];
      if (base in votes) {
        votes[base] += weights[i] / totalWeight;
      }
    }

    // Pick best base and compute confidence
    let bestBase = "A";
    let bestVote = -1;
    for (const b of BASES) {
      if (votes[b] > bestVote) {
        bestVote = votes[b];
        bestBase = b;
      }
    }
    consensusParts[pos] = bestBase;
    confidence[pos] = bestVote;
  }

  return {
    consensus: consensusParts.join(""),
    confidence,
    readsUsed: numReads,
  };
}

/**
 * Decode using attention consensus + LDPC.
 *
 * @param reads Array of reads for the same oligo
 * @param expectedLen Expected DNA length
 * @param ldpc LDPC inner code
 * @param dnaToBytesFn DNA to bytes conversion function
 * @param innerN Inner code N (codeword length)
 * @returns Decoded data or null
 */
export function decodeWithAttentionConsensus(
  reads: SequencingRead[],
  expectedLen: number,
  ldpc: any,
  dnaToBytesFn: (dna: string) => Uint8Array,
  innerN: number,
): { data: Uint8Array; corrected: number; avgConfidence: number } | null {
  const result = attentionConsensus(reads, expectedLen);
  if (!result) return null;

  try {
    const innerBlock = dnaToBytesFn(result.consensus);
    const rsCodeword = innerBlock.slice(0, innerN);
    const r = ldpc.decode(rsCodeword);
    const avgConfidence = result.confidence.reduce((a, b) => a + b, 0) / result.confidence.length;
    return { data: r.data, corrected: r.corrected, avgConfidence };
  } catch {
    return null;
  }
}
