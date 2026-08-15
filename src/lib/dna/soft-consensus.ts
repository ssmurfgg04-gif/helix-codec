/**
 * Soft-Consensus Decoder — Mahoraga-style LLR Fusion
 *
 * Instead of hard-decision consensus (which loses soft information), this module:
 *   1. Aligns all reads for an oligo using Profile HMM (forward-backward)
 *   2. Fuses posteriors across reads via log-product fusion
 *   3. Converts fused posteriors to per-bit LLRs
 *   4. Feeds LLRs to the BP decoder
 *
 * This is the architecture used by Mahoraga (Banal 2026, arXiv:2604.20810)
 * which achieves 1.42x density improvement over DNA-Aeon by retaining
 * per-position posteriors instead of taking hard consensus.
 *
 * Pipeline:
 *   reads[] → Profile HMM → per-read posteriors → log-product fusion
 *   → per-base posteriors → per-bit LLRs → BP decoder → codeword
 *
 * The key insight: hard consensus throws away ~50% of the information
 * (the confidence level). By preserving LLRs, the BP decoder can make
 * much better decisions, especially at low coverage.
 */

import { forwardBackward3, fusePosteriors3, DEFAULT_HMM3_PARAMS } from "./profileHmm3";
import { Base } from "./mapping";
import { LDPCInnerCode } from "./ldpc-codec";
import { SequencingRead } from "./simulate";

const BASES: Base[] = ["A", "C", "G", "T"];

export interface SoftConsensusResult {
  /** Per-bit LLRs for the LDPC codeword (nBits = n * 8). */
  llrs: Float32Array;
  /** Number of reads successfully aligned. */
  readsAligned: number;
  /** Average log-likelihood of alignment (quality metric). */
  avgLogLikelihood: number;
}

/**
 * Compute per-read posteriors using Profile HMM, then fuse across reads.
 *
 * @param reads Array of reads for the same oligo
 * @param referenceDNA The expected DNA sequence (from primer-trimmed read or consensus)
 * @param nBits Number of bits in the LDPC codeword
 * @param useGoldman Whether Goldman mapping was used (affects bit-to-base mapping)
 * @returns Fused LLRs for the BP decoder
 */
export function softConsensus(
  reads: SequencingRead[],
  referenceDNA: string,
  nBits: number,
  useGoldman: boolean = false,
): SoftConsensusResult | null {
  if (reads.length === 0) return null;

  // Step 1: Compute per-read posteriors using Profile HMM
  const allPosteriors: Float32Array[] = [];
  let totalLogLikelihood = 0;
  let alignedCount = 0;

  for (const read of reads) {
    if (!read.sequence || read.sequence.length === 0) continue;

    try {
      const result = forwardBackward3(
        read.sequence,
        referenceDNA,
        read.quality,
        DEFAULT_HMM3_PARAMS,
        15, // v57: increased bandWidth from 10 to 15 for better alignment at low coverage
      );

      if (result.matchPosteriors.length > 0) {
        allPosteriors.push(result.matchPosteriors);
        totalLogLikelihood += result.logLikelihood;
        alignedCount++;
      }
    } catch {
      // HMM alignment failed for this read — skip it
      continue;
    }
  }

  if (alignedCount === 0) return null;

  // Step 2: Fuse posteriors across reads (log-product fusion)
  // Mahoraga uses: fused_posterior = product of per-read posteriors
  // In log domain: log(fused) = sum of log(posteriors)
  // Then normalize to get a proper probability distribution.
  const refLen = referenceDNA.length;
  const fusedPosteriors = fusePosteriors3(allPosteriors);

  if (fusedPosteriors.length === 0) return null;

  // Step 3: Convert per-base posteriors to per-bit LLRs
  // For 2-bit mapping: base b → bits (b >> 1, b & 1)
  //   bit 0 (MSB): 0 if base is A or C, 1 if base is G or T
  //   bit 1 (LSB): 0 if base is A or G, 1 if base is C or T
  //
  // LLR(bit=0) = log(P(bit=0) / P(bit=1))
  //   = log( (P(A)+P(C)) / (P(G)+P(T)) )  for MSB
  //   = log( (P(A)+P(G)) / (P(C)+P(T)) )  for LSB
  //
  // For Goldman mode: the mapping is more complex (trits), so we use a
  // simplified uniform LLR based on the average Q-score.
  const llrs = new Float32Array(nBits);

  if (useGoldman) {
    // Goldman mode: use average Q-score for uniform LLR
    let avgQ = 30;
    if (reads[0]?.quality && reads[0].quality.length > 0) {
      let qSum = 0, qCount = 0;
      for (const read of reads) {
        if (read.quality) {
          for (let i = 0; i < read.quality.length; i++) {
            qSum += read.quality[i];
            qCount++;
          }
        }
      }
      if (qCount > 0) avgQ = qSum / qCount;
    }
    const llrMag = (avgQ / 10) * Math.log(10);
    for (let i = 0; i < nBits; i++) {
      llrs[i] = llrMag; // positive = likely 0
    }
  } else {
    // Direct 2-bit mode: convert per-base posteriors to per-bit LLRs
    // Each base = 2 bits. Base index: A=0(00), C=1(01), G=2(10), T=3(11)
    // For nBits = n * 8, we need n bytes = n * 4 bases.
    const numBases = nBits / 2;

    for (let bitIdx = 0; bitIdx < nBits; bitIdx++) {
      const baseIdx = bitIdx >> 1; // which base this bit belongs to
      const isMSB = (bitIdx & 1) === 0; // MSB or LSB of the base

      if (baseIdx >= refLen || baseIdx * 4 + 3 >= fusedPosteriors.length) {
        // Out of range — use neutral LLR
        llrs[bitIdx] = 0;
        continue;
      }

      // Get posteriors for this base: P(A), P(C), P(G), P(T)
      // v60: fusePosteriors3 returns LINEAR probabilities (not log).
      // The old code did Math.exp() on them again (double-exp bug) which
      // accidentally sharpened the already-peaked (wrong, max-plus) posteriors.
      // Now that profileHmm3 uses correct sum-product, the posteriors are
      // softer and the double-exp would give wrong LLRs. Use them directly.
      const pA = fusedPosteriors[baseIdx * 4 + 0];
      const pC = fusedPosteriors[baseIdx * 4 + 1];
      const pG = fusedPosteriors[baseIdx * 4 + 2];
      const pT = fusedPosteriors[baseIdx * 4 + 3];

      // Normalize (defensive — fusePosteriors3 already normalizes, but
      // numerical drift can occur)
      const total = pA + pC + pG + pT;
      if (total === 0) {
        llrs[bitIdx] = 0;
        continue;
      }

      if (isMSB) {
        // MSB: 0 = {A, C}, 1 = {G, T}
        const p0 = (pA + pC) / total;
        const p1 = (pG + pT) / total;
        llrs[bitIdx] = Math.log((p0 + 1e-10) / (p1 + 1e-10));
      } else {
        // LSB: 0 = {A, G}, 1 = {C, T}
        const p0 = (pA + pG) / total;
        const p1 = (pC + pT) / total;
        llrs[bitIdx] = Math.log((p0 + 1e-10) / (p1 + 1e-10));
      }
    }
  }

  return {
    llrs,
    readsAligned: alignedCount,
    avgLogLikelihood: totalLogLikelihood / alignedCount,
  };
}

/**
 * Fast weighted consensus — Q-score weighted majority vote.
 *
 * This is a simplified soft-consensus that doesn't use HMM alignment.
 * It aligns reads by position (assuming no indels) and takes a Q-score
 * weighted majority vote at each position. The resulting consensus is
 * then decoded with the LDPC hard-decision decoder.
 *
 * Much faster than HMM-based soft consensus (~100x), works well when
 * indels are rare (Illumina < 0.1% indel rate).
 *
 * @param reads All reads for this oligo
 * @param expectedLen Expected DNA length
 * @returns Consensus DNA string, or null if no reads
 */
export function fastWeightedConsensus(
  reads: SequencingRead[],
  expectedLen: number,
): string | null {
  if (reads.length === 0) return null;

  // Align reads by position (no indel handling — assumes perfect alignment)
  // At each position, take Q-score weighted majority vote.
  const parts: string[] = new Array(expectedLen);
  const BASES = ["A", "C", "G", "T"];

  for (let pos = 0; pos < expectedLen; pos++) {
    // Weighted vote: weight = 10^(-Q/10) for error probability
    // Higher Q → higher weight (more confident)
    const votes: Record<string, number> = { A: 0, C: 0, G: 0, T: 0 };

    for (const read of reads) {
      if (pos >= read.sequence.length) continue;
      const base = read.sequence[pos];
      if (!(base in votes)) continue;

      // Weight = confidence. Q-score → probability of correctness.
      // If Q=30, weight = 1 - 10^(-3) = 0.999
      // If Q=10, weight = 1 - 10^(-1) = 0.9
      let weight = 0.9; // default if no Q-score
      if (read.quality && pos < read.quality.length) {
        const q = read.quality[pos];
        weight = 1 - Math.pow(10, -q / 10);
      }
      votes[base] += weight;
    }

    // Pick base with highest weighted vote
    let bestBase = "A";
    let bestWeight = -1;
    for (const b of BASES) {
      if (votes[b] > bestWeight) {
        bestWeight = votes[b];
        bestBase = b;
      }
    }
    parts[pos] = bestBase;
  }

  return parts.join("");
}

/**
 * Decode an oligo using fast weighted consensus + hard-decision LDPC.
 *
 * This is the fast path for soft-consensus: no HMM, no BP, just
 * Q-score weighted majority vote → LDPC hard-decision decode.
 *
 * @param reads All reads for this oligo
 * @param expectedLen Expected DNA length
 * @param ldpc LDPC inner code instance
 * @param dnaToBytesFn Function to convert DNA to bytes
 * @returns Decoded info bytes, or null if failed
 */
export function decodeWithFastConsensus(
  reads: SequencingRead[],
  expectedLen: number,
  ldpc: LDPCInnerCode,
  dnaToBytesFn: (dna: string) => Uint8Array,
  innerN: number,
): { data: Uint8Array; corrected: number } | null {
  const consensus = fastWeightedConsensus(reads, expectedLen);
  if (!consensus) return null;

  try {
    const innerBlock = dnaToBytesFn(consensus);
    const rsCodeword = innerBlock.slice(0, innerN);
    const r = ldpc.decode(rsCodeword);
    return { data: r.data, corrected: r.corrected };
  } catch {
    return null;
  }
}

/**
 * v58: Decode an oligo using Progressive Multiple Sequence Alignment (Kalign-style)
 * + LDPC hard-decision decode.
 *
 * Unlike `fastWeightedConsensus` (which assumes no indels), progressive MSA
 * aligns reads with Needleman-Wunsch profile-profile alignment, placing gaps
 * where insertions/deletions occurred. The column-wise plurality consensus
 * then produces a clean DNA string even when individual reads have many indels.
 *
 * This is the right tool for Nanopore / high-IDS channels where indels
 * dominate. The HMM soft-consensus (decodeWithSoftConsensus) also handles
 * indels but is ~10× slower and can fail to converge on long noisy reads.
 *
 * Architecture (Kalign / Clustal / MAFFT progressive alignment):
 *   1. Compute pairwise k-mer distances between all reads
 *   2. Build UPGMA guide tree
 *   3. Progressively align following guide tree (profile-profile NW)
 *   4. Column-wise plurality vote → consensus DNA
 *   5. DNA → bytes → LDPC hard-decision decode
 *
 * @param reads All reads for this oligo (typically 5-25 reads at 10-25× coverage)
 * @param expectedLen Expected DNA length (length of inner region between primers)
 * @param ldpc LDPC inner code instance
 * @param dnaToBytesFn Function to convert DNA to bytes
 * @param innerN Number of bytes in the LDPC codeword (data + parity, no CRC)
 * @returns Decoded info bytes, or null if MSA / LDPC failed
 *
 * Reference:
 *   - Lassmann & Sonnhammer (2005). "Kalign — an accurate and fast multiple
 *     sequence alignment algorithm." BMC Bioinformatics 6:298.
 *   - Feng & Doolittle (1987). "Progressive sequence alignment."
 */
export function decodeWithProgressiveMSA(
  reads: SequencingRead[],
  expectedLen: number,
  ldpc: LDPCInnerCode,
  dnaToBytesFn: (dna: string) => Uint8Array,
  innerN: number,
  convInner?: { decode: (data: Uint8Array) => { decoded: Uint8Array; corrected: number } } | null,
  convLayout?: { addressBytes: number; convEncodedBytes: number; crcBytes: number } | null,
): { data: Uint8Array; corrected: number; readsAligned: number } | null {
  if (reads.length < 2) return null;

  // Lazy-load progressive MSA to avoid import cycle for callers that don't use it
  const { progressiveMSA, msaConsensus, DEFAULT_MSA_CONFIG } = require("./progressive-msa");

  // Cap the number of reads to align — MSA is O(N²) in reads, and 10-15
  // well-chosen reads give the same consensus quality as 50. Pick the
  // longest reads (highest information content) and a few medium ones.
  let readsToAlign = reads;
  if (reads.length > 15) {
    // Sort by length descending (longest reads have the most anchor points)
    const sorted = [...reads].sort((a, b) => b.sequence.length - a.sequence.length);
    readsToAlign = sorted.slice(0, 15);
  }

  // Trim/pad each read's sequence to roughly expectedLen ± 10% to keep MSA focused
  // on the inner region (post-primer-trim). The decoder already trims primers
  // upstream, so reads[i].sequence should already be inner DNA.
  const trimmedSeqs: string[] = [];
  for (const r of readsToAlign) {
    let s = r.sequence;
    // Allow ±20% length variation (insertions/deletions can change length)
    const maxLen = Math.floor(expectedLen * 1.2);
    const minLen = Math.floor(expectedLen * 0.8);
    if (s.length > maxLen) s = s.slice(0, maxLen);
    else if (s.length < minLen) s = s + "A".repeat(minLen - s.length);
    trimmedSeqs.push(s);
  }

  try {
    // Run progressive MSA (Kalign-style: k-mer distance → UPGMA tree → profile-profile NW)
    const aligned = progressiveMSA(trimmedSeqs, DEFAULT_MSA_CONFIG);
    if (aligned.length === 0) return null;

    // Column-wise plurality consensus (skips gap-heavy columns)
    const consensusResult = msaConsensus(aligned);
    if (!consensusResult.sequence || consensusResult.sequence.length < expectedLen * 0.7) {
      // MSA produced too short a consensus — likely the reads are too noisy
      return null;
    }

    // Pad or trim consensus to exactly expectedLen so dnaToBytes succeeds
    let consensus = consensusResult.sequence;
    if (consensus.length > expectedLen) {
      consensus = consensus.slice(0, expectedLen);
    } else if (consensus.length < expectedLen) {
      consensus = consensus + "A".repeat(expectedLen - consensus.length);
    }

    // DNA → bytes
    const innerBlock = dnaToBytesFn(consensus);

    // v59: If conv inner code is used, run conv Viterbi decode BEFORE LDPC.
    // The MSA consensus has indels corrected, so conv Viterbi should succeed
    // where per-read conv Viterbi failed (the consensus is much cleaner).
    if (convInner && convLayout) {
      try {
        const convBytes = innerBlock.slice(
          convLayout.addressBytes,
          convLayout.addressBytes + convLayout.convEncodedBytes,
        );
        const convResult = convInner.decode(convBytes);
        const rsCodewordDecoded = convResult.decoded;
        // Rebuild innerBlock as [LDPC codeword + CRC]
        const crcFromRead = innerBlock.slice(innerBlock.length - convLayout.crcBytes);
        const rebuilt = new Uint8Array(innerN + convLayout.crcBytes);
        rebuilt.set(rsCodewordDecoded, 0);
        rebuilt.set(crcFromRead, innerN);
        const rsCodeword = rebuilt.slice(0, innerN);
        const r = ldpc.decode(rsCodeword);
        return {
          data: r.data,
          corrected: r.corrected,
          readsAligned: readsToAlign.length,
        };
      } catch {
        // Conv decode failed even on consensus — fall through to direct LDPC
      }
    }

    // Standard path: LDPC hard-decision decode
    const rsCodeword = innerBlock.slice(0, innerN);
    const r = ldpc.decode(rsCodeword);
    return {
      data: r.data,
      corrected: r.corrected,
      readsAligned: readsToAlign.length,
    };
  } catch {
    return null;
  }
}

/**
 * Decode an oligo using soft consensus + BP decoder.
 *
 * @param reads All reads for this oligo
 * @param referenceDNA Expected DNA (e.g., from first read or consensus)
 * @param ldpc LDPC inner code instance
 * @param useGoldman Whether Goldman mapping was used
 * @returns Decoded info bytes, or null if failed
 */
export function decodeWithSoftConsensus(
  reads: SequencingRead[],
  referenceDNA: string,
  ldpc: LDPCInnerCode,
  useGoldman: boolean = false,
): { data: Uint8Array; corrected: number; readsAligned: number } | null {
  const nBits = ldpc.nBits;
  const result = softConsensus(reads, referenceDNA, nBits, useGoldman);
  if (!result) return null;

  try {
    // Convert LLRs to hard-decision bytes for the BP decoder
    const n = ldpc.n;
    const recv = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      let byte = 0;
      for (let bit = 0; bit < 8; bit++) {
        const bitIdx = i * 8 + bit;
        const hardBit = result.llrs[bitIdx] < 0 ? 1 : 0;
        byte |= hardBit << (7 - bit);
      }
      recv[i] = byte;
    }

    // v57: Convert LLRs to Q-scores and pass to BP decoder.
    // Before v57: BP was called with null Q-scores (hard-decision only),
    // discarding the soft information from HMM fusion. This caused BP to
    // fail at 2× coverage because it couldn't distinguish reliable bits
    // (high |LLR|) from uncertain bits (low |LLR|).
    //
    // LLR → Q-score conversion:
    //   P_error = 1 / (1 + exp(|LLR|))
    //   Q = -10 * log10(P_error) = 10 * log10(1 + exp(|LLR|))
    // For large |LLR|, Q ≈ 4.34 * |LLR|.
    //
    // Cap Q at 40 (matching typical Illumina Q-scores) to avoid numerical
    // issues with very large LLRs.
    const numBases = n * 4; // 4 bases per byte (2 bits per base, 4 bases per byte)
    const qScores = new Uint8Array(numBases);
    for (let byteIdx = 0; byteIdx < n; byteIdx++) {
      // Each byte has 4 bases (2 bits each). Average the 2 bits' LLRs per base.
      for (let baseIdx = 0; baseIdx < 4; baseIdx++) {
        const bit0 = byteIdx * 8 + baseIdx * 2;
        const bit1 = byteIdx * 8 + baseIdx * 2 + 1;
        const llr0 = Math.abs(result.llrs[bit0] || 0);
        const llr1 = Math.abs(result.llrs[bit1] || 0);
        const avgLlr = (llr0 + llr1) / 2;
        // Q = 10 * log10(1 + exp(|LLR|)), capped at 40
        const q = Math.min(40, Math.round(10 * Math.log10(1 + Math.exp(avgLlr))));
        qScores[byteIdx * 4 + baseIdx] = q;
      }
    }

    // Run BP decoder with soft info (v57: pass Q-scores derived from HMM LLRs)
    const bpResult = ldpc.decodeBeliefPropagation(recv, qScores, useGoldman, 50);
    return {
      data: bpResult.data,
      corrected: bpResult.corrected,
      readsAligned: result.readsAligned,
    };
  } catch {
    return null;
  }
}
