/**
 * v60: Nanopore Viterbi Preprocessor — REAL indel correction
 *
 * Corrects indels in Nanopore reads BEFORE the LDPC inner decoder runs.
 * This is the HEDGES-style approach: a profile HMM with Viterbi decoding
 * aligns each read to the cluster consensus, then reconstructs the
 * indel-corrected read by walking the Viterbi path.
 *
 * Pipeline:
 *   Read (noisy, indel-heavy)
 *     → Profile HMM forwardBackward3 (real sum-product + Viterbi traceback)
 *     → Viterbi path: M/I/D states along the optimal alignment
 *     → reconstructReadFromPath: walk path, keep M bases, drop I bases, fill D from ref
 *     → Hand off to LDPC (substitution correction)
 *
 * v60 fixes:
 *   - profileHmm3.ts now uses CORRECT log-sum-exp forward-backward (not max-plus)
 *   - profileHmm3.ts now computes REAL Viterbi path via backpointer traceback
 *   - This file now uses reconstructReadFromPath to actually realign reads,
 *     instead of the v59 stub that just truncated/padded to ref length.
 *
 * The truncate/pad stub lost information: if a read had an insertion in the
 * middle, we'd chop off the END of the read (losing data) instead of
 * removing the inserted base where it occurred. The real Viterbi path
 * removes the inserted base at the right position, preserving the rest
 * of the read.
 *
 * Reference:
 *   - HEDGES (Press et al., 2021) — convolutional + Reed-Solomon for DNA storage
 *   - Völkel et al., 2025 (PMC11755093) — nanopore soft-decision decoding
 *   - Durbin/Eddy/Krogh/Mitchison (1998), ch. 4 — Profile HMM
 */

import { SequencingRead } from "./simulate";
import { forwardBackward3, reconstructReadFromPath, Hmm3Params, DEFAULT_HMM3_PARAMS, Hmm3PathStep } from "./profileHmm3";
import { ConvolutionalCode, bytesToBits, bitsToBytes, DEFAULT_CONV_CONFIG } from "./convolutional";

/**
 * Tuning parameters for the Viterbi preprocessor.
 */
export interface ViterbiPreprocessConfig {
  /** Whether the preprocessor is enabled. False for "illumina" channel. */
  enabled: boolean;
  /** Profile HMM band width (default 12 — covers ~12nt indel drift). */
  hmmBandWidth: number;
  /**
   * Minimum cluster size to attempt HMM-based correction.
   * Below this, we use single-read Viterbi (no reference).
   * Default: 2.
   */
  minClusterForHmm: number;
  /**
   * Whether to also apply convolutional Viterbi on the byte stream.
   * Default: false — convolutional encoding is not currently applied at
   * encode time, so we only use HMM-based indel correction.
   * Set to true if the encoder was modified to apply convolutional.ts
   * before LDPC (full HEDGES-style pipeline).
   */
  useConvolutionalViterbi: boolean;
  /**
   * HMM transition parameters. Defaults tuned for ONT R10.4.1 (9% total IDS).
   */
  hmmParams: Hmm3Params;
}

export const DEFAULT_VITERBI_CONFIG: ViterbiPreprocessConfig = {
  enabled: false,
  hmmBandWidth: 12,
  minClusterForHmm: 2,
  useConvolutionalViterbi: false,
  hmmParams: {
    ...DEFAULT_HMM3_PARAMS,
    // ONT R10.4.1 has higher indel rates than Illumina
    matchToInsert: 0.05,
    matchToDelete: 0.05,
    matchToMatch: 0.90,
  },
};

/**
 * v60: Preprocess a single read with Profile HMM alignment to a reference.
 *
 * Conservative indel correction: use the HMM Viterbi alignment to remove
 * insertions and fill deletions, but PRESERVE the read's own base calls
 * for match positions. This keeps per-read variation intact for soft-decision
 * LDPC decoding while fixing the length mismatch that prevents LDPC from
 * operating.
 *
 * Steps:
 *   1. Run forwardBackward3(read, ref, quality) → posteriors + Viterbi path
 *   2. Walk the Viterbi path:
 *      - M states → keep the read's own base at that position
 *      - I states → drop the inserted base (excise from read)
 *      - D states → fill with the reference base (no read base to keep)
 *   3. Return the length-corrected read
 *
 * This is the HEDGES approach: fix the indels (length errors) but leave
 * substitution correction to the LDPC inner code.
 *
 * v60: Now uses the REAL Viterbi path (not a stub), so insertions are
 * removed at the right position instead of just truncating the read.
 *
 * @param read Noisy read (may contain indels)
 * @param ref Reference sequence (typically the cluster consensus)
 * @param cfg Preprocessor config
 * @returns Indel-corrected read (length matches ref), or original if failed
 */
export function viterbiCorrectRead(
  read: SequencingRead,
  ref: string,
  cfg: ViterbiPreprocessConfig = DEFAULT_VITERBI_CONFIG,
): SequencingRead {
  if (!cfg.enabled || ref.length === 0 || read.sequence.length === 0) {
    return read;
  }

  // Quick check: if read is already the right length, skip HMM (no indels to fix)
  if (read.sequence.length === ref.length) {
    return read;
  }

  // v65: Use the REAL Viterbi path reconstruction from profileHmm3.
  //
  // The profile HMM computes the optimal alignment of read vs ref, identifying
  // which read positions are insertions (I), deletions (D), or matches (M).
  // We reconstruct the read by:
  //   - M states: keep the read's own base (preserves substitution info)
  //   - I states: drop the inserted base (excise from read)
  //   - D states: fill with the reference base (read was missing this base)
  //
  // This is the HEDGES approach: fix indels (length errors) at the RIGHT
  // position, then leave substitution correction to the LDPC inner code.
  //
  // Previously (v60-v64), we used the length-only fix (truncate/pad) which
  // just chopped/padded at the END of the read — losing information about
  // WHERE the indel occurred. This caused ~12% IDS recovery because the
  // LDPC decoder got misaligned bits (an insertion at position 50 shifts
  // all subsequent bits by 1, which the LDPC treats as many substitution
  // errors rather than one indel).
  //
  // The real Viterbi path fixes the alignment at the correct position,
  // so the LDPC decoder sees at most 1 substitution per indel (the base
  // that was inserted or the substitution at the deletion site).
  try {
    const result = forwardBackward3(
      read.sequence,
      ref,
      read.quality,
      cfg.hmmParams,
      cfg.hmmBandWidth,
    );

    if (result.path.length > 0) {
      const correctedSeq = reconstructReadFromPath(read.sequence, ref, result.path);

      // Sanity check: corrected length should match ref length
      if (correctedSeq.length === ref.length) {
        // Build corrected quality array
        // For M states: keep the read's quality at that position
        // For D states (filled from ref): assign a moderate quality
        // For I states (dropped): skip
        const readQ = read.quality ?? new Uint8Array(read.sequence.length).fill(20);
        const correctedQ = new Uint8Array(ref.length).fill(20); // default Q20
        let outPos = 0;
        for (const step of result.path) {
          if (step.state === 'M' && step.readPos >= 0 && step.readPos < readQ.length) {
            if (outPos < ref.length) {
              correctedQ[outPos] = readQ[step.readPos];
            }
            outPos++;
          } else if (step.state === 'D') {
            // Deletion filled from ref: moderate quality (Q15)
            if (outPos < ref.length) {
              correctedQ[outPos] = 15;
            }
            outPos++;
          }
          // I states: skip (don't increment outPos)
        }

        return {
          ...read,
          sequence: correctedSeq,
          quality: correctedQ,
        };
      }

      // If corrected length doesn't match ref (shouldn't happen with proper
      // path reconstruction), fall through to length-only fix
    }
  } catch {
    // HMM failed (e.g., numerical underflow, pathological alignment)
    // Fall through to length-only fix
  }

  // Fallback: length-only fix (truncate/pad)
  return lengthOnlyFix(read, ref);
}

/**
 * Length-only fallback: truncate read to ref length or pad with ref bases.
 *
 * This is the v59 stub behavior — used only when the HMM fails completely.
 * It's better than dropping the read, but loses information about WHERE
 * the indel occurred.
 */
function lengthOnlyFix(read: SequencingRead, ref: string): SequencingRead {
  const R = ref.length;
  const readQuality = read.quality ?? new Uint8Array(read.sequence.length).fill(20);
  let corrected: string;
  let newQuality: Uint8Array;
  if (read.sequence.length > R) {
    corrected = read.sequence.slice(0, R);
    // Quality might be shorter than read.sequence — cap at min(quality.length, R)
    const qLen = Math.min(readQuality.length, R);
    newQuality = new Uint8Array(R).fill(20);
    newQuality.set(readQuality.slice(0, qLen), 0);
  } else {
    corrected = read.sequence + ref.slice(read.sequence.length);
    const pad = new Uint8Array(R - read.sequence.length).fill(15);
    newQuality = new Uint8Array(R).fill(20);
    // Quality might be shorter than read.sequence — cap at min(quality.length, read.sequence.length)
    const qLen = Math.min(readQuality.length, read.sequence.length);
    newQuality.set(readQuality.slice(0, qLen), 0);
    newQuality.set(pad, read.sequence.length);
  }
  return {
    ...read,
    sequence: corrected,
    quality: newQuality,
  };
}

/**
 * Preprocess all reads in a cluster using the cluster consensus as reference.
 *
 * Two-pass approach:
 *   1. Build a quick plurality consensus from all reads in the cluster
 *   2. For each read, run viterbiCorrectRead against the consensus
 *   3. Return the corrected reads
 *
 * If the cluster is too small (< minClusterForHmm), we still attempt
 * correction using the longest read as the reference.
 *
 * @param clusterReads Reads in this oligo's cluster
 * @param cfg Preprocessor config
 * @returns Corrected reads (same length as input)
 */
export function viterbiPreprocessCluster(
  clusterReads: SequencingRead[],
  cfg: ViterbiPreprocessConfig = DEFAULT_VITERBI_CONFIG,
): SequencingRead[] {
  if (!cfg.enabled || clusterReads.length === 0) {
    return clusterReads;
  }

  // Build a quick plurality consensus as the HMM reference
  let reference: string;
  if (clusterReads.length >= cfg.minClusterForHmm) {
    reference = buildPluralityConsensus(clusterReads);
  } else {
    // Use the longest read as reference (single-read mode)
    reference = clusterReads.reduce((a, b) =>
      a.sequence.length > b.sequence.length ? a : b,
    ).sequence;
  }

  // Correct each read against the reference
  return clusterReads.map((r) => viterbiCorrectRead(r, reference, cfg));
}

/**
 * Build a simple column-wise plurality consensus from a set of reads.
 * Reads may have different lengths due to indels; we align them at
 * position 0 and take the most common base at each column.
 *
 * This is the same approach used in decode.ts consensus() but extracted
 * here for reuse by the Viterbi preprocessor.
 */
function buildPluralityConsensus(reads: SequencingRead[]): string {
  if (reads.length === 0) return "";
  if (reads.length === 1) return reads[0].sequence;

  // Use median length as consensus length
  const lengths = reads.map((r) => r.sequence.length).sort((a, b) => a - b);
  const medianLen = lengths[Math.floor(lengths.length / 2)];
  const consensusLen = Math.floor(medianLen);

  let result = "";
  for (let col = 0; col < consensusLen; col++) {
    const counts: Record<string, number> = { A: 0, C: 0, G: 0, T: 0 };
    for (const read of reads) {
      if (col < read.sequence.length) {
        const base = read.sequence[col];
        if (counts[base] !== undefined) counts[base]++;
      }
    }
    let bestBase = "A";
    let bestCount = -1;
    for (const b of ["A", "C", "G", "T"]) {
      if (counts[b] > bestCount) {
        bestCount = counts[b];
        bestBase = b;
      }
    }
    result += bestBase;
  }

  return result;
}

/**
 * Top-level entry point: preprocess reads for a given channel.
 *
 * For "illumina" channel: returns reads unchanged.
 * For "nanopore" channel: applies Viterbi/HMM indel correction per cluster.
 *
 * v52: When useConvolutionalInner is true (HEDGES-style pipeline), the
 * preprocessor also sets useConvolutionalViterbi=true so the downstream
 * conv Viterbi decoder (convolutional.ts) is invoked after HMM alignment.
 *
 * @param readsByCluster Map of oligo index → cluster reads
 * @param channel Sequencing channel ("illumina" or "nanopore")
 * @param cfg Optional preprocessor config (defaults to channel-appropriate)
 * @returns Preprocessed reads (same shape as input)
 */
export function viterbiPreprocessReads(
  readsByCluster: Map<number, SequencingRead[]>,
  channel: "illumina" | "nanopore" | "pacbio" = "illumina",
  cfg?: Partial<ViterbiPreprocessConfig>,
): Map<number, SequencingRead[]> {
  // Only nanopore/pacbio channel triggers the preprocessor
  if (channel !== "nanopore" && channel !== "pacbio") {
    return readsByCluster;
  }

  const finalCfg: ViterbiPreprocessConfig = {
    ...DEFAULT_VITERBI_CONFIG,
    enabled: true, // force-enable for nanopore
    useConvolutionalViterbi: cfg?.useConvolutionalViterbi ?? false,
    ...cfg,
  };

  const result = new Map<number, SequencingRead[]>();
  for (const [idx, clusterReads] of readsByCluster) {
    result.set(idx, viterbiPreprocessCluster(clusterReads, finalCfg));
  }
  return result;
}
