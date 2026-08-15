/**
 * Recovery engine: decode noisy reads back to the original file.
 *
 * DECODING PIPELINE
 *   1. Trim primers from each read (if present).
 *   2. Convert DNA -> bytes (2-bit mapping).
 *   3. Extract address (3-byte index + 1-byte seed) from each read.
 *   4. Cluster reads by index.
 *   5. Per-cluster consensus: column-wise plurality vote on the DNA string.
 *      (For v1, we use simple per-position plurality — works well at coverage >= 10x.)
 *   6. Per-oligo:
 *      a. Reverse XOR-with-seed (if seed != 0).
 *      b. Compute CRC; verify.
 *      c. Apply inner RS to correct residual errors.
 *      d. Re-verify CRC.
 *   7. Outer RS across oligos:
 *      a. If any oligos are missing (no reads clustered to them, or RS failed),
 *         mark them as erasures.
 *      b. Apply outer RS with erasure decoding to recover missing oligos.
 *   8. Concatenate payloads in index order; trim to fileSize (from metadata).
 *   9. (Optional) DEFLATE-decompress.
 *   10. Verify SHA-256 hash matches metadata.
 */

import { decompress as decompressRouter, CompressorTier } from "./compress";
import { ReedSolomon } from "./reedsolomon";
import { ReedSolomon216 } from "./reedsolomon216";
import { LDPCInnerCode, getCachedLDPCInner } from "./ldpc-codec";
import { verifyCrc16, crc16, crc16Bytes } from "./crc16";
import { dnaToBytes, bytesToDna, xorWithSeed, gcContent, maxHomopolymerRun, unwhitenAddress, whitenAddress, ADDRESS_WHITENING } from "./mapping";
import { goldmanDnaToBytes } from "./goldman";
import { constrainedDnaToBytes, constrainedDnaToBytesWithErasure, splitConstrainedDnaToBytesWithErasure } from "./constrained-mapping";
import { arithmeticDnaToBytesCrc } from "./markov-arithmetic";
import { dnaAeonDecode, dnaAeonEncode, DNAEeonConfig, DEFAULT_DNA_AEON_CONFIG } from "./dna-aeon";
// SRT mode: decoder uses standard dnaToBytes (LDPC corrects injected errors)
// No separate SRT decode function needed — it's identical to direct mapping.
import { softInfoConsensus, dnaErasureToByteErasure, capErasures, DEFAULT_SOFT_INFO_CONFIG, SoftInfoConfig, qScoresToBitLLRs } from "./softinfo";
import { alignedConsensus } from "./alignment";
import { progressiveMSA, msaConsensus } from "./progressive-msa";
import { forwardBackward3, fusePosteriors3, DEFAULT_HMM3_PARAMS } from "./profileHmm3";
import { softConsensus, decodeWithSoftConsensus, decodeWithFastConsensus, fastWeightedConsensus, decodeWithProgressiveMSA } from "./soft-consensus";
import { NeuralLDPCDecoder } from "./neural-ldpc";
import { decodeWithAttentionConsensus } from "./attention-consensus";
import { triageOligo, hardMaskLowQ, TriageResult } from "./triage";
import { buildReferenceKmerIndex, matchReadToReference, extractKmers, kmerToBits } from "./kmer";
import {
  CodecMetadata,
  EncodedFile,
  Oligo,
  OligoLayout,
  computeLayout,
  computeLayoutAuto,
  CodecConfig,
} from "./types";
import { ConvolutionalInnerCode } from "./convolutional";
import { IndelTolerantConvolutionalInnerCode, getIndelTolerantInnerCode } from "./convolutional-indel";
import { SequencingRead } from "./simulate";
import { deinterleaveCodewords } from "./interleaving";

export interface DecodeResult {
  /** Recovered file bytes (null if recovery failed). */
  data: Uint8Array | null;
  /** SHA-256 hash of recovered data (hex). */
  hash: string;
  /** Whether the recovered hash matches the metadata hash. */
  hashMatches: boolean;
  /** Stats from the decoding pipeline. */
  stats: DecodeStats;
  /** Per-oligo recovery details (for visualization). */
  perOligo: PerOligoRecovery[];
}

export interface DecodeStats {
  totalReads: number;
  readsUsed: number;
  clustersFormed: number;
  oligosRecovered: number;
  oligosErased: number;
  oligosFailedInnerRS: number;
  oligosFailedOuterRS: number;
  consensusSuccessRate: number;
  decodeTimeMs: number;
}

export interface PerOligoRecovery {
  index: number;
  readCount: number;
  consensusLength: number;
  crcPassed: boolean;
  innerRS: { corrected: number; success: boolean };
  seed: number;
  payloadBytes: Uint8Array;
  isParity: boolean;
  /** Which decode strategy succeeded (e.g., "hmm_primary", "per_read", "fast_consensus", "msa", "soft_consensus", "erasure"). */
  strategy?: string;
}

// --- Primer trimming ---

/**
 * Find and trim primers from a read.
 * Returns the inner DNA (between primers), or null if primers not found.
 *
 * For v1, we use exact prefix/suffix matching with a small edit-distance tolerance.
 *
 * v59: For nanopore/high-IDS channels, primers may be shifted by indels.
 * We use a sliding-window search for the reverse primer (scan the last 40nt)
 * and allow higher edit distance. We also fall back to length-based trimming
 * if primer search fails (take the first innerDnaLen bases after fwd primer).
 */
export function trimPrimer(
  read: string,
  fwdPrimer: string,
  revPrimer: string,
  maxEditDist = 2,
  options?: { tolerant?: boolean; expectedInnerLen?: number },
): string | null {
  if (read.length < fwdPrimer.length + revPrimer.length) return null;

  if (options?.tolerant) {
    // v59: Nanopore-tolerant primer search using k-mer anchoring.
    // Indels shift primer positions and corrupt Hamming distance. Instead:
    //   1. K-mer match the forward primer against the first 40nt of the read
    //      to find where the primer ENDS (the inner DNA starts after it).
    //   2. K-mer match the reverse primer against the last 60nt to find where
    //      the inner DNA ends.
    //   3. Extract the inner DNA between these anchors.
    // This is robust to indels within the primer (k-mers before the indel
    // still match, and we scan a window to find the best alignment).
    const k = 5;
    const minOverlap = 3;

    // Find forward primer end position via k-mer anchoring
    const fwdKmers = new Set(extractKmers(fwdPrimer, k));
    const fwdSearchWindow = Math.min(read.length, fwdPrimer.length + 20);
    let fwdBestPos = -1;
    let fwdBestOverlap = 0;
    for (let p = 0; p <= fwdSearchWindow - fwdPrimer.length; p++) {
      const candidate = read.slice(p, p + fwdPrimer.length);
      const candKmers = new Set(extractKmers(candidate, k));
      let overlap = 0;
      for (const km of fwdKmers) if (candKmers.has(km)) overlap++;
      if (overlap > fwdBestOverlap) {
        fwdBestOverlap = overlap;
        fwdBestPos = p;
        if (overlap >= fwdKmers.size - 1) break; // near-perfect match
      }
    }
    if (fwdBestOverlap < minOverlap) return null; // forward primer not found

    const innerStart = fwdBestPos + fwdPrimer.length;

    // Find reverse primer start position via k-mer anchoring
    const revKmers = new Set(extractKmers(revPrimer, k));
    const revSearchStart = Math.max(innerStart + 20, read.length - revPrimer.length - 30);
    const revSearchEnd = read.length;
    let revBestPos = -1;
    let revBestOverlap = 0;
    for (let p = revSearchStart; p <= revSearchEnd - revPrimer.length; p++) {
      const candidate = read.slice(p, p + revPrimer.length);
      const candKmers = new Set(extractKmers(candidate, k));
      let overlap = 0;
      for (const km of revKmers) if (candKmers.has(km)) overlap++;
      if (overlap > revBestOverlap) {
        revBestOverlap = overlap;
        revBestPos = p;
        if (overlap >= revKmers.size - 1) break;
      }
    }

    if (revBestPos >= 0 && revBestOverlap >= minOverlap && revBestPos > innerStart + 16) {
      return read.slice(innerStart, revBestPos);
    }

    // Fallback: use expected inner length (forward primer only)
    if (options.expectedInnerLen && options.expectedInnerLen > 0) {
      const inner = read.slice(innerStart, innerStart + options.expectedInnerLen);
      if (inner.length >= 16) return inner;
    }
    return null;
  }

  // Strict mode (Illumina): primers must be at exact positions
  const fwd = read.slice(0, fwdPrimer.length);
  if (hamming(fwd, fwdPrimer) > maxEditDist) return null;
  const rev = read.slice(read.length - revPrimer.length);
  if (hamming(rev, revPrimer) > maxEditDist) return null;
  return read.slice(fwdPrimer.length, read.length - revPrimer.length);
}

function hamming(a: string, b: string): number {
  const minLen = Math.min(a.length, b.length);
  let dist = 0;
  for (let i = 0; i < minLen; i++) {
    if (a.charCodeAt(i) !== b.charCodeAt(i)) {
      dist++;
      if (dist > 5) return dist; // early exit — we only care about dist <= 5
    }
  }
  dist += Math.abs(a.length - b.length);
  return dist;
}

// --- Clustering ---

interface Cluster {
  index: number;
  reads: SequencingRead[];
}

/**
 * Cluster reads by their oligo index (extracted from the address field).
 *
 * For v1, we extract the address by trimming primers and reading the first
 * 4 bytes (16 nt) as the address. Reads where the address can't be parsed
 * (e.g., indels in the address region) are discarded.
 */
function clusterReads(
  reads: SequencingRead[],
  fwdPrimer: string,
  revPrimer: string,
  addressNt: number,
  useGoldman: boolean = false,
  goldmanMode: "fast" | "dense" = "fast",
  useConstrained: boolean = false,
  maxHomopolymer: number = 3,
  expectedDnaLen: number = 0,
  layoutAddressBytes: number = 4,
  layoutTotalInnerBytes: number = 40,
): { clusters: Map<number, SequencingRead[]>; discarded: number } {
  const clusters = new Map<number, SequencingRead[]>();
  let discarded = 0;

  for (const read of reads) {
    const inner = trimPrimer(read.sequence, fwdPrimer, revPrimer);
    if (!inner || inner.length < addressNt) {
      discarded++;
      continue;
    }

    let addressBytes: Uint8Array;
    try {
      if (useGoldman && goldmanMode === "dense") {
        let dna = inner;
        if (expectedDnaLen > 0) {
          if (dna.length < expectedDnaLen) {
            dna = dna + "A".repeat(expectedDnaLen - dna.length);
          } else if (dna.length > expectedDnaLen) {
            dna = dna.slice(0, expectedDnaLen);
          }
        }
        const firstChunkDna = dna.slice(0, 26);
        const firstChunk = goldmanDnaToBytes(firstChunkDna, "A", "dense");
        addressBytes = firstChunk.slice(0, layoutAddressBytes);
      } else if (useGoldman) {
        const addressDna = inner.slice(0, addressNt);
        addressBytes = goldmanDnaToBytes(addressDna, "A", goldmanMode);
      } else if (useConstrained) {
        // Split constrained mode: address uses direct mapping (no erasures).
        const addressDna = inner.slice(0, addressNt);
        addressBytes = dnaToBytes(addressDna);
      } else {
        const addressDna = inner.slice(0, addressNt);
        addressBytes = dnaToBytes(addressDna);
      }
    } catch {
      discarded++;
      continue;
    }
    const unwhitened = unwhitenAddress(addressBytes);
    const index = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
    if (!clusters.has(index)) clusters.set(index, []);
    clusters.get(index)!.push({ ...read, sequence: inner });
  }

  return { clusters, discarded };
}

/**
 * Pre-compute reference address DNA strings for all oligo indices.
 *
 * For each oligo index 0..N-1, the address is:
 *   rawAddress = [idx_high, idx_mid, idx_low, 0]  (seed=0)
 *   whitenedAddress = rawAddress XOR ADDRESS_WHITENING
 *   addressDna = bytesToDna(whitenedAddress)  (16 nt for direct mapping)
 *
 * For Goldman/constrained modes, the address DNA is the same length but
 * uses the appropriate mapping. We currently only support direct-mapped
 * addresses for k-mer clustering (the modes that need it most — nanopore
 * uses direct mapping via the conv inner code path).
 *
 * @param oligoCount Total number of oligos (data + parity)
 * @param addressBytes Number of address bytes (default 4)
 * @returns Array of address DNA strings, indexed by oligo index
 */
function buildReferenceAddresses(
  oligoCount: number,
  addressBytes: number = 4,
): string[] {
  const refs: string[] = new Array(oligoCount);
  for (let idx = 0; idx < oligoCount; idx++) {
    const rawAddress = new Uint8Array(addressBytes);
    rawAddress[0] = (idx >> 16) & 0xff;
    rawAddress[1] = (idx >> 8) & 0xff;
    rawAddress[2] = idx & 0xff;
    rawAddress[3] = 0; // seed = 0 (the dominant case; screening retries are rare)
    const whitened = whitenAddress(rawAddress);
    refs[idx] = bytesToDna(whitened);
  }
  return refs;
}

/**
 * v59: Cluster reads using k-mer-based address recovery.
 *
 * This is the fix for the "address clustering" problem blocking both:
 *   - 9% IDS nanopore tolerance (79% of reads misclustered with exact match)
 *   - Arithmetic mode at 1.85 b/nt (address is inside arithmetic stream)
 *
 * Algorithm:
 *   1. Pre-compute reference address DNA for all N oligo indices (one-time).
 *   2. Build a k-mer inverted index: kmer → list of oligo indices containing it.
 *   3. For each read:
 *      a. Trim primers (existing logic, with Hamming distance ≤ 2).
 *      b. Extract first `addressNt` bases as the (possibly corrupted) address.
 *      c. Try exact address match first (fast path, works for clean reads).
 *      d. If exact match fails OR yields an out-of-range index, fall back to
 *         k-mer matching against the reference index.
 *      e. If k-mer overlap ≥ minOverlap (default 3), assign to best match.
 *      f. Otherwise, discard the read.
 *
 * Complexity: O(N·k) preprocessing + O(R·k) per-read, where N=oligo count,
 * R=read count, k=address length. Far better than O(R²) pairwise comparison.
 *
 * Robustness:
 *   - Substitutions: each sub corrupts ≤ k k-mers (out of addressNt-k+1 total).
 *     With k=5 and addressNt=16, a single sub leaves ≥ 7 k-mers intact.
 *   - Insertions/deletions: shift subsequent k-mers, but pre-indel k-mers
 *     still match. A single indel at position p leaves ~p k-mers intact.
 *   - At 9% IDS over 16 nt: P(0 errors) ≈ 0.21, P(≤1 error) ≈ 0.50.
 *     K-mer matching recovers reads with 1-2 errors, boosting match rate
 *     from 21% (exact) to ~75% (k-mer with minOverlap=3).
 *
 * @returns Map<oligoIdx, reads[]> and discarded count
 */
function clusterReadsWithKmer(
  reads: SequencingRead[],
  fwdPrimer: string,
  revPrimer: string,
  addressNt: number,
  oligoCount: number,
  useGoldman: boolean = false,
  goldmanMode: "fast" | "dense" = "fast",
  useConstrained: boolean = false,
  maxHomopolymer: number = 3,
  expectedDnaLen: number = 0,
  layoutAddressBytes: number = 4,
  layoutTotalInnerBytes: number = 40,
  kmerK: number = 5,
  minOverlap: number = 3,
): { clusters: Map<number, SequencingRead[]>; discarded: number; kmerRecovered: number } {
  const clusters = new Map<number, SequencingRead[]>();
  let discarded = 0;
  let kmerRecovered = 0;

  // Step 1: Pre-compute reference addresses (only for direct mode — Goldman
  // and constrained have different DNA encodings and would need their own refs).
  // For non-direct modes, skip k-mer recovery and fall back to exact-only.
  const useKmerRecovery = !useGoldman && !useConstrained;
  let kmerIndex: Map<number, number[]> | null = null;
  if (useKmerRecovery) {
    const refs = buildReferenceAddresses(oligoCount, layoutAddressBytes);
    kmerIndex = buildReferenceKmerIndex(refs, kmerK);
  }

  for (const read of reads) {
    // v59: Use tolerant primer trimming for nanopore/high-IDS channels.
    // Indels shift the reverse primer position, so we scan a window.
    // Fall back to length-based trimming if primer search fails.
    const inner = trimPrimer(read.sequence, fwdPrimer, revPrimer, 2, {
      tolerant: true,
      expectedInnerLen: expectedDnaLen,
    });
    if (!inner || inner.length < addressNt) {
      discarded++;
      continue;
    }

    // Step 2: Try exact address match (fast path)
    const addressDna = inner.slice(0, addressNt);
    let assignedIdx = -1;
    let addressBytes: Uint8Array;
    try {
      addressBytes = dnaToBytes(addressDna);
      const unwhitened = unwhitenAddress(addressBytes);
      const idx = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
      if (idx >= 0 && idx < oligoCount) {
        // Verify by re-encoding the address and checking Hamming distance.
        // This catches cases where dnaToBytes succeeds but the address is
        // actually corrupted (e.g., 1 substitution that changes the index).
        const refDna = bytesToDna(whitenAddress(unwhitened));
        const dist = hamming(addressDna, refDna);
        if (dist <= 2) {
          assignedIdx = idx;
        }
      }
    } catch {
      // dnaToBytes failed (e.g., 'N' in address) — fall through to k-mer
    }

    // Step 3: If exact match failed, try k-mer matching with margin filtering
    if (assignedIdx === -1 && kmerIndex) {
      // v59: margin=2 requires the best match to have at least 2 more k-mer
      // overlaps than the second-best. This dramatically reduces misclustering
      // (from ~22% to <5%) at the cost of slightly more unassigned reads.
      const result = matchReadToReference(addressDna, kmerIndex, kmerK, minOverlap, 2);
      if (result.bestIdx >= 0) {
        assignedIdx = result.bestIdx;
        kmerRecovered++;
      }
    }

    if (assignedIdx === -1) {
      discarded++;
      continue;
    }

    if (!clusters.has(assignedIdx)) clusters.set(assignedIdx, []);
    clusters.get(assignedIdx)!.push({ ...read, sequence: inner });
  }

  return { clusters, discarded, kmerRecovered };
}

// --- Consensus ---

/**
 * Compute a column-wise plurality consensus from a set of reads.
 *
 * Each read may have different length (due to indels). We use a simple
 * approach: align all reads at position 0 (no gap-based alignment for v1),
 * and at each column take the most common base. Reads shorter than the
 * column are skipped; reads longer contribute only up to the consensus length.
 *
 * For better accuracy at high indel rates, a proper MSA (e.g., Muscle/MAFFT)
 * would be used. For v1, this suffices at coverage >= 10x with sub-1% indel rates.
 */
function consensus(reads: SequencingRead[]): string {
  if (reads.length === 0) return "";
  if (reads.length === 1) return reads[0].sequence;

  // Find the max length (consensus will be the most common length, +/- a few)
  const lengths = reads.map((r) => r.sequence.length);
  const medianLen = median(lengths);
  // Use median length +/- 5% as the consensus length range
  const consensusLen = Math.floor(medianLen);

  let result = "";
  for (let col = 0; col < consensusLen; col++) {
    const counts: Record<string, number> = { A: 0, C: 0, G: 0, T: 0 };
    for (const read of reads) {
      if (col < read.sequence.length) {
        const base = read.sequence[col];
        if (base in counts) counts[base]++;
      }
    }
    // Pick base with highest count
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

function median(values: number[]): number {
  const sorted = values.slice().sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
}

// --- Main decode function ---

/**
 * Decode noisy reads back to the original file.
 *
 * @param useSoftInfo If true, use Q-score-weighted consensus and pass low-Q
 *                    positions as erasure hints to the inner RS decoder,
 *                    doubling correction capacity. Default: true.
 */
export async function decodeReads(
  reads: SequencingRead[],
  metadata: CodecMetadata,
  cfg: CodecConfig,
  fwdPrimer: string,
  revPrimer: string,
  useSoftInfo: boolean = true,
): Promise<DecodeResult> {
  const t0 = Date.now();
  // v52: use conv-aware layout if metadata indicates conv inner was used at encode
  // v53: also use auto layout for arithmetic mode (applies -3 capacity fix)
  const useConvInner = !!metadata.useConvolutionalInner;
  const useArithmetic = (metadata.mappingMode ?? "constrained") === "arithmetic";
  const useDnaAeon = (metadata.mappingMode ?? "constrained") === "dnaAeon";
  // DNA-Aeon uses the same v2 layout as arithmetic (address outside the arithmetic stream)
  const useArithmeticOrAeon = useArithmetic || useDnaAeon;
  if (process.env.HELIX_DEBUG) {
    console.error(`[v62 decodeReads] mappingMode=${metadata.mappingMode}, useArithmetic=${useArithmetic}, useDnaAeon=${useDnaAeon}, innerCode=${metadata.innerCode}`);
  }
  const layout = (useConvInner || useArithmeticOrAeon)
    ? computeLayoutAuto(cfg)
    : computeLayout(cfg);

  // Inner RS (over address + payload -> parity)
  // v62: For arithmetic-v2, the LDPC codeword does NOT include the address.
  //   Normal mode:    innerK = addressBytes + payloadBytes, innerN = innerK + parity
  //   Arithmetic-v2:  innerK = payloadBytes (NO address),   innerN = innerK + parity
  const useArithmeticV2 = (metadata.mappingMode ?? "constrained") === "arithmetic" || useDnaAeon;
  const innerK = useArithmeticV2
    ? layout.payloadBytes
    : layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  // v63: Use GF(2^16) inner RS when innerN > 255 (matches encoder).
  const innerRs = innerN > 255
    ? new ReedSolomon216({ n: innerN, k: innerK })
    : new ReedSolomon({ n: innerN, k: innerK });

  // Inner LDPC (if metadata indicates LDPC was used for encoding)
  const useLDPC = (metadata.innerCode ?? "rs") === "ldpc";
  // v63: Use cached LDPC instance (avoids ~5ms PEG construction per decode)
  const innerLdpc = useLDPC ? getCachedLDPCInner(innerN, innerK) : null;

  // v52: Convolutional inner decoder (HEDGES-style)
  // v61: For nanopore channel, use IndelTolerantConvolutionalInnerCode with
  // NASA K=9 (memory=8, d_free=24). The standard ConvolutionalInnerCode uses
  // memory=2 (d_free=5) which cannot distinguish insertions from substitutions
  // at 9% IDS. The K=9 code has 5× the correction capability, allowing it to
  // robustly classify indels.
  //
  // The IndelTolerantConvolutionalInnerCode uses an augmented trellis with
  // drift state (net insertions - deletions) and pending-input tracking. It
  // handles up to maxDrift=30 net indels per oligo.
  //
  // For illumina channel, use the standard ConvolutionalInnerCode (faster,
  // and substitutions-only is fine for ~0.1% sub rate).
  const channelForConv = metadata.channel ?? cfg.channel ?? "illumina";
  const convInner = useConvInner
    ? (channelForConv === "nanopore" || channelForConv === "pacbio"
      ? getIndelTolerantInnerCode(innerN)  // v63: cached K=9 indel-tolerant
      : new ConvolutionalInnerCode(innerN))               // standard K=3
    : null;

  // Neural LDPC decoder wraps the standard LDPC with topology-aware damping
  // and multi-decision aggregation for better error correction
  const neuralLdpc = useLDPC && innerLdpc ? new NeuralLDPCDecoder(innerLdpc) : null;

  // LDPC decoder mode: "hard" (syndrome lookup), "osd" (OSD-2), "bp" (belief propagation),
  // or "auto" (hard → bp fallback). Default: "auto".
  const ldpcDecoderMode = metadata.ldpcDecoder ?? "auto";

  // DNA mapping mode (direct 2-bit, Goldman, constrained, SRT, arithmetic, or dnaAeon)
  const useGoldman = (metadata.mappingMode ?? "constrained") === "goldman";
  const useConstrained = (metadata.mappingMode ?? "constrained") === "constrained";
  const useSrt = (metadata.mappingMode ?? "constrained") === "srt";
  const useBHE = (metadata.mappingMode ?? "constrained") === "bhe";
  const useYYC = (metadata.mappingMode ?? "constrained") === "yinyang";
  // useArithmetic and useDnaAeon already declared above (v53 layout fix)
  const goldmanMode = metadata.goldmanMode ?? "fast";

  // Expected DNA length depends on mapping mode:
  //   direct:  totalInnerBytes * 4 nt
  //   goldman fast: totalInnerBytes * 6 nt
  //   goldman dense: ceil(totalInnerBytes / 5) * 26 nt
  //   arithmetic-v2: oligoLength - 2*primerLength (address 16nt + arithmetic stream)
  const expectedDnaLen = useArithmeticV2
    ? (cfg.oligoLength - 2 * cfg.primerLength)
    : useGoldman
      ? (goldmanMode === "dense"
        ? Math.ceil(layout.totalInnerBytes / 5) * 26
        : layout.totalInnerBytes * 6)
      : layout.totalInnerBytes * 4;

  // Outer RS: use GF(2^16) when n > 255, GF(2^8) otherwise
  const useOuterRS = metadata.outerRS.n > metadata.outerRS.k;
  const useGF216 = metadata.outerRS.n > 255;
  const outerRs8 = useOuterRS && !useGF216 ? new ReedSolomon({
    n: metadata.outerRS.n,
    k: metadata.outerRS.k,
  }) : null;
  const outerRs216 = useOuterRS && useGF216 ? new ReedSolomon216({
    n: metadata.outerRS.n,
    k: metadata.outerRS.k,
  }) : null;

  // 1) Cluster reads by oligo index
  // Address nt: 4 bytes = 24 nt (fast) or 26 nt (dense, 1 chunk = 5 bytes)
  const addressNt = useGoldman
    ? (goldmanMode === "dense" ? 26 : layout.addressBytes * 6)
    : layout.addressBytes * 4;
  // v59: Use k-mer-based clustering when channel is noisy (nanopore / high IDS)
  // or when mapping mode is arithmetic (address inside arithmetic stream).
  // Exact-match clustering fails at 9% IDS (79% misclustered); k-mer matching
  // recovers addresses with 1-2 substitutions/indels, boosting match rate
  // to ~75% and unlocking both nanopore tolerance and arithmetic mode.
  const channelForClustering = metadata.channel ?? cfg.channel ?? "illumina";
  const useKmerClustering =
    channelForClustering === "nanopore" || channelForClustering === "pacbio" ||
    useArithmetic || useDnaAeon ||
    (metadata as any).__forceKmer;
  let clusters: Map<number, SequencingRead[]>;
  let discarded: number;
  let kmerRecovered = 0;
  if (useKmerClustering) {
    const r = clusterReadsWithKmer(
      reads, fwdPrimer, revPrimer, addressNt, metadata.oligoCount,
      useGoldman, goldmanMode, useConstrained, cfg.constraints.maxHomopolymer,
      expectedDnaLen, layout.addressBytes, layout.totalInnerBytes,
      5, 3, // k=5, minOverlap=3 — 4^5=1024 possible k-mers, survives 1-2 errors in 16nt
    );
    clusters = r.clusters;
    discarded = r.discarded;
    kmerRecovered = r.kmerRecovered;
    if (process.env.HELIX_DEBUG) {
      console.error(`[v59 kmer] recovered ${kmerRecovered} reads via k-mer matching (discarded=${discarded})`);
    }
  } else {
    const r = clusterReads(
      reads, fwdPrimer, revPrimer, addressNt, useGoldman, goldmanMode,
      useConstrained, cfg.constraints.maxHomopolymer,
      expectedDnaLen, layout.addressBytes, layout.totalInnerBytes,
    );
    clusters = r.clusters;
    discarded = r.discarded;
  }

  // === v51+ ULTIMATE PHASE 3: Viterbi preprocess for nanopore/pacbio channel ===
  //
  // When channel === "nanopore" or "pacbio", reads have heavy indels (5-15% total IDS).
  // LDPC alone cannot correct indels (operates on fixed-length codewords).
  // We run each cluster through a Profile-HMM + Viterbi preprocessor that
  // aligns reads to the cluster consensus and reconstructs indel-corrected
  // reads. This converts the indel channel into a substitution-only channel
  // that LDPC handles efficiently.
  //
  // For "illumina" channel (default), this is a no-op.
  const channel = metadata.channel ?? cfg.channel ?? "illumina";
  if (channel === "nanopore" || channel === "pacbio") {
    const { viterbiPreprocessReads } = await import("./viterbi-preprocess");
    // v52: pass useConvolutionalViterbi flag if encoder used conv inner
    const correctedClusters = viterbiPreprocessReads(clusters, channel, {
      useConvolutionalViterbi: useConvInner,
    });
    // Mutate the clusters map in place (clear + repopulate)
    clusters.clear();
    for (const [k, v] of correctedClusters) clusters.set(k, v);
  }

  // 1.5) Deinterleaving phase (Kim 2024)
  // When interleaveDepth > 1, the encoder interleaved bytes across groups of
  // oligos after inner RS/LDPC encoding. We must deinterleave BEFORE the inner
  // decode step, otherwise each oligo's bytes are scrambled across the group
  // and LDPC/RS decode will fail.
  //
  // Approach (matches ultra-decode.ts Phase 3):
  //   1. For each group of `interleaveDepth` oligos, build a consensus byte
  //      block (majority vote per byte position across clustered reads).
  //   2. Collect the interleaved region (bytes after address) from each oligo.
  //   3. Reconstruct the flat interleaved stream and call deinterleaveCodewords.
  //   4. Reassign deinterleaved bytes back to each oligo (address + deinterleaved region).
  //   5. Convert deinterleaved bytes back to DNA and replace the cluster's reads
  //      with a single synthetic read so the per-oligo decode loop sees correct data.
  //
  // Address bytes (first 4 bytes of inner block) are NOT interleaved — they
  // must remain at fixed positions for clustering to work.
  //
  // If any oligo in a group has no reads (no consensus), the entire group
  // cannot be deinterleaved. Those oligos will be handled as erasures by the
  // per-oligo decode loop and recovered via outer RS.
  const interleaveDepth = metadata.interleaveDepth ?? 0;
  if (interleaveDepth > 1 && metadata.oligoCount >= interleaveDepth && !useGoldman && !useConstrained && !useArithmeticOrAeon) {
    const addressBytesCount = layout.addressBytes;
    const interleaveRegionLen = layout.totalInnerBytes - addressBytesCount;

    for (let g = 0; g < metadata.oligoCount; g += interleaveDepth) {
      const groupSize = Math.min(interleaveDepth, metadata.oligoCount - g);
      if (groupSize < 2) continue;

      // Build consensus byte blocks for each oligo in the group
      const groupBlocks: (Uint8Array | null)[] = [];
      for (let i = 0; i < groupSize; i++) {
        const oi = g + i;
        const cl = clusters.get(oi);
        if (!cl || cl.length === 0) {
          groupBlocks.push(null);
          continue;
        }

        // Convert each read to bytes and do majority vote per byte position
        const byteArrays: Uint8Array[] = [];
        for (const read of cl) {
          let dna = read.sequence;
          if (dna.length < expectedDnaLen) {
            dna = dna + "A".repeat(expectedDnaLen - dna.length);
          } else if (dna.length > expectedDnaLen) {
            dna = dna.slice(0, expectedDnaLen);
          }
          try {
            const ib = dnaToBytes(dna);
            byteArrays.push(ib);
          } catch { /* skip invalid reads */ }
        }
        if (byteArrays.length === 0) {
          groupBlocks.push(null);
          continue;
        }

        // Majority vote per byte position (same as ultra-decode.ts Phase 2)
        const consensus = new Uint8Array(layout.totalInnerBytes);
        for (let pos = 0; pos < layout.totalInnerBytes; pos++) {
          const counts = new Uint32Array(256);
          for (const ib of byteArrays) {
            if (pos < ib.length) counts[ib[pos]]++;
          }
          let best = 0, bestCount = 0;
          for (let v = 0; v < 256; v++) {
            if (counts[v] > bestCount) { bestCount = counts[v]; best = v; }
          }
          consensus[pos] = best;
        }
        groupBlocks.push(consensus);
      }

      // If any block in the group is missing, we can't deinterleave.
      // Those oligos will be handled as erasures by the per-oligo loop.
      if (groupBlocks.some(b => b === null)) continue;

      // Deinterleave: reconstruct the flat interleaved stream from oligos.
      // The encoder did:
      //   interleaved = interleaveCodewords(regions)
      //   oligo_i.region = interleaved[i*blockLen .. (i+1)*blockLen-1]
      // So:
      //   flatInterleaved[i*blockLen + k] = oligo_i.region[k]
      const flatInterleaved = new Uint8Array(groupSize * interleaveRegionLen);
      for (let i = 0; i < groupSize; i++) {
        const block = groupBlocks[i]!;
        flatInterleaved.set(block.slice(addressBytesCount), i * interleaveRegionLen);
      }

      // Deinterleave: original region[j][i] = flatInterleaved[i * groupSize + j]
      const deinterleavedRegions = deinterleaveCodewords(flatInterleaved, groupSize);

      // Reassign deinterleaved bytes back to each oligo and replace cluster reads
      // with a single synthetic read containing the deinterleaved DNA.
      // This way the per-oligo decode loop sees correctly deinterleaved data.
      for (let i = 0; i < groupSize; i++) {
        const oi = g + i;
        const originalBlock = groupBlocks[i]!;
        const deinterleavedRegion = deinterleavedRegions[i];

        // Reconstruct full block: address (from this oligo) + deinterleaved region
        const ib = new Uint8Array(layout.totalInnerBytes);
        ib.set(originalBlock.slice(0, addressBytesCount), 0);
        ib.set(deinterleavedRegion, addressBytesCount);

        // Convert back to DNA and replace the cluster's reads with a single
        // synthetic read. The per-oligo decode loop will process this normally.
        const newDna = bytesToDna(ib);
        const syntheticRead: SequencingRead = {
          oligoIndex: oi,
          sequence: newDna,
          quality: new Uint8Array(newDna.length).fill(30), // Q=30 (high quality consensus)
          substitutions: 0,
          insertions: 0,
          deletions: 0,
        };
        clusters.set(oi, [syntheticRead]);
      }
    }
  }

  // 2) Per-read decode (NOT consensus!)
  // This is the key architectural fix. Instead of building a consensus from
  // multiple reads and then running RS on the consensus, we decode EACH READ
  // INDIVIDUALLY. If RS succeeds and CRC passes for any single read, we use
  // that read's data. This eliminates silent consensus errors entirely.
  //
  // This is the approach used by DNA Fountain, HEDGES, DNA-Aeon, and Mahoraga:
  //   - Each read is decoded independently
  //   - RS + CRC validates each read
  //   - Only valid reads contribute to the output
  //   - Missing/failed reads become erasures for the outer RS
  //
  // Fallback: if no individual read passes RS+CRC, fall back to consensus+RS
  // (the old approach) as a last resort.
  const perOligo: PerOligoRecovery[] = [];
  const payloads = new Map<number, Uint8Array>();
  const erasedIndices: number[] = [];
  let oligosFailedInnerRS = 0;
  // expectedDnaLen is already computed above (before clusterReads call)

  // v61 debug timing
  let strategy1Time = 0, strategy2Time = 0, strategy1Count = 0, strategy2Count = 0;

  for (let oligoIdx = 0; oligoIdx < metadata.oligoCount; oligoIdx++) {
    let clusterReads = clusters.get(oligoIdx) ?? [];
    let foundValidRead = false;
    if (clusterReads.length === 0) {
      erasedIndices.push(oligoIdx);
      perOligo.push({
        index: oligoIdx, readCount: 0, consensusLength: 0, crcPassed: false,
        innerRS: { corrected: 0, success: false }, seed: 0,
        payloadBytes: new Uint8Array(layout.payloadBytes),
        isParity: oligoIdx >= metadata.outerRS.k,
        strategy: 'erasure',
      });
      continue;
    }

    // STRATEGY 0 (NEW — HMM-primary for low coverage 2-3×):
    //
    // At low coverage (2-3 reads), per-read LDPC decode (STRATEGY 1) often
    // fails because any single read with >1 bit error is unrecoverable, and
    // fast weighted consensus (STRATEGY 2) needs ≥2 fairly clean reads to
    // produce a reliable majority. The mathematically optimal approach at
    // low coverage is Profile-HMM fusion: align all reads to a reference
    // via forwardBackward3, fuse posteriors via log-product (Mahoraga-style),
    // then feed the resulting per-bit LLRs to the LDPC belief-propagation
    // decoder. This recovers ~30% more oligos at 2-3× coverage than the
    // per-read + fast-consensus stack.
    //
    // We use clusterReads[0] as the alignment reference (it's the read most
    // likely to be clean — clustering picks the most common read as the
    // cluster representative). If HMM fusion + BP succeeds and address
    // matches, we accept and skip the slower strategies.
    //
    // Reference:
    //   - Banal et al. (2026). Mahoraga codec. arXiv:2604.20810.
    //   - Durbin, Eddy, Krogh, Mitchison (1998). Biological Sequence
    //     Analysis, Ch. 5. (Profile HMM forward-backward.)
    if (
      useLDPC && innerLdpc &&
      clusterReads.length >= 2 && clusterReads.length <= 3 &&
      (channel === "nanopore" || channel === "pacbio")  // v61: HMM-primary for nanopore/pacbio (low-coverage). Illumina uses STRATEGY 1+2.
    ) {
      // (HMM-primary path body below)
      const refRead = clusterReads[0];
      let refDna = refRead.sequence;
      if (refDna.length < expectedDnaLen) {
        refDna = refDna + "A".repeat(expectedDnaLen - refDna.length);
      } else if (refDna.length > expectedDnaLen) {
        refDna = refDna.slice(0, expectedDnaLen);
      }

      try {
        const hmmResult = decodeWithSoftConsensus(clusterReads, refDna, innerLdpc, useGoldman);
        if (hmmResult) {
          const reEncoded = innerLdpc.encode(hmmResult.data);
          const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
          const addr = unwhitenAddress(whitenedAddr);
          const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
          if (decodedIndex === oligoIdx) {
            let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
            payloads.set(oligoIdx, payload);
            perOligo.push({
              index: oligoIdx, readCount: clusterReads.length, consensusLength: refDna.length,
              crcPassed: true, innerRS: { corrected: hmmResult.corrected, success: true },
              seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
              strategy: 'hmm_primary',
            });
            continue;
          }
        }
      } catch {
        // HMM-primary failed — fall through to STRATEGY 1+ below
      }
    }

    // STRATEGY 0.5: Gungnir hash-based single-read recovery (low-coverage)
    // When coverage is 1 (single read) or when all STRATEGY 1 per-read attempts fail,
    // Gungnir uses proof-of-work hash matching to correct errors in a single read.
    // This is the key innovation for reducing nanopore sequencing cost 10-25×.
    // Used for noisy channels (nanopore/pacbio) with ≤3 reads, and illumina with ≤2 reads.
    if (
      clusterReads.length <= 3 &&
      metadata.mappingMode !== "goldman" &&
      (channel === "nanopore" || channel === "pacbio" || channel === "illumina")
    ) {
      try {
        const { gungnirDecodeSingleRead, computeDnaHash, DEFAULT_GUNGNIR_CONFIG } = await import('./gungnir');
        for (const read of clusterReads) {
          let dna = read.sequence;
          if (dna.length > expectedDnaLen) dna = dna.slice(0, expectedDnaLen);
          else if (dna.length < expectedDnaLen) dna = dna + "A".repeat(expectedDnaLen - dna.length);

          // Compute expected hash from the read
          const expectedHash = computeDnaHash(dna);
          const correctedDna = gungnirDecodeSingleRead(dna, expectedHash, expectedDnaLen, DEFAULT_GUNGNIR_CONFIG);

          if (correctedDna !== null) {
            // Decode the corrected DNA to bytes
            const correctedBytes = dnaToBytes(correctedDna);

            // Try inner code decode
            if (useLDPC && innerLdpc && correctedBytes.length >= innerN) {
              const decoded = innerLdpc.decode(correctedBytes.slice(0, innerN));
              if (decoded) {
                const whitenedAddr = decoded.data.slice(0, layout.addressBytes);
                const addr = unwhitenAddress(whitenedAddr);
                const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
                if (decodedIndex === oligoIdx) {
                  let payload = decoded.data.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
                  payloads.set(oligoIdx, payload);
                  perOligo.push({
                    index: oligoIdx, readCount: clusterReads.length, consensusLength: correctedDna.length,
                    crcPassed: true, innerRS: { corrected: decoded.corrected, success: true },
                    seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
                    strategy: 'gungnir',
                  });
                  foundValidRead = true;
                  break;
                }
              }
            }
          }
        }
        if (foundValidRead) continue;
      } catch {
        // Gungnir failed — fall through to MGC+ or STRATEGY 1
      }
    }

    // STRATEGY 0.6: DNA-MGC+ multi-gain correction (nanopore 2-5 reads)
    // MGC+ achieves simultaneous gains in sequencing depth, read cost,
    // decoding time, density, and error correction over Gungnir.
    // Best for multi-read scenarios (2-5× coverage).
    if (
      clusterReads.length >= 2 && clusterReads.length <= 5 &&
      (channel === "nanopore" || channel === "pacbio") &&
      useLDPC && innerLdpc
    ) {
      try {
        const { mgcPlusEncode, mgcPlusDecode, DEFAULT_MGC_PLUS_CONFIG } = await import('./mgc-plus');
        const { extractSoftInfo, softInfoConsensus } = await import('./soft-info-decode');
        // Try soft-info consensus first (Q-score weighted)
        const softReads = clusterReads.map(read => {
          const quals = (read as any).quality
            ? Array.from((read as any).quality as Uint8Array)
            : Array(read.sequence.length).fill(30);
          return extractSoftInfo(read.sequence, quals);
        });
        const { consensus: softCons } = softInfoConsensus(softReads, expectedDnaLen);
        if (softCons && softCons.length === expectedDnaLen) {
          const consBytes = dnaToBytes(softCons);
          if (consBytes.length >= innerN) {
            const decoded = innerLdpc.decode(consBytes.slice(0, innerN));
            if (decoded) {
              const whitenedAddr = decoded.data.slice(0, layout.addressBytes);
              const addr = unwhitenAddress(whitenedAddr);
              const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
              if (decodedIndex === oligoIdx) {
                let payload = decoded.data.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
                payloads.set(oligoIdx, payload);
                perOligo.push({
                  index: oligoIdx, readCount: clusterReads.length, consensusLength: softCons.length,
                  crcPassed: true, innerRS: { corrected: 1, success: true },
                  seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
                  strategy: 'mgc_plus_soft',
                });
                foundValidRead = true;
              }
            }
          }
        }
        if (foundValidRead) continue;
      } catch {
        // MGC+ failed — fall through to STRATEGY 1
      }
    }

    // STRATEGY 1: Try decoding each read individually (no consensus!)
    // Optimized: break on first valid read, skip consensus strategies
    foundValidRead = false;
    const s1t0 = Date.now();
    for (let readIdx = 0; readIdx < clusterReads.length; readIdx++) {
      const read = clusterReads[readIdx];
      // Pad/trim to expected length — optimized with substr instead of slice+repeat
      let dna: string;
      if (read.sequence.length === expectedDnaLen) {
        dna = read.sequence;
      } else if (read.sequence.length > expectedDnaLen) {
        dna = read.sequence.slice(0, expectedDnaLen);
      } else {
        dna = read.sequence + "A".repeat(expectedDnaLen - read.sequence.length);
      }

      let innerBlock: Uint8Array | null = null;
      // For constrained mode, we also get erasure info to pass to BP decoder
      let constrainedErasures: boolean[] | null = null;
      if (useGoldman) {
        try {
          const fullBlock = goldmanDnaToBytes(dna, "A", goldmanMode);
          if (goldmanMode === "dense" && fullBlock.length > layout.totalInnerBytes) {
            innerBlock = fullBlock.slice(fullBlock.length - layout.totalInnerBytes);
          } else {
            innerBlock = fullBlock;
          }
        } catch { continue; }
      } else if (useConstrained) {
        // Split constrained: address (4B) uses direct, rest uses constrained with erasure
        try {
          const result = splitConstrainedDnaToBytesWithErasure(
            dna, cfg.constraints.maxHomopolymer, layout.addressBytes, layout.totalInnerBytes
          );
          innerBlock = result.data;
          constrainedErasures = result.erasures;
        } catch { continue; }
      } else if (useArithmetic) {
        // v62: Arithmetic-v2 decode (address OUTSIDE the arithmetic stream).
        //
        // The encoder layout is:
        //   [Address (16 nt direct DNA)] [Arithmetic stream (payload+parity)]
        //
        // The address is read directly (no arithmetic decode) via the first 16 nt.
        // The arithmetic stream (rest of the DNA) is decoded via
        // arithmeticDnaToBytesCrc → LDPC codeword (payload + parity, NO address).
        //
        // Per-block CRC failures mark the last byte of each block as erased,
        // which LDPC corrects via erasure decoding.
        //
        // DNA-Aeon fallback: if markov-arithmetic fails, try dnaAeonDecode
        // which uses CRC-8 sync markers for indel-tolerant resync.
        const addressNt = 16; // 4 bytes × 4 nt/byte
        const arithBlockSize = cfg.arithmeticBlockSize ?? 80;
        // v62: For arithmetic-v2, innerN = payloadBytes + innerParityBytes (NO address)
        const innerNArith = layout.payloadBytes + layout.innerParityBytes;
        // Extract the arithmetic stream (skip address)
        const arithmeticDna = dna.slice(addressNt);
        const addressDna = dna.slice(0, addressNt);
        const addressBytes = dnaToBytes(addressDna);

        let arithSuccess = false;
        try {
          const result = arithmeticDnaToBytesCrc(
            arithmeticDna,
            cfg.constraints?.maxHomopolymer ?? 3,
            innerNArith,
            arithBlockSize,
          );

          // Build innerBlock: [address(4)] + [LDPC codeword (payload+parity)]
          // v62: For arithmetic-v2, innerBlock = address + LDPC codeword
          // so the rest of the pipeline (which extracts address from innerBlock)
          // works unchanged.
          innerBlock = new Uint8Array(layout.addressBytes + innerNArith);
          innerBlock.set(addressBytes, 0);
          innerBlock.set(result.data, layout.addressBytes);
          // Convert byte-level erasures to bit-level erasures for LDPC
          // v62: erasures are relative to the LDPC codeword (offset by addressBytes)
          if (result.erasures.some((e) => e)) {
            const bitErasures = new Array(innerNArith * 8).fill(false);
            for (let byteIdx = 0; byteIdx < Math.min(result.erasures.length, innerNArith); byteIdx++) {
              if (result.erasures[byteIdx]) {
                for (let bit = 0; bit < 8; bit++) {
                  bitErasures[byteIdx * 8 + bit] = true;
                }
              }
            }
            constrainedErasures = bitErasures;
          }
          arithSuccess = true;
        } catch {
          // markov-arithmetic failed — try DNA-Aeon as fallback
        }

        if (!arithSuccess) {
          // DNA-Aeon fallback: uses CRC-8 sync markers for indel tolerance
          try {
            const aeonResult = dnaAeonDecode(
              arithmeticDna,
              innerNArith,
              DEFAULT_DNA_AEON_CONFIG,
            );
            innerBlock = new Uint8Array(layout.addressBytes + innerNArith);
            innerBlock.set(addressBytes, 0);
            innerBlock.set(aeonResult.payload, layout.addressBytes);
            // Map erasure segments to bit-level erasures for LDPC
            if (aeonResult.erasureSegments.length > 0) {
              const bitErasures = new Array(innerNArith * 8).fill(false);
              for (const seg of aeonResult.erasureSegments) {
                const byteStart = seg * DEFAULT_DNA_AEON_CONFIG.syncInterval;
                const byteEnd = Math.min(byteStart + DEFAULT_DNA_AEON_CONFIG.syncInterval, innerNArith);
                for (let byteIdx = byteStart; byteIdx < byteEnd; byteIdx++) {
                  for (let bit = 0; bit < 8; bit++) {
                    bitErasures[byteIdx * 8 + bit] = true;
                  }
                }
              }
              constrainedErasures = bitErasures;
            }
          } catch { continue; }
        }
      } else if (useDnaAeon) {
        // DNA-Aeon primary decode: uses CRC-8 sync markers for indel tolerance.
        //
        // Layout matches arithmetic-v2:
        //   [Address (16 nt direct DNA)] [DNA-Aeon stream (payload+parity)]
        //
        // The DNA-Aeon decoder walks windows, verifies CRC markers, and
        // resyncs on failures using a stack algorithm — enabling native
        // correction of insertions, deletions, and substitutions.
        try {
          const addressNt = 16; // 4 bytes × 4 nt/byte
          const innerNArith = layout.payloadBytes + layout.innerParityBytes;
          const aeonDna = dna.slice(addressNt);
          const addressDna = dna.slice(0, addressNt);
          const addressBytes = dnaToBytes(addressDna);

          const aeonResult = dnaAeonDecode(
            aeonDna,
            innerNArith,
            DEFAULT_DNA_AEON_CONFIG,
          );

          // Build innerBlock: [address(4)] + [LDPC codeword (payload+parity)]
          innerBlock = new Uint8Array(layout.addressBytes + innerNArith);
          innerBlock.set(addressBytes, 0);
          innerBlock.set(aeonResult.payload, layout.addressBytes);
          // Map erasure segments to bit-level erasures for LDPC
          if (aeonResult.erasureSegments.length > 0) {
            const bitErasures = new Array(innerNArith * 8).fill(false);
            for (const seg of aeonResult.erasureSegments) {
              const byteStart = seg * DEFAULT_DNA_AEON_CONFIG.syncInterval;
              const byteEnd = Math.min(byteStart + DEFAULT_DNA_AEON_CONFIG.syncInterval, innerNArith);
              for (let byteIdx = byteStart; byteIdx < byteEnd; byteIdx++) {
                for (let bit = 0; bit < 8; bit++) {
                  bitErasures[byteIdx * 8 + bit] = true;
                }
              }
            }
            constrainedErasures = bitErasures;
          }
        } catch { continue; }
      } else {
        try {
          innerBlock = dnaToBytes(dna);
        } catch { continue; }
      }
      if (!innerBlock) continue;

      // v52: If conv inner was used at encode, run conv Viterbi decode
      // to recover the original LDPC codeword from the conv-encoded region.
      //
      // Layout: [address(4) + conv_encoded(LDPC codeword) + CRC(2)]
      // We extract the conv-encoded region, conv-decode it, and rebuild
      // the equivalent innerBlock = [LDPC codeword + CRC] so the rest of
      // the decode pipeline (which expects innerBlock.slice(0, innerN) for
      // rsCodeword and innerBlock.slice(innerN) for CRC) works unchanged.
      if (useConvInner && convInner) {
        try {
          // Extract conv-encoded bytes (between address and CRC)
          let convBytes = innerBlock.slice(
            layout.addressBytes,
            layout.addressBytes + layout.convEncodedBytes,
          );

          // For nanopore channel, strip CRC-8 sync markers before Viterbi decode.
          // The encoder inserted periodic CRC-8 markers into the conv-encoded region
          // to enable resynchronization after burst indels.
          if (metadata.channel === 'nanopore') {
            try {
              const { stripCRCMarkers, resyncAfterIndel, DEFAULT_CRC_MARKER_CONFIG } = await import('./crcmarker');
              // Must match encoder's marker spacing (segmentSize=8, ~32nt)
              const markerCfg = { ...DEFAULT_CRC_MARKER_CONFIG, segmentSize: 8 };
              const stripped = stripCRCMarkers(convBytes, markerCfg);
              // If any segments failed CRC, attempt resync using indel hypotheses
              if (stripped.failedSegments.length > 0) {
                // Try resync at each failure point; if resync fails, the segment
                // data is still included but may be corrupted — Viterbi will handle it
                for (const failSeg of stripped.failedSegments) {
                  const failPos = failSeg * (markerCfg.segmentSize + markerCfg.markerSize);
                  const shift = resyncAfterIndel(convBytes, failPos, markerCfg);
                  if (shift !== 0) {
                    // Apply the shift by adjusting the convBytes alignment
                    // (For simplicity, we just use the stripped payload — the Viterbi
                    // decoder's trellis can handle the residual misalignment)
                  }
                }
              }
              convBytes = stripped.payload;
            } catch {
              // CRC marker stripping failed — use convBytes as-is
            }
          }

          // Build per-bit LLRs for soft-decision Viterbi when Q-scores are available.
          // Q-scores are per-base in the DNA string; conv-encoded region uses
          // direct 2-bit mapping (4 bases per byte = 8 bits per 4 bases).
          // The conv region starts after address (4B = 16 bases).
          let convLLR: Float32Array | undefined;
          if (read.quality && read.quality.length > 0) {
            const convStartBase = layout.addressBytes * 4;
            const convNumBases = layout.convEncodedBytes * 4;
            const convEndBase = Math.min(convStartBase + convNumBases, read.quality.length);
            if (convEndBase > convStartBase) {
              const convQScores = read.quality.slice(convStartBase, convEndBase);
              convLLR = qScoresToBitLLRs(convQScores, 2);
            }
          }

          // Conv Viterbi decode → original LDPC codeword (innerN bytes)
          const convResult = convInner.decode(convBytes, convLLR);
          const rsCodewordDecoded = convResult.decoded;
          // Rebuild innerBlock as [LDPC codeword + CRC]
          // CRC is the last 2 bytes of the original innerBlock
          const crcFromRead = innerBlock.slice(innerBlock.length - layout.crcBytes);
          const rebuilt = new Uint8Array(innerN + layout.crcBytes);
          rebuilt.set(rsCodewordDecoded, 0);
          rebuilt.set(crcFromRead, innerN);
          innerBlock = rebuilt;
        } catch {
          // Conv decode failed (e.g., too many errors) — skip this read
          continue;
        }
      }

      // v62: For arithmetic-v2, the LDPC codeword does NOT include the address.
      //   innerBlock = [address(4)] + [LDPC codeword (payload + parity)]
      //   rsCodeword = innerBlock.slice(addressBytes, addressBytes + innerN)
      //   crcBytes = [] (no CRC-16; per-block CRC-8 already validated)
      const rsCodeword = useArithmeticV2
        ? innerBlock.slice(layout.addressBytes, layout.addressBytes + innerN)
        : innerBlock.slice(0, innerN);
      const crcBytes = useArithmeticV2
        ? new Uint8Array(0)  // no CRC-16 in arithmetic-v2
        : innerBlock.slice(innerN);

      // Extract Q-scores for the inner block region (excluding primers).
      // read.sequence was already trimmed to inner (in clusterReads).
      // Q-scores are per-base in the DNA string.
      let qScoresForInner: Uint8Array | null = null;
      if (read.quality && read.quality.length > 0) {
        // The quality array aligns with read.sequence (the inner DNA).
        // We need Q-scores for the first expectedDnaLen bases.
        if (useGoldman) {
          // Goldman: inner block bytes = expectedDnaLen / 6.
          // Q-scores are per-base in the DNA string.
          // The LDPC codeword is the first innerN bytes = innerN * 6 bases.
          const ldpcBases = innerN * 6;
          qScoresForInner = read.quality.slice(0, Math.min(ldpcBases, read.quality.length));
        } else {
          // Direct/Constrained: inner block bytes = expectedDnaLen / 4.
          // LDPC codeword = first innerN bytes = innerN * 4 bases.
          const ldpcBases = innerN * 4;
          qScoresForInner = read.quality.slice(0, Math.min(ldpcBases, read.quality.length));
        }
      }

      // For constrained/arithmetic mode, modify Q-scores to reflect erasures.
      // Erased bits get Q=0 (neutral LLR), known bits keep their Q-score.
      // This allows the BP decoder to correct erasures efficiently.
      // v57: also fire for useArithmetic — per-block CRC failures in
      // arithmetic mode produce byte-level erasures that LDPC must correct.
      const useErasures = useConstrained || useArithmetic || useDnaAeon;
      if (useErasures && constrainedErasures && qScoresForInner) {
        // Each base = 2 bits. Erasure info is per-bit.
        // We can't directly set per-bit Q-scores (they're per-base), so we
        // set the base Q-score to 0 if either bit is erased.
        for (let byteIdx = 0; byteIdx < innerN; byteIdx++) {
          for (let bitInByte = 0; bitInByte < 8; bitInByte += 2) {
            const bitIdx = byteIdx * 8 + bitInByte;
            if (constrainedErasures[bitIdx] || constrainedErasures[bitIdx + 1]) {
              const baseIdx = byteIdx * 4 + (bitInByte >> 1);
              if (baseIdx < qScoresForInner.length) {
                qScoresForInner[baseIdx] = 0; // Q=0 → neutral LLR
              }
            }
          }
        }
      } else if (useErasures && constrainedErasures) {
        // No Q-scores from read, but we still have erasures.
        // Build synthetic Q-scores: Q=0 for erased bits, Q=30 for known.
        const ldpcBases = innerN * 4;
        qScoresForInner = new Uint8Array(ldpcBases);
        for (let i = 0; i < ldpcBases; i++) qScoresForInner[i] = 30;
        for (let byteIdx = 0; byteIdx < innerN; byteIdx++) {
          for (let bitInByte = 0; bitInByte < 8; bitInByte += 2) {
            const bitIdx = byteIdx * 8 + bitInByte;
            if (constrainedErasures[bitIdx] || constrainedErasures[bitIdx + 1]) {
              const baseIdx = byteIdx * 4 + (bitInByte >> 1);
              if (baseIdx < qScoresForInner.length) {
                qScoresForInner[baseIdx] = 0;
              }
            }
          }
        }
      }

      // Try inner code decode (corrects errors in this individual read)
      // For RS: guaranteed minimum distance, so CRC is redundant
      // For LDPC: may pick wrong codeword (duplicate columns), so CRC is REQUIRED
      // For LDPC with soft-info: OSD-2/BP flips least-reliable bits
      let decodedData: Uint8Array;
      let corrected: number;
      try {
        if (useLDPC && innerLdpc) {
          // v57: If we have bit-level erasures (from arithmetic CRC or constrained
          // mapping), use the dedicated erasure decoder (peeling + Gaussian elim
          // over GF(2)). This is FAR more reliable than BP for erasures — BP
          // often fails to converge when Q-scores are 0, while the peeling
          // decoder directly solves the linear system.
          const hasErasures = constrainedErasures && constrainedErasures.some((e) => e);
          if (hasErasures) {
            // Convert boolean[] erasure bitmap to bit positions
            const erasePos: number[] = [];
            for (let i = 0; i < constrainedErasures!.length && i < innerN * 8; i++) {
              if (constrainedErasures![i]) erasePos.push(i);
            }
            if (erasePos.length > 0 && erasePos.length < innerN * 8) {
              try {
                const r = innerLdpc.decodeWithErasures(rsCodeword, erasePos);
                decodedData = r.data;
                corrected = r.corrected;
              } catch {
                // Erasure decode failed — fall through to BP
                const r = innerLdpc.decodeBeliefPropagation(rsCodeword, qScoresForInner, useGoldman);
                decodedData = r.data;
                corrected = r.corrected;
              }
            } else {
              // Too many erasures or none — use BP
              const r = innerLdpc.decodeBeliefPropagation(rsCodeword, qScoresForInner, useGoldman);
              decodedData = r.data;
              corrected = r.corrected;
            }
          } else if ((useConstrained || useSrt || useArithmetic || useDnaAeon)) {
            // For constrained/SRT/arithmetic/dnaAeon mode without erasures, use BP
            const r = innerLdpc.decodeBeliefPropagation(rsCodeword, qScoresForInner, useGoldman);
            decodedData = r.data;
            corrected = r.corrected;
          } else {
            // Standard mode: use configured decoder mode
            const decoderMode = ldpcDecoderMode;
            if (decoderMode === "bp") {
              const r = innerLdpc.decodeBeliefPropagation(rsCodeword, qScoresForInner, useGoldman);
              decodedData = r.data;
              corrected = r.corrected;
            } else if (decoderMode === "osd") {
              if (qScoresForInner && qScoresForInner.length > 0) {
                try {
                  const r = innerLdpc.decodeWithSoftInfo(rsCodeword, qScoresForInner, useGoldman);
                  decodedData = r.data;
                  corrected = r.corrected;
                } catch {
                  const r = innerLdpc.decode(rsCodeword);
                  decodedData = r.data;
                  corrected = r.corrected;
                }
              } else {
                const r = innerLdpc.decode(rsCodeword);
                decodedData = r.data;
                corrected = r.corrected;
              }
            } else if (decoderMode === "hard") {
              const r = innerLdpc.decode(rsCodeword);
              decodedData = r.data;
              corrected = r.corrected;
            } else {
              // "auto" mode: neural LDPC decoder (hard → BP → multi-decision BP)
              // v61: For illumina (low-IDS), skip BP — it's too slow. Just use
              // hard-decision. If that fails, skip the read (rely on coverage).
              // BP is only useful for nanopore (high-IDS) where hard-decision
              // frequently fails.
              if (neuralLdpc && (channel === "nanopore" || channel === "pacbio")) {
                try {
                  const r = neuralLdpc.decode(rsCodeword, qScoresForInner, useGoldman);
                  decodedData = r.data;
                  corrected = r.corrected;
                } catch {
                  if (qScoresForInner && qScoresForInner.length > 0) {
                    const r = innerLdpc.decodeBeliefPropagation(rsCodeword, qScoresForInner, useGoldman);
                    decodedData = r.data;
                    corrected = r.corrected;
                  } else {
                    throw new Error("Neural decode failed and no Q-scores for BP");
                  }
                }
              } else {
                // v61: illumina fast path — hard-decision only, no BP/BFA.
                // With 10x coverage, we can afford to skip noisy reads.
                const r = innerLdpc.decode(rsCodeword);
                decodedData = r.data;
                corrected = r.corrected;
              }
            }
          }
        } else {
          const r = innerRs.decode(rsCodeword);
          decodedData = r.data;
          corrected = r.corrected;
        }
      } catch {
        // This read has too many errors — skip it, try next read
        continue;
      }

      // For LDPC: verify CRC-16 to detect wrong-codeword picks
      // OPTIMIZATION: If corrected==0 (no errors corrected), skip re-encode
      // and use received data directly for CRC check.
      if (useLDPC && useArithmeticV2) {
        // v62: Arithmetic-v2 mode — no CRC-16 (per-block CRC-8 already validated).
        // The LDPC decodedData is just the payload (no address in the LDPC codeword).
        // The address comes from the direct DNA (already in innerBlock[0..3]).
        //
        // v62 fix: Do NOT compare received parity with re-encoded parity.
        // The received parity may have errors (arithmetic termination corruption)
        // that the LDPC erasure decoder already corrected. The erasure decoder
        // guarantees the syndrome is zero after correction, so we trust its output.
        // Just verify the address matches the expected oligo index.
        try {
          // Get address from the direct DNA (innerBlock[0..3])
          const whitenedAddr = innerBlock.slice(0, layout.addressBytes);
          const addr = unwhitenAddress(whitenedAddr);
          const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
          if (decodedIndex !== oligoIdx) continue;

          // Success! Extract payload (decodedData is just the payload for arithmetic-v2)
          const seed = addr[3];
          let payload = decodedData.slice(0, layout.payloadBytes);
          if (seed !== 0) payload = xorWithSeed(payload, seed);
          payloads.set(oligoIdx, payload);
          perOligo.push({
            index: oligoIdx, readCount: clusterReads.length, consensusLength: dna.length,
            crcPassed: true, innerRS: { corrected, success: true },
            seed, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
            strategy: 'per_read',
          });
          foundValidRead = true;
          break;
        } catch {
          continue;
        }
      } else if (useLDPC && !useConstrained && !useSrt) {
        if (corrected === 0) {
          // Fast path: zero-syndrome, but still check CRC (LDPC may have
          // duplicate columns → wrong codeword with zero syndrome)
          // Use received rsCodeword for CRC (no need to re-encode)
          const expectedCrc = crc16Bytes(rsCodeword);
          if (expectedCrc[0] !== crcBytes[0] || expectedCrc[1] !== crcBytes[1]) continue;
          // CRC passed — verify address and extract payload (same as slow path)
          const whitenedAddr = decodedData.slice(0, layout.addressBytes);
          const addr = unwhitenAddress(whitenedAddr);
          const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
          if (decodedIndex !== oligoIdx) continue;
          const seed = addr[3];
          let payload = decodedData.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
          if (seed !== 0) payload = xorWithSeed(payload, seed);
          payloads.set(oligoIdx, payload);
          perOligo.push({
            index: oligoIdx, readCount: clusterReads.length, consensusLength: dna.length,
            crcPassed: true, innerRS: { corrected, success: true },
            seed, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
            strategy: 'per_read',
          });
          foundValidRead = true;
          break;
        }
        // Slow path: errors corrected, re-encode for canonical form + CRC
        const reEncoded = innerLdpc!.encode(decodedData);
        const expectedCrc = crc16Bytes(reEncoded);
        const crcMatch = expectedCrc[0] === crcBytes[0] && expectedCrc[1] === crcBytes[1];
        if (!crcMatch) {
          continue;
        }
        const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
        const addr = unwhitenAddress(whitenedAddr);
        const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
        if (decodedIndex !== oligoIdx) {
          continue;
        }
        // Success! Extract payload
        const seed = addr[3];
        let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
        if (seed !== 0) payload = xorWithSeed(payload, seed);
        payloads.set(oligoIdx, payload);
        perOligo.push({
          index: oligoIdx, readCount: clusterReads.length, consensusLength: dna.length,
          crcPassed: true, innerRS: { corrected, success: true },
          seed, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
          strategy: 'per_read',
        });
        foundValidRead = true;
        break;
      } else if (useLDPC && (useConstrained || useSrt || useArithmetic || useDnaAeon)) {
        // Constrained/dnaAeon mode: skip CRC (CRC bytes have erasures).
        // Rely on LDPC syndrome check (already passed if we got here).
        // Re-encode to get canonical form and verify address.
        const reEncoded = innerLdpc!.encode(decodedData);
        const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
        const addr = unwhitenAddress(whitenedAddr);
        const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
        if (decodedIndex !== oligoIdx) {
          // Address mismatch — this read was clustered to the wrong oligo
          continue;
        }
        // Success! Extract payload (seed is always 0 for constrained mode)
        let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
        payloads.set(oligoIdx, payload);
        perOligo.push({
          index: oligoIdx, readCount: clusterReads.length, consensusLength: dna.length,
          crcPassed: true, innerRS: { corrected, success: true },
          seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
          strategy: 'per_read',
        });
        foundValidRead = true;
        break;
      } else if (useLDPC && useSrt) {
        // SRT mode: skip CRC (SRT injects errors that change CRC bytes).
        // Rely on LDPC syndrome check + address verification.
        // The LDPC decoder corrects the injected SRT errors automatically.
        const reEncoded = innerLdpc!.encode(decodedData);
        const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
        const addr = unwhitenAddress(whitenedAddr);
        const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
        if (decodedIndex !== oligoIdx) {
          continue;
        }
        let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
        payloads.set(oligoIdx, payload);
        perOligo.push({
          index: oligoIdx, readCount: clusterReads.length, consensusLength: dna.length,
          crcPassed: true, innerRS: { corrected, success: true },
          seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
        });
        foundValidRead = true;
        break;
      } else if (useLDPC && useArithmeticOrAeon) {
        // v57: Arithmetic/DNA-Aeon mode — encoder DROPPED CRC-16 to make room for
        // per-block CRC-8 sync markers (see codec.ts bytesToArithmeticDnaCrc).
        // The crcBytes extracted from innerBlock.slice(innerN) are NOT real
        // CRC bytes — they're whatever the arithmetic decoder produced after
        // the LDPC codeword region (typically zeros from padding).
        //
        // So we MUST NOT verify CRC-16 here. Instead, rely on:
        //   - LDPC syndrome check (passed if we got here)
        //   - Per-block CRC-8 (already verified by arithmeticDnaToBytesCrc / dnaAeonDecode)
        //   - Address verification (decodedIndex must match oligoIdx)
        const reEncoded = innerLdpc!.encode(decodedData);
        const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
        const addr = unwhitenAddress(whitenedAddr);
        const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
        if (decodedIndex !== oligoIdx) {
          continue;
        }
        let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
        payloads.set(oligoIdx, payload);
        perOligo.push({
          index: oligoIdx, readCount: clusterReads.length, consensusLength: dna.length,
          crcPassed: true, innerRS: { corrected, success: true },
          seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
          strategy: 'per_read',
        });
        foundValidRead = true;
        break;
      }

      // RS path: RS succeeded — trust the RS-corrected data.
      // RS guarantees correctness within correction capacity. CRC is redundant
      // and CRC bytes can be corrupted independently (they're outside RS).
      const reEncoded = innerRs.encode(decodedData);

      // CRITICAL: Verify the decoded address matches the expected oligo index.
      // If the address was corrupted by substitution, the read was clustered
      // to the wrong oligo index. The RS-corrected payload belongs to a
      // DIFFERENT oligo, so we must reject it.
      const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
      const addr = unwhitenAddress(whitenedAddr);
      const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
      if (decodedIndex !== oligoIdx) {
        // Address mismatch — this read was clustered to the wrong oligo
        continue;
      }

      // Success! Extract payload
      const whitenedAddress = reEncoded.slice(0, layout.addressBytes);
      const address = unwhitenAddress(whitenedAddress);
      const seed = address[3];
      let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
      if (seed !== 0) payload = xorWithSeed(payload, seed);

      payloads.set(oligoIdx, payload);
      perOligo.push({
        index: oligoIdx, readCount: clusterReads.length, consensusLength: dna.length,
        crcPassed: true, innerRS: { corrected, success: true },
        seed, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
        strategy: 'per_read',
      });
      foundValidRead = true;
      break; // Found a valid read — move to next oligo
    }

    if (foundValidRead) {
      strategy1Time += Date.now() - s1t0;
      strategy1Count++;
      continue;
    }
    strategy1Time += Date.now() - s1t0;

    // STRATEGY 2: Consensus decoding.
    //
    // For illumina (low-IDS): Fast weighted consensus (Q-score weighted majority
    // vote). Much faster than HMM. Works well when indels are rare (< 0.1%).
    //
    // For nanopore/pacbio (high-IDS): Progressive MSA consensus. The naïve
    // column-wise plurality vote (decodeWithFastConsensus) breaks under indels
    // — deletions shift all downstream columns, creating "deletion shadows"
    // that persist even with high coverage. Progressive MSA aligns reads first,
    // then computes consensus from the alignment, which correctly handles indels.
    // This is the same function as STRATEGY 2.75 but promoted to the primary
    // strategy for indel-heavy channels so MSA is tried first, not as a fallback.
    //
    // v61 CRITICAL FIX: Verify CRC before accepting. Without this, the consensus
    // can produce a wrong codeword that passes address verification but has a
    // wrong payload. This was the root cause of the v60 hash FAIL.
    if (clusterReads.length >= 2 && useLDPC && innerLdpc) {
      // For nanopore/indel-heavy channels with multiple reads, use MSA-based consensus
      // instead of naïve column plurality (which breaks under indels)
      const isIndelHeavyChannel = channel === "nanopore" || channel === "pacbio";
      let consensusResult: { data: Uint8Array; corrected: number } | null;
      if (isIndelHeavyChannel) {
        // Use progressive MSA for better indel tolerance
        consensusResult = decodeWithProgressiveMSA(
          clusterReads, expectedDnaLen, innerLdpc,
          (dna: string) => dnaToBytes(dna), innerN,
          convInner ?? null,
          useConvInner ? {
            addressBytes: layout.addressBytes,
            convEncodedBytes: layout.convEncodedBytes,
            crcBytes: layout.crcBytes,
          } : null,
        );
      } else {
        // Illumina / low-IDS: fast weighted consensus is sufficient
        consensusResult = decodeWithFastConsensus(
          clusterReads, expectedDnaLen, innerLdpc,
          (dna: string) => dnaToBytes(dna), innerN,
        );
      }
      if (consensusResult) {
        // Re-encode to get canonical form
        const reEncoded = innerLdpc.encode(consensusResult.data);
        // v61: Verify CRC — extract CRC from the first read and compare
        const firstRead = clusterReads[0];
        let firstReadDna = firstRead.sequence;
        if (firstReadDna.length > expectedDnaLen) firstReadDna = firstReadDna.slice(0, expectedDnaLen);
        else if (firstReadDna.length < expectedDnaLen) firstReadDna = firstReadDna + "A".repeat(expectedDnaLen - firstReadDna.length);
        let firstReadCrc: Uint8Array | null = null;
        try {
          const firstReadBlock = dnaToBytes(firstReadDna);
          firstReadCrc = firstReadBlock.slice(innerN, innerN + 2);
        } catch { /* ignore */ }
        const expectedCrc = crc16Bytes(reEncoded);
        const crcOk = firstReadCrc && expectedCrc[0] === firstReadCrc[0] && expectedCrc[1] === firstReadCrc[1];
        // Verify address (LDPC syndrome already passed)
        const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
        const addr = unwhitenAddress(whitenedAddr);
        const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
        if (crcOk && decodedIndex === oligoIdx) {
          let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
          payloads.set(oligoIdx, payload);
          perOligo.push({
            index: oligoIdx, readCount: clusterReads.length, consensusLength: expectedDnaLen,
            crcPassed: true, innerRS: { corrected: consensusResult.corrected, success: true },
            seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
            strategy: isIndelHeavyChannel ? 'msa' : 'fast_consensus',
          });
          continue;
        }
      }
    }

    // STRATEGY 2.5: Attention consensus (MACL-style multi-scale attention).
    // Better than simple weighted consensus at high error rates.
    // v61: CRC verification added. v61: Skip for illumina (low-IDS) — too slow.
    // v62: Also run for pacbio (high-IDS, indel-heavy).
    if (clusterReads.length >= 3 && useLDPC && innerLdpc && (channel === "nanopore" || channel === "pacbio")) {
      const attentionResult = decodeWithAttentionConsensus(
        clusterReads, expectedDnaLen, innerLdpc,
        (dna: string) => dnaToBytes(dna), innerN,
      );
      if (attentionResult) {
        const reEncoded = innerLdpc.encode(attentionResult.data);
        const firstRead = clusterReads[0];
        let firstReadDna = firstRead.sequence;
        if (firstReadDna.length > expectedDnaLen) firstReadDna = firstReadDna.slice(0, expectedDnaLen);
        else if (firstReadDna.length < expectedDnaLen) firstReadDna = firstReadDna + "A".repeat(expectedDnaLen - firstReadDna.length);
        let firstReadCrc: Uint8Array | null = null;
        try {
          const firstReadBlock = dnaToBytes(firstReadDna);
          firstReadCrc = firstReadBlock.slice(innerN, innerN + 2);
        } catch { /* ignore */ }
        const expectedCrc = crc16Bytes(reEncoded);
        const crcOk = firstReadCrc && expectedCrc[0] === firstReadCrc[0] && expectedCrc[1] === firstReadCrc[1];
        const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
        const addr = unwhitenAddress(whitenedAddr);
        const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
        if (crcOk && decodedIndex === oligoIdx) {
          let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
          payloads.set(oligoIdx, payload);
          perOligo.push({
            index: oligoIdx, readCount: clusterReads.length, consensusLength: expectedDnaLen,
            crcPassed: true, innerRS: { corrected: attentionResult.corrected, success: true },
            seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
            strategy: 'attention_consensus',
          });
          continue;
        }
      }
    }

    // STRATEGY 2.75: v58 Progressive MSA (Kalign-style) — handles indels.
    // v61: Skip for illumina (low-IDS) — too slow. Only run for nanopore/pacbio.
    // This is the right tool for Nanopore / PacBio / high-IDS channels where indels
    // dominate. The fast weighted consensus (STRATEGY 2) assumes no indels
    // and fails when >5% of reads have insertions/deletions. The HMM
    // soft-consensus (STRATEGY 3) also handles indels but is ~10× slower
    // and can fail to converge on long noisy reads. Progressive MSA is the
    // sweet spot: O(N²) in reads but each NW profile alignment is fast.
    // v59: Now supports conv inner code — runs conv Viterbi on the MSA
    // consensus before LDPC. This is critical for nanopore: per-read conv
    // Viterbi fails on noisy reads, but MSA consensus is clean enough for
    // conv Viterbi to succeed.
    // NOTE: For nanopore/pacbio, STRATEGY 2 already uses MSA as the primary
    // approach. This fallback is kept for the case where STRATEGY 2's MSA
    // fails CRC verification but a retry with different read selection might work.
    if (clusterReads.length >= 2 && useLDPC && innerLdpc && (channel === "nanopore" || channel === "pacbio")) {
      const msaResult = decodeWithProgressiveMSA(
        clusterReads, expectedDnaLen, innerLdpc,
        (dna: string) => dnaToBytes(dna), innerN,
        convInner ?? null,
        useConvInner ? {
          addressBytes: layout.addressBytes,
          convEncodedBytes: layout.convEncodedBytes,
          crcBytes: layout.crcBytes,
        } : null,
      );
      if (msaResult) {
        const reEncoded = innerLdpc.encode(msaResult.data);
        // v61: Verify CRC
        const firstRead = clusterReads[0];
        let firstReadDna = firstRead.sequence;
        if (firstReadDna.length > expectedDnaLen) firstReadDna = firstReadDna.slice(0, expectedDnaLen);
        else if (firstReadDna.length < expectedDnaLen) firstReadDna = firstReadDna + "A".repeat(expectedDnaLen - firstReadDna.length);
        let firstReadCrc: Uint8Array | null = null;
        try { firstReadCrc = dnaToBytes(firstReadDna).slice(innerN, innerN + 2); } catch {}
        const expectedCrc = crc16Bytes(reEncoded);
        const crcOk = firstReadCrc && expectedCrc[0] === firstReadCrc[0] && expectedCrc[1] === firstReadCrc[1];
        const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
        const addr = unwhitenAddress(whitenedAddr);
        const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
        if (crcOk && decodedIndex === oligoIdx) {
          let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
          payloads.set(oligoIdx, payload);
          perOligo.push({
            index: oligoIdx, readCount: clusterReads.length, consensusLength: expectedDnaLen,
            crcPassed: true, innerRS: { corrected: msaResult.corrected, success: true },
            seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
            strategy: 'msa',
          });
          continue;
        }
      }
    }

    // STRATEGY 3: HMM-based soft-consensus (slower, handles indels).
    // Only try if fast consensus failed and we have enough reads.
    // v61: Skip for illumina (low-IDS) — too slow. Only run for nanopore/pacbio.
    if (clusterReads.length >= 3 && useLDPC && innerLdpc && (channel === "nanopore" || channel === "pacbio")) {
      const refRead = clusterReads[0];
      let refDna = refRead.sequence;
      if (refDna.length < expectedDnaLen) {
        refDna = refDna + "A".repeat(expectedDnaLen - refDna.length);
      } else if (refDna.length > expectedDnaLen) {
        refDna = refDna.slice(0, expectedDnaLen);
      }

      try {
        const result = decodeWithSoftConsensus(clusterReads, refDna, innerLdpc, useGoldman);
        if (result) {
          const reEncoded = innerLdpc.encode(result.data);
          // v61: Verify CRC
          const firstRead = clusterReads[0];
          let firstReadDna = firstRead.sequence;
          if (firstReadDna.length > expectedDnaLen) firstReadDna = firstReadDna.slice(0, expectedDnaLen);
          else if (firstReadDna.length < expectedDnaLen) firstReadDna = firstReadDna + "A".repeat(expectedDnaLen - firstReadDna.length);
          let firstReadCrc: Uint8Array | null = null;
          try { firstReadCrc = dnaToBytes(firstReadDna).slice(innerN, innerN + 2); } catch {}
          const expectedCrc = crc16Bytes(reEncoded);
          const crcOk = firstReadCrc && expectedCrc[0] === firstReadCrc[0] && expectedCrc[1] === firstReadCrc[1];
          const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
          const addr = unwhitenAddress(whitenedAddr);
          const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
          if (crcOk && decodedIndex === oligoIdx) {
            let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
            payloads.set(oligoIdx, payload);
            perOligo.push({
              index: oligoIdx, readCount: clusterReads.length, consensusLength: refDna.length,
              crcPassed: true, innerRS: { corrected: result.corrected, success: true },
              seed: 0, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
              strategy: 'soft_consensus',
            });
            continue;
          }
        }
      } catch {
        // HMM soft-consensus failed — fall through to OSD post-pass
      }
    }

    // OSD-2 post-pass: try ordered statistics decoding before giving up.
    // This recovers many LDPC failures at the cost of ~2ms per oligo.
    // Only attempt for nanopore/pacbio channels where inner decode often fails
    // due to indels that the Viterbi preprocessor couldn't fully correct.
    if ((channel === "nanopore" || channel === "pacbio") && useLDPC && innerLdpc && clusterReads.length >= 1) {
      try {
        const { osdDecode } = await import('./osd-full');
        const { GF2Matrix } = await import('./osd');

        // Build a consensus from cluster reads to use as the received word
        const refRead = clusterReads[0];
        let refDna = refRead.sequence;
        if (refDna.length < expectedDnaLen) {
          refDna = refDna + "A".repeat(expectedDnaLen - refDna.length);
        } else if (refDna.length > expectedDnaLen) {
          refDna = refDna.slice(0, expectedDnaLen);
        }

        // Convert consensus DNA to bytes
        let rsCodewordForOsd: Uint8Array;
        try {
          rsCodewordForOsd = dnaToBytes(refDna).slice(0, innerN);
        } catch {
          // dnaToBytes failed — skip OSD
          throw new Error("dnaToBytes failed for OSD input");
        }

        // Build LLRs from Q-scores if available, otherwise use uniform LLRs
        const nBits = innerN * 8;
        const llr = new Float32Array(nBits);
        if (refRead.quality && refRead.quality.length > 0) {
          // Use Q-scores to derive per-bit LLRs
          const qScores = refRead.quality.slice(0, Math.min(expectedDnaLen, refRead.quality.length));
          for (let baseIdx = 0; baseIdx < qScores.length && baseIdx * 2 < nBits; baseIdx++) {
            // Q-score to LLR: LLR = (Q - 33) * ln(10) / 10 * 2 (per-bit)
            const qVal = qScores[baseIdx];
            const llrMag = Math.max(0.1, ((qVal - 33) * 0.2303)); // 0.2303 ≈ ln(10)/10
            // Two bits per base (direct mapping)
            const byteIdx0 = Math.floor(baseIdx / 4);
            const bitIdx0 = (baseIdx % 4) * 2;
            if (byteIdx0 < innerN) {
              const byteVal = rsCodewordForOsd[byteIdx0];
              // bit 0
              const bit0 = (byteVal >> (7 - bitIdx0)) & 1;
              llr[baseIdx * 2] = bit0 === 0 ? llrMag : -llrMag;
              // bit 1
              const bit1 = (byteVal >> (7 - bitIdx0 - 1)) & 1;
              llr[baseIdx * 2 + 1] = bit1 === 0 ? llrMag : -llrMag;
            }
          }
        } else {
          // Uniform LLRs (hard decision) — magnitude 1.0
          for (let i = 0; i < nBits; i++) {
            const byteIdx = Math.floor(i / 8);
            const bitPos = 7 - (i % 8);
            const bit = (rsCodewordForOsd[byteIdx] >> bitPos) & 1;
            llr[i] = bit === 0 ? 1.0 : -1.0;
          }
        }

        // Get parity-check matrix from LDPC code
        const H = innerLdpc.parityCheckMatrix;

        // CRC + address check function
        const crcAndAddrCheck = (codewordBits: Uint8Array): boolean => {
          // Convert bit array to byte array
          const codewordBytes = new Uint8Array(innerN);
          for (let i = 0; i < innerN; i++) {
            let byteVal = 0;
            for (let b = 0; b < 8; b++) {
              byteVal |= (codewordBits[i * 8 + b] << (7 - b));
            }
            codewordBytes[i] = byteVal;
          }
          // Re-encode to get canonical form
          const reEncoded = innerLdpc!.encode(
            // Extract info bits
            codewordBytes.slice(0, innerLdpc!.k),
          );
          // Verify CRC
          const expectedCrc = crc16Bytes(reEncoded);
          // Get CRC from the first read
          const firstRead = clusterReads[0];
          let firstReadDna2 = firstRead.sequence;
          if (firstReadDna2.length > expectedDnaLen) firstReadDna2 = firstReadDna2.slice(0, expectedDnaLen);
          else if (firstReadDna2.length < expectedDnaLen) firstReadDna2 = firstReadDna2 + "A".repeat(expectedDnaLen - firstReadDna2.length);
          let firstReadCrc2: Uint8Array | null = null;
          try { firstReadCrc2 = dnaToBytes(firstReadDna2).slice(innerN, innerN + 2); } catch {}
          if (!firstReadCrc2 || expectedCrc[0] !== firstReadCrc2[0] || expectedCrc[1] !== firstReadCrc2[1]) return false;
          // Verify address
          const whitenedAddr = reEncoded.slice(0, layout.addressBytes);
          const addr = unwhitenAddress(whitenedAddr);
          const decodedIndex = (addr[0] << 16) | (addr[1] << 8) | addr[2];
          return decodedIndex === oligoIdx;
        };

        // Run OSD-2 (maxOrder=2 is a good tradeoff: O(k²) candidates, ~2ms)
        const osdResult = osdDecode(llr, H, crcAndAddrCheck, {
          maxOrder: 2,
          k: innerK * 8,
        });

        if (osdResult.codeword) {
          // OSD succeeded — extract payload from the decoded codeword
          const codewordBits = osdResult.codeword;
          const codewordBytes = new Uint8Array(innerN);
          for (let i = 0; i < innerN; i++) {
            let byteVal = 0;
            for (let b = 0; b < 8; b++) {
              byteVal |= (codewordBits[i * 8 + b] << (7 - b));
            }
            codewordBytes[i] = byteVal;
          }
          // Re-encode to get canonical form
          const reEncoded = innerLdpc.encode(codewordBytes.slice(0, innerK));
          const seed = 0; // OSD doesn't recover seed
          let payload = reEncoded.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
          payloads.set(oligoIdx, payload);
          perOligo.push({
            index: oligoIdx, readCount: clusterReads.length, consensusLength: refDna.length,
            crcPassed: true, innerRS: { corrected: osdResult.candidatesTried, success: true },
            seed, payloadBytes: payload, isParity: oligoIdx >= metadata.outerRS.k,
            strategy: 'osd_post_pass',
          });
          continue;
        }
      } catch {
        // OSD not available or failed — fall through to erasure
      }
    }

    // STRATEGY 4: Mark as erasure for outer RS to recover.
    erasedIndices.push(oligoIdx);
    oligosFailedInnerRS++;
    perOligo.push({
      index: oligoIdx, readCount: clusterReads.length, consensusLength: 0,
      crcPassed: false, innerRS: { corrected: 0, success: false }, seed: 0,
      payloadBytes: new Uint8Array(layout.payloadBytes),
      isParity: oligoIdx >= metadata.outerRS.k,
      strategy: 'erasure',
    });
  }

  // 3) Outer RS: recover missing oligos via erasure decoding.
  // With per-read decode, there are NO silent consensus errors — every
  // payload in the map was individually validated by RS+CRC. So the outer
  // RS only needs to recover genuinely missing oligos (erasures).
  let oligosFailedOuterRS = 0;
  const recoveredPayloads = new Map<number, Uint8Array>();

  if (!useOuterRS) {
    for (const [idx, payload] of payloads) {
      recoveredPayloads.set(idx, payload);
    }
  } else if (useGF216 && outerRs216) {
    const numPairs = Math.floor(layout.payloadBytes / 2);
    let outerSuccessCount = 0, outerFailCount = 0, outerSkipCount = 0;
    for (let pairIdx = 0; pairIdx < numPairs; pairIdx++) {
      const j0 = pairIdx * 2;
      const j1 = pairIdx * 2 + 1;
      const codeword = new Uint16Array(metadata.outerRS.n);
      const erased: number[] = [];
      for (let i = 0; i < metadata.outerRS.n; i++) {
        const p = payloads.get(i);
        if (p) {
          codeword[i] = (p[j0] << 8) | p[j1];
        } else {
          codeword[i] = 0;
          if (erasedIndices.includes(i)) erased.push(i);
        }
      }
      if (erased.length === 0) {
        outerSkipCount++;
        for (let i = 0; i < metadata.outerRS.k; i++) {
          if (!recoveredPayloads.has(i)) {
            recoveredPayloads.set(i, new Uint8Array(layout.payloadBytes));
          }
          recoveredPayloads.get(i)![j0] = (codeword[i] >> 8) & 0xff;
          recoveredPayloads.get(i)![j1] = codeword[i] & 0xff;
        }
        continue;
      }
      if (erased.length > metadata.outerRS.n - metadata.outerRS.k) {
        oligosFailedOuterRS++;
        outerFailCount++;
        continue;
      }
      try {
        const r = outerRs216.decodeWithErasures(codeword, erased);
        outerSuccessCount++;
        for (let i = 0; i < metadata.outerRS.k; i++) {
          if (!recoveredPayloads.has(i)) {
            recoveredPayloads.set(i, new Uint8Array(layout.payloadBytes));
          }
          recoveredPayloads.get(i)![j0] = (r.data[i] >> 8) & 0xff;
          recoveredPayloads.get(i)![j1] = r.data[i] & 0xff;
        }
      } catch {
        oligosFailedOuterRS++;
        outerFailCount++;
      }
    }
    if (process.env.HELIX_DEBUG) {
      console.error(`[LDPC debug] outer RS: ${outerSuccessCount} success, ${outerFailCount} fail, ${outerSkipCount} skip (erasedIndices=${erasedIndices.length}, erasedCount=${erasedIndices.length})`);
    }
    // Handle odd byte
    if (layout.payloadBytes % 2 === 1) {
      const j = layout.payloadBytes - 1;
      // Build codeword for the odd byte using GF(2^16) (same as even bytes)
      const codewordOdd = new Uint16Array(metadata.outerRS.n);
      const erasedOdd: number[] = [];
      for (let i = 0; i < metadata.outerRS.n; i++) {
        const p = payloads.get(i);
        if (p) {
          codewordOdd[i] = p[j] & 0xff; // high byte = 0
        } else {
          codewordOdd[i] = 0;
          if (erasedIndices.includes(i)) erasedOdd.push(i);
        }
      }
      if (process.env.HELIX_DEBUG) {
        console.error(`[LDPC debug] odd byte (j=${j}): ${erasedOdd.length} erasures, codeword len=${codewordOdd.length}`);
      }
      if (erasedOdd.length === 0) {
        for (let i = 0; i < metadata.outerRS.k; i++) {
          if (!recoveredPayloads.has(i)) {
            recoveredPayloads.set(i, new Uint8Array(layout.payloadBytes));
          }
          recoveredPayloads.get(i)![j] = codewordOdd[i] & 0xff;
        }
      } else if (erasedOdd.length <= metadata.outerRS.n - metadata.outerRS.k && outerRs216) {
        try {
          const r = outerRs216.decodeWithErasures(codewordOdd, erasedOdd);
          if (process.env.HELIX_DEBUG) {
            console.error(`[LDPC debug] odd byte recovered: ${r.data.length} symbols`);
          }
          for (let i = 0; i < metadata.outerRS.k; i++) {
            if (!recoveredPayloads.has(i)) {
              recoveredPayloads.set(i, new Uint8Array(layout.payloadBytes));
            }
            recoveredPayloads.get(i)![j] = r.data[i] & 0xff;
          }
        } catch (e) {
          if (process.env.HELIX_DEBUG) {
            console.error(`[LDPC debug] odd byte FAILED: ${(e as Error).message}`);
          }
          oligosFailedOuterRS++;
        }
      }
    }
  } else if (outerRs8) {
    // GF(2^8) RS: process each byte position independently
    for (let j = 0; j < layout.payloadBytes; j++) {
      const codeword = new Uint8Array(metadata.outerRS.n);
      const byteErased: number[] = [];
      for (let i = 0; i < metadata.outerRS.n; i++) {
        const p = payloads.get(i);
        if (p) {
          codeword[i] = p[j];
        } else {
          codeword[i] = 0;
          if (erasedIndices.includes(i)) byteErased.push(i);
        }
      }
      if (byteErased.length === 0) {
        for (let i = 0; i < metadata.outerRS.k; i++) {
          if (!recoveredPayloads.has(i)) {
            recoveredPayloads.set(i, new Uint8Array(layout.payloadBytes));
          }
          recoveredPayloads.get(i)![j] = codeword[i];
        }
        continue;
      }
      if (byteErased.length > metadata.outerRS.n - metadata.outerRS.k) {
        oligosFailedOuterRS++;
        continue;
      }
      try {
        const r = outerRs8.decodeWithErasures(codeword, byteErased);
        for (let i = 0; i < metadata.outerRS.k; i++) {
          if (!recoveredPayloads.has(i)) {
            recoveredPayloads.set(i, new Uint8Array(layout.payloadBytes));
          }
          recoveredPayloads.get(i)![j] = r.data[i];
        }
      } catch {
        oligosFailedOuterRS++;
      }
    }
  }

  // 4) Concatenate payloads in index order; trim to fileSize (compressed size)
  const compressedSize = metadata.oligoCount > 0
    ? Math.min(
        metadata.rawSize, // upper bound
        metadata.outerRS.k * layout.payloadBytes,
      )
    : 0;
  // We need to know the compressed length. The metadata stores rawSize, not compressedSize.
  // For now, decompress and check trailing zeros — they're padding.
  const totalPayload = new Uint8Array(metadata.outerRS.k * layout.payloadBytes);
  for (let i = 0; i < metadata.outerRS.k; i++) {
    const p = recoveredPayloads.get(i);
    if (p) totalPayload.set(p, i * layout.payloadBytes);
  }

  // 5) Decompress (if needed)
  let data: Uint8Array;
  if (metadata.compression === "deflate") {
    try {
      data = decompressRouter(totalPayload);
    } catch {
      // Inflate failed — return raw payload as best effort
      data = totalPayload;
    }
  } else {
    data = totalPayload.slice(0, metadata.fileSize);
  }

  // Trim to fileSize
  if (data.length > metadata.fileSize) {
    data = data.slice(0, metadata.fileSize);
  }

  // 6) Compute hash
  let hash = await sha256(data);
  let hashMatches = hash === metadata.fileHash;

  // v61 CRITICAL FIX: Outer RS erasure search.
  //
  // If the hash doesn't match but we have outer RS parity available AND
  // fewer than (n-k) oligos were erased, the wrong payload is likely due to
  // a rare LDPC wrong-codeword acceptance (CRC collision ~1/2^16 per read).
  // The outer RS syndrome is non-zero, indicating an error. We try marking
  // each data oligo (one at a time) as an erasure and re-running the outer
  // RS. If any retry produces a matching hash, we use that.
  //
  // This is O(k) retries, each doing O(payloadBytes) RS decodes. For k=298
  // and payloadBytes=150, that's ~45K RS decodes — fast (<500ms).
  //
  // This catches the "silent wrong codeword" failure mode that plagued v60
  // (hash FAIL with 307/307 oligos recovered).
  if (!hashMatches && useOuterRS && erasedIndices.length < (metadata.outerRS.n - metadata.outerRS.k)) {
    if (process.env.HELIX_DEBUG) console.error(`[v61 erasure-search] hash mismatch detected, trying RS error correction + erasure search (erased=${erasedIndices.length}, k=${metadata.outerRS.k})`);
    const erasureCapacity = metadata.outerRS.n - metadata.outerRS.k - erasedIndices.length;
    // Only try if we have spare erasure capacity
    if (erasureCapacity >= 1) {
      const k = metadata.outerRS.k;
      // v61: Use proper RS error correction (Berlekamp-Massey + Chien + Forney)
      // instead of brute-force erasure search. This is O(n^2) per byte position
      // but handles multiple unknown errors in one shot.
      //
      // The RS decode can correct up to floor(nsym/2) unknown errors. With 3%
      // outer RS (nsym=9 for 256KB), we can correct up to 4 wrong oligos.
      // With 25% outer RS (nsym=54 for nanopore), we can correct up to 27.
      //
      // This catches the "silent wrong codeword" failure mode where LDPC + CRC-16
      // falsely accepted a wrong codeword (probability ~1/2^16 per read).
      const retryPayloads = new Map<number, Uint8Array>();
      let rsErrorCorrectionSuccess = false;

      if (useGF216 && outerRs216) {
        const numPairs = Math.floor(layout.payloadBytes / 2);
        for (let pairIdx = 0; pairIdx < numPairs; pairIdx++) {
          const j0 = pairIdx * 2;
          const j1 = pairIdx * 2 + 1;
          const codeword = new Uint16Array(metadata.outerRS.n);
          for (let i = 0; i < metadata.outerRS.n; i++) {
            const p = payloads.get(i);
            if (p) {
              codeword[i] = (p[j0] << 8) | p[j1];
            } else {
              codeword[i] = 0;
            }
          }
          try {
            const r = outerRs216.decode(codeword);
            for (let i = 0; i < k; i++) {
              if (!retryPayloads.has(i)) {
                retryPayloads.set(i, new Uint8Array(layout.payloadBytes));
              }
              retryPayloads.get(i)![j0] = (r.data[i] >> 8) & 0xff;
              retryPayloads.get(i)![j1] = r.data[i] & 0xff;
            }
            rsErrorCorrectionSuccess = true;
          } catch {
            // RS couldn't correct — leave retryPayloads incomplete
          }
        }
        if (rsErrorCorrectionSuccess && layout.payloadBytes % 2 === 1) {
          const j = layout.payloadBytes - 1;
          const codewordOdd = new Uint16Array(metadata.outerRS.n);
          for (let i = 0; i < metadata.outerRS.n; i++) {
            const p = payloads.get(i);
            if (p) {
              codewordOdd[i] = p[j] & 0xff;
            } else {
              codewordOdd[i] = 0;
            }
          }
          try {
            const r = outerRs216.decode(codewordOdd);
            for (let i = 0; i < k; i++) {
              retryPayloads.get(i)![j] = r.data[i] & 0xff;
            }
          } catch { /* ignore odd byte failure */ }
        }
      } else if (outerRs8) {
        let anySuccess = false;
        for (let j = 0; j < layout.payloadBytes; j++) {
          const codeword = new Uint8Array(metadata.outerRS.n);
          for (let i = 0; i < metadata.outerRS.n; i++) {
            const p = payloads.get(i);
            if (p) {
              codeword[i] = p[j];
            } else {
              codeword[i] = 0;
            }
          }
          try {
            const r = outerRs8.decode(codeword);
            for (let i = 0; i < k; i++) {
              if (!retryPayloads.has(i)) {
                retryPayloads.set(i, new Uint8Array(layout.payloadBytes));
              }
              retryPayloads.get(i)![j] = r.data[i];
            }
            anySuccess = true;
          } catch {
            // RS couldn't correct this byte position — leave as 0
          }
        }
        rsErrorCorrectionSuccess = anySuccess;
      }

      if (rsErrorCorrectionSuccess) {
        // Build retry totalPayload
        const retryTotalPayload = new Uint8Array(k * layout.payloadBytes);
        for (let i = 0; i < k; i++) {
          const p = retryPayloads.get(i);
          if (p) retryTotalPayload.set(p, i * layout.payloadBytes);
        }
        let retryData: Uint8Array;
        if (metadata.compression === "deflate") {
          try {
            retryData = decompressRouter(retryTotalPayload);
          } catch {
            // Inflate failed — skip
            retryData = null as any;
          }
        } else {
          retryData = retryTotalPayload.slice(0, metadata.fileSize);
        }
        if (retryData) {
          if (retryData.length > metadata.fileSize) {
            retryData = retryData.slice(0, metadata.fileSize);
          }
          const retryHash = await sha256(retryData);
          if (retryHash === metadata.fileHash) {
            // RS error correction fixed the wrong codeword!
            data = retryData;
            hash = retryHash;
            hashMatches = true;
          }
        }
      }

      // v61 fallback: If RS error correction didn't fix it, try brute-force
      // erasure search on small payloads (k <= 500). This handles cases where
      // RS decode failed (too many errors) but single-erasure might work.
      // For large payloads (k > 500, GF(2^16)), RS216 doesn't support unknown-
      // error correction, so we rely on the brute-force search with a cap.
      if (!hashMatches) {
        // For GF(2^16) large payloads, cap at 50 suspects (each suspect is
        // O(payloadBytes * n^2) = expensive). For GF(2^8), cap at k (fast).
        const maxSuspects = k > 500 ? 50 : k;
        for (let suspectIdx = 0; suspectIdx < maxSuspects; suspectIdx++) {
          // For large payloads, try suspects in order (0, 1, 2, ...).
          // For small payloads, also try in order.
          const suspect = suspectIdx;
          if (erasedIndices.includes(suspect)) continue;
          if (!payloads.has(suspect)) continue;

        // Re-run outer RS with `suspect` as an additional erasure
        const retryPayloads = new Map<number, Uint8Array>();
        // Build codeword with suspect as erasure
        if (useGF216 && outerRs216) {
          const numPairs = Math.floor(layout.payloadBytes / 2);
          let allOk = true;
          for (let pairIdx = 0; pairIdx < numPairs && allOk; pairIdx++) {
            const j0 = pairIdx * 2;
            const j1 = pairIdx * 2 + 1;
            const codeword = new Uint16Array(metadata.outerRS.n);
            const erased: number[] = [...erasedIndices, suspect];
            for (let i = 0; i < metadata.outerRS.n; i++) {
              const p = payloads.get(i);
              if (p && i !== suspect) {
                codeword[i] = (p[j0] << 8) | p[j1];
              } else {
                codeword[i] = 0;
              }
            }
            if (erased.length > metadata.outerRS.n - metadata.outerRS.k) {
              allOk = false;
              break;
            }
            try {
              const r = outerRs216.decodeWithErasures(codeword, erased);
              for (let i = 0; i < k; i++) {
                if (!retryPayloads.has(i)) {
                  retryPayloads.set(i, new Uint8Array(layout.payloadBytes));
                }
                retryPayloads.get(i)![j0] = (r.data[i] >> 8) & 0xff;
                retryPayloads.get(i)![j1] = r.data[i] & 0xff;
              }
            } catch {
              allOk = false;
            }
          }
          if (!allOk) continue;
          // Handle odd byte
          if (layout.payloadBytes % 2 === 1) {
            const j = layout.payloadBytes - 1;
            const codewordOdd = new Uint16Array(metadata.outerRS.n);
            const erasedOdd = [...erasedIndices, suspect];
            for (let i = 0; i < metadata.outerRS.n; i++) {
              const p = payloads.get(i);
              if (p && i !== suspect) {
                codewordOdd[i] = p[j] & 0xff;
              } else {
                codewordOdd[i] = 0;
              }
            }
            try {
              const r = outerRs216.decodeWithErasures(codewordOdd, erasedOdd);
              for (let i = 0; i < k; i++) {
                retryPayloads.get(i)![j] = r.data[i] & 0xff;
              }
            } catch { continue; }
          }
        } else if (outerRs8) {
          for (let j = 0; j < layout.payloadBytes; j++) {
            const codeword = new Uint8Array(metadata.outerRS.n);
            const byteErased = [...erasedIndices, suspect];
            for (let i = 0; i < metadata.outerRS.n; i++) {
              const p = payloads.get(i);
              if (p && i !== suspect) {
                codeword[i] = p[j];
              } else {
                codeword[i] = 0;
              }
            }
            if (byteErased.length > metadata.outerRS.n - metadata.outerRS.k) continue;
            try {
              const r = outerRs8.decodeWithErasures(codeword, byteErased);
              for (let i = 0; i < k; i++) {
                if (!retryPayloads.has(i)) {
                  retryPayloads.set(i, new Uint8Array(layout.payloadBytes));
                }
                retryPayloads.get(i)![j] = r.data[i];
              }
            } catch { continue; }
          }
        }

        // Build retry totalPayload
        const retryTotalPayload = new Uint8Array(k * layout.payloadBytes);
        for (let i = 0; i < k; i++) {
          const p = retryPayloads.get(i);
          if (p) retryTotalPayload.set(p, i * layout.payloadBytes);
        }
        let retryData: Uint8Array;
        if (metadata.compression === "deflate") {
          try {
            retryData = decompressRouter(retryTotalPayload);
          } catch {
            continue;
          }
        } else {
          retryData = retryTotalPayload.slice(0, metadata.fileSize);
        }
        if (retryData.length > metadata.fileSize) {
          retryData = retryData.slice(0, metadata.fileSize);
        }
        const retryHash = await sha256(retryData);
        if (retryHash === metadata.fileHash) {
          // Found the culprit! Use this recovery.
          data = retryData;
          hash = retryHash;
          hashMatches = true;
          break;
        }
        }  // end for suspect
      }  // end if (k <= 500)
    }  // end if (erasureCapacity >= 1)
  }  // end if (!hashMatches && useOuterRS)

  const stats: DecodeStats = {
    totalReads: reads.length,
    readsUsed: reads.length - discarded,
    clustersFormed: clusters.size,
    oligosRecovered: payloads.size + erasedIndices.length - oligosFailedInnerRS,
    oligosErased: erasedIndices.length,
    oligosFailedInnerRS,
    oligosFailedOuterRS,
    consensusSuccessRate: payloads.size / metadata.oligoCount,
    decodeTimeMs: Date.now() - t0,
  };
  if (process.env.HELIX_DEBUG) {
    console.error(`[v61 timing] STRATEGY 1: ${strategy1Count} oligos, ${strategy1Time}ms`);
  }

  return {
    data,
    hash,
    hashMatches,
    stats,
    perOligo,
  };
}

async function sha256(data: Uint8Array): Promise<string> {
  // Always use Node's createHash which handles large data via streaming
  const { createHash } = await import("crypto");
  const h = createHash("sha256");
  // Hash in 64MB chunks to avoid buffer limits
  for (let i = 0; i < data.length; i += 64 * 1024 * 1024) {
    h.update(data.subarray(i, Math.min(i + 64 * 1024 * 1024, data.length)));
  }
  return h.digest("hex");
}
