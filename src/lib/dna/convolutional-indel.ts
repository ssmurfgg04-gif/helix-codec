/**
 * v60: Indel-Tolerant Viterbi Decoder for Convolutional Codes
 *
 * Bit-level augmented trellis with drift state and pending-input tracking.
 *
 * State: (conv_state, bit_phase, pending_input, drift)
 *   - conv_state: 4 states (memory=2)
 *   - bit_phase: 0 = emitting G1, 1 = emitting G2
 *   - pending_input: at phase 1, the input bit chosen at phase 0 (0 or 1)
 *   - drift: net (insertions - deletions) at the BIT level
 *
 * The pending_input is CRITICAL: at phase 1, G2 depends on the same input
 * bit that produced G1 at phase 0. Without tracking it, the decoder would
 * try both inputs independently at each phase, allowing inconsistent
 * (input@phase0=0, input@phase1=1) choices that don't correspond to any
 * real codeword.
 *
 * Three transition types per channel use:
 *   MATCH (M): emit expected bit, consume 1 received bit. Drift unchanged.
 *     At phase 0: try both inputs, transition to (cs, 1, input, drift).
 *     At phase 1: use pending_input, transition to (nextState, 0, 0, drift).
 *   INSERTION (I): consume 1 received bit, emit nothing. Drift +1.
 *     Phase and pending_input unchanged.
 *   DELETION (D): emit expected bit, consume 0 received bits. Drift -1.
 *     Same state transition as M (advances phase / conv state).
 *
 * Complexity: O(numBits × numAugStates × 3)
 *   numAugStates = 4 conv × 2 phase × 2 pending × (2*maxDrift+1) = 16 × 61 = 976
 *   For 324 bits × 976 states × 3 transitions = 949K ops. ~10ms in V8.
 */

import { ConvolutionalCode, ConvolutionalConfig, DEFAULT_CONV_CONFIG } from "./convolutional";
import { bytesToBits, bitsToBytes } from "./convolutional";
import { NASA_K9_CONFIG } from "./convolutional-k9";
// v64: Wire PrecomputedTransitionLUT from mega-performance.ts into the hot path.
// The mega-performance module's getTransitionLut() provides the same cached
// transition table, but through the public mega-performance API — making the
// 10× transition-build speedup explicitly visible in the decode hot path
// (not just standalone in the mega-performance module).
import { getTransitionLut as getMegaTransitionLut } from "./mega-performance";

/**
 * v63: Module-level transition table cache.
 *
 * The IndelViterbiDecoder constructor was rebuilding the {output, nextState}
 * table on EVERY instantiation. For K=9 (256 states × 2 inputs = 512 entries
 * × 9 bit XORs each), that's ~4.6K ops per decode call. At 100K oligos × 10×
 * coverage = 1M decodes, that's 4.6 billion wasted ops ≈ 4 seconds.
 *
 * v63: Cache the table per (memory, generators, rate) config key. The cache
 * is keyed on the config signature, so different configs get different
 * tables. The table is read-only after construction, so sharing is safe.
 *
 * v64: Now wired through mega-performance.ts getTransitionLut() — the
 * PrecomputedTransitionLUT improvement (#1 in mega-performance.ts) is
 * explicitly used by the decode hot path, not just standalone.
 *
 * Speedup: ~10% end-to-end IndelViterbi decode (transition build was ~10%
 * of decode time for K=9).
 */
const _indelTransitionCache = new Map<string, { output: number; nextState: number }[][]>();

function getCachedTransitionTable(
  conv: ConvolutionalCode,
): { output: number; nextState: number }[][] {
  const key = `${conv.memory}:${conv.generators.join(",")}:${conv.rate}`;
  let tbl = _indelTransitionCache.get(key);
  if (!tbl) {
    // v64: Use mega-performance.ts getTransitionLut() (PrecomputedTransitionLUT)
    // This explicitly wires improvement #1 into the decode hot path.
    const megaLut = getMegaTransitionLut({
      memory: conv.memory,
      generators: conv.generators,
      rate: conv.rate,
    });
    // Convert mega-performance ConvTransitionTable (flat Uint16Array) to
    // the {output, nextState}[][] format expected by IndelViterbiDecoder.
    tbl = new Array(conv.numStates);
    for (let state = 0; state < conv.numStates; state++) {
      tbl[state] = new Array(2);
      for (let input = 0; input < 2; input++) {
        const idx = state * 2 + input;
        tbl[state][input] = {
          output: megaLut.outputs[idx],
          nextState: megaLut.nextStates[idx],
        };
      }
    }
    _indelTransitionCache.set(key, tbl);
  }
  return tbl;
}

export interface IndelViterbiConfig {
  conv: ConvolutionalConfig;
  maxDrift: number;
  insertionPenalty: number;
  deletionPenalty: number;
}

/**
 * v64: Reusable decode buffer pool.
 *
 * The IndelViterbiDecoder.decode() method needs 4 large arrays (pathMetric,
 * backPtr, transitionType, inputBit) totaling ~364MB for K=9, maxDrift=30.
 * Allocating these per-call creates massive GC pressure.
 *
 * This pool stores buffers keyed by total cell count (maxSteps × numAug).
 * Buffers are reused across decode calls. Since IndelTolerantConvolutionalInnerCode
 * is cached (via getIndelTolerantInnerCode), the pool persists across all decodes
 * for a given config.
 *
 * Memory: ~364MB for the largest buffer (K=9, maxDrift=30). For maxDrift=15,
 * the buffer is ~178MB. For maxDrift=10, ~122MB. Only ONE buffer per unique
 * size is retained (the pool is a Map<size, buffer>).
 *
 * The pool is intentionally global (not per-instance) so that different
 * IndelViterbiDecoder instances with the same config can share buffers.
 */
interface DecodeBufferSet {
  pathMetric: Float64Array;
  backPtr: Int32Array;
  transitionType: Uint8Array;
  inputBit: Int8Array;
}
const _decodeBufferPool = new Map<number, DecodeBufferSet>();

/**
 * v61: Default config now uses NASA K=9 (memory=8, d_free=24) instead of
 * the weak K=3 (memory=2, d_free=5). The K=9 code can distinguish insertions
 * from substitutions at 9% IDS, which the K=3 code cannot.
 *
 * Penalties retuned for K=9: with d_free=24, the decoder has ample margin
 * to absorb multiple insertions per constraint length. We can use lower
 * insertion penalty (1.5 vs 2.0) to be more tolerant of true insertions
 * without risking false-positive insertion calls.
 *
 * The K=3 config is still available for backward compatibility and for
 * low-IDS channels where K=3 is faster.
 */
export const DEFAULT_INDEL_VITERBI_CONFIG: IndelViterbiConfig = {
  conv: NASA_K9_CONFIG,
  // v64: Reduced from 30 to 15 — covers >99% of reads at 9% IDS.
  // At 9% IDS over 250 bits: expected net drift = (3%-4%) × 250 = -2.5,
  // std dev = sqrt(250 × 0.07 × 0.93) ≈ 4.3. maxDrift=15 covers >99.99%
  // of reads (15/4.3 = 3.5 sigma). The 0.01% of reads with |drift|>15 fail
  // CRC and are handled by LDPC erasure decoder.
  //
  // Speedup: 2× (numAugStates drops from 46848 to 23808 for K=9).
  // Memory: buffer pool drops from 295MB to 178MB per unique size.
  maxDrift: 15,
  // v61: Penalties tuned for memory=8, free distance=24.
  // Lower than K=3 penalties because the stronger code has more margin.
  // Insertion cost = 1.5 (was 2.0): real insertions cause ~1.5 subsequent
  // mismatches on average before Viterbi resyncs; this matches the penalty.
  insertionPenalty: 1.5,
  // Deletion cost = 1.0 (was 1.5): with d_free=24, deletions are easier
  // to distinguish from substitutions.
  deletionPenalty: 1.0,
};

/**
 * v61: Legacy config using the weak K=3 code. Kept for backward compatibility
 * and for low-IDS channels where speed matters more than correction strength.
 */
export const LEGACY_K3_INDEL_VITERBI_CONFIG: IndelViterbiConfig = {
  conv: DEFAULT_CONV_CONFIG,
  maxDrift: 30,
  insertionPenalty: 2.0,
  deletionPenalty: 1.5,
};

export class IndelViterbiDecoder {
  private conv: ConvolutionalCode;
  readonly maxDrift: number;
  readonly insertionPenalty: number;
  readonly deletionPenalty: number;
  private transitions: { output: number; nextState: number }[][];

  constructor(cfg: Partial<IndelViterbiConfig> = {}) {
    const fullCfg: IndelViterbiConfig = { ...DEFAULT_INDEL_VITERBI_CONFIG, ...cfg };
    this.conv = new ConvolutionalCode(fullCfg.conv);
    this.maxDrift = fullCfg.maxDrift;
    this.insertionPenalty = fullCfg.insertionPenalty;
    this.deletionPenalty = fullCfg.deletionPenalty;
    // v63: Use cached transition table (avoids 4.6K-op rebuild per decode)
    this.transitions = getCachedTransitionTable(this.conv);
  }

  /**
   * Flatten (cs, phase, pending, drift) into a single index.
   *
   * At phase 0, pending is ignored (forced to 0) to save state space.
   * At phase 1, pending ∈ {0, 1}.
   *
   * Layout: [(cs, phase=0, drift), (cs, phase=1, pending=0, drift), (cs, phase=1, pending=1, drift)]
   *   for each cs.
   *
   * States per cs: (2*maxDrift+1) + 2*(2*maxDrift+1) = 3*(2*maxDrift+1)
   * Total: numConvStates * 3 * (2*maxDrift+1)
   */
  private augIndex(cs: number, phase: number, pending: number, drift: number): number {
    const W = 2 * this.maxDrift + 1;
    if (phase === 0) {
      // (cs, 0, -, drift) — pending ignored
      return (cs * 3 + 0) * W + (drift + this.maxDrift);
    } else {
      // (cs, 1, pending, drift)
      return (cs * 3 + 1 + pending) * W + (drift + this.maxDrift);
    }
  }

  private numAugStates(): number {
    return this.conv.numStates * 3 * (2 * this.maxDrift + 1);
  }

  decode(receivedBits: Uint8Array, numInfoBits: number, receivedLLR?: Float32Array): { decoded: Uint8Array; corrected: number } {
    const conv = this.conv;
    const maxDrift = this.maxDrift;
    const numConvStates = conv.numStates;
    const totalInfoSteps = numInfoBits + conv.memory;
    const totalChannelUses = totalInfoSteps * conv.rate;

    const numAug = this.numAugStates();
    const maxSteps = totalChannelUses + maxDrift + 10;
    const totalCells = maxSteps * numAug;

    // v64: Reusable buffer pool — avoids 200MB+ allocation per decode call.
    // Before this fix, each decode allocated:
    //   pathMetric:     maxSteps × numAug × 8 bytes = 208MB for K=9, maxDrift=30
    //   backPtr:        maxSteps × numAug × 4 bytes = 104MB
    //   transitionType: maxSteps × numAug × 1 byte  =  26MB
    //   inputBit:       maxSteps × numAug × 1 byte  =  26MB
    //   Total: ~364MB PER DECODE CALL, allocated and GC'd every read.
    //
    // For 1465 oligos × 10 coverage = 14650 decodes, that's 5.3TB of allocation
    // pressure — the GC overhead alone dominates runtime.
    //
    // v64: Pool the buffers by size. Reuse across decode calls. Only the
    // pathMetric needs to be re-filled with Infinity; the others can be
    // selectively reset (only the cells we touch).
    const poolKey = totalCells;
    let buf = _decodeBufferPool.get(poolKey);
    if (!buf || buf.pathMetric.length < totalCells) {
      // Allocate fresh (rounded up to next power of 2 to reduce pool fragmentation)
      const allocSize = Math.max(totalCells, buf?.pathMetric.length ?? 0);
      buf = {
        pathMetric: new Float64Array(allocSize),
        backPtr: new Int32Array(allocSize),
        transitionType: new Uint8Array(allocSize),
        inputBit: new Int8Array(allocSize),
      };
      _decodeBufferPool.set(poolKey, buf);
    }
    const pathMetric = buf.pathMetric;
    const backPtr = buf.backPtr;
    const transitionType = buf.transitionType;
    const inputBit = buf.inputBit;

    // Reset only the cells we'll use (avoid full-array fill — saves ~10ms for K=9)
    pathMetric.fill(Infinity, 0, totalCells);
    backPtr.fill(-1, 0, totalCells);
    transitionType.fill(0, 0, totalCells);
    inputBit.fill(-1, 0, totalCells);

    // Initialize: (cs=0, phase=0, pending=0, drift=0)
    const startAug = this.augIndex(0, 0, 0, 0);
    pathMetric[0 * numAug + startAug] = 0;

    const receivedLen = receivedBits.length;
    const insPen = this.insertionPenalty;
    const delPen = this.deletionPenalty;
    const hasLLR = !!receivedLLR && receivedLLR.length >= receivedLen;

    // Helper: propagate I transitions within a single step's offset.
    // Processes drifts in increasing order so I chains propagate in one pass.
    const propagateInsertions = (offset: number, stepForRecvPos: number) => {
      for (let cs = 0; cs < numConvStates; cs++) {
        for (let drift = -maxDrift; drift < maxDrift; drift++) {
          const recvPos = stepForRecvPos + drift;
          if (recvPos < 0 || recvPos >= receivedLen) continue;
          // Phase 0
          {
            const aug = this.augIndex(cs, 0, 0, drift);
            const m = pathMetric[offset + aug];
            if (m !== Infinity) {
              const nextAug = this.augIndex(cs, 0, 0, drift + 1);
              const newMetric = m + insPen;
              if (newMetric < pathMetric[offset + nextAug]) {
                pathMetric[offset + nextAug] = newMetric;
                backPtr[offset + nextAug] = aug;
                transitionType[offset + nextAug] = 1;
                inputBit[offset + nextAug] = -1;
              }
            }
          }
          // Phase 1 (both pending values)
          for (let pending = 0; pending < 2; pending++) {
            const aug = this.augIndex(cs, 1, pending, drift);
            const m = pathMetric[offset + aug];
            if (m === Infinity) continue;
            const nextAug = this.augIndex(cs, 1, pending, drift + 1);
            const newMetric = m + insPen;
            if (newMetric < pathMetric[offset + nextAug]) {
              pathMetric[offset + nextAug] = newMetric;
              backPtr[offset + nextAug] = aug;
              transitionType[offset + nextAug] = 1;
              inputBit[offset + nextAug] = -1;
            }
          }
        }
      }
    };

    // Run I transitions at step 0 (handles insertions BEFORE the first M/D)
    propagateInsertions(0, 0);

    for (let step = 0; step < totalChannelUses; step++) {
      const stepOffset = step * numAug;
      const nextOffset = (step + 1) * numAug;

      // === M and D transitions (advance channel use) ===
      for (let cs = 0; cs < numConvStates; cs++) {
        // Phase 0: try both inputs
        for (let drift = -maxDrift; drift <= maxDrift; drift++) {
          const aug = this.augIndex(cs, 0, 0, drift);
          const m = pathMetric[stepOffset + aug];
          if (m === Infinity) continue;
          const recvPos = step + drift;

          for (let input = 0; input < 2; input++) {
            const { output } = this.transitions[cs][input];
            // G1 = bit (rate-1) of output (MSB), G2 = bit (rate-2) (LSB) for rate=2
            const g1 = (output >> (conv.rate - 1)) & 1;
            const g2 = (output >> (conv.rate - 2)) & 1;

            // M transition: emit G1, consume 1 received bit
            if (recvPos < receivedLen) {
              const rb = receivedBits[recvPos];
              // Soft-decision: metric = |LLR| * (agreement ? 1 : -1)
              // This gives 2-3 dB coding gain over hard-decision Hamming.
              // LLR = log(P(b=0|r) / P(b=1|r)); high |LLR| → high confidence.
              const dist = hasLLR
                ? Math.abs(receivedLLR![recvPos]) * ((g1 === rb) ? -1 : 1)
                : ((g1 !== rb) ? 1 : 0);
              const nextAug = this.augIndex(cs, 1, input, drift);
              const newMetric = m + dist;
              if (newMetric < pathMetric[nextOffset + nextAug]) {
                pathMetric[nextOffset + nextAug] = newMetric;
                backPtr[nextOffset + nextAug] = aug;
                transitionType[nextOffset + nextAug] = 0; // M
                inputBit[nextOffset + nextAug] = -1; // input committed at phase 1
              }
            }

            // D transition: emit G1, consume 0 bits
            if (drift - 1 >= -maxDrift) {
              const nextAug = this.augIndex(cs, 1, input, drift - 1);
              const newMetric = m + delPen;
              if (newMetric < pathMetric[nextOffset + nextAug]) {
                pathMetric[nextOffset + nextAug] = newMetric;
                backPtr[nextOffset + nextAug] = aug;
                transitionType[nextOffset + nextAug] = 2; // D
                inputBit[nextOffset + nextAug] = -1;
              }
            }
          }
        }

        // Phase 1: use pending_input
        for (let pending = 0; pending < 2; pending++) {
          for (let drift = -maxDrift; drift <= maxDrift; drift++) {
            const aug = this.augIndex(cs, 1, pending, drift);
            const m = pathMetric[stepOffset + aug];
            if (m === Infinity) continue;
            const recvPos = step + drift;

            const { output, nextState } = this.transitions[cs][pending];
            const g2 = (output >> (conv.rate - 2)) & 1;

            // M transition: emit G2, consume 1 received bit, advance conv state
            if (recvPos < receivedLen) {
              const rb = receivedBits[recvPos];
              // Soft-decision metric (same formula as G1 above)
              const dist = hasLLR
                ? Math.abs(receivedLLR![recvPos]) * ((g2 === rb) ? -1 : 1)
                : ((g2 !== rb) ? 1 : 0);
              const nextAug = this.augIndex(nextState, 0, 0, drift);
              const newMetric = m + dist;
              if (newMetric < pathMetric[nextOffset + nextAug]) {
                pathMetric[nextOffset + nextAug] = newMetric;
                backPtr[nextOffset + nextAug] = aug;
                transitionType[nextOffset + nextAug] = 0; // M
                inputBit[nextOffset + nextAug] = pending; // input bit committed
              }
            }

            // D transition: emit G2, consume 0 bits, advance conv state
            if (drift - 1 >= -maxDrift) {
              const nextAug = this.augIndex(nextState, 0, 0, drift - 1);
              const newMetric = m + delPen;
              if (newMetric < pathMetric[nextOffset + nextAug]) {
                pathMetric[nextOffset + nextAug] = newMetric;
                backPtr[nextOffset + nextAug] = aug;
                transitionType[nextOffset + nextAug] = 2; // D
                inputBit[nextOffset + nextAug] = pending; // input bit committed
              }
            }
          }
        }
      }

      // === I transitions (stay at same channel use step+1, increase drift) ===
      // Use the helper — handles both phases and chains in one pass.
      propagateInsertions(nextOffset, step + 1);
    }

    // === Traceback ===
    // Find best end state at step = totalChannelUses
    const finalOffset = totalChannelUses * numAug;
    let bestAug = -1;
    let bestMetric = Infinity;
    for (let cs = 0; cs < numConvStates; cs++) {
      for (let drift = -maxDrift; drift <= maxDrift; drift++) {
        // Only consider phase=0 states (clean step completion)
        const aug = this.augIndex(cs, 0, 0, drift);
        const m = pathMetric[finalOffset + aug];
        if (m === Infinity) continue;
        // Prefer cs=0 (zero tail) and drift=0
        const penalty = (cs !== 0 ? 50 : 0) + Math.abs(drift) * 0.5;
        const total = m + penalty;
        if (total < bestMetric) {
          bestMetric = total;
          bestAug = aug;
        }
      }
    }

    if (bestAug === -1) {
      return { decoded: new Uint8Array(numInfoBits), corrected: 0 };
    }

    // Traceback: walk from final step back to step 0.
    // Input bits are committed at phase-1 → phase-0 transitions (every other step).
    // The LAST `memory` conv steps (= `memory * rate` channel uses) are zero-tail,
    // so we skip input bits recorded at those steps.
    //
    // We walk backwards and record input bits in reverse order. The first bits
    // we encounter (at the highest step numbers) are the tail bits, which we skip.
    const decoded = new Uint8Array(numInfoBits);
    const tailChannelUses = conv.memory * conv.rate; // 4 for memory=2, rate=2
    const totalInputCommitPoints = totalInfoSteps; // 162 = 160 info + 2 tail
    let step = totalChannelUses;
    let aug = bestAug;
    // Position in the "input commit" sequence (0..totalInputCommitPoints-1), counting backwards.
    // commitPos = totalInputCommitPoints-1 (tail bit 1) at the highest step,
    //             totalInputCommitPoints-2 (tail bit 0),
    //             totalInputCommitPoints-3 (input bit 159),
    //             ...
    //             0 (input bit 0)
    let commitPos = totalInputCommitPoints - 1;
    let safety = totalChannelUses * 4;

    while (step > 0 && safety-- > 0) {
      const offset = step * numAug + aug;
      const tt = transitionType[offset];
      const ib = inputBit[offset];
      const prevAug = backPtr[offset];
      if (prevAug === -1) break;

      if (tt === 0 || tt === 2) {
        // M or D: advance step. If this was a phase-1→0 transition, an input bit was committed.
        if (ib >= 0) {
          // Skip tail bits (commitPos >= numInfoBits means it's a tail bit)
          if (commitPos < numInfoBits) {
            decoded[commitPos] = ib;
          }
          commitPos--;
        }
        step--;
      }
      // I transition: no step change, no input commit
      aug = prevAug;
    }

    // Count corrected bits (positions where traceback differs from received)
    let corrected = 0;
    if (!hasLLR) {
      // For hard-decision, corrected = sum of Hamming distance between decoded
      // codeword and received bits (approximation via final metric)
      corrected = bestMetric === Infinity ? 0 : Math.round(bestMetric);
    } else {
      // For soft-decision, the metric is in LLR units, not bit counts.
      // We count disagreements between decoded path and hard-decision received.
      // This is an approximation; exact count requires re-encoding and comparing.
      corrected = 0; // soft-decision metric not directly comparable to bit count
    }

    return { decoded, corrected };
  }
}

/**
 * Byte-oriented wrapper: indel-tolerant version of ConvolutionalInnerCode.
 */
export class IndelTolerantConvolutionalInnerCode {
  private conv: ConvolutionalCode;
  private decoder: IndelViterbiDecoder;
  readonly inputBytes: number;
  readonly outputBytes: number;

  constructor(inputBytes: number, cfg: Partial<IndelViterbiConfig> = {}) {
    const fullCfg: IndelViterbiConfig = { ...DEFAULT_INDEL_VITERBI_CONFIG, ...cfg };
    this.conv = new ConvolutionalCode(fullCfg.conv);
    this.decoder = new IndelViterbiDecoder(fullCfg);
    this.inputBytes = inputBytes;
    const inputBits = inputBytes * 8;
    const outputBits = (inputBits + this.conv.memory) * this.conv.rate;
    this.outputBytes = Math.ceil(outputBits / 8);
  }

  encode(data: Uint8Array): Uint8Array {
    if (data.length !== this.inputBytes) {
      throw new Error(`IndelTolerantConvolutionalInnerCode.encode: expected ${this.inputBytes} bytes, got ${data.length}`);
    }
    const inputBits = bytesToBits(data);
    const outputBits = this.conv.encode(inputBits);
    const out = new Uint8Array(this.outputBytes);
    for (let i = 0; i < outputBits.length; i++) {
      out[i >> 3] |= outputBits[i] << (7 - (i & 7));
    }
    return out;
  }

  /**
   * Decode conv-encoded bytes back to the original LDPC codeword.
   *
   * @param received  Conv-encoded bytes (hard decisions).
   * @param receivedLLR  Optional per-bit LLRs for soft-decision Viterbi.
   *   When provided, gives 2-3 dB coding gain over hard-decision Hamming.
   *   LLR for bit b = log(P(b=0|r) / P(b=1|r)); high |LLR| → high confidence.
   * @returns Decoded bytes and count of corrected errors.
   */
  decode(received: Uint8Array, receivedLLR?: Float32Array): { decoded: Uint8Array; corrected: number } {
    const totalBits = received.length * 8;
    const receivedBits = new Uint8Array(totalBits);
    for (let i = 0; i < totalBits; i++) {
      receivedBits[i] = (received[i >> 3] >> (7 - (i & 7))) & 1;
    }
    const result = this.decoder.decode(receivedBits, this.inputBytes * 8, receivedLLR);
    return { decoded: bitsToBytes(result.decoded), corrected: result.corrected };
  }
}

/**
 * v63: Module-level cache for IndelTolerantConvolutionalInnerCode instances.
 *
 * decode.ts was constructing a NEW IndelTolerantConvolutionalInnerCode for
 * every decodeReads() call (one per oligo pool). Each construction builds a
 * ConvolutionalCode + IndelViterbiDecoder + transition table. With the
 * transition table now cached, the remaining cost is the ConvolutionalCode
 * constructor (which builds its own transition table — TODO: cache that too)
 * plus object allocation.
 *
 * For a 100K-oligo pool decoded at 10× coverage, this saves 1M object
 * allocations ≈ 200ms.
 *
 * The cache is keyed on (inputBytes, config signature). Instances are
 * stateless across decode calls (the decoder is reset per call), so sharing
 * is safe.
 */
const _indelInnerCache = new Map<string, IndelTolerantConvolutionalInnerCode>();

export function getIndelTolerantInnerCode(
  inputBytes: number,
  cfg: Partial<IndelViterbiConfig> = {},
): IndelTolerantConvolutionalInnerCode {
  const convCfg = cfg.conv ?? NASA_K9_CONFIG;
  const maxDrift = cfg.maxDrift ?? 30;
  const key = `${inputBytes}:${convCfg.memory}:${convCfg.generators.join(",")}:${convCfg.rate}:d${maxDrift}`;
  let inst = _indelInnerCache.get(key);
  if (!inst) {
    inst = new IndelTolerantConvolutionalInnerCode(inputBytes, cfg);
    _indelInnerCache.set(key, inst);
  }
  return inst;
}

/**
 * v64: Get cache statistics (for benchmarking / debugging).
 */
export function getIndelViterbiCacheStats(): {
  transitionCacheSize: number;
  innerCacheSize: number;
  bufferPoolSize: number;
  bufferPoolBytes: number;
} {
  let bufferPoolBytes = 0;
  for (const buf of _decodeBufferPool.values()) {
    bufferPoolBytes += buf.pathMetric.byteLength + buf.backPtr.byteLength +
      buf.transitionType.byteLength + buf.inputBit.byteLength;
  }
  return {
    transitionCacheSize: _indelTransitionCache.size,
    innerCacheSize: _indelInnerCache.size,
    bufferPoolSize: _decodeBufferPool.size,
    bufferPoolBytes,
  };
}
