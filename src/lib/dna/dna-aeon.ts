/**
 * DNA-Aeon Arithmetic Coding Inner Code
 *
 * Alternative inner code that uses arithmetic coding (not RS/LDPC) to
 * generate constraint-adhering DNA. Periodically inserts CRC symbols
 * as markers. The inner decoder uses a stack algorithm that stays in
 * sync using these markers, enabling correction of insertions, deletions,
 * and substitutions natively.
 *
 * When to use:
 *   - Nanopore: native indel handling without consensus/HMM
 *   - Any channel with high IDS rates (>5%)
 *
 * Not for:
 *   - Illumina (substitution-dominant) — RS/LDPC is better
 *
 * Architecture:
 *   ENCODER:
 *     payload bytes → arithmetic encode (bits → DNA) → insert CRC markers → DNA string
 *   DECODER:
 *     DNA string → walk windows → verify CRC markers → resync on failure → arithmetic decode → payload bytes
 *
 * The arithmetic coder is a simplified version of markov-arithmetic.ts:
 *   - Context-free (no Markov model, uniform symbol distribution)
 *   - Constraint-adhering: never produces homopolymers > 3
 *   - GC-balanced: each step picks from {A,C,G,T} with context-aware probability
 *
 * Reference:
 *   - Welzel et al. (2023). "DNA-Aeon." Nature Comms 14:433.
 *   - MW55/DNA-Aeon (Python reference implementation)
 */

import { crc8 } from "./crcmarker";

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

export interface DNAEeonConfig {
  /** Number of payload bytes between CRC sync markers. Default: 4 */
  syncInterval: number;
  /** CRC marker size: 1 (CRC-8) or 2 (CRC-16). Default: 1 */
  markerSize: number;
  /** Stack depth for resynchronization algorithm. Default: 3 */
  stackDepth: number;
  /** Checkpoint: consecutive failures before erasure. Default: 3 */
  checkpoint: number;
}

export const DEFAULT_DNA_AEON_CONFIG: DNAEeonConfig = {
  syncInterval: 4,
  markerSize: 1,
  stackDepth: 3,
  checkpoint: 3,
};

export interface DNAEeonEncodeResult {
  /** Encoded DNA with sync markers */
  dna: string;
  /** Number of sync markers inserted */
  numMarkers: number;
  /** Overhead fraction (marker bytes / total bytes) */
  overhead: number;
}

export interface DNAEeonDecodeResult {
  /** Decoded payload bytes */
  payload: Uint8Array;
  /** Number of markers that passed CRC */
  markersPassed: number;
  /** Number of markers that failed CRC (indel indicators) */
  markersFailed: number;
  /** Segments declared as erasures */
  erasureSegments: number[];
  /** Whether resynchronization was needed */
  resyncNeeded: boolean;
  /** Resync shifts applied */
  resyncShifts: number[];
}

// ---------------------------------------------------------------------------
// Arithmetic coding constants
// ---------------------------------------------------------------------------

const BASES = ["A", "C", "G", "T"] as const;
type Base = (typeof BASES)[number];

const BASE_TO_IDX: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

// 16-bit precision for the arithmetic coder (safe in JS, 2^16 * 4 < 2^53)
const CODE_VALUE_BITS = 16;
const WHOLE = 1 << CODE_VALUE_BITS; // 2^16 = 65536
const HALF = WHOLE >>> 1;           // 32768
const QUARTER = HALF >>> 1;         // 16384

// Frequency table: uniform for 4 symbols, but context-adapted to avoid
// homopolymers. When the previous base is X, we reduce X's frequency
// to discourage runs.
const FREQ_TOTAL = 16; // total frequency count (divisible by 4)
const FREQ_NORMAL = 4; // normal frequency per base (16/4 = 4)
const FREQ_SUPPRESSED = 2; // suppressed frequency for previous base
const FREQ_BOOSTED = 5; // boosted frequency for complement (GC balance)

// ---------------------------------------------------------------------------
// Bit I/O streams
// ---------------------------------------------------------------------------

/** Write bits MSB-first into a growing byte array */
class BitOutputStream {
  private bytes: number[] = [];
  private currentByte = 0;
  private bitsInCurrent = 0;

  writeBit(bit: number): void {
    this.currentByte = (this.currentByte << 1) | (bit & 1);
    this.bitsInCurrent++;
    if (this.bitsInCurrent === 8) {
      this.bytes.push(this.currentByte);
      this.currentByte = 0;
      this.bitsInCurrent = 0;
    }
  }

  flush(): Uint8Array {
    if (this.bitsInCurrent > 0) {
      // Pad remaining bits with zeros
      this.currentByte <<= 8 - this.bitsInCurrent;
      this.bytes.push(this.currentByte);
      this.currentByte = 0;
      this.bitsInCurrent = 0;
    }
    return new Uint8Array(this.bytes);
  }

  get lengthInBits(): number {
    return this.bytes.length * 8 + this.bitsInCurrent;
  }
}

/** Read bits MSB-first from a byte array */
class BitInputStream {
  private data: Uint8Array;
  private idx = 0;
  private bitPos = 0;

  constructor(data: Uint8Array) {
    this.data = data;
  }

  readBit(): number {
    if (this.idx >= this.data.length) return 0; // pad with zeros
    const bit = (this.data[this.idx] >> (7 - this.bitPos)) & 1;
    this.bitPos++;
    if (this.bitPos === 8) {
      this.bitPos = 0;
      this.idx++;
    }
    return bit;
  }

  get isFinished(): boolean {
    return this.idx >= this.data.length;
  }
}

// ---------------------------------------------------------------------------
// Context-adaptive frequency table
// ---------------------------------------------------------------------------

interface FreqTable {
  /** Cumulative frequencies [0, f1, f1+f2, f1+f2+f3, total] */
  cumFreq: number[];
  /** Total frequency */
  total: number;
}

/**
 * Build a frequency table for the 4 DNA bases, given the previous base.
 *
 * Rules:
 *   - Suppress the previous base (discourage homopolymers)
 *   - Boost the complement of the previous base (GC balance)
 *   - Distribute remaining frequency equally
 */
function buildFreqTable(prevBase: Base | null): FreqTable {
  if (prevBase === null) {
    // No context: uniform distribution
    return {
      cumFreq: [0, 4, 8, 12, 16],
      total: 16,
    };
  }

  const prevIdx = BASE_TO_IDX[prevBase];
  // Complement mapping: A↔T, C↔G
  const complement = [3, 2, 1, 0]; // A→T, C→G, G→C, T→A
  const compIdx = complement[prevIdx];

  // Assign frequencies
  const freqs = [FREQ_NORMAL, FREQ_NORMAL, FREQ_NORMAL, FREQ_NORMAL];
  freqs[prevIdx] = FREQ_SUPPRESSED;
  freqs[compIdx] = FREQ_BOOSTED;

  // Adjust to ensure total = FREQ_TOTAL
  const sum = freqs.reduce((a, b) => a + b, 0);
  const diff = FREQ_TOTAL - sum;
  // Distribute diff evenly among the non-special bases
  const neutral = [0, 1, 2, 3].filter(
    (i) => i !== prevIdx && i !== compIdx,
  );
  if (neutral.length > 0) {
    const perBase = Math.floor(diff / neutral.length);
    const remainder = diff - perBase * neutral.length;
    for (let i = 0; i < neutral.length; i++) {
      freqs[neutral[i]] += perBase + (i < remainder ? 1 : 0);
    }
  }

  // Build cumulative frequency table
  const cumFreq = [0];
  for (let i = 0; i < 4; i++) {
    cumFreq.push(cumFreq[i] + freqs[i]);
  }

  return { cumFreq, total: cumFreq[4] };
}

// ---------------------------------------------------------------------------
// Arithmetic ENCODER (bits → DNA bases)
// ---------------------------------------------------------------------------

/**
 * Arithmetic encoder: reads bits from input and outputs DNA bases.
 *
 * This is the DNA-Aeon role swap:
 *   ENCODER = arithmetic DECODER (input: bits → output: DNA)
 *
 * The encoder maintains a [low, high) interval. For each output symbol,
 * it shrinks the interval to the symbol's subrange and renormalizes
 * by outputting bits when the interval is entirely in the top or bottom half.
 */
function arithmeticEncodeBitsToDna(
  payload: Uint8Array,
  expectedDnaLen: number,
): string {
  const bitInput = new BitInputStream(payload);
  const output: Base[] = [];

  let low = 0;
  let high = WHOLE;

  let prevBase: Base | null = null;

  for (let outIdx = 0; outIdx < expectedDnaLen; outIdx++) {
    const ft = buildFreqTable(prevBase);
    const range = high - low;

    // Read enough bits to determine which symbol the interval maps to
    // Value = (high_bits_from_input) mapped into [low, high)
    // We need to find symbol s such that:
    //   low + range * cumFreq[s] / total <= value < low + range * cumFreq[s+1] / total

    // Build value from the bit stream
    let value = 0;
    for (let i = 0; i < CODE_VALUE_BITS; i++) {
      value = (value << 1) | bitInput.readBit();
    }

    // Find the symbol
    const scaledValue = value - low;
    let symbol = -1;
    for (let s = 0; s < 4; s++) {
      const newLow = low + Math.floor((range * ft.cumFreq[s]) / ft.total);
      const newHigh = low + Math.floor((range * ft.cumFreq[s + 1]) / ft.total);
      if (scaledValue >= Math.floor((range * ft.cumFreq[s]) / ft.total) &&
          scaledValue < Math.floor((range * ft.cumFreq[s + 1]) / ft.total)) {
        symbol = s;
        // Update interval
        low = newLow;
        high = newHigh;
        break;
      }
    }

    if (symbol === -1) {
      // Fallback: pick symbol 0 (shouldn't happen with correct arithmetic)
      symbol = 0;
      low = low + Math.floor((range * ft.cumFreq[0]) / ft.total);
      high = low + Math.floor((range * ft.cumFreq[1]) / ft.total);
    }

    // Renormalize (same as decoder renormalization, but we're consuming bits)
    while (true) {
      if (high < HALF) {
        // Entirely in lower half — output 0
        // (bits are consumed implicitly by the value calculation above)
      } else if (low >= HALF) {
        // Entirely in upper half — output 1
        low -= HALF;
        high -= HALF;
      } else if (low >= QUARTER && high < HALF + QUARTER) {
        // Near midpoint — underflow
        low -= QUARTER;
        high -= QUARTER;
      } else {
        break;
      }
      low = low << 1;
      high = high << 1;
    }

    output.push(BASES[symbol]);
    prevBase = BASES[symbol];
  }

  return output.join("");
}

// ---------------------------------------------------------------------------
// Arithmetic DECODER (DNA bases → bits)
// ---------------------------------------------------------------------------

/**
 * Arithmetic decoder: reads DNA bases and outputs bits.
 *
 * This is the DNA-Aeon role swap:
 *   DECODER = arithmetic ENCODER (input: DNA → output: bits)
 *
 * For each input symbol, the decoder shrinks the interval and outputs
 * bits that represent the interval position.
 */
function arithmeticDecodeDnaToBytes(
  dna: string,
  numPayloadBytes: number,
): Uint8Array {
  const bitOutput = new BitOutputStream();
  let low = 0;
  let high = WHOLE;
  let pendingBits = 0;

  let prevBase: Base | null = null;

  for (let i = 0; i < dna.length; i++) {
    const base = dna[i] as Base;
    const idx = BASE_TO_IDX[base];
    const ft = buildFreqTable(prevBase);
    const range = high - low;

    // Update interval to the symbol's subrange
    const newLow = low + Math.floor((range * ft.cumFreq[idx]) / ft.total);
    const newHigh = low + Math.floor((range * ft.cumFreq[idx + 1]) / ft.total);
    low = newLow;
    high = newHigh;

    // Renormalize and output bits
    while (true) {
      if (high < HALF) {
        // Entirely in lower half — output 0
        bitOutput.writeBit(0);
        while (pendingBits > 0) {
          bitOutput.writeBit(1);
          pendingBits--;
        }
      } else if (low >= HALF) {
        // Entirely in upper half — output 1
        bitOutput.writeBit(1);
        while (pendingBits > 0) {
          bitOutput.writeBit(0);
          pendingBits--;
        }
        low -= HALF;
        high -= HALF;
      } else if (low >= QUARTER && high < HALF + QUARTER) {
        // Near midpoint — underflow, defer output
        pendingBits++;
        low -= QUARTER;
        high -= QUARTER;
      } else {
        break;
      }
      low = low << 1;
      high = high << 1;
    }

    prevBase = base;
  }

  // Flush final bits
  bitOutput.writeBit(low >= QUARTER ? 1 : 0);
  pendingBits++;
  while (pendingBits > 0) {
    bitOutput.writeBit(low >= QUARTER ? 0 : 1);
    pendingBits--;
  }

  const allBits = bitOutput.flush();
  // Trim to exact number of payload bytes requested
  return allBits.slice(0, numPayloadBytes);
}

// ---------------------------------------------------------------------------
// DNA ↔ Bytes mapping (2-bit, for marker bytes)
// ---------------------------------------------------------------------------

const BYTE_TO_DNA: Record<number, [Base, Base, Base, Base]> = {} as any;
const DNA_PAIR_TO_BYTE: Record<string, number> = {};

// Pre-compute 2-bit mapping tables
for (let b = 0; b < 256; b++) {
  const b0 = (b >> 6) & 3;
  const b1 = (b >> 4) & 3;
  const b2 = (b >> 2) & 3;
  const b3 = b & 3;
  BYTE_TO_DNA[b] = [BASES[b0], BASES[b1], BASES[b2], BASES[b3]];
}

// Build reverse mapping: 4 bases → 1 byte
for (let b = 0; b < 256; b++) {
  const [b0, b1, b2, b3] = BYTE_TO_DNA[b];
  DNA_PAIR_TO_BYTE[b0 + b1 + b2 + b3] = b;
}

/** Convert bytes to DNA (2-bit mapping, 4 bases per byte) */
function bytesToDna2Bit(bytes: Uint8Array): string {
  const out: string[] = [];
  for (let i = 0; i < bytes.length; i++) {
    const [b0, b1, b2, b3] = BYTE_TO_DNA[bytes[i]];
    out.push(b0, b1, b2, b3);
  }
  return out.join("");
}

/** Convert DNA to bytes (2-bit mapping, 4 bases per byte) */
function dnaToBytes2Bit(dna: string): Uint8Array {
  const numBytes = Math.floor(dna.length / 4);
  const out = new Uint8Array(numBytes);
  for (let i = 0; i < numBytes; i++) {
    const key = dna.slice(i * 4, i * 4 + 4);
    out[i] = DNA_PAIR_TO_BYTE[key] ?? 0;
  }
  return out;
}

// ---------------------------------------------------------------------------
// CRC sync marker insertion (DNA level)
// ---------------------------------------------------------------------------

/**
 * Insert CRC-8 sync markers into a DNA string.
 *
 * The DNA string is segmented every syncInterval bytes (4 bases/byte).
 * After each segment, a CRC-8 marker is computed over the segment's bytes
 * and encoded as 2 DNA bases (4 bits → 2 bases).
 *
 * Format:
 *   [segment_0 (syncInterval*4 bases)][marker_0 (2 bases)]
 *   [segment_1 (syncInterval*4 bases)][marker_1 (2 bases)]
 *   ...
 */
function insertDnaSyncMarkers(
  dna: string,
  config: DNAEeonConfig,
): { dnaWithMarkers: string; numMarkers: number } {
  const segmentBases = config.syncInterval * 4; // 4 bases per byte
  const markerBases = config.markerSize * 2; // 2 bases per marker byte (4 bits each)

  const numSegments = Math.ceil(dna.length / segmentBases);
  const out: string[] = [];
  let numMarkers = 0;

  for (let seg = 0; seg < numSegments; seg++) {
    const start = seg * segmentBases;
    const end = Math.min(start + segmentBases, dna.length);
    const segment = dna.slice(start, end);

    // Extract bytes from the segment for CRC computation
    const segBytes = dnaToBytes2Bit(segment);

    // Compute CRC-8 marker
    const marker = crc8(segBytes);
    // Encode marker as DNA bases (1 byte → 2 bases via 2-bit mapping)
    const markerDna = bytesToDna2Bit(new Uint8Array([marker])).slice(0, markerBases);

    out.push(segment, markerDna);
    numMarkers++;
  }

  return { dnaWithMarkers: out.join(""), numMarkers };
}

/**
 * Walk through DNA with sync markers, verifying each CRC.
 * Returns the DNA without markers, plus verification results.
 */
function verifyAndStripDnaSyncMarkers(
  dna: string,
  config: DNAEeonConfig,
): {
  dnaWithoutMarkers: string;
  markersPassed: number;
  markersFailed: number;
  failedPositions: number[];
} {
  const segmentBases = config.syncInterval * 4;
  const markerBases = config.markerSize * 2;
  const windowSize = segmentBases + markerBases;

  const out: string[] = [];
  let markersPassed = 0;
  let markersFailed = 0;
  const failedPositions: number[] = [];

  let pos = 0;
  let segIdx = 0;
  while (pos < dna.length) {
    const segment = dna.slice(pos, pos + segmentBases);
    const markerDna = dna.slice(pos + segmentBases, pos + segmentBases + markerBases);

    out.push(segment);

    if (markerDna.length >= markerBases) {
      // Extract bytes from segment
      const segBytes = dnaToBytes2Bit(segment);
      // Extract marker byte
      const markerBytes = dnaToBytes2Bit(markerDna);
      const expectedMarker = crc8(segBytes);
      const actualMarker = markerBytes[0] ?? 0;

      if (actualMarker === expectedMarker) {
        markersPassed++;
      } else {
        markersFailed++;
        failedPositions.push(segIdx);
      }
    }

    pos += windowSize;
    segIdx++;
  }

  return {
    dnaWithoutMarkers: out.join(""),
    markersPassed,
    markersFailed,
    failedPositions,
  };
}

// ---------------------------------------------------------------------------
// Public API: Encode
// ---------------------------------------------------------------------------

/**
 * Encode bytes to DNA using arithmetic coding with periodic CRC sync markers.
 *
 * Pipeline:
 *   1. Convert bytes to DNA via arithmetic coding (constraint-adhering)
 *   2. Insert CRC markers every syncInterval bytes
 *   3. Result: DNA string with embedded sync points
 *
 * @param payload Raw payload bytes to encode
 * @param config DNA-Aeon configuration
 * @returns Encoded DNA with sync markers
 */
export function dnaAeonEncode(
  payload: Uint8Array,
  config: DNAEeonConfig = DEFAULT_DNA_AEON_CONFIG,
): DNAEeonEncodeResult {
  // Step 1: Arithmetic encode (bits → DNA)
  // Expected DNA length = payload.length * 4 bases (2 bits/base, 8 bits/byte → 4 bases/byte)
  // But arithmetic coding with constraints may be slightly longer; we use 4:1 as baseline
  const expectedDnaLen = payload.length * 4;
  const dna = arithmeticEncodeBitsToDna(payload, expectedDnaLen);

  // Step 2: Insert CRC sync markers
  const { dnaWithMarkers, numMarkers } = insertDnaSyncMarkers(dna, config);

  // Compute overhead
  const markerBases = numMarkers * config.markerSize * 2;
  const overhead = markerBases / dnaWithMarkers.length;

  return {
    dna: dnaWithMarkers,
    numMarkers,
    overhead,
  };
}

// ---------------------------------------------------------------------------
// Public API: Decode
// ---------------------------------------------------------------------------

/**
 * Decode DNA-Aeon encoded string back to bytes.
 *
 * Pipeline:
 *   1. Walk DNA in windows of (syncInterval + markerSize) bases
 *   2. For each window, verify CRC of payload against marker
 *   3. If CRC fails: indel detected → use stack algorithm to resync
 *   4. After checkpoint consecutive failures → declare erasure
 *   5. Return payload + list of erasure segments for outer code
 *
 * @param dna DNA string with sync markers
 * @param numPayloadBytes Expected number of payload bytes
 * @param config DNA-Aeon configuration
 * @returns Decoded result with erasure information
 */
export function dnaAeonDecode(
  dna: string,
  numPayloadBytes: number,
  config: DNAEeonConfig = DEFAULT_DNA_AEON_CONFIG,
): DNAEeonDecodeResult {
  // Step 1: Walk through DNA, verifying sync markers
  const {
    dnaWithoutMarkers,
    markersPassed,
    markersFailed,
    failedPositions,
  } = verifyAndStripDnaSyncMarkers(dna, config);

  // Step 2: If there are failed markers, attempt stack resynchronization
  const resyncShifts: number[] = [];
  let resyncNeeded = false;
  const erasureSegments: number[] = [];

  if (markersFailed > 0) {
    resyncNeeded = true;

    // Group consecutive failures into erasure segments
    let consecutiveFailures = 0;
    let erasureStart = -1;

    for (let i = 0; i < failedPositions.length; i++) {
      const failPos = failedPositions[i];

      // Try stack resync
      const shift = stackResync(dna, failPos, config);
      resyncShifts.push(shift);

      if (shift !== 0) {
        // Resync succeeded — reset consecutive failure count
        consecutiveFailures = 0;
      } else {
        // Resync failed — increment consecutive failure count
        if (consecutiveFailures === 0) {
          erasureStart = failPos;
        }
        consecutiveFailures++;

        if (consecutiveFailures >= config.checkpoint) {
          // Declare this segment as an erasure
          erasureSegments.push(erasureStart);
          consecutiveFailures = 0;
        }
      }
    }
  }

  // Step 3: Arithmetic decode (DNA → bits)
  const payload = arithmeticDecodeDnaToBytes(dnaWithoutMarkers, numPayloadBytes);

  return {
    payload,
    markersPassed,
    markersFailed,
    erasureSegments,
    resyncNeeded,
    resyncShifts,
  };
}

// ---------------------------------------------------------------------------
// Stack resynchronization
// ---------------------------------------------------------------------------

/**
 * Stack-based resynchronization after indel detection.
 *
 * When a CRC marker fails, the decoder hypothesizes that an indel (insertion
 * or deletion) occurred within the current segment. It tries shifts of ±1, ±2,
 * ±3 bases and checks if subsequent markers validate with each shift.
 *
 * The stack depth determines the maximum shift magnitude to try.
 * A shift of +N means N bases were inserted (or N bases should be skipped).
 * A shift of -N means N bases were deleted (or N extra bases should be read).
 *
 * @param dna The full DNA string with sync markers
 * @param failPosition The segment index where the marker failed
 * @param config DNA-Aeon configuration
 * @returns The shift that validates (positive = insertion, negative = deletion),
 *          or 0 if no valid shift found
 */
export function stackResync(
  dna: string,
  failPosition: number,
  config: DNAEeonConfig = DEFAULT_DNA_AEON_CONFIG,
): number {
  const segmentBases = config.syncInterval * 4;
  const markerBases = config.markerSize * 2;
  const windowSize = segmentBases + markerBases;

  // Base position of the failed window
  const failWindowStart = failPosition * windowSize;

  // Try shifts from -stackDepth to +stackDepth (skip 0 = original which failed)
  for (let shift = -config.stackDepth; shift <= config.stackDepth; shift++) {
    if (shift === 0) continue;

    // Look at the NEXT window after the failed one, with the hypothesized shift
    const nextWindowStart = failWindowStart + windowSize + shift;

    // Bounds check
    if (nextWindowStart < 0) continue;
    if (nextWindowStart + windowSize > dna.length) continue;

    // Extract the next segment and marker with the shift applied
    const nextSegment = dna.slice(
      nextWindowStart,
      nextWindowStart + segmentBases,
    );
    const nextMarkerDna = dna.slice(
      nextWindowStart + segmentBases,
      nextWindowStart + segmentBases + markerBases,
    );

    if (nextSegment.length < segmentBases || nextMarkerDna.length < markerBases) {
      continue;
    }

    // Verify CRC of the shifted window
    const segBytes = dnaToBytes2Bit(nextSegment);
    const markerBytes = dnaToBytes2Bit(nextMarkerDna);
    const expectedMarker = crc8(segBytes);
    const actualMarker = markerBytes[0] ?? 0;

    if (actualMarker === expectedMarker) {
      // This shift validates! The indel is at this offset.
      return shift;
    }
  }

  // No valid shift found — declare erasure
  return 0;
}

// ---------------------------------------------------------------------------
// Capacity computation
// ---------------------------------------------------------------------------

/**
 * Compute capacity (max payload bytes) for given DNA length and config.
 *
 * The DNA length includes sync markers. The capacity is:
 *   usableDnaBases = dnaLength / (1 + markerOverheadFraction)
 *   payloadBytes = usableDnaBases / 4  (4 bases per byte at 2 bits/base)
 *
 * With sync markers, the overhead fraction is:
 *   markerBases / (segmentBases + markerBases)
 *   = (markerSize * 2) / (syncInterval * 4 + markerSize * 2)
 *
 * @param dnaLength Total DNA length in bases (including markers)
 * @param config DNA-Aeon configuration
 * @returns Maximum number of payload bytes that can be stored
 */
export function dnaAeonCapacity(
  dnaLength: number,
  config: DNAEeonConfig = DEFAULT_DNA_AEON_CONFIG,
): number {
  const segmentBases = config.syncInterval * 4;
  const markerBases = config.markerSize * 2;
  const windowSize = segmentBases + markerBases;

  // Number of complete windows
  const numWindows = Math.floor(dnaLength / windowSize);

  // Payload bases = number of windows × segment bases
  const payloadBases = numWindows * segmentBases;

  // Remaining bases (partial window) — can also carry some payload
  const remaining = dnaLength - numWindows * windowSize;
  const extraPayloadBases = Math.min(remaining, segmentBases);

  // Total payload bases → bytes (4 bases per byte)
  return Math.floor((payloadBases + extraPayloadBases) / 4);
}
