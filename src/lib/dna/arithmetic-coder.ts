/**
 * Binary arithmetic coder for DNA-specific compression.
 *
 * This is the KEY differentiator between real DNA compressors (NAF, AGC, DeepGeCo,
 * MBGC2, JARVIS3) and DEFLATE-based approximations. Arithmetic coding encodes
 * symbols with fractional bits, achieving compression ratios closer to the
 * entropy limit than Huffman coding (which is limited to integer bits per symbol).
 *
 * Implementation uses 32-bit integer arithmetic (standard for production coders):
 *   - WHOLE  = 0x10000 (65536)
 *   - HALF   = 0x8000  (32768)
 *   - QUARTER = 0x4000  (16384)
 *   - PRECISION = 16 bits
 *
 * The adaptive frequency model updates counts after each symbol, allowing
 * the coder to track changing statistics in the input stream — critical for
 * DNA residual streams where local composition varies.
 *
 * References:
 *   - Witten, Neal, Cleary (1987). "Arithmetic Coding for Data Compression."
 *   - Sayood (2017). "Introduction to Data Compression." 5th ed.
 *   - Varshney et al. (2024). "A universal nucleotide archive format." arXiv.
 */

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/** Full range of the code value (2^16). */
const WHOLE = 0x10000; // 65536
/** Half of the full range. */
const HALF = 0x8000;   // 32768
/** Quarter of the full range. */
const QUARTER = 0x4000; // 16384

// ---------------------------------------------------------------------------
// Frequency Table
// ---------------------------------------------------------------------------

/**
 * Immutable frequency table for arithmetic coding.
 * Stores cumulative frequencies for efficient range lookups.
 */
export interface FrequencyTable {
  /** Total count of all symbols. */
  total: number;
  /** Cumulative frequency: cumFreq[i] = sum of freq[0..i-1]. cumFreq[0] = 0. */
  cumFreq: Uint32Array;
  /** Number of distinct symbols. */
  numSymbols: number;
}

/**
 * Create a FrequencyTable from raw frequency counts.
 *
 * @param freq Raw frequency counts (freq[i] = count of symbol i)
 * @returns FrequencyTable with cumulative frequencies
 */
export function createFrequencyTable(freq: Uint32Array): FrequencyTable {
  const numSymbols = freq.length;
  const cumFreq = new Uint32Array(numSymbols + 1);
  cumFreq[0] = 0;
  for (let i = 0; i < numSymbols; i++) {
    cumFreq[i + 1] = cumFreq[i] + freq[i];
  }
  return { total: cumFreq[numSymbols], cumFreq, numSymbols };
}

// ---------------------------------------------------------------------------
// Adaptive Frequency Model
// ---------------------------------------------------------------------------

/**
 * Adaptive order-0 frequency model.
 *
 * Updates symbol counts after each observation, allowing the arithmetic
 * coder to track changing statistics. Uses increment rescaling to prevent
 * counts from overflowing.
 *
 * For DNA compression, the 4-symbol alphabet {A,C,G,T} benefits greatly
 * from adaptive modeling because base frequencies vary across genomic regions.
 */
export class AdaptiveModel {
  /** Number of distinct symbols. */
  readonly numSymbols: number;
  /** Current frequency counts. */
  freq: Uint32Array;
  /** Cumulative frequencies (recomputed as needed). */
  cumFreq: Uint32Array;
  /** Total count. */
  total: number;
  /** Maximum count before rescaling. */
  private maxCount: number;

  /**
   * @param numSymbols Number of distinct symbols (e.g., 4 for DNA)
   * @param initialValue Initial count per symbol (default: 1 for Laplace smoothing)
   * @param maxCount Maximum count before rescaling (default: 8192)
   */
  constructor(numSymbols: number, initialValue: number = 1, maxCount: number = 8192) {
    this.numSymbols = numSymbols;
    this.maxCount = maxCount;
    this.freq = new Uint32Array(numSymbols).fill(initialValue);
    this.cumFreq = new Uint32Array(numSymbols + 1);
    this.total = numSymbols * initialValue;
    this.recomputeCumFreq();
  }

  /** Recompute cumulative frequencies from raw counts. */
  private recomputeCumFreq(): void {
    this.cumFreq[0] = 0;
    for (let i = 0; i < this.numSymbols; i++) {
      this.cumFreq[i + 1] = this.cumFreq[i] + this.freq[i];
    }
    this.total = this.cumFreq[this.numSymbols];
  }

  /** Get the current frequency table. */
  getTable(): FrequencyTable {
    return { total: this.total, cumFreq: this.cumFreq, numSymbols: this.numSymbols };
  }

  /** Update the model after observing a symbol. */
  update(symbol: number): void {
    this.freq[symbol]++;
    this.total++;

    // Rescale if any count exceeds maxCount
    if (this.total > this.maxCount) {
      this.total = 0;
      for (let i = 0; i < this.numSymbols; i++) {
        this.freq[i] = Math.max(1, this.freq[i] >>> 1); // halve, min 1
        this.total += this.freq[i];
      }
    }

    this.recomputeCumFreq();
  }
}

// ---------------------------------------------------------------------------
// Adaptive Order-k Model
// ---------------------------------------------------------------------------

/**
 * Adaptive order-k context model for DNA sequences.
 *
 * Maintains separate frequency tables for each context (previous k symbols),
 * allowing the arithmetic coder to exploit higher-order dependencies.
 * For DNA with order-2, this gives 16 contexts (4^2), capturing
 * dinucleotide frequencies like CG suppression in vertebrate genomes.
 */
export class AdaptiveContextModel {
  /** Context order. */
  readonly order: number;
  /** Number of symbols (4 for DNA). */
  readonly numSymbols: number;
  /** Models indexed by context. */
  models: Map<number, AdaptiveModel>;
  /** Current context (encoded as base-4 number). */
  context: number;
  /** Number of possible contexts. */
  readonly numContexts: number;

  constructor(order: number, numSymbols: number = 4) {
    this.order = order;
    this.numSymbols = numSymbols;
    this.numContexts = numSymbols ** order;
    this.models = new Map();
    this.context = 0;

    // Pre-create all models
    for (let i = 0; i < this.numContexts; i++) {
      this.models.set(i, new AdaptiveModel(numSymbols));
    }
  }

  /** Get the frequency table for the current context. */
  getTable(): FrequencyTable {
    return (this.models.get(this.context) ?? this.models.get(0)!).getTable();
  }

  /** Update the model after observing a symbol. */
  update(symbol: number): void {
    const model = this.models.get(this.context);
    if (model) model.update(symbol);

    // Update context: shift and add current symbol
    this.context = ((this.context * this.numSymbols) + symbol) % this.numContexts;
  }

  /** Reset context to initial state. */
  reset(): void {
    this.context = 0;
  }
}

// ---------------------------------------------------------------------------
// Arithmetic Encoder
// ---------------------------------------------------------------------------

/**
 * Binary arithmetic encoder.
 *
 * Encodes symbols to a compressed bit stream using frequency tables.
 * Each symbol narrows the current interval [low, low+range) based on
 * its cumulative frequency. When the interval shrinks below half,
 * bits are output to the stream.
 *
 * Usage:
 *   const enc = new ArithmeticEncoder();
 *   const model = new AdaptiveModel(4);
 *   for (const sym of symbols) {
 *     enc.encode(sym, model.getTable());
 *     model.update(sym);
 *   }
 *   const compressed = enc.finish();
 */
export class ArithmeticEncoder {
  /** Lower bound of current interval. */
  private low: number = 0;
  /** Range of current interval (high - low). */
  private range: number = WHOLE;
  /** Pending bits to output (for the E3 mapping / underflow count). */
  private pendingBits: number = 0;
  /** Output byte buffer. */
  private buffer: number[] = [];
  /** Current byte being assembled (bits accumulated from MSB). */
  private currentByte: number = 0;
  /** Number of bits written to currentByte. */
  private bitsInByte: number = 0;

  /**
   * Encode a single symbol.
   *
   * @param symbol Symbol index (0-based)
   * @param table Frequency table for the current context
   */
  encode(symbol: number, table: FrequencyTable): void {
    const { total, cumFreq } = table;

    // Narrow the interval
    const range = this.range;
    const newLow = this.low + Math.floor((range * cumFreq[symbol]) / total);
    const newRange = Math.floor((range * cumFreq[symbol + 1]) / total) - Math.floor((range * cumFreq[symbol]) / total);

    this.low = newLow;
    this.range = Math.max(1, newRange); // Ensure range >= 1

    // Renormalize: output bits while interval is in upper or lower half
    while (true) {
      if (this.low + this.range <= HALF) {
        // Interval is entirely in lower half — output 0
        this.outputBitPlusPending(0);
      } else if (this.low >= HALF) {
        // Interval is entirely in upper half — output 1
        this.outputBitPlusPending(1);
        this.low -= HALF;
      } else if (this.low >= QUARTER && this.low + this.range <= 3 * QUARTER) {
        // E3 mapping: interval straddles midpoint
        this.pendingBits++;
        this.low -= QUARTER;
      } else {
        // Cannot normalize further
        break;
      }

      this.low *= 2;
      this.range *= 2;

      // Prevent overflow
      if (this.range > WHOLE) this.range = WHOLE;
    }
  }

  /**
   * Output a bit plus all pending bits (for E3 mapping).
   */
  private outputBitPlusPending(bit: number): void {
    this.emitBit(bit);
    for (let i = 0; i < this.pendingBits; i++) {
      this.emitBit(bit ^ 1);
    }
    this.pendingBits = 0;
  }

  /**
   * Emit a single bit to the output stream.
   */
  private emitBit(bit: number): void {
    this.currentByte = (this.currentByte << 1) | bit;
    this.bitsInByte++;
    if (this.bitsInByte === 8) {
      this.buffer.push(this.currentByte);
      this.currentByte = 0;
      this.bitsInByte = 0;
    }
  }

  /**
   * Finish encoding and return the compressed data.
   *
   * Flushes the pending bits and the final state of the interval.
   *
   * @returns Compressed bytes (Uint8Array)
   */
  finish(): Uint8Array {
    // Output enough bits to uniquely identify the final interval
    this.pendingBits++;
    if (this.low < QUARTER) {
      this.outputBitPlusPending(0);
    } else {
      this.outputBitPlusPending(1);
    }

    // Pad remaining bits in the current byte
    if (this.bitsInByte > 0) {
      this.currentByte <<= (8 - this.bitsInByte);
      this.buffer.push(this.currentByte);
    }

    return new Uint8Array(this.buffer);
  }
}

// ---------------------------------------------------------------------------
// Arithmetic Decoder
// ---------------------------------------------------------------------------

/**
 * Binary arithmetic decoder.
 *
 * Decodes symbols from a compressed bit stream using frequency tables.
 * Maintains a code value read from the stream and narrows it using
 * the same interval arithmetic as the encoder.
 *
 * Usage:
 *   const dec = new ArithmeticDecoder(compressed);
 *   const model = new AdaptiveModel(4);
 *   for (let i = 0; i < numSymbols; i++) {
 *     const sym = dec.decode(model.getTable());
 *     model.update(sym);
 *   }
 */
export class ArithmeticDecoder {
  /** Lower bound of current interval. */
  private low: number = 0;
  /** Range of current interval. */
  private range: number = WHOLE;
  /** Code value read from the stream. */
  private code: number = 0;
  /** Input bytes. */
  private input: Uint8Array;
  /** Current bit position in the input. */
  private bitPos: number = 0;

  /**
   * @param data Compressed bytes from ArithmeticEncoder.finish()
   */
  constructor(data: Uint8Array) {
    this.input = data;

    // Read the initial code value (first PRECISION bits)
    for (let i = 0; i < 16; i++) {
      this.code = (this.code << 1) | this.readBit();
    }
  }

  /**
   * Read a single bit from the input stream.
   */
  private readBit(): number {
    const byteIdx = this.bitPos >>> 3;
    const bitIdx = 7 - (this.bitPos & 7); // MSB first
    this.bitPos++;
    if (byteIdx < this.input.length) {
      return (this.input[byteIdx] >>> bitIdx) & 1;
    }
    return 0; // Past end of input — return 0 (standard practice)
  }

  /**
   * Decode a single symbol.
   *
   * @param table Frequency table for the current context
   * @returns Decoded symbol index (0-based)
   */
  decode(table: FrequencyTable): number {
    const { total, cumFreq, numSymbols } = table;

    // Find the symbol whose cumulative frequency range contains the scaled code value
    const scaledValue = Math.floor(((this.code - this.low + 1) * total - 1) / this.range);

    // Binary search for the symbol
    let lo = 0;
    let hi = numSymbols - 1;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (cumFreq[mid + 1] <= scaledValue) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    const symbol = lo;

    // Narrow the interval (same as encoder)
    const newLow = this.low + Math.floor((this.range * cumFreq[symbol]) / total);
    const newRange = Math.floor((this.range * cumFreq[symbol + 1]) / total) - Math.floor((this.range * cumFreq[symbol]) / total);

    this.low = newLow;
    this.range = Math.max(1, newRange);

    // Renormalize
    while (true) {
      if (this.low + this.range <= HALF) {
        // Lower half — no output
      } else if (this.low >= HALF) {
        // Upper half
        this.low -= HALF;
        this.code -= HALF;
      } else if (this.low >= QUARTER && this.low + this.range <= 3 * QUARTER) {
        // E3 mapping
        this.low -= QUARTER;
        this.code -= QUARTER;
      } else {
        break;
      }

      this.low *= 2;
      this.range *= 2;
      this.code = (this.code << 1) | this.readBit();

      if (this.range > WHOLE) this.range = WHOLE;
    }

    return symbol;
  }
}

// ---------------------------------------------------------------------------
// Convenience: Encode/Decode with adaptive order-0 model
// ---------------------------------------------------------------------------

/**
 * Encode an array of symbols (0-based integers) using arithmetic coding
 * with an adaptive order-0 model.
 *
 * @param symbols Array of symbol indices
 * @param numSymbols Number of distinct symbols (e.g., 4 for DNA)
 * @returns Compressed bytes
 */
export function arithmeticEncode(symbols: Uint8Array, numSymbols: number): Uint8Array {
  const enc = new ArithmeticEncoder();
  const model = new AdaptiveModel(numSymbols);

  for (let i = 0; i < symbols.length; i++) {
    enc.encode(symbols[i], model.getTable());
    model.update(symbols[i]);
  }

  return enc.finish();
}

/**
 * Decode compressed data back to symbols using arithmetic coding
 * with an adaptive order-0 model.
 *
 * @param data Compressed bytes from arithmeticEncode
 * @param numSymbols Number of distinct symbols
 * @param count Number of symbols to decode
 * @returns Decoded symbol indices
 */
export function arithmeticDecode(data: Uint8Array, numSymbols: number, count: number): Uint8Array {
  const dec = new ArithmeticDecoder(data);
  const model = new AdaptiveModel(numSymbols);
  const symbols = new Uint8Array(count);

  for (let i = 0; i < count; i++) {
    symbols[i] = dec.decode(model.getTable());
    model.update(symbols[i]);
  }

  return symbols;
}

// ---------------------------------------------------------------------------
// AdaptiveFrequencyModel — compat wrapper for dna-compress-real.ts
// ---------------------------------------------------------------------------

/**
 * Adaptive frequency model with order-k context and serialization.
 *
 * This is a compatibility wrapper used by dna-compress-real.ts. It combines
 * AdaptiveModel (order-0) and AdaptiveContextModel (order-k) into a single
 * class with serialize/deserialize support.
 */
export class AdaptiveFrequencyModel {
  readonly numSymbols: number;
  readonly order: number;
  private inner: AdaptiveModel | AdaptiveContextModel;

  constructor(numSymbols: number, order: number = 0) {
    this.numSymbols = numSymbols;
    this.order = order;
    if (order === 0) {
      this.inner = new AdaptiveModel(numSymbols);
    } else {
      this.inner = new AdaptiveContextModel(order, numSymbols);
    }
  }

  /** Get the current frequency table. */
  getFrequencyTable(): FrequencyTable {
    return this.inner.getTable();
  }

  /** Update the model after observing a symbol. */
  update(symbol: number): void {
    this.inner.update(symbol);
  }

  /**
   * Serialize the model state to bytes.
   * Format: [order(1)] [numSymbols(2 LE)] [freq0(4 LE)] [freq1(4 LE)] ...
   * For order-k: also includes all context tables.
   */
  serialize(): Uint8Array {
    if (this.order === 0) {
      const m = this.inner as AdaptiveModel;
      // [order(1)] [numSymbols(2 LE)] [freq0(4 LE)] ... [freqN(4 LE)]
      const out = new Uint8Array(3 + m.numSymbols * 4);
      out[0] = this.order;
      const dv = new DataView(out.buffer, out.byteOffset, out.byteLength);
      dv.setUint16(1, m.numSymbols, true);
      for (let i = 0; i < m.numSymbols; i++) {
        dv.setUint32(3 + i * 4, m.freq[i], true);
      }
      return out;
    } else {
      const m = this.inner as AdaptiveContextModel;
      // Serialize all context tables
      // [order(1)] [numSymbols(2 LE)] [numContexts(4 LE)] [ctx0_freq0(4 LE)] ... per context
      const numCtx = m.numContexts;
      const ctxSize = m.numSymbols * 4;
      const out = new Uint8Array(7 + numCtx * ctxSize);
      out[0] = this.order;
      const dv = new DataView(out.buffer, out.byteOffset, out.byteLength);
      dv.setUint16(1, m.numSymbols, true);
      dv.setUint32(3, numCtx, true);
      let off = 7;
      for (let ctx = 0; ctx < numCtx; ctx++) {
        const model = m.models.get(ctx)!;
        for (let i = 0; i < m.numSymbols; i++) {
          dv.setUint32(off, model.freq[i], true);
          off += 4;
        }
      }
      return out;
    }
  }

  /**
   * Deserialize a model from bytes.
   */
  static deserialize(data: Uint8Array, offset: number = 0): { model: AdaptiveFrequencyModel; bytesRead: number } {
    const order = data[offset];
    const dv = new DataView(data.buffer, data.byteOffset, data.byteLength);
    const numSymbols = dv.getUint16(offset + 1, true);

    const result = new AdaptiveFrequencyModel(numSymbols, order);

    if (order === 0) {
      const m = result.inner as AdaptiveModel;
      for (let i = 0; i < numSymbols; i++) {
        m.freq[i] = dv.getUint32(offset + 3 + i * 4, true);
      }
      m.total = 0;
      for (let i = 0; i < numSymbols; i++) m.total += m.freq[i];
      // Recompute cumFreq
      m.cumFreq[0] = 0;
      for (let i = 0; i < numSymbols; i++) m.cumFreq[i + 1] = m.cumFreq[i] + m.freq[i];
      return { model: result, bytesRead: 3 + numSymbols * 4 };
    } else {
      const m = result.inner as AdaptiveContextModel;
      const numCtx = dv.getUint32(offset + 3, true);
      let off = offset + 7;
      for (let ctx = 0; ctx < numCtx; ctx++) {
        const model = m.models.get(ctx)!;
        model.total = 0;
        for (let i = 0; i < numSymbols; i++) {
          model.freq[i] = dv.getUint32(off, true);
          model.total += model.freq[i];
          off += 4;
        }
        model.cumFreq[0] = 0;
        for (let i = 0; i < numSymbols; i++) model.cumFreq[i + 1] = model.cumFreq[i] + model.freq[i];
      }
      return { model: result, bytesRead: 7 + numCtx * numSymbols * 4 };
    }
  }
}

/**
 * Encode with an adaptive order-k context model.
 *
 * @param symbols Array of symbol indices
 * @param numSymbols Number of distinct symbols
 * @param order Context order (0 = no context, 1 = previous symbol, 2 = dinucleotide)
 * @returns Compressed bytes (includes header with order and count)
 */
export function arithmeticEncodeContext(symbols: Uint8Array, numSymbols: number, order: number): Uint8Array {
  const enc = new ArithmeticEncoder();
  const model = new AdaptiveContextModel(order, numSymbols);

  for (let i = 0; i < symbols.length; i++) {
    enc.encode(symbols[i], model.getTable());
    model.update(symbols[i]);
  }

  const compressed = enc.finish();

  // Header: [order(1)] [numSymbols(1)] [count(4 LE)] [compressed...]
  const out = new Uint8Array(6 + compressed.length);
  out[0] = order;
  out[1] = numSymbols;
  const view = new DataView(out.buffer, out.byteOffset, out.byteLength);
  view.setUint32(2, symbols.length, true);
  out.set(compressed, 6);

  return out;
}

/**
 * Decode with an adaptive order-k context model.
 *
 * @param data Compressed bytes from arithmeticEncodeContext
 * @returns Decoded symbol indices
 */
export function arithmeticDecodeContext(data: Uint8Array): Uint8Array {
  const order = data[0];
  const numSymbols = data[1];
  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
  const count = view.getUint32(2, true);
  const compressed = data.slice(6);

  const dec = new ArithmeticDecoder(compressed);
  const model = new AdaptiveContextModel(order, numSymbols);
  const symbols = new Uint8Array(count);

  for (let i = 0; i < count; i++) {
    symbols[i] = dec.decode(model.getTable());
    model.update(symbols[i]);
  }

  return symbols;
}
