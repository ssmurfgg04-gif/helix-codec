/**
 * True Arithmetic Coding for Constrained DNA — DNA-Aeon Architecture
 *
 * KEY INSIGHT: DNA-Aeon SWAPS the encoder/decoder roles:
 *   - ENCODER = arithmetic DECODER (input: bits → output: DNA)
 *   - DECODER = arithmetic ENCODER (input: DNA → output: bits)
 *
 * The arithmetic decoder reads bits from the input and outputs symbols (DNA bases).
 * The arithmetic encoder reads symbols (DNA bases) and outputs bits.
 *
 * This is FULLY REVERSIBLE with ZERO ambiguity because both sides track
 * the same interval state.
 *
 * Implementation based on Nayuki's reference arithmetic coder:
 *   https://github.com/nayuki/Reference-arithmetic-coding
 *
 * Reference:
 *   - Welzel et al. (2023). Nature Comms 14:628 (DNA-Aeon).
 *   - Witten, Neal, Cleary (1987). CACM 30:6.
 *   - Nayuki (2020). Reference arithmetic coding.
 */

import { Base, gcContent, maxHomopolymerRun } from "./mapping";

const BASES: Base[] = ["A", "C", "G", "T"];
const BASE_TO_IDX: Record<Base, number> = { A: 0, C: 1, G: 2, T: 3 };

// 24-bit precision (safe for JavaScript — 2^24 * 4 = 2^26 < 2^53)
const NUM_STATE_BITS = 24;
const STATE_MASK = 0xFFFFFF;
const HALF_STATE = 0x800000;
const QUARTER_STATE = 0x400000;

/**
 * Bit input stream — reads bits from a byte array MSB-first.
 */
class BitInputStream {
  private bits: number[];
  private idx = 0;

  constructor(data: Uint8Array) {
    this.bits = [];
    for (let i = 0; i < data.length; i++) {
      for (let b = 7; b >= 0; b--) {
        this.bits.push((data[i] >> b) & 1);
      }
    }
  }

  readBit(): number {
    if (this.idx < this.bits.length) return this.bits[this.idx++];
    return 0; // pad with zeros
  }

  hasBits(): boolean {
    return this.idx < this.bits.length;
  }
}

/**
 * Bit output stream — collects bits and converts to bytes.
 */
class BitOutputStream {
  private bits: number[] = [];

  writeBit(bit: number) {
    this.bits.push(bit & 1);
  }

  toBytes(numBytes?: number): Uint8Array {
    const n = numBytes ?? Math.ceil(this.bits.length / 8);
    const out = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      let byte = 0;
      for (let b = 0; b < 8; b++) {
        byte = (byte << 1) | (this.bits[i * 8 + b] ?? 0);
      }
      out[i] = byte;
    }
    return out;
  }
}

/**
 * Arithmetic DECODER — reads bits, outputs symbols.
 * This is used as the ENCODER in DNA-Aeon (bits → DNA).
 *
 * The decoder maintains a "code" register filled from the input bitstream.
 * At each step, it determines which symbol the code corresponds to,
 * narrows the interval, and renormalizes.
 */
class ArithmeticDecoder {
  private low = 0;
  private high = STATE_MASK;
  private code = 0;
  private input: BitInputStream;
  private pendingOutput: Base[] = [];

  constructor(input: BitInputStream) {
    this.input = input;
    // Fill code with initial bits — mask to 24 bits (>>> 0 only truncates to 32)
    for (let i = 0; i < NUM_STATE_BITS; i++) {
      this.code = (((this.code << 1) | this.input.readBit()) & STATE_MASK) >>> 0;
    }
  }

  /**
   * Read a symbol from the arithmetic stream.
   * The symbol is selected from `numSymbols` equally-probable options.
   * Returns the symbol index (0 to numSymbols-1).
   */
  read(numSymbols: number): number {
    const range = (this.high - this.low + 1) >>> 0;
    const symbolWidth = Math.floor(range / numSymbols);

    // Determine which symbol the code corresponds to
    const valueOffset = ((this.code - this.low) >>> 0);
    let symbolIdx = Math.floor(valueOffset / symbolWidth);
    if (symbolIdx >= numSymbols) symbolIdx = numSymbols - 1;

    // Narrow the interval
    this.high = (this.low + symbolWidth * (symbolIdx + 1) - 1) >>> 0;
    this.low = (this.low + symbolWidth * symbolIdx) >>> 0;

    // Renormalize — standard Nayuki reference arithmetic coder conditions.
    // v56 fix: The previous XOR-based condition
    //   `(this.high ^ HALF_STATE) > HALF_STATE || ((this.high ^ HALF_STATE) === 0)`
    // was WRONG. It evaluates to TRUE when high === HALF_STATE (due to the `=== 0`
    // check), but the standard decoder should NOT take the first branch when
    // high === HALF. Instead, it should fall through to the underflow branch
    // (low >= QUARTER && high < 3*QUARTER) or break.
    //
    // When high === HALF_STATE and the buggy code took the first branch,
    // it shifted high = (HALF << 1) | 1 = 0x1000001 (25 bits, > STATE_MASK),
    // corrupting the interval state. This caused the encoder to emit wrong DNA
    // symbols, which the decoder then couldn't recover — manifesting as
    // "termination corruption" of the last 1-3 bytes of each block.
    //
    // The fix: use the straightforward `this.high < HALF_STATE` comparison,
    // exactly matching the encoder (ArithmeticEncoder.write) and the Rust
    // decoder (arithmetic_dna_to_bytes in lib.rs).
    while (true) {
      if (this.high < HALF_STATE) {
        // high < HALF — encoder output 0, decoder shifts in next bit
      } else if (this.low >= HALF_STATE) {
        // low >= HALF — encoder output 1
        this.code = (this.code - HALF_STATE) >>> 0;
        this.low = (this.low - HALF_STATE) >>> 0;
        this.high = (this.high - HALF_STATE) >>> 0;
      } else if (this.low >= QUARTER_STATE && this.high < (QUARTER_STATE * 3)) {
        // underflow — encoder incremented pending bits
        this.code = (this.code - QUARTER_STATE) >>> 0;
        this.low = (this.low - QUARTER_STATE) >>> 0;
        this.high = (this.high - QUARTER_STATE) >>> 0;
      } else {
        break;
      }
      this.low = ((this.low << 1) & STATE_MASK) >>> 0;
      this.high = (((this.high << 1) | 1) & STATE_MASK) >>> 0;
      this.code = (((this.code << 1) | this.input.readBit()) & STATE_MASK) >>> 0;
    }

    return symbolIdx;
  }
}

/**
 * Arithmetic ENCODER — reads symbols, outputs bits.
 * This is used as the DECODER in DNA-Aeon (DNA → bits).
 *
 * The encoder maintains an interval [low, high). At each step, it narrows
 * the interval based on the input symbol, and renormalizes by outputting bits.
 */
class ArithmeticEncoder {
  private low = 0;
  private high = STATE_MASK;
  private pendingBits = 0;
  private output: BitOutputStream;

  constructor(output: BitOutputStream) {
    this.output = output;
  }

  /**
   * Write a symbol to the arithmetic stream.
   * The symbol is selected from `numSymbols` equally-probable options.
   */
  write(symbolIdx: number, numSymbols: number) {
    const range = (this.high - this.low + 1) >>> 0;
    const symbolWidth = Math.floor(range / numSymbols);

    this.high = (this.low + symbolWidth * (symbolIdx + 1) - 1) >>> 0;
    this.low = (this.low + symbolWidth * symbolIdx) >>> 0;

    // Renormalize
    while (true) {
      if (this.high < HALF_STATE) {
        this.outputBit(0);
      } else if (this.low >= HALF_STATE) {
        this.outputBit(1);
        this.low = (this.low - HALF_STATE) >>> 0;
        this.high = (this.high - HALF_STATE) >>> 0;
      } else if (this.low >= QUARTER_STATE && this.high < (QUARTER_STATE * 3)) {
        this.pendingBits++;
        this.low = (this.low - QUARTER_STATE) >>> 0;
        this.high = (this.high - QUARTER_STATE) >>> 0;
      } else {
        break;
      }
      this.low = ((this.low << 1) & STATE_MASK) >>> 0;
      this.high = (((this.high << 1) | 1) & STATE_MASK) >>> 0;
    }
  }

  private outputBit(bit: number) {
    this.output.writeBit(bit);
    while (this.pendingBits > 0) {
      this.output.writeBit(1 - bit);
      this.pendingBits--;
    }
  }

  finish() {
    this.pendingBits++;
    if (this.low < QUARTER_STATE) {
      this.outputBit(0);
    } else {
      this.outputBit(1);
    }
  }
}

/**
 * Encode bytes to DNA using true arithmetic coding.
 *
 * ENCODER = ArithmeticDecoder (reads bits, outputs DNA symbols).
 * At each position, the allowed symbols depend on the homopolymer state.
 *
 * @param data Input bytes
 * @param maxHomopolymer Maximum allowed homopolymer run (default 3)
 * @param targetLen Target DNA length (pad if needed)
 * @returns DNA string with homopolymer runs ≤ maxHomopolymer
 */
export function bytesToArithmeticDna(
  data: Uint8Array,
  maxHomopolymer: number = 3,
  targetLen?: number,
): string {
  const input = new BitInputStream(data);
  const numBases = targetLen ?? Math.ceil(data.length * 8 / 1.9);
  const decoder = new ArithmeticDecoder(input);

  const result: Base[] = [];
  let prev: Base = "A";
  let runLength = 0;

  for (let pos = 0; pos < numBases; pos++) {
    // Determine allowed symbols
    let allowed: Base[];
    if (runLength >= maxHomopolymer) {
      allowed = BASES.filter((b) => b !== prev);
    } else {
      allowed = [...BASES];
    }
    const numAllowed = allowed.length;

    // Read a symbol from the arithmetic decoder
    const symbolIdx = decoder.read(numAllowed);
    const base = allowed[Math.min(symbolIdx, numAllowed - 1)];

    result.push(base);

    // Update homopolymer state
    if (base === prev) {
      runLength++;
    } else {
      runLength = 1;
      prev = base;
    }
  }

  let dna = result.join("");

  // Pad to target length
  if (targetLen && dna.length < targetLen) {
    const parts = dna.split("");
    while (parts.length < targetLen) {
      const lastBase = parts[parts.length - 1] as Base;
      const allowed2 = BASES.filter((b) => b !== lastBase);
      parts.push(allowed2[0]);
    }
    dna = parts.join("");
  }
  if (targetLen && dna.length > targetLen) {
    dna = dna.slice(0, targetLen);
  }

  return dna;
}

/**
 * v57 capacity constant: effective arithmetic coding rate in bits/nt.
 *
 * Theoretical max is 2.0 bits/nt (4 equally-likely bases), but homopolymer
 * constraints reduce this. v56 used 1.85 bits/nt to leave slack for the
 * arithmetic coder's `finish()` termination — but this was too conservative:
 * with per-block CRC-8 overhead (1 byte per block), the LDPC codeword
 * didn't fit for typical oligo sizes (700nt → 156B LDPC vs 150B capacity).
 *
 * v57: bumped to 1.95 bits/nt. The termination corruption (1 byte per block)
 * is detected by the per-block CRC-8 and corrected by the LDPC erasure
 * decoder (peeling + Gaussian elimination over GF(2)). This gives ~5% more
 * capacity, enough to fit the LDPC codeword with margin.
 *
 * The trade-off: ~1 byte per block is marked as erased at decode time,
 * which the LDPC erasure decoder corrects using the inner parity bytes
 * (typically 4 bytes → can correct 4 byte erasures).
 */
export const ARITH_CAPACITY_RATE = 1.95;

/**
 * Block-wise arithmetic encode: split data into chunks, encode each chunk into
 * a DNA block of `blockSize` nucleotides with FRESH arithmetic state.
 *
 * This confines error propagation to a single block (matching the Rust
 * `arithmetic_dna_to_bytes_blocked` decoder). With block_size=20, one DNA
 * error corrupts at most 5 bytes, not the entire stream.
 *
 * v56: Each block encodes floor(blockSize * 1.85 / 8) bytes (was blockSize/4).
 * The 1.85 bits/nt rate leaves ~7.5% slack for arithmetic termination.
 *
 * @param data Input bytes
 * @param maxHomopolymer Max homopolymer run (default 3)
 * @param targetLen Target DNA length (must be multiple of blockSize)
 * @param blockSize Block size in nucleotides (default 20). 0 = no blocking.
 * @returns DNA string with homopolymer runs ≤ maxHomopolymer
 */
export function bytesToArithmeticDnaBlocked(
  data: Uint8Array,
  maxHomopolymer: number = 3,
  targetLen?: number,
  blockSize: number = 20,
): string {
  if (blockSize === 0) {
    return bytesToArithmeticDna(data, maxHomopolymer, targetLen);
  }

  // v56: Use 1.85 bits/nt capacity (was 2.0 = blockSize/4).
  // The slack absorbs the arithmetic `finish()` termination overhead.
  const bytesPerBlock = Math.max(1, Math.floor((blockSize * ARITH_CAPACITY_RATE) / 8));
  const numBlocks = Math.ceil(data.length / bytesPerBlock);
  const totalDnaLen = numBlocks * blockSize;

  // Adjust targetLen to be a multiple of blockSize
  const adjustedTarget = targetLen ?? totalDnaLen;
  const paddedTarget = Math.ceil(adjustedTarget / blockSize) * blockSize;

  const result: Base[] = [];

  for (let block = 0; block < numBlocks; block++) {
    const blockStart = block * bytesPerBlock;
    const blockEnd = Math.min(blockStart + bytesPerBlock, data.length);
    const blockData = data.slice(blockStart, blockEnd);

    // Encode this block with FRESH arithmetic state
    // (new BitInputStream + new ArithmeticDecoder = reset state)
    const input = new BitInputStream(blockData);
    const decoder = new ArithmeticDecoder(input);

    let prev: Base = "A";
    let runLength = 0;

    for (let pos = 0; pos < blockSize; pos++) {
      let allowed: Base[];
      if (runLength >= maxHomopolymer) {
        allowed = BASES.filter((b) => b !== prev);
      } else {
        allowed = [...BASES];
      }
      const numAllowed = allowed.length;

      const symbolIdx = decoder.read(numAllowed);
      const base = allowed[Math.min(symbolIdx, numAllowed - 1)];

      result.push(base);

      if (base === prev) {
        runLength++;
      } else {
        runLength = 1;
        prev = base;
      }
    }
  }

  let dna = result.join("");

  // Pad to target length (pad with alternating bases to avoid homopolymers)
  if (paddedTarget > dna.length) {
    const parts = dna.split("");
    while (parts.length < paddedTarget) {
      const lastBase = parts[parts.length - 1] as Base;
      const allowed2 = BASES.filter((b) => b !== lastBase);
      parts.push(allowed2[0]);
    }
    dna = parts.join("");
  }
  if (paddedTarget < dna.length) {
    dna = dna.slice(0, paddedTarget);
  }

  return dna;
}

/**
 * CRC-8 (poly 0x07, init 0x00) — matches Rust crc8().
 * Used for per-block CRC sync markers in DNA-Aeon-style error confinement.
 */
function crc8(data: Uint8Array): number {
  let crc = 0;
  for (const b of data) {
    crc ^= b;
    for (let i = 0; i < 8; i++) {
      if (crc & 0x80) crc = ((crc << 1) ^ 0x07) & 0xff;
      else crc = (crc << 1) & 0xff;
    }
  }
  return crc & 0xff;
}

/**
 * Block-wise arithmetic encode WITH per-block CRC-8 sync markers.
 *
 * DNA-Aeon approach: each block of `blockSize` nucleotides encodes
 * `bytesPerBlock` data bytes + 1 CRC-8 byte. At decode, each block is
 * verified independently. Blocks with failed CRC are marked as erased.
 *
 * This confines error propagation to single blocks and allows the decoder
 * to detect which blocks are corrupted (rather than blindly trusting
 * corrupted data that would fail LDPC).
 *
 * @param data Input bytes (the LDPC codeword + CRC-16, total = totalInnerBytes)
 * @param maxHomopolymer Max homopolymer run
 * @param targetLen Target DNA length (total, across all blocks)
 * @param blockSize Block size in nucleotides (default 20). 0 = no blocking.
 * @returns DNA string with homopolymer runs ≤ maxHomopolymer
 */
export function bytesToArithmeticDnaCrc(
  data: Uint8Array,
  maxHomopolymer: number = 3,
  targetLen?: number,
  blockSize: number = 20,
): string {
  if (blockSize === 0) {
    return bytesToArithmeticDna(data, maxHomopolymer, targetLen);
  }

  // v56: Use 1.85 bits/nt capacity (was 2.0 = blockSize/4).
  // The slack from 1.85 (vs 2.0) absorbs the arithmetic `finish()` termination
  // overhead, so we only need 1 byte for CRC (no separate padding bytes).
  //   bytesPerBlockTotal = floor(blockSize * 1.85 / 8)  (total bytes decoded per block)
  //   bytesPerBlockData = bytesPerBlockTotal - 1        (data bytes, excl 1 CRC byte)
  //
  // Before v56: bytesPerBlockTotal = blockSize/4, bytesPerBlockData = total - 4
  // (1 CRC + 3 padding). The 3 padding bytes were a workaround for the
  // termination bug — now that we use 1.85 bits/nt with built-in slack,
  // we don't need the padding. This recovers ~3 bytes per block of capacity.
  const bytesPerBlockTotal = Math.max(2, Math.floor((blockSize * ARITH_CAPACITY_RATE) / 8));
  const bytesPerBlockData = bytesPerBlockTotal - 1; // 1 byte for CRC-8, no padding needed
  if (bytesPerBlockData <= 0) {
    return bytesToArithmeticDna(data, maxHomopolymer, targetLen);
  }

  const numBlocks = Math.ceil(data.length / bytesPerBlockData);
  const totalDnaLen = numBlocks * blockSize;

  const adjustedTarget = targetLen ?? totalDnaLen;
  const paddedTarget = Math.ceil(adjustedTarget / blockSize) * blockSize;

  const result: Base[] = [];

  for (let block = 0; block < numBlocks; block++) {
    const blockStart = block * bytesPerBlockData;
    const blockEnd = Math.min(blockStart + bytesPerBlockData, data.length);
    const blockData = data.slice(blockStart, blockEnd);

    // Append CRC-8 of the data bytes
    const crc = crc8(blockData);
    const blockWithCrc = new Uint8Array(bytesPerBlockTotal);
    // Put CRC FIRST so arithmetic termination corruption affects padding, not CRC
    blockWithCrc[0] = crc;
    blockWithCrc.set(blockData, 1);
    // Pad to bytesPerBlockData if needed (last block may be shorter)
    for (let i = 1 + blockData.length; i < bytesPerBlockTotal; i++) blockWithCrc[i] = 0;

    // Encode this block with FRESH arithmetic state
    const input = new BitInputStream(blockWithCrc);
    const decoder = new ArithmeticDecoder(input);

    let prev: Base = "A";
    let runLength = 0;

    for (let pos = 0; pos < blockSize; pos++) {
      let allowed: Base[];
      if (runLength >= maxHomopolymer) {
        allowed = BASES.filter((b) => b !== prev);
      } else {
        allowed = [...BASES];
      }
      const numAllowed = allowed.length;

      const symbolIdx = decoder.read(numAllowed);
      const base = allowed[Math.min(symbolIdx, numAllowed - 1)];

      result.push(base);

      if (base === prev) {
        runLength++;
      } else {
        runLength = 1;
        prev = base;
      }
    }
  }

  let dna = result.join("");

  // Pad to target length
  if (paddedTarget > dna.length) {
    const parts = dna.split("");
    while (parts.length < paddedTarget) {
      const lastBase = parts[parts.length - 1] as Base;
      const allowed2 = BASES.filter((b) => b !== lastBase);
      parts.push(allowed2[0]);
    }
    dna = parts.join("");
  }
  if (paddedTarget < dna.length) {
    dna = dna.slice(0, paddedTarget);
  }

  return dna;
}

/**
 * Decode DNA (arithmetic coded) back to bytes.
 *
 * DECODER = ArithmeticEncoder (reads DNA symbols, outputs bits).
 * At each position, the allowed symbols depend on the homopolymer state.
 *
 * @param dna DNA string
 * @param maxHomopolymer Maximum allowed homopolymer run (must match encoder)
 * @param expectedBytes Expected number of output bytes
 * @returns Original bytes (exact, no ambiguity)
 */
export function arithmeticDnaToBytes(
  dna: string,
  maxHomopolymer: number = 3,
  expectedBytes?: number,
): Uint8Array {
  const output = new BitOutputStream();
  const encoder = new ArithmeticEncoder(output);

  let prev: Base = "A";
  let runLength = 0;

  for (let i = 0; i < dna.length; i++) {
    const base = dna[i] as Base;

    // Determine allowed symbols (same as encoder)
    let allowed: Base[];
    if (runLength >= maxHomopolymer) {
      allowed = BASES.filter((b) => b !== prev);
    } else {
      allowed = [...BASES];
    }
    const numAllowed = allowed.length;

    // Find the symbol index
    const symbolIdx = allowed.indexOf(base);
    if (symbolIdx < 0) continue; // invalid base

    // Write the symbol to the arithmetic encoder
    encoder.write(symbolIdx, numAllowed);

    // Update homopolymer state
    if (base === prev) {
      runLength++;
    } else {
      runLength = 1;
      prev = base;
    }
  }

  encoder.finish();

  return output.toBytes(expectedBytes);
}

export function satisfiesHomopolymer(dna: string, maxHomopolymer: number = 3): boolean {
  return maxHomopolymerRun(dna) <= maxHomopolymer;
}

/**
 * v57: Block-wise arithmetic DECODE with per-block CRC-8 verification.
 *
 * JS mirror of the Rust `arithmetic_dna_to_bytes_crc` (lib.rs ~L2017).
 * Used by the JS low-coverage decode path (decode.ts) when
 * `metadata.mappingMode === "arithmetic"` AND avg cluster size is below
 * the lowCoverageTrigger. Previously the JS path used direct `dnaToBytes`
 * for arithmetic mode — which produced garbage bytes (arithmetic-coded DNA
 * is NOT direct 2-bit decodable), so the 2× coverage test for arithmetic
 * mode failed at hash verification.
 *
 * Algorithm (matches Rust `arithmetic_dna_to_bytes_crc` exactly):
 *   1. Compute bytesPerBlockTotal = max(2, floor(blockSize * 1.85 / 8))
 *      bytesPerBlockData  = bytesPerBlockTotal - 1   (1 byte for CRC-8)
 *   2. For each block of `blockSize` nucleotides:
 *      a. Run arithmeticDnaToBytes(blockDna, maxHomopolymer, bytesPerBlockTotal)
 *         → recovers [CRC, data...]
 *      b. Recompute CRC-8 of data bytes
 *      c. If CRC matches → copy data bytes into result
 *      d. If CRC fails → copy data bytes anyway, but mark LAST 1 byte as erased
 *         (termination corruption zone — v56 1.85 bits/nt slack confines
 *         corruption to ≤1 byte)
 *   3. Return (result, erasures) where erasures is per-byte boolean[]
 *
 * @param dna DNA string (arithmetic-coded with per-block CRC-8 sync markers)
 * @param maxHomopolymer Max homopolymer run (must match encoder)
 * @param expectedBytes Total expected decoded bytes
 * @param blockSize Block size in nucleotides
 * @returns { data: Uint8Array, erasures: boolean[] } — erasures[i]=true means byte i is unreliable
 */
export function arithmeticDnaToBytesCrc(
  dna: string,
  maxHomopolymer: number = 3,
  expectedBytes: number,
  blockSize: number = 20,
): { data: Uint8Array; erasures: boolean[] } {
  if (blockSize === 0 || expectedBytes === 0) {
    const bytes = arithmeticDnaToBytes(dna, maxHomopolymer, expectedBytes);
    return { data: bytes, erasures: new Array(bytes.length).fill(false) };
  }

  // v56: 1.85 bits/nt capacity (matches encoder + Rust decoder)
  const bytesPerBlockTotal = Math.max(2, Math.floor((blockSize * ARITH_CAPACITY_RATE) / 8));
  const bytesPerBlockData = bytesPerBlockTotal - 1; // 1 byte for CRC-8
  if (bytesPerBlockData <= 0) {
    const bytes = arithmeticDnaToBytes(dna, maxHomopolymer, expectedBytes);
    return { data: bytes, erasures: new Array(bytes.length).fill(false) };
  }

  const numBlocks = Math.ceil(expectedBytes / bytesPerBlockData);

  const result = new Uint8Array(expectedBytes);
  const erasures = new Array(expectedBytes).fill(false);
  let byteOffset = 0;

  for (let block = 0; block < numBlocks; block++) {
    const blockStart = block * blockSize;
    const blockEnd = Math.min(blockStart + blockSize, dna.length);
    if (blockStart >= dna.length) break;
    const blockDna = dna.slice(blockStart, blockEnd);

    // v57 fix: ALWAYS decode bytesPerBlockData bytes per block (matching encoder).
    // The encoder always encodes bytesPerBlockData bytes per block, padding the
    // last block with zeros if data.length isn't a multiple of bytesPerBlockData.
    // The CRC is computed over ALL bytesPerBlockData bytes (including padding).
    //
    // Before v57: the decoder used `min(bytesPerBlockData, remaining)` for the
    // last block. This caused CRC mismatch because:
    //   - Encoder CRC: computed over bytesPerBlockData bytes (79 for bs=330)
    //   - Decoder CRC: computed over `remaining` bytes (e.g., 77 for last block)
    //   → CRC always failed for the last block → false erasures → LDPC confused
    //
    // After v57: the decoder always decodes bytesPerBlockData bytes and computes
    // CRC over all of them. Only the real data bytes (up to `remaining`) are
    // copied to the result; the rest (padding) are discarded.
    const dataBytesThisBlock = bytesPerBlockData;
    const remaining = expectedBytes - byteOffset;
    const realBytesThisBlock = Math.min(bytesPerBlockData, remaining);

    // Decode block: returns bytesPerBlockTotal bytes [CRC, data...]
    const blockBytes = arithmeticDnaToBytes(blockDna, maxHomopolymer, bytesPerBlockTotal);

    if (blockBytes.length >= bytesPerBlockTotal) {
      const crcStored = blockBytes[0];
      // Compute CRC-8 over ALL dataBytesThisBlock bytes (matches encoder)
      let crc = 0;
      for (let i = 0; i < dataBytesThisBlock; i++) {
        crc ^= blockBytes[1 + i];
        for (let bit = 0; bit < 8; bit++) {
          if (crc & 0x80) crc = ((crc << 1) ^ 0x07) & 0xff;
          else crc = (crc << 1) & 0xff;
        }
      }
      const crcComputed = crc & 0xff;

      // Copy only the REAL data bytes (not padding) into result
      for (let i = 0; i < realBytesThisBlock; i++) {
        result[byteOffset + i] = blockBytes[1 + i];
      }

      if (crcStored !== crcComputed) {
        // CRC FAIL — mark last 1 byte as erased (termination corruption zone)
        const eraseCount = Math.min(1, realBytesThisBlock);
        for (let i = 0; i < eraseCount; i++) {
          erasures[byteOffset + realBytesThisBlock - 1 - i] = true;
        }
      }
    } else {
      // Short block — mark as erased
      for (let i = 0; i < realBytesThisBlock; i++) {
        erasures[byteOffset + i] = true;
      }
    }

    byteOffset += realBytesThisBlock;
    if (byteOffset >= expectedBytes) break;
  }

  return { data: result, erasures };
}
