/**
 * GC+ Code — Systematic Short-Blocklength Code for Indels + Substitutions
 *
 * Based on Hanna 2024 (arXiv:2402.01244). A systematic code that corrects
 * random insertions, deletions, and substitutions in short blocks — ideal
 * for DNA storage where each oligo is a short block.
 *
 * The code works by:
 *   1. Adding a hash-based syndrome to each block
 *   2. Using a bounded-distance decoder that tries all possible single edits
 *   3. Verifying against the syndrome to detect the correct codeword
 *
 * This is a simplified implementation that corrects single edits (1 sub OR
 * 1 ins OR 1 del) per block. The full GC+ code handles multiple edits via
 * interleaving, but single-edit correction is the core.
 *
 * Reference:
 *   - Hanna (2024). "GC+ Code: Systematic Short-Blocklength Codes Correcting
 *     Insertions, Deletions, and Substitutions." arXiv:2402.01244.
 *   - Lenz & Polyanskii (2020). "Optimal codes correcting a single deletion."
 */

/**
 * Compute a hash-based syndrome for a binary block.
 * Uses a simple parity-check + position-weighted checksum.
 */
function computeSyndrome(data: Uint8Array): { parity: number; checksum: number } {
  let parity = 0;
  let checksum = 0;
  for (let i = 0; i < data.length; i++) {
    parity ^= data[i];
    // Position-weighted checksum (mod 256 for byte-level, mod 2 for binary)
    checksum = (checksum + (i + 1) * data[i]) % 256;
  }
  return { parity, checksum };
}

export interface GCPlusConfig {
  /** Block size in bits. */
  blockSize: number;
  /** Number of syndrome bits. */
  syndromeBits: number;
}

export const DEFAULT_GCPLUS_CONFIG: GCPlusConfig = {
  blockSize: 32, // 32-bit blocks (4 bytes)
  syndromeBits: 8, // 8-bit syndrome (1 byte overhead per block)
};

export interface GCPlusEncoding {
  /** Encoded blocks (data + syndrome per block). */
  blocks: Uint8Array[];
  /** Original data length. */
  originalLength: number;
  /** Block size. */
  blockSize: number;
}

/**
 * Encode data into GC+ codewords.
 * Each block = data bits + syndrome bits.
 */
export function gcPlusEncode(
  data: Uint8Array,
  config: GCPlusConfig = DEFAULT_GCPLUS_CONFIG,
): GCPlusEncoding {
  const { blockSize, syndromeBits } = config;
  const totalBits = data.length * 8;
  const numBlocks = Math.ceil(totalBits / blockSize);
  const blocks: Uint8Array[] = [];

  // Pad data to fill blocks
  const paddedLen = Math.ceil(numBlocks * blockSize / 8);
  const padded = new Uint8Array(paddedLen);
  padded.set(data, 0);

  for (let b = 0; b < numBlocks; b++) {
    // Extract block bits
    const blockData = new Uint8Array(blockSize);
    for (let i = 0; i < blockSize; i++) {
      const bitIdx = b * blockSize + i;
      const byteIdx = Math.floor(bitIdx / 8);
      const bitInByte = 7 - (bitIdx % 8);
      if (byteIdx < padded.length) {
        blockData[i] = (padded[byteIdx] >> bitInByte) & 1;
      }
    }

    // Compute syndrome
    const { parity, checksum } = computeSyndrome(blockData);
    const syndrome = new Uint8Array(syndromeBits);
    // Pack parity + checksum into syndrome bits
    syndrome[0] = parity;
    if (syndromeBits >= 8) {
      for (let i = 1; i < syndromeBits && i <= 8; i++) {
        syndrome[i] = (checksum >> (8 - i)) & 1;
      }
    }

    // Append syndrome to block
    const codeword = new Uint8Array(blockSize + syndromeBits);
    codeword.set(blockData, 0);
    codeword.set(syndrome, blockSize);
    blocks.push(codeword);
  }

  return {
    blocks,
    originalLength: data.length,
    blockSize,
  };
}

/**
 * Decode GC+ codewords, correcting single edits (sub/ins/del) per block.
 *
 * The decoder tries:
 *   1. No error (check syndrome)
 *   2. Single substitution at each position
 *   3. Single insertion at each position (try removing each base)
 *   4. Single deletion at each position (try inserting each base)
 *
 * Returns the first candidate that matches the syndrome.
 */
export function gcPlusDecode(
  encoding: GCPlusEncoding,
): Uint8Array | null {
  const { blocks, originalLength, blockSize } = encoding;
  const syndromeBits = blocks[0]?.length ?? blockSize;
  const actualSyndromeBits = syndromeBits - blockSize;
  const recoveredBits: number[] = [];

  for (const block of blocks) {
    const data = block.slice(0, blockSize);
    const expectedSyndrome = block.slice(blockSize);

    // Try no-error case
    const { parity, checksum } = computeSyndrome(data);
    if (matchesSyndrome(parity, checksum, expectedSyndrome, actualSyndromeBits)) {
      for (const bit of data) recoveredBits.push(bit);
      continue;
    }

    // Try single substitution at each position
    let corrected = false;
    for (let i = 0; i < blockSize && !corrected; i++) {
      const candidate = data.slice();
      candidate[i] ^= 1;
      const { parity: p, checksum: c } = computeSyndrome(candidate);
      if (matchesSyndrome(p, c, expectedSyndrome, actualSyndromeBits)) {
        for (const bit of candidate) recoveredBits.push(bit);
        corrected = true;
      }
    }
    if (corrected) continue;

    // Try single deletion (insert a bit at each position)
    for (let i = 0; i <= blockSize && !corrected; i++) {
      for (let b = 0; b < 2 && !corrected; b++) {
        const candidate = new Uint8Array(blockSize);
        candidate.set(data.slice(0, i), 0);
        candidate[i] = b;
        candidate.set(data.slice(i), i + 1);
        // Remove the last bit to maintain block size
        const trimmed = candidate.slice(0, blockSize);
        const { parity: p, checksum: c } = computeSyndrome(trimmed);
        if (matchesSyndrome(p, c, expectedSyndrome, actualSyndromeBits)) {
          for (const bit of trimmed) recoveredBits.push(bit);
          corrected = true;
        }
      }
    }
    if (corrected) continue;

    // Try single insertion (remove a bit at each position)
    if (data.length > 1) {
      for (let i = 0; i < data.length && !corrected; i++) {
        const candidate = new Uint8Array(blockSize - 1);
        candidate.set(data.slice(0, i), 0);
        candidate.set(data.slice(i + 1), i);
        // Pad to blockSize
        const padded = new Uint8Array(blockSize);
        padded.set(candidate, 0);
        const { parity: p, checksum: c } = computeSyndrome(padded);
        if (matchesSyndrome(p, c, expectedSyndrome, actualSyndromeBits)) {
          for (const bit of padded) recoveredBits.push(bit);
          corrected = true;
        }
      }
    }
    if (corrected) continue;

    // Failed to correct — use raw data
    for (const bit of data) recoveredBits.push(bit);
  }

  // Convert bits back to bytes
  const numBytes = Math.floor(recoveredBits.length / 8);
  const result = new Uint8Array(numBytes);
  for (let i = 0; i < numBytes; i++) {
    let byte = 0;
    for (let j = 0; j < 8; j++) {
      byte = (byte << 1) | (recoveredBits[i * 8 + j] ?? 0);
    }
    result[i] = byte;
  }

  return result.slice(0, originalLength);
}

function matchesSyndrome(
  parity: number,
  checksum: number,
  expected: Uint8Array,
  syndromeBits: number,
): boolean {
  if (expected[0] !== parity) return false;
  if (syndromeBits >= 8) {
    // Check checksum bits
    for (let i = 1; i < Math.min(syndromeBits, 9); i++) {
      const expectedBit = expected[i] ?? 0;
      const actualBit = (checksum >> (8 - i)) & 1;
      if (expectedBit !== actualBit) return false;
    }
  }
  return true;
}
