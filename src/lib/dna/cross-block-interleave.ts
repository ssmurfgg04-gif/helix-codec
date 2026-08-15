/**
 * v51+ Cross-Block Oligo Address Interleaving
 *
 * The existing `interleaving.ts` handles intra-block codeword interleaving:
 * it spreads errors WITHIN a single block of N oligos. But a burst error
 * in a localized region of the DNA pool (e.g., a thermal hot-spot during
 * PCR amplification, or a pipetting dropout) can still wipe out a whole
 * contiguous run of oligo addresses — taking out an entire block.
 *
 * This module provides cross-block address interleaving: oligo addresses
 * are NOT assigned sequentially per block. Instead, they are permuted so
 * that any contiguous run of K physical oligos touches at most
 * ceil(K / interleaveStride) different logical blocks.
 *
 * Example (interleaveStride=4, 4 blocks of 8 oligos each):
 *   Without interleaving (sequential):
 *     Physical positions 0-7  → Block 0 (oligos 0-7)
 *     Physical positions 8-15 → Block 1 (oligos 8-15)
 *     A burst error at positions 5-12 takes out 4 oligos of block 0
 *     and 4 oligos of block 1 — but block 0 loses 50% of its oligos,
 *     exceeding the outer RS correction capacity.
 *
 *   With cross-block interleaving (stride=4):
 *     Physical pos 0 → Block 0, oligo 0
 *     Physical pos 1 → Block 1, oligo 0
 *     Physical pos 2 → Block 2, oligo 0
 *     Physical pos 3 → Block 3, oligo 0
 *     Physical pos 4 → Block 0, oligo 1
 *     ...
 *     A burst error at positions 5-12 takes out 2 oligos from each of
 *     the 4 blocks — every block loses only 25% of its oligos, well
 *     within outer RS correction capacity.
 *
 * This is the same principle as RAID-5 striping, applied to DNA pools.
 *
 * Reference: Kim et al., IEEE TNB 2024 (PubMed 38512749) — showed that
 * inter-oligo LDPC + interleaving cuts required oligo reads by 26-38%.
 */

/**
 * Configuration for cross-block interleaving.
 */
export interface CrossBlockInterleaveConfig {
  /** Number of logical blocks. */
  numBlocks: number;
  /** Oligos per block. */
  oligosPerBlock: number;
  /** Interleave stride (default = numBlocks for full round-robin). */
  stride: number;
}

/**
 * Map a logical (block, oligo-within-block) pair to a physical oligo index.
 *
 * With stride=N, physical index = blockIndex + oligoIndex * stride.
 * This places the i-th oligo of every block N positions apart, so a burst
 * of K consecutive physical positions affects ceil(K/N) blocks at most
 * once each.
 *
 * @param blockIndex Logical block index [0, numBlocks)
 * @param oligoIndex Oligo index within block [0, oligosPerBlock)
 * @param cfg Interleave configuration
 * @returns Physical oligo index [0, numBlocks * oligosPerBlock)
 */
export function logicalToPhysical(
  blockIndex: number,
  oligoIndex: number,
  cfg: CrossBlockInterleaveConfig,
): number {
  if (blockIndex < 0 || blockIndex >= cfg.numBlocks) {
    throw new Error(`blockIndex ${blockIndex} out of range [0, ${cfg.numBlocks})`);
  }
  if (oligoIndex < 0 || oligoIndex >= cfg.oligosPerBlock) {
    throw new Error(`oligoIndex ${oligoIndex} out of range [0, ${cfg.oligosPerBlock})`);
  }
  return blockIndex + oligoIndex * cfg.stride;
}

/**
 * Inverse mapping: physical oligo index → logical (block, oligo-within-block).
 *
 * @param physicalIndex Physical oligo index [0, numBlocks * oligosPerBlock)
 * @param cfg Interleave configuration
 * @returns [blockIndex, oligoIndex]
 */
export function physicalToLogical(
  physicalIndex: number,
  cfg: CrossBlockInterleaveConfig,
): [number, number] {
  const blockIndex = physicalIndex % cfg.stride;
  const oligoIndex = Math.floor(physicalIndex / cfg.stride);
  if (blockIndex >= cfg.numBlocks || oligoIndex >= cfg.oligosPerBlock) {
    throw new Error(`physicalIndex ${physicalIndex} out of range`);
  }
  return [blockIndex, oligoIndex];
}

/**
 * Generate the full physical-to-logical permutation for a given config.
 *
 * Returns an array where perm[physicalIndex] = logicalIndex
 * (logicalIndex = blockIndex * oligosPerBlock + oligoIndex).
 *
 * Use this to remap oligo addresses at encode time:
 *   oligo[physicalIndex].address = perm[physicalIndex]
 *
 * And at decode time, invert:
 *   logicalIndex = perm.indexOf(physicalIndex)
 *
 * For large pools, use physicalToLogical() directly instead of building
 * the full permutation array.
 */
export function buildPermutation(cfg: CrossBlockInterleaveConfig): Uint32Array {
  const total = cfg.numBlocks * cfg.oligosPerBlock;
  const perm = new Uint32Array(total);
  for (let phys = 0; phys < total; phys++) {
    const [blockIdx, oligoIdx] = physicalToLogical(phys, cfg);
    perm[phys] = blockIdx * cfg.oligosPerBlock + oligoIdx;
  }
  return perm;
}

/**
 * Compute the burst tolerance of a given interleave configuration.
 *
 * @param cfg Interleave configuration
 * @param outerParityRatio Outer RS parity ratio (e.g., 0.10 for 10%)
 * @returns Maximum burst length (in oligos) that can be fully corrected
 */
export function burstTolerance(
  cfg: CrossBlockInterleaveConfig,
  outerParityRatio: number,
): number {
  // Each block can tolerate up to floor(oligosPerBlock * outerParityRatio)
  // erasures via outer RS.
  const erasuresPerBlock = Math.floor(cfg.oligosPerBlock * outerParityRatio);
  // With stride=N, a burst of length B affects ceil(B/N) oligos per block.
  // For full correction, ceil(B/N) <= erasuresPerBlock → B <= N * erasuresPerBlock.
  return cfg.stride * erasuresPerBlock;
}

/**
 * Estimate the optimal stride for a given configuration.
 *
 * The optimal stride balances:
 *   - Larger stride → better burst tolerance
 *   - Smaller stride → less fragmentation, faster encode/decode
 *
 * A common choice is stride = numBlocks (full round-robin), but this can
 * be reduced if numBlocks is very large.
 *
 * @param numBlocks Number of logical blocks
 * @param oligosPerBlock Oligos per block
 * @param outerParityRatio Outer RS parity ratio
 * @param targetBurstTolerance Desired burst tolerance (in oligos)
 * @returns Recommended stride
 */
export function optimalStride(
  numBlocks: number,
  oligosPerBlock: number,
  outerParityRatio: number,
  targetBurstTolerance: number,
): number {
  // Required: stride * floor(oligosPerBlock * outerParityRatio) >= targetBurstTolerance
  const erasuresPerBlock = Math.max(1, Math.floor(oligosPerBlock * outerParityRatio));
  const minStride = Math.ceil(targetBurstTolerance / erasuresPerBlock);
  // Cap at numBlocks (full round-robin is the max useful stride)
  return Math.min(numBlocks, Math.max(1, minStride));
}
