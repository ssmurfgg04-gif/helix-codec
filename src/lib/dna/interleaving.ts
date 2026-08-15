/**
 * Interleaving for LDPC codewords (Kim 2024 approach).
 *
 * Problem: burst errors in DNA sequencing (e.g., from homopolymer runs or
 * synthesis artifacts) corrupt consecutive bytes in a single oligo. A single
 * burst of 4+ byte errors exceeds LDPC's correction capacity (~4 single-bit
 * errors per codeword).
 *
 * Solution: interleave the LDPC codewords across multiple oligos. After
 * interleaving, a burst error in one oligo affects 1 byte in each of N
 * different LDPC codewords — each codeword sees only 1 error, easily corrected.
 *
 * Architecture:
 *   Encode:
 *     1. Group oligos into blocks of `interleaveDepth` (e.g., 4)
 *     2. For each block, collect the LDPC codewords from all oligos
 *     3. Interleave: byte i of codeword j → position (i * depth + j) in the
 *        interleaved stream
 *     4. Map the interleaved stream to DNA (each oligo gets interleaved bytes)
 *
 *   Decode:
 *     1. For each block, collect DNA from all oligos → de-interleave
 *     2. Each LDPC codeword is reconstructed from bytes across N oligos
 *     3. LDPC decode each codeword (each has at most 1 burst error → 1 byte)
 *
 * This is transparent to the existing LDPC/CRC/RS pipeline — it just changes
 * how bytes are assigned to oligos. The oligo's address still identifies which
 * block and position it belongs to.
 *
 * Reference:
 *   - Kim et al. (2024). "Design of DNA Storage Coding Scheme With LDPC."
 *   - Interleaving is standard in digital communications (e.g., DVB-T, WiFi).
 */

/**
 * Interleave bytes from multiple codewords.
 *
 * Given `depth` codewords of `len` bytes each, produce a single interleaved
 * buffer of `depth * len` bytes where:
 *   interleaved[i * depth + j] = codeword[j][i]
 *
 * A burst error of length B in the interleaved buffer affects at most
 * ceil(B / depth) bytes in any single codeword.
 *
 * @param codewords Array of `depth` Uint8Arrays, each of length `len`
 * @returns Interleaved buffer of length `depth * len`
 */
export function interleaveCodewords(codewords: Uint8Array[]): Uint8Array {
  if (codewords.length === 0) return new Uint8Array(0);
  const depth = codewords.length;
  const len = codewords[0].length;
  const result = new Uint8Array(depth * len);

  for (let i = 0; i < len; i++) {
    for (let j = 0; j < depth; j++) {
      result[i * depth + j] = codewords[j][i];
    }
  }
  return result;
}

/**
 * De-interleave a buffer back into `depth` codewords.
 *
 * Given an interleaved buffer of `depth * len` bytes, produce `depth`
 * codewords of `len` bytes each where:
 *   codeword[j][i] = interleaved[i * depth + j]
 *
 * @param interleaved The interleaved buffer
 * @param depth Number of codewords to extract
 * @returns Array of `depth` Uint8Arrays, each of length `interleaved.length / depth`
 */
export function deinterleaveCodewords(
  interleaved: Uint8Array,
  depth: number,
): Uint8Array[] {
  if (depth === 0) return [];
  const len = Math.floor(interleaved.length / depth);
  const codewords: Uint8Array[] = [];

  for (let j = 0; j < depth; j++) {
    const cw = new Uint8Array(len);
    for (let i = 0; i < len; i++) {
      cw[i] = interleaved[i * depth + j];
    }
    codewords.push(cw);
  }
  return codewords;
}

/**
 * Interleave a single codeword into a group.
 *
 * Given a group of `depth` existing codewords (some may be placeholder/empty),
 * insert the new codeword at position `pos` and return the interleaved buffer.
 *
 * This is used during encode: as each oligo's LDPC codeword is generated,
 * it's added to its group. When the group is full (depth oligos), the
 * interleaved buffer is mapped to DNA.
 *
 * @param groupCodewords Array of `depth` codewords (some may be null)
 * @param newCodeword The codeword to insert
 * @param pos Position in the group (0 to depth-1)
 * @returns The interleaved buffer (only valid when all positions are filled)
 */
export function interleaveIntoGroup(
  groupCodewords: (Uint8Array | null)[],
  newCodeword: Uint8Array,
  pos: number,
): void {
  groupCodewords[pos] = newCodeword;
}

/**
 * Check if a group is complete (all positions filled).
 */
export function isGroupComplete(groupCodewords: (Uint8Array | null)[]): boolean {
  return groupCodewords.every((cw) => cw !== null);
}

/**
 * Get the interleave group index and position for a given oligo index.
 *
 * @param oligoIdx Global oligo index
 * @param depth Interleave depth (number of oligos per group)
 * @returns { groupIdx, pos } where groupIdx is the group number and pos is
 *          the position within the group (0 to depth-1)
 */
export function getInterleavePosition(
  oligoIdx: number,
  depth: number,
): { groupIdx: number; pos: number } {
  return {
    groupIdx: Math.floor(oligoIdx / depth),
    pos: oligoIdx % depth,
  };
}

/**
 * Burst error analysis: given a burst of length B and interleave depth D,
 * compute the maximum number of errors in any single codeword.
 *
 * @param burstLength Length of the burst error in bytes
 * @param depth Interleave depth
 * @returns Maximum errors per codeword
 */
export function maxErrorsPerCodeword(
  burstLength: number,
  depth: number,
): number {
  return Math.ceil(burstLength / depth);
}

/**
 * Compute the minimum interleave depth needed to correct a burst of length B
 * given that LDPC can correct E errors per codeword.
 *
 * @param burstLength Expected burst error length
 * @param ldpcCapacity LDPC error correction capacity (e.g., 4 for our code)
 * @returns Minimum interleave depth
 */
export function minDepthForBurst(
  burstLength: number,
  ldpcCapacity: number,
): number {
  return Math.ceil(burstLength / ldpcCapacity);
}
