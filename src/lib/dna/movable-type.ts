/**
 * 
 * ⚠️ EXPERIMENTAL / ROADMAP — This module is a scaffolding prototype.
 * It is NOT wired into the production encode/decode pipeline.
 * See CHANGELOG.md for v59 status and remaining work.
 *
 * Helix Movable-Type Encoder — Templated Ligation Synthesis
 *
 * Based on Wang et al., Adv. Sci. 2024 (10.1002/advs.202411354).
 * "DNA Movable Type" uses pre-fabricated double-stranded "type blocks"
 * that are enzymatically ligated in one step to build the final oligo.
 *
 * Economic advantage: 1-2 orders of magnitude write-cost reduction vs
 * de novo column synthesis because writing = assembly of pre-made blocks
 * rather than base-by-base phosphoramidite addition.
 *
 * Helix's angle: The codec's constrained coding is directly portable to
 * "constrain to the movable-type block library". The allowed sequences
 * are exactly the movable-type block concatenations.
 *
 * Block library economics:
 *   - 256 different 8-mer blocks → encodes log₂(256) = 8 bits per 8 nt = 1.0 bits/nt
 *   - 1024 different 10-mer blocks → encodes log₂(1024) ≈ 10 bits per 10 nt = 1.0 bits/nt
 *   - With assembly-decoding gains on top (composite DNA letters), can reach 1.5-2.0 bits/nt
 *
 * This module:
 *   1. Defines a block library (pre-fabricated DNA blocks)
 *   2. Encodes data bytes → block indices → ligated DNA sequence
 *   3. Decodes DNA sequence → block indices → data bytes
 *   4. Computes economic metrics (cost per MB)
 */

import { createHash } from "crypto";

export interface BlockLibrary {
  /** Block length (nt) — typically 8-12 */
  blockLength: number;
  /** Number of unique blocks in the library */
  librarySize: number;
  /** The blocks themselves (DNA strings) */
  blocks: string[];
  /** Bits encoded per block */
  bitsPerBlock: number;
  /** Theoretical density (bits/nt) */
  theoreticalDensity: number;
  /** Synthesis cost of the library (USD, one-time) */
  libraryCost: number;
  /** Ligation cost per oligo (USD) */
  ligationCostPerOligo: number;
}

/**
 * Generate a block library of `librarySize` unique blocks, each `blockLength` nt.
 * Blocks are designed to:
 *   - Have balanced GC content (40-60%)
 *   - Have no homopolymers ≥ 4
 *   - Have no restriction sites (EcoRI, BamHI, etc.)
 *   - Be ligatable (no self-complementarity at the ligation junction)
 */
export function generateBlockLibrary(
  blockLength: number = 10,
  librarySize: number = 1024,
): BlockLibrary {
  const blocks: string[] = [];
  const bases = ["A", "C", "G", "T"];
  const restrictionSites = ["GAATTC", "GGATCC", "AAGCTT", "CTGCAG"]; // EcoRI, BamHI, HindIII, PstI

  let attempts = 0;
  while (blocks.length < librarySize && attempts < librarySize * 100) {
    attempts++;
    let block = "";
    let prev = "";
    let runLength = 0;

    for (let i = 0; i < blockLength; i++) {
      // Choose a base that doesn't create a homopolymer ≥ 4
      let available = bases.filter(b => !(b === prev && runLength >= 3));
      const base = available[Math.floor(Math.random() * available.length)];

      if (base === prev) runLength++;
      else { runLength = 1; prev = base; }
      block += base;
    }

    // Check GC content
    const gc = (block.match(/[GC]/g) || []).length;
    if (gc < blockLength * 0.4 || gc > blockLength * 0.6) continue;

    // Check restriction sites
    if (restrictionSites.some(site => block.includes(site))) continue;

    // Check self-complementarity at junction (last 4 + first 4)
    const junction = block.slice(-4) + block.slice(0, 4);
    const complement = junction.replace(/A/g, "t").replace(/T/g, "a").replace(/C/g, "g").replace(/G/g, "c").toUpperCase();
    if (junction === complement.split("").reverse().join("")) continue; // self-complementary

    // Check uniqueness
    if (blocks.includes(block)) continue;

    blocks.push(block);
  }

  const bitsPerBlock = Math.log2(blocks.length);
  const theoreticalDensity = bitsPerBlock / blockLength;

  // Cost estimates (based on Wang et al. 2024)
  // Library synthesis: ~$122.20 for 470 OD of movable types
  // Ligation cost: ~$0.23 per MB
  const libraryCost = (librarySize / 470) * 122.20;
  const ligationCostPerOligo = 0.23 / 1e6 * blockLength * 10; // rough estimate

  return {
    blockLength,
    librarySize: blocks.length,
    blocks,
    bitsPerBlock,
    theoreticalDensity,
    libraryCost,
    ligationCostPerOligo,
  };
}

/**
 * Encode data bytes → block indices → ligated DNA sequence.
 *
 * The data is split into chunks of `bitsPerBlock` bits. Each chunk is
 * mapped to a block index, and the blocks are concatenated to form the
 * final DNA sequence.
 */
export function encodeWithMovableType(
  data: Uint8Array,
  library: BlockLibrary,
): { sequence: string; blockIndices: number[] } {
  const bitsPerBlock = Math.floor(library.bitsPerBlock);
  const numBlocks = Math.ceil((data.length * 8) / bitsPerBlock);

  // Convert data to bit array
  const bits: number[] = [];
  for (const byte of data) {
    for (let b = 7; b >= 0; b--) {
      bits.push((byte >> b) & 1);
    }
  }

  // Pad bits to a multiple of bitsPerBlock
  while (bits.length < numBlocks * bitsPerBlock) {
    bits.push(0);
  }

  // Map each chunk of bits to a block index
  const blockIndices: number[] = [];
  for (let i = 0; i < numBlocks; i++) {
    let index = 0;
    for (let b = 0; b < bitsPerBlock; b++) {
      index = (index << 1) | bits[i * bitsPerBlock + b];
    }
    // Ensure index is within library size
    index = index % library.librarySize;
    blockIndices.push(index);
  }

  // Concatenate blocks
  const sequence = blockIndices.map(idx => library.blocks[idx]).join("");

  return { sequence, blockIndices };
}

/**
 * Decode a DNA sequence → block indices → data bytes.
 *
 * The sequence is split into chunks of `blockLength` nt. Each chunk is
 * looked up in the library to find the block index, which gives the
 * original bits.
 */
export function decodeWithMovableType(
  sequence: string,
  library: BlockLibrary,
  expectedBytes: number,
): Uint8Array {
  const bitsPerBlock = Math.floor(library.bitsPerBlock);
  const numBlocks = Math.floor(sequence.length / library.blockLength);

  // Build reverse lookup: block sequence → index
  const blockToIndex = new Map<string, number>();
  for (let i = 0; i < library.blocks.length; i++) {
    blockToIndex.set(library.blocks[i], i);
  }

  // Decode each block
  const bits: number[] = [];
  for (let i = 0; i < numBlocks; i++) {
    const blockSeq = sequence.slice(i * library.blockLength, (i + 1) * library.blockLength);
    const index = blockToIndex.get(blockSeq);

    if (index === undefined) {
      // Block not found — could be a synthesis/sequencing error
      // Try nearest-neighbor matching (simplified — real impl would use HMM)
      let bestIdx = 0;
      let bestDist = Infinity;
      for (let j = 0; j < library.blocks.length; j++) {
        const dist = hammingDistance(blockSeq, library.blocks[j]);
        if (dist < bestDist) {
          bestDist = dist;
          bestIdx = j;
        }
      }
      // Extract bits from best match
      for (let b = bitsPerBlock - 1; b >= 0; b--) {
        bits.push((bestIdx >> b) & 1);
      }
    } else {
      // Extract bits from index
      for (let b = bitsPerBlock - 1; b >= 0; b--) {
        bits.push((index >> b) & 1);
      }
    }
  }

  // Pack bits into bytes
  const result = new Uint8Array(expectedBytes);
  for (let i = 0; i < expectedBytes; i++) {
    let byte = 0;
    for (let b = 0; b < 8; b++) {
      byte = (byte << 1) | (bits[i * 8 + b] ?? 0);
    }
    result[i] = byte;
  }

  return result;
}

/**
 * Compute Hamming distance between two equal-length strings.
 */
function hammingDistance(a: string, b: string): number {
  let dist = 0;
  for (let i = 0; i < Math.min(a.length, b.length); i++) {
    if (a[i] !== b[i]) dist++;
  }
  return dist + Math.abs(a.length - b.length);
}

/**
 * Compute economic metrics for a movable-type encoding.
 *
 * Compares cost per MB for:
 *   1. Column-based synthesis (retail): ~$100,000/MB
 *   2. Array-based synthesis (Twist): ~$10-100/MB
 *   3. Movable-type (this approach): ~$0.23/MB + library amortization
 */
export function computeEconomics(
  library: BlockLibrary,
  dataPayloadMB: number,
): {
  columnCost: number;
  arrayCost: number;
  movableTypeCost: number;
  savingsVsColumn: number;
  savingsVsArray: number;
  amortizedLibraryCost: number;
  breakEvenMB: number;
} {
  const columnCostPerMB = 100000;
  const arrayCostPerMB = 50; // midpoint of $10-100
  const ligationCostPerMB = 0.23;

  const columnCost = dataPayloadMB * columnCostPerMB;
  const arrayCost = dataPayloadMB * arrayCostPerMB;
  const ligationCost = dataPayloadMB * ligationCostPerMB;
  const amortizedLibraryCost = library.libraryCost / dataPayloadMB;
  const movableTypeCost = ligationCost + amortizedLibraryCost * Math.min(dataPayloadMB, 1);

  const savingsVsColumn = (1 - movableTypeCost / columnCost) * 100;
  const savingsVsArray = (1 - movableTypeCost / arrayCost) * 100;
  const breakEvenMB = library.libraryCost / (arrayCostPerMB - ligationCostPerMB);

  return {
    columnCost,
    arrayCost,
    movableTypeCost,
    savingsVsColumn,
    savingsVsArray,
    amortizedLibraryCost: amortizedLibraryCost * Math.min(dataPayloadMB, 1),
    breakEvenMB,
  };
}

/**
 * Serialize a block library to JSON (for storage/sharing).
 */
export function serializeLibrary(library: BlockLibrary): string {
  return JSON.stringify(library);
}

/**
 * Deserialize a block library from JSON.
 */
export function deserializeLibrary(json: string): BlockLibrary {
  return JSON.parse(json) as BlockLibrary;
}

/**
 * Compute the library ID (deterministic hash of all blocks).
 */
export function getLibraryId(library: BlockLibrary): string {
  return createHash("sha256")
    .update(library.blocks.join(""))
    .digest("hex")
    .slice(0, 16);
}
