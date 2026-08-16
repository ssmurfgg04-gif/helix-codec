/** DNA-MT Archive — Molecular Tape ligation recipe archive mode */

import { blake3 } from '@noble/hashes/blake3.js';

// ---------------------------------------------------------------------------
// Interfaces
// ---------------------------------------------------------------------------

/** A single block in the DNA-MT archive */
export interface DnaMTBlock {
  /** Which MTs to ligate, in order (indices into the MT library) */
  mtIndices: number[];
  /** Linker sequence and orientation for ligation */
  ligationRecipe: string;
  /** Payload bytes encoded in this block */
  payloadBytes: number;
}

/** DNA-MT Archive format */
export interface DnaMTArchive {
  /** Format version */
  version: 1;
  /** BLAKE3 hash of the MT library for content addressing */
  mtLibraryHash: Uint8Array;
  /** Archive blocks */
  blocks: DnaMTBlock[];
  /** Total original payload size in bytes */
  totalPayloadSize: number;
  /** Number of MT library entries */
  mtLibrarySize: number;
}

/** MT Library entry */
export interface MTLibraryEntry {
  /** Index in the library */
  index: number;
  /** DNA sequence of this MT */
  sequence: string;
  /** Length in nucleotides */
  length: number;
  /** GC content */
  gc: number;
  /** Max homopolymer run */
  maxHomopolymer: number;
}

// ---------------------------------------------------------------------------
// Default MT Library generation
// ---------------------------------------------------------------------------

/**
 * Generate 256 well-designed MT entries (8-bit index space).
 * Each entry is a 30-nt DNA oligo that satisfies GC 40-60% and max homopolymer 3.
 * Uses a deterministic PRNG (xorshift32) for reproducibility.
 */
function generateDefaultMTLibrary(): MTLibraryEntry[] {
  const entries: MTLibraryEntry[] = [];
  const bases = ['A', 'C', 'G', 'T'];
  let state = 0x12345678;

  for (let i = 0; i < 256; i++) {
    let seq = '';
    let prev = '';
    let gcCount = 0;
    let runLength = 1;

    // Generate 30-nt sequence with GC 40-60% and no homopolymer > 3
    while (seq.length < 30) {
      // xorshift32 PRNG step
      state ^= state << 13;
      state ^= state >>> 17;
      state ^= state << 5;
      state = state >>> 0;

      // Bias toward GC balance
      const needGC = gcCount < seq.length * 0.4;
      const needAT = gcCount > seq.length * 0.6;

      let base: string;
      if (needGC) {
        base = bases[1 + (state % 2)]; // C or G
      } else if (needAT) {
        base = state % 2 === 0 ? 'A' : 'T';
      } else {
        base = bases[state % 4];
      }

      // Enforce max homopolymer run of 3
      if (base === prev && runLength >= 3) {
        // Would exceed homopolymer limit, pick a different base
        continue;
      }

      seq += base;
      if (base === prev) {
        runLength++;
      } else {
        runLength = 1;
      }
      prev = base;
      if (base === 'G' || base === 'C') gcCount++;
    }

    // Compute actual max homopolymer run
    let maxH = 1;
    let curH = 1;
    for (let j = 1; j < seq.length; j++) {
      if (seq[j] === seq[j - 1]) {
        curH++;
        if (curH > maxH) maxH = curH;
      } else {
        curH = 1;
      }
    }

    entries.push({
      index: i,
      sequence: seq,
      length: seq.length,
      gc: gcCount / seq.length,
      maxHomopolymer: maxH,
    });
  }

  return entries;
}

/** Default MT Library — a set of well-designed oligos for DNA storage */
export const DEFAULT_MT_LIBRARY: MTLibraryEntry[] = generateDefaultMTLibrary();

// ---------------------------------------------------------------------------
// BLAKE3 hash helper
// ---------------------------------------------------------------------------

/**
 * Compute BLAKE3 hash of an MT library for content-addressing.
 *
 * The hash is computed over the canonical serialization of the library:
 * for each entry, concatenate the index (uint32 LE), the sequence as UTF-8,
 * and a null terminator.
 *
 * @param library MT library entries
 * @returns 32-byte BLAKE3 hash
 */
export function computeMTLibraryHash(library: MTLibraryEntry[]): Uint8Array {
  // Build canonical byte representation of the library
  let totalLen = 0;
  const encodedSeqs: Uint8Array[] = [];

  for (const entry of library) {
    // 4 bytes for index + sequence UTF-8 + 1 null byte
    const seqBytes = new TextEncoder().encode(entry.sequence);
    encodedSeqs.push(seqBytes);
    totalLen += 4 + seqBytes.length + 1;
  }

  const buf = new Uint8Array(totalLen);
  const view = new DataView(buf.buffer);
  let offset = 0;

  for (let i = 0; i < library.length; i++) {
    view.setUint32(offset, library[i].index, true);
    offset += 4;
    buf.set(encodedSeqs[i], offset);
    offset += encodedSeqs[i].length;
    buf[offset++] = 0; // null terminator
  }

  return blake3(buf);
}

// ---------------------------------------------------------------------------
// Build MT Library
// ---------------------------------------------------------------------------

/**
 * Build an MT library from a list of DNA oligos.
 *
 * Each oligo is validated for DNA characters and assigned an index.
 * GC content and max homopolymer are computed.
 *
 * @param oligos Array of DNA sequences
 * @returns MT library entries
 * @throws Error if any sequence contains invalid DNA characters
 */
export function buildMTLibrary(oligos: string[]): MTLibraryEntry[] {
  const validBases = new Set(['A', 'C', 'G', 'T']);
  const entries: MTLibraryEntry[] = [];

  for (let i = 0; i < oligos.length; i++) {
    const seq = oligos[i].toUpperCase();

    // Validate
    for (let j = 0; j < seq.length; j++) {
      if (!validBases.has(seq[j])) {
        throw new Error(
          `Invalid DNA character '${seq[j]}' at position ${j} in oligo ${i}: "${oligos[i]}"`,
        );
      }
    }

    // Compute GC content
    let gcCount = 0;
    for (let j = 0; j < seq.length; j++) {
      if (seq[j] === 'G' || seq[j] === 'C') gcCount++;
    }

    // Compute max homopolymer run
    let maxH = 1;
    let curH = 1;
    for (let j = 1; j < seq.length; j++) {
      if (seq[j] === seq[j - 1]) {
        curH++;
        if (curH > maxH) maxH = curH;
      } else {
        curH = 1;
      }
    }

    entries.push({
      index: i,
      sequence: seq,
      length: seq.length,
      gc: gcCount / seq.length,
      maxHomopolymer: maxH,
    });
  }

  return entries;
}

// ---------------------------------------------------------------------------
// Encode to MT Archive
// ---------------------------------------------------------------------------

/**
 * Encode a list of DNA oligos into a DNA-MT archive format.
 *
 * The encoding strategy:
 * - Each oligo is broken into segments that match entries in the MT library
 * - For each segment, we find the best-matching MT library entry
 * - The ligation recipe encodes linker sequences between MT segments
 * - Payload bytes are derived from the original oligo data
 *
 * @param oligos Array of DNA sequences to encode
 * @param payload Original payload bytes that were encoded into these oligos
 * @param library Optional MT library (defaults to DEFAULT_MT_LIBRARY)
 * @returns DNA-MT archive
 */
export function encodeToMTArchive(
  oligos: string[],
  payload: Uint8Array,
  library: MTLibraryEntry[] = DEFAULT_MT_LIBRARY,
): DnaMTArchive {
  const blocks: DnaMTBlock[] = [];
  const libraryHash = computeMTLibraryHash(library);

  // Build a lookup map from sequence to index for fast matching
  const seqToIndex = new Map<string, number>();
  for (const entry of library) {
    seqToIndex.set(entry.sequence, entry.index);
  }

  let totalPayloadBytes = 0;

  for (let oligoIdx = 0; oligoIdx < oligos.length; oligoIdx++) {
    const oligo = oligos[oligoIdx].toUpperCase();
    const mtIndices: number[] = [];
    const linkerParts: string[] = [];

    let pos = 0;
    let matchedLen = 0;

    while (pos < oligo.length) {
      // Try to find the longest MT entry that matches starting at pos
      let bestMatch = -1;
      let bestLen = 0;

      for (const entry of library) {
        if (oligo.startsWith(entry.sequence, pos) && entry.length > bestLen) {
          bestMatch = entry.index;
          bestLen = entry.length;
        }
      }

      if (bestMatch >= 0) {
        mtIndices.push(bestMatch);
        matchedLen += bestLen;

        // Add ligation linker between MT segments
        if (mtIndices.length > 1) {
          // Linker format: "L<from_idx>:<to_idx>" meaning ligate entry from_idx
          // to entry to_idx with standard BsaI linker
          linkerParts.push(
            `L${mtIndices[mtIndices.length - 2]}:${bestMatch}`,
          );
        }

        pos += bestLen;
      } else {
        // No exact MT match — find the closest MT by Hamming distance
        // for the remaining sequence, using a sliding window approach
        let closest = 0;
        let closestDist = Infinity;
        let closestLen = 0;

        for (const entry of library) {
          if (entry.length > oligo.length - pos) continue;
          let dist = 0;
          for (let j = 0; j < entry.length; j++) {
            if (oligo[pos + j] !== entry.sequence[j]) dist++;
          }
          if (dist < closestDist || (dist === closestDist && entry.length > closestLen)) {
            closestDist = dist;
            closest = entry.index;
            closestLen = entry.length;
          }
        }

        // Store the residue (unmatched bases) as part of the ligation recipe
        const residue = oligo.slice(pos, pos + closestLen > oligo.length - pos
          ? oligo.length - pos
          : closestLen);

        mtIndices.push(closest);
        matchedLen += closestLen;

        if (mtIndices.length > 1) {
          linkerParts.push(
            `L${mtIndices[mtIndices.length - 2]}:${closest}:${residue}`,
          );
        } else {
          linkerParts.push(`R0:${closest}:${residue}`);
        }

        pos += closestLen || 1; // advance at least 1 base
      }
    }

    // Calculate payload bytes for this oligo (distribute evenly)
    const oligoPayloadBytes = Math.ceil(payload.length / oligos.length);
    totalPayloadBytes += oligoPayloadBytes;

    blocks.push({
      mtIndices: mtIndices,
      ligationRecipe: linkerParts.join('|') || 'identity',
      payloadBytes: oligoPayloadBytes,
    });
  }

  // Adjust total to match actual payload
  totalPayloadBytes = payload.length;

  return {
    version: 1,
    mtLibraryHash: libraryHash,
    blocks,
    totalPayloadSize: payload.length,
    mtLibrarySize: library.length,
  };
}

// ---------------------------------------------------------------------------
// Decode from MT Archive
// ---------------------------------------------------------------------------

/**
 * Decode a DNA-MT archive back to a list of DNA oligos.
 *
 * Reconstructs each oligo by looking up MT indices in the library
 * and applying the ligation recipe to insert linkers.
 *
 * @param archive DNA-MT archive to decode
 * @param library MT library entries (must match the archive's library hash)
 * @returns Array of reconstructed DNA sequences
 * @throws Error if the library hash doesn't match the archive
 */
export function decodeFromMTArchive(
  archive: DnaMTArchive,
  library: MTLibraryEntry[],
): string[] {
  // Verify library hash matches
  const libraryHash = computeMTLibraryHash(library);
  if (libraryHash.length !== archive.mtLibraryHash.length) {
    throw new Error('MT library hash length mismatch');
  }
  for (let i = 0; i < libraryHash.length; i++) {
    if (libraryHash[i] !== archive.mtLibraryHash[i]) {
      throw new Error(
        'MT library hash mismatch: the provided library does not match the archive',
      );
    }
  }

  const oligos: string[] = [];

  // Build index lookup
  const indexToEntry = new Map<number, MTLibraryEntry>();
  for (const entry of library) {
    indexToEntry.set(entry.index, entry);
  }

  for (const block of archive.blocks) {
    const parts: string[] = [];

    // Reconstruct from MT indices
    for (const idx of block.mtIndices) {
      const entry = indexToEntry.get(idx);
      if (!entry) {
        throw new Error(`MT index ${idx} not found in library`);
      }
      parts.push(entry.sequence);
    }

    // Apply ligation recipe to handle linkers and residues
    if (block.ligationRecipe !== 'identity' && block.ligationRecipe.length > 0) {
      const recipes = block.ligationRecipe.split('|');
      for (const recipe of recipes) {
        // Recipes with residue: "L<from>:<to>:<residue>" or "R0:<to>:<residue>"
        // The residue contains original bases that weren't exactly matched
        const segments = recipe.split(':');
        if (segments.length >= 3 && segments[0].startsWith('R')) {
          // Initial residue: replace first MT match with residue
          // R0:<mtIdx>:<residue>
          const residue = segments.slice(2).join(':');
          if (parts.length > 0) {
            parts[0] = residue;
          }
        } else if (segments.length >= 3 && segments[0].startsWith('L')) {
          // Linker with residue: replace the last matched MT with residue
          // L<from>:<to>:<residue>
          const residue = segments.slice(2).join(':');
          if (parts.length > 1) {
            parts[parts.length - 1] = residue;
          }
        }
        // Simple linker "L<from>:<to>" — just concat, already handled by parts
      }
    }

    oligos.push(parts.join(''));
  }

  return oligos;
}

// ---------------------------------------------------------------------------
// Binary serialization
// ---------------------------------------------------------------------------

/** DNA-MT archive magic bytes: ".dmt" = 0x2E 0x64 0x6D 0x74 */
const DMT_MAGIC = new Uint8Array([0x2E, 0x64, 0x6D, 0x74]);

/** DNA-MT archive format version */
const DMT_VERSION = 1;

/**
 * Serialize a DNA-MT archive to binary format.
 *
 * Layout:
 *   ┌──────────────────────────────────────────────────────────────────┐
 *   │ Header                                                          │
 *   │   magic          4 bytes  ".dmt"                                │
 *   │   version        2 bytes  uint16 (1)                            │
 *   │   mtLibrarySize  4 bytes  uint32                                │
 *   │   mtLibraryHash 32 bytes  BLAKE3                                │
 *   │   totalPayload   8 bytes  uint64                                │
 *   │   numBlocks      4 bytes  uint32                                │
 *   ├──────────────────────────────────────────────────────────────────┤
 *   │ For each block:                                                 │
 *   │   numMTIndices   2 bytes  uint16                                │
 *   │   mtIndices      numMTIndices × 1 byte  (uint8 indices)        │
 *   │   recipeLen      2 bytes  uint16                                │
 *   │   ligationRecipe recipeLen bytes  UTF-8                         │
 *   │   payloadBytes   4 bytes  uint32                                │
 *   └──────────────────────────────────────────────────────────────────┘
 *
 * @param archive DNA-MT archive to serialize
 * @returns Binary representation
 */
export function serializeMTArchive(archive: DnaMTArchive): Uint8Array {
  // Calculate total size
  // Header: 4 + 2 + 4 + 32 + 8 + 4 = 54 bytes
  const HEADER_SIZE = 54;
  let bodySize = 0;

  const encodedRecipes: Uint8Array[] = [];
  for (const block of archive.blocks) {
    const recipeBytes = new TextEncoder().encode(block.ligationRecipe);
    encodedRecipes.push(recipeBytes);
    // 2 (numMTIndices) + numMTIndices (1 byte each) + 2 (recipeLen) + recipeBytes + 4 (payloadBytes)
    bodySize += 2 + block.mtIndices.length + 2 + recipeBytes.length + 4;
  }

  const totalSize = HEADER_SIZE + bodySize;
  const out = new Uint8Array(totalSize);
  const view = new DataView(out.buffer);
  let pos = 0;

  // Header
  out.set(DMT_MAGIC, pos); pos += 4;
  view.setUint16(pos, DMT_VERSION, true); pos += 2;
  view.setUint32(pos, archive.mtLibrarySize, true); pos += 4;
  out.set(archive.mtLibraryHash, pos); pos += 32;
  // totalPayloadSize as uint64
  view.setUint32(pos, archive.totalPayloadSize & 0xFFFFFFFF, true); pos += 4;
  view.setUint32(pos, Math.floor(archive.totalPayloadSize / 0x100000000), true); pos += 4;
  view.setUint32(pos, archive.blocks.length, true); pos += 4;

  // Body
  for (let i = 0; i < archive.blocks.length; i++) {
    const block = archive.blocks[i];
    const recipeBytes = encodedRecipes[i];

    // numMTIndices (uint16)
    view.setUint16(pos, block.mtIndices.length, true); pos += 2;
    // mtIndices (1 byte each — supports libraries up to 256 entries)
    for (const idx of block.mtIndices) {
      out[pos++] = idx & 0xFF;
    }
    // recipeLen (uint16)
    view.setUint16(pos, recipeBytes.length, true); pos += 2;
    // ligationRecipe (UTF-8)
    out.set(recipeBytes, pos); pos += recipeBytes.length;
    // payloadBytes (uint32)
    view.setUint32(pos, block.payloadBytes, true); pos += 4;
  }

  return out;
}

/**
 * Deserialize a DNA-MT archive from binary data.
 *
 * @param data Binary representation (from serializeMTArchive)
 * @returns Parsed DNA-MT archive
 * @throws Error if the magic bytes or version are invalid
 */
export function deserializeMTArchive(data: Uint8Array): DnaMTArchive {
  const HEADER_SIZE = 54;
  if (data.length < HEADER_SIZE) {
    throw new Error(
      `Data too short for DNA-MT header: ${data.length} < ${HEADER_SIZE}`,
    );
  }

  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
  let pos = 0;

  // Magic
  if (
    data[0] !== DMT_MAGIC[0] || data[1] !== DMT_MAGIC[1] ||
    data[2] !== DMT_MAGIC[2] || data[3] !== DMT_MAGIC[3]
  ) {
    throw new Error(
      `Invalid DNA-MT magic: expected [${DMT_MAGIC}], got [${data[0]}, ${data[1]}, ${data[2]}, ${data[3]}]`,
    );
  }
  pos += 4;

  // Version
  const version = view.getUint16(pos, true); pos += 2;
  if (version !== DMT_VERSION) {
    throw new Error(`Unsupported DNA-MT version: ${version} (expected ${DMT_VERSION})`);
  }

  // mtLibrarySize
  const mtLibrarySize = view.getUint32(pos, true); pos += 4;

  // mtLibraryHash (32 bytes)
  const mtLibraryHash = data.slice(pos, pos + 32); pos += 32;

  // totalPayloadSize (uint64)
  const payloadLow = view.getUint32(pos, true); pos += 4;
  const payloadHigh = view.getUint32(pos, true); pos += 4;
  const totalPayloadSize = payloadLow + payloadHigh * 0x100000000;

  // numBlocks
  const numBlocks = view.getUint32(pos, true); pos += 4;

  // Parse blocks
  const blocks: DnaMTBlock[] = [];
  for (let i = 0; i < numBlocks; i++) {
    // numMTIndices (uint16)
    const numMTIndices = view.getUint16(pos, true); pos += 2;

    // mtIndices
    const mtIndices: number[] = [];
    for (let j = 0; j < numMTIndices; j++) {
      mtIndices.push(data[pos++]);
    }

    // recipeLen (uint16)
    const recipeLen = view.getUint16(pos, true); pos += 2;

    // ligationRecipe (UTF-8)
    const recipeBytes = data.slice(pos, pos + recipeLen); pos += recipeLen;
    const ligationRecipe = new TextDecoder().decode(recipeBytes);

    // payloadBytes (uint32)
    const payloadBytes = view.getUint32(pos, true); pos += 4;

    blocks.push({ mtIndices, ligationRecipe, payloadBytes });
  }

  return {
    version: 1,
    mtLibraryHash,
    blocks,
    totalPayloadSize,
    mtLibrarySize,
  };
}
