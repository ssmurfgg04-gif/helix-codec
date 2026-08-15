/**
 * O(1) seek index for .hlx archives.
 *
 * The block index maps logical byte offsets in the uncompressed payload to
 * physical block positions in the compressed archive. This enables random
 * access: given a logical byte offset, binary search finds the containing
 * block in O(log N), then a single seek reads that block.
 *
 * Binary format (compact, little-endian):
 *   - numEntries: uint32
 *   - For each entry:
 *       - offset:        uint64 (byte offset of block start in uncompressed stream)
 *       - compressedPos: uint64 (byte offset of block start in compressed archive)
 *       - uncompressedSize: uint32
 *       - crc32:         uint32
 *     = 24 bytes per entry
 *
 * The index is stored as the footer of the .hlx archive, enabling single-pass
 * write (append index at end) and single-pass read (seek to end, read index,
 * then random-access any block).
 *
 * Reference:
 *   - BGZF (Li et al. 2009): similar indexing scheme for BAM files.
 *   - Tabix (Li 2011): interval-based index for genomic coordinates.
 */

/** A single entry in the block index. */
export interface IndexEntry {
  /** Byte offset of this block's start in the uncompressed stream. */
  offset: number;
  /** Uncompressed size of this block in bytes. */
  uncompressedSize: number;
  /** CRC-32 checksum of the uncompressed block data. */
  crc32: number;
}

/** Block index for O(log N) random access into a compressed archive. */
export interface BlockIndex {
  /** Sorted array of index entries (by offset, ascending). */
  entries: IndexEntry[];
}

// ---------------------------------------------------------------------------
// Build
// ---------------------------------------------------------------------------

/**
 * Build a block index from an array of block metadata.
 *
 * The input blocks must be sorted by offset. The offset of each block
 * is the cumulative position in the uncompressed stream.
 *
 * @param blocks Array of block metadata (offset, uncompressedSize, crc32)
 * @returns BlockIndex structure
 *
 * @example
 *   buildBlockIndex([
 *     { offset: 0,    uncompressedSize: 65536, crc32: 0xDEADBEEF },
 *     { offset: 65536, uncompressedSize: 65536, crc32: 0xCAFEBABE },
 *   ])
 */
export function buildBlockIndex(
  blocks: { offset: number; uncompressedSize: number; crc32: number }[],
): BlockIndex {
  const entries: IndexEntry[] = blocks.map((b) => ({
    offset: b.offset,
    uncompressedSize: b.uncompressedSize,
    crc32: b.crc32,
  }));

  // Verify sorted order
  for (let i = 1; i < entries.length; i++) {
    if (entries[i].offset <= entries[i - 1].offset) {
      throw new Error(
        `Block index entries must be sorted by offset, but entry ${i} ` +
        `(offset=${entries[i].offset}) <= entry ${i - 1} (offset=${entries[i - 1].offset})`,
      );
    }
  }

  return { entries };
}

// ---------------------------------------------------------------------------
// Find (binary search)
// ---------------------------------------------------------------------------

/**
 * Find the block containing a given byte offset in the uncompressed stream.
 *
 * Uses binary search for O(log N) lookup. Returns the block index and the
 * offset within that block where the requested byte starts.
 *
 * @param index Block index (entries sorted by offset)
 * @param byteOffset Logical byte offset in the uncompressed stream
 * @returns Object with blockIndex (0-based) and offsetInBlock (bytes from
 *          start of that block). If byteOffset is past the end of the indexed
 *          data, returns { blockIndex: -1, offsetInBlock: 0 }.
 *
 * @example
 *   const idx = buildBlockIndex([
 *     { offset: 0, uncompressedSize: 100, crc32: 0 },
 *     { offset: 100, uncompressedSize: 200, crc32: 0 },
 *   ]);
 *   findBlock(idx, 50)   // { blockIndex: 0, offsetInBlock: 50 }
 *   findBlock(idx, 150)  // { blockIndex: 1, offsetInBlock: 50 }
 *   findBlock(idx, 400)  // { blockIndex: -1, offsetInBlock: 0 } (past end)
 */
export function findBlock(
  index: BlockIndex,
  byteOffset: number,
): { blockIndex: number; offsetInBlock: number } {
  const entries = index.entries;
  if (entries.length === 0 || byteOffset < 0) {
    return { blockIndex: -1, offsetInBlock: 0 };
  }

  // Check if offset is past the last block
  const last = entries[entries.length - 1];
  if (byteOffset >= last.offset + last.uncompressedSize) {
    return { blockIndex: -1, offsetInBlock: 0 };
  }

  // Binary search: find the largest entry where entry.offset <= byteOffset
  let lo = 0;
  let hi = entries.length - 1;
  while (lo < hi) {
    // Round up mid to avoid infinite loop when lo = hi - 1
    const mid = lo + ((hi - lo + 1) >> 1);
    if (entries[mid].offset <= byteOffset) {
      lo = mid;
    } else {
      hi = mid - 1;
    }
  }

  return {
    blockIndex: lo,
    offsetInBlock: byteOffset - entries[lo].offset,
  };
}

// ---------------------------------------------------------------------------
// Serialize / Deserialize
// ---------------------------------------------------------------------------

/**
 * Serialize a block index to a compact binary format.
 *
 * Format (little-endian):
 *   numEntries (uint32)
 *   For each entry (24 bytes):
 *     offset          (uint64 — stored as two uint32: low, high)
 *     uncompressedSize (uint32)
 *     crc32           (uint32)
 *
 * Total size: 4 + numEntries * 24 bytes.
 *
 * @param index BlockIndex to serialize
 * @returns Compact binary representation
 */
export function serializeIndex(index: BlockIndex): Uint8Array {
  const n = index.entries.length;
  const buf = new Uint8Array(4 + n * 24);
  const view = new DataView(buf.buffer);

  // numEntries
  view.setUint32(0, n, true);

  for (let i = 0; i < n; i++) {
    const entry = index.entries[i];
    const base = 4 + i * 24;

    // offset as uint64 (low 32 bits, high 32 bits)
    view.setUint32(base, entry.offset & 0xFFFFFFFF, true);
    view.setUint32(base + 4, Math.floor(entry.offset / 0x100000000), true);

    // uncompressedSize
    view.setUint32(base + 8, entry.uncompressedSize, true);

    // crc32
    view.setUint32(base + 12, entry.crc32, true);
  }

  return buf;
}

/**
 * Deserialize a block index from its compact binary format.
 *
 * @param data Binary data produced by serializeIndex
 * @returns BlockIndex structure
 *
 * @throws Error if the data is too short or malformed
 */
export function deserializeIndex(data: Uint8Array): BlockIndex {
  if (data.length < 4) {
    throw new Error('Index data too short: need at least 4 bytes for numEntries');
  }

  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
  const n = view.getUint32(0, true);

  if (data.length < 4 + n * 24) {
    throw new Error(
      `Index data too short: expected ${4 + n * 24} bytes for ${n} entries, ` +
      `got ${data.length}`,
    );
  }

  const entries: IndexEntry[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const base = 4 + i * 24;

    // Reconstruct uint64 offset from two uint32 halves
    const offsetLow = view.getUint32(base, true);
    const offsetHigh = view.getUint32(base + 4, true);
    const offset = offsetLow + offsetHigh * 0x100000000;

    const uncompressedSize = view.getUint32(base + 8, true);
    const crc32 = view.getUint32(base + 12, true);

    entries[i] = { offset, uncompressedSize, crc32 };
  }

  return { entries };
}
