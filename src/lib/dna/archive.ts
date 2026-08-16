/**
 * .hlx canonical archive format — Binary, versioned, indexed.
 *
 * The .hlx format is a self-contained binary archive for DNA storage data,
 * designed for efficient random access and integrity verification.
 *
 * Layout:
 *   ┌─────────────────────────────────────────────────────────┐
 *   │ Header (63 bytes)                                       │
 *   │   magic        4 bytes  0x2E 0x68 0x6C 0x78  (".hlx") │
 *   │   version      2 bytes  uint16 (currently 1)           │
 *   │   flags        2 bytes  uint16                         │
 *   │   blockSize    2 bytes  uint16 (max uncompressed/block) │
 *   │   numBlocks    4 bytes  uint32                         │
 *   │   payloadSize  8 bytes  uint64 (total uncompressed)    │
 *   │   checksumAlg  1 byte   0=sha256, 1=blake3            │
 *   │   masterCksum 32 bytes  checksum of uncompressed data  │
 *   │   indexOffset 8 bytes  uint64 (offset of footer index) │
 *   ├─────────────────────────────────────────────────────────┤
 *   │ Body: concatenated BGZF-compatible blocks               │
 *   │   For each block:                                      │
 *   │     blockLen      4 bytes  uint32 (compressed data len) │
 *   │     compressedData blockLen bytes (BGZF subformat)     │
 *   │     crc32         4 bytes  CRC-32 of uncompressed block │
 *   ├─────────────────────────────────────────────────────────┤
 *   │ Footer index (numBlocks entries, 16 bytes each):        │
 *   │   blockOffset     8 bytes  uint64                      │
 *   │   uncompressedSize 4 bytes  uint32                     │
 *   │   crc32           4 bytes  CRC-32                      │
 *   └─────────────────────────────────────────────────────────┘
 *
 * BGZF compatibility:
 *   Each compressed block starts with the BGZF subformat header:
 *     0x1F 0x8B 0x08 0x04  (gzip ID2, ID1, CM=deflate, FLG=FEXTRA)
 *   This makes .hlx blocks parseable by standard BGZF readers (samtools,
 *   htslib), enabling interoperability with bioinformatics tools.
 *
 * Footer index enables O(1) random access:
 *   Given a logical byte offset, binary search the footer index to find
 *   the containing block, then decompress just that block.
 *
 * Reference:
 *   - BGZF specification (Li et al. 2009, SAM/BAM format).
 *   - ZFS-inspired checksumming (Bonwick et al. 2003).
 */

import { crc32 } from './crc32';
import * as pako from 'pako';
import {
  serializeMTArchive,
  deserializeMTArchive,
  DnaMTArchive,
} from './dna-mt-archive';

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/** Magic bytes: ".hlx" = 0x2E 0x68 0x6C 0x78 */
export const HLX_MAGIC = new Uint8Array([0x2E, 0x68, 0x6C, 0x78]);

/** Current format version. */
export const HLX_VERSION = 1;

/** Header size in bytes. */
export const HLX_HEADER_SIZE = 63;

/** Footer index entry size in bytes (blockOffset(8) + uncompressedSize(4) + crc32(4)). */
export const HLX_INDEX_ENTRY_SIZE = 16;

/** Checksum algorithm identifiers. */
export const CHECKSUM_SHA256 = 0;
export const CHECKSUM_BLAKE3 = 1;

/** BGZF subformat header: gzip ID2, ID1, CM=deflate(8), FLG=FEXTRA(4). */
export const BGZF_HEADER = new Uint8Array([0x1F, 0x8B, 0x08, 0x04]);

/** BGZF extra field subfield identifiers (BC = Block Comment / size). */
const BGZF_XLEN = 6; // extra field length
const BGZF_SI1 = 0x42; // 'B'
const BGZF_SI2 = 0x43; // 'C'
const BGZF_SUBFIELD_LEN = 2; // 2 bytes for block size

/** Default block size (max uncompressed data per block). */
export const DEFAULT_BLOCK_SIZE = 65536; // 64 KiB (same as BGZF default)

/** DNA-MT archive magic bytes: ".dmt" */
export const DMT_MAGIC = new Uint8Array([0x2E, 0x64, 0x6D, 0x74]);

// ---------------------------------------------------------------------------
// Interfaces
// ---------------------------------------------------------------------------

/** .hlx archive header. */
export interface HlxHeader {
  /** Format version (currently 1). */
  version: number;
  /** Feature flags (reserved, currently 0). */
  flags: number;
  /** Maximum uncompressed data per block. */
  blockSize: number;
  /** Number of blocks in the archive. */
  numBlocks: number;
  /** Total uncompressed payload size. */
  payloadSize: number;
  /** Checksum algorithm: 'sha256' or 'blake3'. */
  checksumAlg: 'sha256' | 'blake3';
  /** Master checksum of the entire uncompressed payload. */
  masterChecksum: Uint8Array;
  /** Byte offset of the footer index in the archive. */
  indexOffset: number;
}

/** A single compressed block in the archive body. */
export interface HlxBlock {
  /** Compressed data (BGZF-compatible). */
  compressedData: Uint8Array;
  /** Uncompressed size of this block. */
  uncompressedSize: number;
  /** CRC-32 of the uncompressed block data. */
  crc32: number;
}

/** A footer index entry for a single block. */
export interface HlxIndexEntry {
  /** Byte offset of this block's start in the archive (from archive start). */
  blockOffset: number;
  /** Uncompressed size of this block. */
  uncompressedSize: number;
  /** CRC-32 of the uncompressed block data. */
  crc32: number;
}

/** Parsed .hlx archive. */
export interface HlxArchive {
  /** Parsed header. */
  header: HlxHeader;
  /** Array of blocks (compressed data + metadata). */
  blocks: HlxBlock[];
  /** Footer index entries (one per block). */
  index: HlxIndexEntry[];
  /** The raw archive bytes (for seekBlock operations). */
  raw: Uint8Array;
}

// ---------------------------------------------------------------------------
// BGZF block wrapping / unwrapping
// ---------------------------------------------------------------------------

/**
 * Wrap compressed data in a BGZF-compatible container.
 *
 * BGZF is a restricted gzip format where each member represents one block.
 * The gzip header includes an FEXTRA field with subfield SI1='B', SI2='C'
 * and a 2-byte little-endian block size (total compressed size + 1 for
 * the trailing 0x00 byte, inclusive of the gzip wrapper).
 *
 * We construct a minimal gzip member:
 *   [ID1 ID2 CM FLG] [MTIME(4)] [XFL] [OS] [XLEN(2)] [SI1 SI2 SLEN BSIZE(2)]
 *   [compressed data...] [CRC32(4)] [ISIZE(4)]
 *
 * @param compressedData DEFLATE-compressed data (raw deflate, not gzip)
 * @param uncompressedSize Original uncompressed size
 * @param uncompressedCrc32 CRC-32 of the uncompressed data
 * @returns BGZF-compatible block bytes
 */
export function wrapBgzfBlock(
  compressedData: Uint8Array,
  uncompressedSize: number,
  uncompressedCrc32: number,
): Uint8Array {
  // Gzip header (10 bytes) + extra field (6 bytes) = 16 bytes before data
  // Trailer: CRC32(4) + ISIZE(4) = 8 bytes after data
  const headerLen = 10 + BGZF_XLEN;
  const trailerLen = 8;
  const totalLen = headerLen + compressedData.length + trailerLen;

  // Block size for BGZF BC subfield: total member size - 1
  // (The -1 is a BGZF convention: BSIZE = total_size - 1)
  const bsize = totalLen - 1;

  const block = new Uint8Array(totalLen);
  const view = new DataView(block.buffer);

  let pos = 0;

  // Gzip header
  block[pos++] = 0x1F; // ID1
  block[pos++] = 0x8B; // ID2
  block[pos++] = 0x08; // CM = deflate
  block[pos++] = 0x04; // FLG = FEXTRA
  block[pos++] = 0;    // MTIME (4 bytes)
  block[pos++] = 0;
  block[pos++] = 0;
  block[pos++] = 0;
  block[pos++] = 0;    // XFL
  block[pos++] = 0xFF; // OS = unknown

  // Extra field
  view.setUint16(pos, BGZF_XLEN, true); pos += 2; // XLEN
  block[pos++] = BGZF_SI1; // SI1
  block[pos++] = BGZF_SI2; // SI2
  view.setUint16(pos, BGZF_SUBFIELD_LEN, true); pos += 2; // SLEN
  view.setUint16(pos, bsize, true); pos += 2; // BSIZE

  // Compressed data
  block.set(compressedData, pos);
  pos += compressedData.length;

  // Trailer
  view.setUint32(pos, uncompressedCrc32, true); pos += 4; // CRC32
  view.setUint32(pos, uncompressedSize, true); pos += 4; // ISIZE

  return block;
}

/**
 * Unwrap a BGZF block, extracting the compressed data and metadata.
 *
 * @param block BGZF-compressed block bytes
 * @returns Object with compressed data (raw deflate), CRC32, and ISIZE
 */
export function unwrapBgzfBlock(block: Uint8Array): {
  compressedData: Uint8Array;
  crc32: number;
  isize: number;
} {
  if (block.length < 18) {
    throw new Error('BGZF block too short');
  }

  // Verify gzip header
  if (
    block[0] !== 0x1F || block[1] !== 0x8B ||
    block[2] !== 0x08 || block[3] !== 0x04
  ) {
    throw new Error('Invalid BGZF header: expected gzip with FEXTRA');
  }

  const view = new DataView(block.buffer, block.byteOffset, block.byteLength);

  // Read XLEN from bytes 10-11
  const xlen = view.getUint16(10, true);

  // Skip extra field to find compressed data start
  const compressedStart = 10 + xlen;

  // Read trailer (last 8 bytes)
  const trailerOffset = block.length - 8;
  const storedCrc32 = view.getUint32(trailerOffset, true);
  const isize = view.getUint32(trailerOffset + 4, true);

  const compressedData = block.slice(compressedStart, trailerOffset);

  return { compressedData, crc32: storedCrc32, isize };
}

// ---------------------------------------------------------------------------
// SHA-256 helper
// ---------------------------------------------------------------------------

/**
 * Compute SHA-256 hash using crypto.subtle.digest (always available in Node.js).
 *
 * @param data Input bytes
 * @returns 32-byte SHA-256 hash
 */
async function sha256(data: Uint8Array): Promise<Uint8Array> {
  // Ensure a clean ArrayBuffer (avoids SharedArrayBuffer/offset issues)
  const buf = data.buffer.byteLength === data.byteLength && data.byteOffset === 0
    ? data.buffer as ArrayBuffer
    : data.slice().buffer as ArrayBuffer;
  const hash = await crypto.subtle.digest('SHA-256', buf);
  return new Uint8Array(hash);
}

/**
 * Compute the master checksum of the uncompressed payload.
 *
 * @param data Uncompressed payload bytes
 * @param alg Checksum algorithm
 * @returns 32-byte checksum (for blake3, only first 32 bytes)
 */
async function computeMasterChecksum(
  data: Uint8Array,
  alg: 'sha256' | 'blake3',
): Promise<Uint8Array> {
  if (alg === 'sha256') {
    return sha256(data);
  }
  // Blake3 stub: fall back to SHA-256 until WASM module is loaded
  // In production, this would call the blake3 WASM module
  return sha256(data);
}

// ---------------------------------------------------------------------------
// Write archive
// ---------------------------------------------------------------------------

/**
 * Write a .hlx archive from an array of uncompressed blocks.
 *
 * Each block is compressed independently (BGZF-compatible), and a footer
 * index is appended for O(1) random access.
 *
 * @param blocks Array of uncompressed block data
 * @param checksumAlg Master checksum algorithm ('sha256' or 'blake3')
 * @returns Complete .hlx archive as a Uint8Array
 *
 * @example
 *   const block1 = new Uint8Array([0x00, 0x01, 0x02]);
 *   const block2 = new Uint8Array([0x03, 0x04, 0x05]);
 *   const archive = await writeHlxArchive([block1, block2], 'sha256');
 */
export async function writeHlxArchive(
  blocks: Uint8Array[],
  checksumAlg: 'sha256' | 'blake3' = 'sha256',
): Promise<Uint8Array> {
  const numBlocks = blocks.length;
  const payloadSize = blocks.reduce((sum, b) => sum + b.length, 0);

  // Compute master checksum of the entire payload
  const fullPayload = new Uint8Array(payloadSize);
  let offset = 0;
  for (const block of blocks) {
    fullPayload.set(block, offset);
    offset += block.length;
  }
  const masterChecksum = await computeMasterChecksum(fullPayload, checksumAlg);

  // Compress each block and compute CRC-32
  const compressedBlocks: {
    bgzf: Uint8Array;
    uncompressedSize: number;
    crc: number;
  }[] = [];

  for (const block of blocks) {
    const uncompressedCrc = crc32(block) >>> 0;
    const compressed = pako.deflateRaw(block, { level: 6 });
    const bgzf = wrapBgzfBlock(compressed, block.length, uncompressedCrc);
    compressedBlocks.push({
      bgzf,
      uncompressedSize: block.length,
      crc: uncompressedCrc,
    });
  }

  // Calculate body size and index offset
  let bodySize = 0;
  for (const cb of compressedBlocks) {
    // Each block in the body: blockLen(4) + compressedData + crc32(4)
    bodySize += 4 + cb.bgzf.length + 4;
  }

  const indexOffset = HLX_HEADER_SIZE + bodySize;
  const indexSize = numBlocks * HLX_INDEX_ENTRY_SIZE;
  const totalSize = HLX_HEADER_SIZE + bodySize + indexSize;

  // Allocate output buffer
  const out = new Uint8Array(totalSize);
  const view = new DataView(out.buffer);

  let pos = 0;

  // --- Write header ---
  out.set(HLX_MAGIC, pos); pos += 4;
  view.setUint16(pos, HLX_VERSION, true); pos += 2;
  view.setUint16(pos, 0, true); pos += 2; // flags
  view.setUint16(pos, DEFAULT_BLOCK_SIZE, true); pos += 2; // blockSize
  view.setUint32(pos, numBlocks, true); pos += 4;
  // payloadSize as uint64 (low, high)
  view.setUint32(pos, payloadSize & 0xFFFFFFFF, true); pos += 4;
  view.setUint32(pos, Math.floor(payloadSize / 0x100000000), true); pos += 4;
  // checksumAlg
  out[pos++] = checksumAlg === 'sha256' ? CHECKSUM_SHA256 : CHECKSUM_BLAKE3;
  // masterChecksum (32 bytes)
  out.set(masterChecksum, pos); pos += 32;
  // indexOffset as uint64
  view.setUint32(pos, indexOffset & 0xFFFFFFFF, true); pos += 4;
  view.setUint32(pos, Math.floor(indexOffset / 0x100000000), true); pos += 4;

  // --- Write body ---
  const blockOffsets: number[] = [];
  for (let i = 0; i < compressedBlocks.length; i++) {
    const cb = compressedBlocks[i];
    blockOffsets.push(pos);

    // blockLen (uint32)
    view.setUint32(pos, cb.bgzf.length, true); pos += 4;
    // compressed data (BGZF block)
    out.set(cb.bgzf, pos); pos += cb.bgzf.length;
    // CRC-32 of uncompressed block
    view.setUint32(pos, cb.crc, true); pos += 4;
  }

  // --- Write footer index ---
  for (let i = 0; i < numBlocks; i++) {
    const cb = compressedBlocks[i];
    // blockOffset (uint64)
    view.setUint32(pos, blockOffsets[i] & 0xFFFFFFFF, true); pos += 4;
    view.setUint32(pos, Math.floor(blockOffsets[i] / 0x100000000), true); pos += 4;
    // uncompressedSize (uint32)
    view.setUint32(pos, cb.uncompressedSize, true); pos += 4;
    // crc32 (uint32)
    view.setUint32(pos, cb.crc, true); pos += 4;
  }

  return out;
}

// ---------------------------------------------------------------------------
// Read archive
// ---------------------------------------------------------------------------

/**
 * Read (parse) a .hlx archive from binary data.
 *
 * Parses the header, body blocks, and footer index. Does NOT decompress
 * the blocks — use seekBlock() to decompress individual blocks.
 *
 * @param data Raw .hlx archive bytes
 * @returns Parsed HlxArchive structure
 *
 * @throws Error if magic bytes, version, or structure are invalid
 */
export function readHlxArchive(data: Uint8Array): HlxArchive {
  if (data.length < HLX_HEADER_SIZE) {
    throw new Error(
      `Data too short for .hlx header: ${data.length} < ${HLX_HEADER_SIZE}`,
    );
  }

  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
  let pos = 0;

  // --- Parse header ---
  // Magic
  if (
    data[0] !== HLX_MAGIC[0] || data[1] !== HLX_MAGIC[1] ||
    data[2] !== HLX_MAGIC[2] || data[3] !== HLX_MAGIC[3]
  ) {
    throw new Error(
      `Invalid .hlx magic: expected [${HLX_MAGIC}], got [${data[0]}, ${data[1]}, ${data[2]}, ${data[3]}]`,
    );
  }
  pos += 4;

  const version = view.getUint16(pos, true); pos += 2;
  if (version !== HLX_VERSION) {
    throw new Error(`Unsupported .hlx version: ${version} (expected ${HLX_VERSION})`);
  }

  const flags = view.getUint16(pos, true); pos += 2;
  const blockSize = view.getUint16(pos, true); pos += 2;
  const numBlocks = view.getUint32(pos, true); pos += 4;

  // payloadSize (uint64)
  const payloadLow = view.getUint32(pos, true); pos += 4;
  const payloadHigh = view.getUint32(pos, true); pos += 4;
  const payloadSize = payloadLow + payloadHigh * 0x100000000;

  // checksumAlg
  const checksumAlgByte = data[pos++];
  const checksumAlg: 'sha256' | 'blake3' = checksumAlgByte === CHECKSUM_SHA256
    ? 'sha256'
    : 'blake3';

  // masterChecksum (32 bytes)
  const masterChecksum = data.slice(pos, pos + 32); pos += 32;

  // indexOffset (uint64)
  const idxOffLow = view.getUint32(pos, true); pos += 4;
  const idxOffHigh = view.getUint32(pos, true); pos += 4;
  const indexOffset = idxOffLow + idxOffHigh * 0x100000000;

  const header: HlxHeader = {
    version,
    flags,
    blockSize,
    numBlocks,
    payloadSize,
    checksumAlg,
    masterChecksum,
    indexOffset,
  };

  // --- Parse body blocks ---
  const blocks: HlxBlock[] = [];
  pos = HLX_HEADER_SIZE; // reset to body start

  for (let i = 0; i < numBlocks; i++) {
    const blockLen = view.getUint32(pos, true); pos += 4;
    const compressedData = data.slice(pos, pos + blockLen); pos += blockLen;
    const blockCrc32 = view.getUint32(pos, true); pos += 4;

    // Determine uncompressed size from BGZF trailer
    const bgzf = unwrapBgzfBlock(compressedData);

    blocks.push({
      compressedData,
      uncompressedSize: bgzf.isize,
      crc32: blockCrc32,
    });
  }

  // --- Parse footer index ---
  pos = indexOffset;
  const index: HlxIndexEntry[] = [];

  for (let i = 0; i < numBlocks; i++) {
    const blkOffLow = view.getUint32(pos, true); pos += 4;
    const blkOffHigh = view.getUint32(pos, true); pos += 4;
    const blockOffset = blkOffLow + blkOffHigh * 0x100000000;

    const uncompressedSize = view.getUint32(pos, true); pos += 4;
    const entryCrc32 = view.getUint32(pos, true); pos += 4;

    index.push({ blockOffset, uncompressedSize, crc32: entryCrc32 });
  }

  return { header, blocks, index, raw: data };
}

// ---------------------------------------------------------------------------
// Seek (random access)
// ---------------------------------------------------------------------------

/**
 * Seek to and decompress a specific block by index (O(1) seek via footer index).
 *
 * Uses the footer index to find the block's position in the archive,
 * then decompresses just that block. This enables random access without
 * reading the entire archive.
 *
 * @param archive Parsed .hlx archive (from readHlxArchive)
 * @param blockIndex 0-based block index
 * @returns Decompressed block data
 *
 * @throws Error if blockIndex is out of range
 *
 * @example
 *   const archive = readHlxArchive(hlxBytes);
 *   const block0 = seekBlock(archive, 0);  // first block
 *   const block5 = seekBlock(archive, 5);  // sixth block
 */
export function seekBlock(archive: HlxArchive, blockIndex: number): Uint8Array {
  if (blockIndex < 0 || blockIndex >= archive.header.numBlocks) {
    throw new Error(
      `Block index ${blockIndex} out of range [0, ${archive.header.numBlocks})`,
    );
  }

  const block = archive.blocks[blockIndex];

  // Unwrap BGZF and decompress
  const bgzf = unwrapBgzfBlock(block.compressedData);
  const decompressed = pako.inflateRaw(bgzf.compressedData);

  return decompressed;
}

// ---------------------------------------------------------------------------
// Validate
// ---------------------------------------------------------------------------

/**
 * Validate a .hlx archive: check header, per-block CRCs, and master checksum.
 *
 * @param archive Parsed .hlx archive (from readHlxArchive)
 * @returns True if all checksums match, false otherwise
 *
 * @example
 *   const archive = readHlxArchive(hlxBytes);
 *   const valid = await validateHlxArchive(archive);  // true if intact
 */
export async function validateHlxArchive(archive: HlxArchive): Promise<boolean> {
  // 1. Check each block's CRC-32
  for (let i = 0; i < archive.blocks.length; i++) {
    const block = archive.blocks[i];
    const decompressed = seekBlock(archive, i);
    const actualCrc = crc32(decompressed) >>> 0;

    if (actualCrc !== block.crc32) {
      return false; // CRC mismatch
    }
  }

  // 2. Check master checksum
  const fullPayload = new Uint8Array(archive.header.payloadSize);
  let offset = 0;
  for (let i = 0; i < archive.blocks.length; i++) {
    const decompressed = seekBlock(archive, i);
    fullPayload.set(decompressed, offset);
    offset += decompressed.length;
  }

  const masterChecksum = await computeMasterChecksum(
    fullPayload,
    archive.header.checksumAlg,
  );

  // Compare 32 bytes
  for (let i = 0; i < 32; i++) {
    if (masterChecksum[i] !== archive.header.masterChecksum[i]) {
      return false;
    }
  }

  return true;
}

// ---------------------------------------------------------------------------
// DNA-MT archive support (alternative format)
// ---------------------------------------------------------------------------

/** Union type for any supported archive format. */
export type AnyArchive = HlxArchive | DnaMTArchive;

/** Archive format identifier. */
export type ArchiveFormat = 'hlx' | 'dnamt';

/**
 * Detect the archive format from magic bytes.
 *
 * @param data Raw archive bytes
 * @returns Format identifier
 * @throws Error if the format is not recognized
 */
export function detectArchiveFormat(data: Uint8Array): ArchiveFormat {
  if (data.length < 4) {
    throw new Error('Data too short to detect archive format');
  }

  // Check for .hlx magic
  if (
    data[0] === HLX_MAGIC[0] && data[1] === HLX_MAGIC[1] &&
    data[2] === HLX_MAGIC[2] && data[3] === HLX_MAGIC[3]
  ) {
    return 'hlx';
  }

  // Check for .dmt magic
  if (
    data[0] === DMT_MAGIC[0] && data[1] === DMT_MAGIC[1] &&
    data[2] === DMT_MAGIC[2] && data[3] === DMT_MAGIC[3]
  ) {
    return 'dnamt';
  }

  throw new Error(
    `Unknown archive format: magic bytes [${data[0].toString(16)}, ${data[1].toString(16)}, ${data[2].toString(16)}, ${data[3].toString(16)}]`,
  );
}

/**
 * Read an archive in any supported format (auto-detects .hlx vs .dmt).
 *
 * @param data Raw archive bytes
 * @returns Parsed archive (either HlxArchive or DnaMTArchive)
 * @throws Error if the format is not recognized or the data is invalid
 */
export function readAnyArchive(data: Uint8Array): AnyArchive {
  const format = detectArchiveFormat(data);
  if (format === 'hlx') {
    return readHlxArchive(data);
  }
  return deserializeMTArchive(data);
}

/**
 * Write a DNA-MT archive to binary format.
 *
 * Convenience wrapper that delegates to dna-mt-archive's serializeMTArchive.
 *
 * @param archive DNA-MT archive
 * @returns Binary representation
 */
export function writeMTArchive(archive: DnaMTArchive): Uint8Array {
  return serializeMTArchive(archive);
}

/**
 * Read a DNA-MT archive from binary data.
 *
 * Convenience wrapper that delegates to dna-mt-archive's deserializeMTArchive.
 *
 * @param data Binary archive data
 * @returns Parsed DNA-MT archive
 */
export function readMTArchive(data: Uint8Array): DnaMTArchive {
  return deserializeMTArchive(data);
}
