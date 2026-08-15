/**
 * Streaming encode/decode for large files.
 *
 * For files larger than available RAM, process in chunks:
 *   - Streaming encode: read file in chunks → encode each chunk → yield oligos
 *   - Streaming decode: receive reads in batches → decode → yield recovered chunks
 *
 * The stream processes one outer RS block at a time (typically 100-1000 oligos),
 * keeping peak RAM bounded to O(block_size) rather than O(file_size).
 *
 * Peak memory:
 *   - Encode: O(chunkSize) — one chunk in memory at a time
 *   - Decode: O(chunkSize * coverage) — stores reads per chunk for clustering
 */

import { encodeFile, EncodeResult } from './codec';
import { decodeReads, DecodeResult } from './decode';
import { CodecConfig, CodecMetadata, DEFAULT_CONFIG, EncodedFile } from './types';
import { SequencingRead, SimulationResult, simulate, PRESET_CLEAN } from './simulate';

// ---------------------------------------------------------------------------
// Interfaces
// ---------------------------------------------------------------------------

/** A chunk of oligos produced by streaming encode. */
export interface StreamEncodeChunk {
  /** 0-based chunk index. */
  chunkIndex: number;
  /** Byte offset in the original file. */
  offset: number;
  /** Number of bytes in this chunk. */
  length: number;
  /** Encoded file for this chunk (metadata + oligos). */
  encoded: EncodedFile;
  /** SHA-256 hash of the chunk data (hex). */
  hash: string;
}

/** A chunk of recovered data produced by streaming decode. */
export interface StreamDecodeChunk {
  /** 0-based chunk index. */
  chunkIndex: number;
  /** Recovered data for this chunk. */
  data: Uint8Array;
  /** SHA-256 hash of the recovered data (hex). */
  hash: string;
  /** Whether the recovered hash matches the expected hash from metadata. */
  hashMatches: boolean;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/**
 * Create an async iterable that splits a large buffer into fixed-size chunks.
 *
 * @param data   The full data buffer to split.
 * @param chunkSize  Maximum bytes per chunk (last chunk may be smaller).
 * @yields Uint8Array views (subarrays) of the original buffer.
 */
export async function* createChunkIterator(
  data: Uint8Array,
  chunkSize: number,
): AsyncIterable<Uint8Array> {
  if (chunkSize <= 0) {
    throw new Error(`chunkSize must be > 0, got ${chunkSize}`);
  }
  let offset = 0;
  while (offset < data.length) {
    const end = Math.min(offset + chunkSize, data.length);
    yield data.subarray(offset, end);
    offset = end;
  }
}

/**
 * Create an async iterable from a Node.js ReadableStream or web ReadableStream.
 * Enables true streaming encode/decode from files without loading entire file into RAM.
 * Peak memory is O(chunkSize).
 *
 * @param stream  A ReadableStream (web) or Node.js readable stream.
 * @param chunkSize  Maximum bytes per buffer hint (default: 64KB).
 * @yields Uint8Array chunks read from the stream.
 */
export async function* createStreamIterator(
  stream: AsyncIterable<Uint8Array> | ReadableStream<Uint8Array>,
  chunkSize: number = 65536,
): AsyncIterable<Uint8Array> {
  if (chunkSize <= 0) {
    throw new Error(`chunkSize must be > 0, got ${chunkSize}`);
  }

  // Web ReadableStream
  if (typeof ReadableStream !== 'undefined' && stream instanceof ReadableStream) {
    const reader = (stream as ReadableStream<Uint8Array>).getReader();
    try {
      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        if (value && value.length > 0) yield value;
      }
    } finally {
      reader.releaseLock();
    }
    return;
  }

  // Async iterable (Node.js readable stream)
  for await (const chunk of stream as AsyncIterable<Uint8Array>) {
    if (chunk && chunk.length > 0) {
      yield chunk instanceof Uint8Array ? chunk : new Uint8Array(chunk as unknown as ArrayBuffer);
    }
  }
}

/**
 * Compute SHA-256 hash of a Uint8Array, returning hex string.
 * Works in both Node.js and browser environments.
 */
async function sha256Hex(data: Uint8Array): Promise<string> {
  if (typeof crypto !== 'undefined' && crypto.subtle) {
    const hash = await crypto.subtle.digest('SHA-256', data.buffer as ArrayBuffer);
    return Array.from(new Uint8Array(hash))
      .map((b) => b.toString(16).padStart(2, '0'))
      .join('');
  }
  // Node.js fallback
  const { createHash } = await import('crypto');
  return createHash('sha256').update(Buffer.from(data)).digest('hex');
}

// ---------------------------------------------------------------------------
// Streaming Encode
// ---------------------------------------------------------------------------

/**
 * Streaming encode: read data in chunks, encode each independently, yield results.
 *
 * Each chunk is encoded as a self-contained outer RS block. This keeps peak
 * memory at O(chunkSize) rather than O(fileSize).
 *
 * @param data   Async iterable of data chunks (e.g., from createChunkIterator or a file stream).
 * @param cfg    Codec configuration (applied to each chunk).
 * @param meta   File metadata (fileName, contentType).
 * @yields StreamEncodeChunk for each encoded chunk.
 */
export async function* streamEncode(
  data: AsyncIterable<Uint8Array>,
  cfg: CodecConfig = DEFAULT_CONFIG,
  meta: { fileName: string; contentType: string },
): AsyncGenerator<StreamEncodeChunk> {
  let chunkIndex = 0;
  let byteOffset = 0;

  for await (const chunk of data) {
    // Skip zero-length chunks (can occur at file boundaries).
    if (chunk.length === 0) {
      chunkIndex++;
      continue;
    }

    const result: EncodeResult = await encodeFile(chunk, cfg, {
      fileName: `${meta.fileName}.chunk_${chunkIndex.toString().padStart(6, '0')}`,
      contentType: meta.contentType,
    });

    const hash = await sha256Hex(chunk);

    yield {
      chunkIndex,
      offset: byteOffset,
      length: chunk.length,
      encoded: result.encoded,
      hash,
    };

    byteOffset += chunk.length;
    chunkIndex++;
  }
}

// ---------------------------------------------------------------------------
// Streaming Decode
// ---------------------------------------------------------------------------

/**
 * Streaming decode: receive reads in batches, decode each chunk, yield recovered data.
 *
 * Reads are grouped by chunk index (extracted from metadata). Each chunk is
 * decoded independently. Peak memory is O(chunkSize * coverage) because we
 * buffer reads for the current chunk before decoding.
 *
 * @param reads       Async iterable of sequencing reads.
 * @param metadata    Codec metadata for the current chunk.
 * @param cfg         Codec configuration (must match encode config).
 * @param fwdPrimer   Forward primer used during encoding.
 * @param revPrimer   Reverse primer used during encoding.
 * @yields StreamDecodeChunk for each decoded chunk.
 */
export async function* streamDecode(
  reads: AsyncIterable<SequencingRead>,
  metadata: CodecMetadata,
  cfg: CodecConfig = DEFAULT_CONFIG,
  fwdPrimer: string,
  revPrimer: string,
): AsyncGenerator<StreamDecodeChunk> {
  // Buffer reads for the current chunk, then decode.
  // Since reads arrive as an async iterable, we collect them all
  // for the current chunk before decoding (the caller controls batching).
  const bufferedReads: SequencingRead[] = [];
  let currentChunkIndex = 0;

  for await (const read of reads) {
    // We use a simple strategy: collect all reads, then decode once the
    // iterable is exhausted. For true multi-chunk streaming, the caller
    // should split reads by chunk index and call this function per chunk.
    bufferedReads.push(read);
  }

  // Decode the buffered reads for this chunk.
  if (bufferedReads.length > 0) {
    const result: DecodeResult = await decodeReads(
      bufferedReads,
      metadata,
      cfg,
      fwdPrimer,
      revPrimer,
      true, // useSoftInfo
    );

    const recoveredData = result.data ?? new Uint8Array(0);
    const hash = await sha256Hex(recoveredData);

    yield {
      chunkIndex: currentChunkIndex,
      data: recoveredData,
      hash,
      hashMatches: result.hashMatches,
    };
  }
}

// ---------------------------------------------------------------------------
// Convenience: full streaming roundtrip
// ---------------------------------------------------------------------------

/**
 * Perform a full streaming roundtrip: encode in chunks → simulate clean reads → decode.
 *
 * Useful for testing and benchmarking. Not for production (use real sequencing reads).
 *
 * @param data       Full data to encode/decode.
 * @param cfg        Codec configuration.
 * @param chunkSize  Chunk size in bytes.
 * @param meta       File metadata.
 * @returns Concatenated recovered data.
 */
export async function streamingRoundtrip(
  data: Uint8Array,
  cfg: CodecConfig = DEFAULT_CONFIG,
  chunkSize: number = 64 * 1024,
  meta: { fileName: string; contentType: string } = {
    fileName: 'stream-rp',
    contentType: 'application/octet-stream',
  },
): Promise<Uint8Array> {
  // Encode phase: collect all encode chunks.
  const encodeChunks: StreamEncodeChunk[] = [];
  const chunkIter = createChunkIterator(data, chunkSize);
  for await (const chunk of streamEncode(chunkIter, cfg, meta)) {
    encodeChunks.push(chunk);
  }

  // Decode phase: for each chunk, simulate clean reads and decode.
  const decodedParts: Uint8Array[] = [];
  for (let i = 0; i < encodeChunks.length; i++) {
    const ec = encodeChunks[i];
    const simResult: SimulationResult = simulate(ec.encoded.oligos, PRESET_CLEAN);

    // Create async iterable from reads.
    const readsIter = async function* (): AsyncGenerator<SequencingRead> {
      for (const read of simResult.reads) {
        yield read;
      }
    };

    for await (const dc of streamDecode(
      readsIter(),
      ec.encoded.metadata,
      cfg,
      ec.encoded.forwardPrimer,
      ec.encoded.reversePrimer,
    )) {
      decodedParts.push(dc.data);
    }
  }

  // Concatenate.
  const totalLen = decodedParts.reduce((s, p) => s + p.length, 0);
  const result = new Uint8Array(totalLen);
  let offset = 0;
  for (const part of decodedParts) {
    result.set(part, offset);
    offset += part.length;
  }
  return result;
}
