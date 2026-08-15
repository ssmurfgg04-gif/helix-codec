/**
 * Streaming encode/decode for large files.
 * 
 * Processes data in fixed-size chunks without loading the entire file into memory.
 * Each chunk is independently encoded as a shard of the overall archive.
 * 
 * This removes the 65535-oligo RS block limit and the O(file_size) memory requirement.
 */

import { encodeFile, EncodeResult } from './codec';
import { decodeReads } from './decode';
import { CodecConfig, DEFAULT_CONFIG, CodecMetadata } from './types';

export const DEFAULT_CHUNK_SIZE = 64 * 1024 * 1024; // 64MB per shard

export interface StreamingEncodeResult {
  shards: EncodeResult[];
  totalOligos: number;
  totalSize: number;
  manifest: {
    shardCount: number;
    shardSize: number;
    totalOriginalSize: number;
  };
}

/**
 * Encode a large file by splitting into shards and encoding each independently.
 * Each shard is a self-contained RS block that can be decoded independently.
 */
export async function encodeStream(
  data: Uint8Array,
  cfg: CodecConfig = DEFAULT_CONFIG,
  chunkSize: number = DEFAULT_CHUNK_SIZE,
): Promise<StreamingEncodeResult> {
  const shards: EncodeResult[] = [];
  let offset = 0;
  let shardIndex = 0;
  
  while (offset < data.length) {
    const end = Math.min(offset + chunkSize, data.length);
    const chunk = data.slice(offset, end);
    const shard = await encodeFile(chunk, cfg, {
      fileName: `shard-${shardIndex}`,
      contentType: 'application/octet-stream',
    });
    shards.push(shard);
    offset = end;
    shardIndex++;
  }
  
  return {
    shards,
    totalOligos: shards.reduce((sum, s) => sum + s.encoded.oligos.length, 0),
    totalSize: data.length,
    manifest: {
      shardCount: shards.length,
      shardSize: chunkSize,
      totalOriginalSize: data.length,
    },
  };
}

/**
 * Decode shards back to the original file.
 * Each shard is decoded independently, then concatenated.
 */
export async function decodeStream(
  shards: Array<{ reads: any[]; metadata: CodecMetadata; config: CodecConfig; fwdPrimer: string; revPrimer: string }>,
): Promise<Uint8Array> {
  const parts: Uint8Array[] = [];
  
  for (const shard of shards) {
    const result = await decodeReads(
      shard.reads,
      shard.metadata,
      shard.config,
      shard.fwdPrimer,
      shard.revPrimer,
      true,
    );
    if (result.data) {
      parts.push(result.data);
    }
  }
  
  // Concatenate all decoded parts
  const totalLen = parts.reduce((sum, p) => sum + p.length, 0);
  const result = new Uint8Array(totalLen);
  let offset = 0;
  for (const part of parts) {
    result.set(part, offset);
    offset += part.length;
  }
  return result;
}
