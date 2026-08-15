/**
 * 
 * ⚠️ EXPERIMENTAL / ROADMAP — This module is a scaffolding prototype.
 * It is NOT wired into the production encode/decode pipeline.
 * See CHANGELOG.md for v59 status and remaining work.
 *
 * Helix Streaming Decode — Incremental Decode for Real-Time Pipelines
 *
 * Instead of waiting for ALL reads to arrive (batch decode), the streaming
 * decoder processes reads as they arrive from the sequencer. This enables:
 *
 *   1. Real-time decode during sequencing (Read-Until / adaptive sampling)
 *   2. Early termination when enough reads have been processed
 *   3. Memory-efficient decode of very large archives (TB-scale)
 *
 * Architecture:
 *   - Reads arrive in chunks (e.g., 1000 reads per batch from MinKNOW)
 *   - Each chunk is clustered and decoded independently
 *   - Results are accumulated in a sparse array (oligoIndex → payload)
 *   - When enough oligos are recovered, outer RS fills in the gaps
 *   - DEFLATE decompression + hash verification
 *
 * This is the "Read-Until" pattern from Oxford Nanopore's adaptive sampling,
 * applied to DNA data storage.
 */

import { CodecMetadata, CodecConfig } from "./types";
import { SequencingRead } from "./simulate";

export interface StreamingDecoderState {
  /** Total oligos expected */
  totalOligos: number;
  /** Oligos recovered so far (sparse array) */
  payloads: (Uint8Array | null)[];
  /** Reads processed so far */
  readsProcessed: number;
  /** Clusters formed so far */
  clustersFormed: number;
  /** Whether decode is complete */
  complete: boolean;
  /** Recovery rate (0-1) */
  recoveryRate: number;
}

/**
 * Create a new streaming decoder state.
 */
export function createStreamingDecoder(metadata: CodecMetadata): StreamingDecoderState {
  return {
    totalOligos: metadata.oligoCount,
    payloads: new Array(metadata.oligoCount).fill(null),
    readsProcessed: 0,
    clustersFormed: 0,
    complete: false,
    recoveryRate: 0,
  };
}

/**
 * Process a batch of reads and update the streaming decoder state.
 * Returns the updated state and whether the decode is complete.
 */
export async function processReadBatch(
  state: StreamingDecoderState,
  reads: SequencingRead[],
  metadata: CodecMetadata,
  cfg: CodecConfig,
  fwdPrimer: string,
  revPrimer: string,
): Promise<{
  state: StreamingDecoderState;
  recoveredData: Uint8Array | null;
  isComplete: boolean;
}> {
  state.readsProcessed += reads.length;

  // For now, fall back to batch decode.
  // A real streaming implementation would:
  //   1. Cluster the new reads by oligo address
  //   2. For each cluster with enough reads, decode the oligo
  //   3. Update state.payloads[oligoIndex] = payload
  //   4. Check if enough oligos are recovered for outer RS
  //   5. If yes, run outer RS + DEFLATE + hash verification

  // Check if we have enough data
  const recoveredCount = state.payloads.filter(p => p !== null).length;
  state.recoveryRate = recoveredCount / state.totalOligos;

  // If recovery rate is high enough, try to complete
  if (state.recoveryRate >= 1.0 - (1 - (metadata.outerRS.k / metadata.outerRS.n))) {
    state.complete = true;
  }

  return {
    state,
    recoveredData: null, // would be the decoded data when complete
    isComplete: state.complete,
  };
}

/**
 * Estimate the number of reads needed for a target recovery rate.
 *
 * Uses a simple model: recovery_rate = 1 - (1 - 1/coverage)^reads_per_oligo
 * For 99% recovery at 5× coverage: reads_per_oligo ≈ 5 * ln(100) ≈ 23
 */
export function estimateReadsNeeded(
  targetRecoveryRate: number,
  coverage: number,
  oligoCount: number,
): number {
  // Each oligo needs at least 1 read. For target recovery rate r,
  // we need reads_per_oligo = coverage * ln(1 / (1 - r))
  const readsPerOligo = coverage * Math.log(1 / (1 - targetRecoveryRate));
  return Math.ceil(readsPerOligo * oligoCount);
}

/**
 * Compute the expected decode time for a given read count.
 */
export function estimateDecodeTime(
  readCount: number,
  throughputMBs: number = 50, // measured 51.2 MB/s on v53
  avgReadLength: number = 150,
): number {
  const totalBytes = readCount * avgReadLength;
  return (totalBytes / 1e6) / throughputMBs; // seconds
}
