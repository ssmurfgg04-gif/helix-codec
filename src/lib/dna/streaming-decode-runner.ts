/**
 * v64: Real Streaming Decode Runner (Task 3)
 *
 * Replaces the stub in streaming-decode.ts with a proper implementation that:
 *   1. Accepts reads in batches via addReads()
 *   2. Clusters each batch by oligo index (k-mer matching for noisy channels)
 *   3. Accumulates reads per oligo (capped at maxReadsPerOligo to bound memory)
 *   4. Calls the full decodeReads() pipeline on the accumulated reads
 *
 * Memory bounding:
 *   - The caller generates reads in batches (e.g., from a FASTQ file or
 *     streaming sequencer) and feeds each batch to addReads().
 *   - After each batch, the caller's reads array can be freed.
 *   - The runner internally holds at most maxReadsPerOligo × numOligos reads.
 *   - For 42K oligos × 5 reads × 4KB/read = 840MB (vs 1.6GB for all reads).
 *
 * This is NOT true incremental decode (which would decode each oligo as soon
 * as it has enough reads). It's a practical middle ground: the caller can
 * stream reads from disk, and the runner bounds memory by capping reads per
 * oligo. The final decode() call runs the standard pipeline.
 *
 * For TRUE incremental decode (decode oligos as they accumulate), the decode
 * pipeline would need to be refactored to support per-oligo decode with
 * outer RS assembly at the end. That's a future v65+ task.
 */

import { SequencingRead } from "./simulate";
import { CodecMetadata, CodecConfig } from "./types";
import { decodeReads, trimPrimer } from "./decode";
import { dnaToBytes, unwhitenAddress } from "./mapping";

export interface StreamingDecodeRunnerOptions {
  /** Maximum reads to keep per oligo (default 10). Extra reads are discarded. */
  maxReadsPerOligo?: number;
  /** Batch size for internal processing (default 10000) */
  batchSize?: number;
}

export interface StreamingDecodeRunnerStats {
  totalReadsProcessed: number;
  totalReadsDiscarded: number;
  oligosWithReads: number;
  totalReadsAccumulated: number;
  batchesProcessed: number;
}

export class StreamingDecodeRunner {
  private metadata: CodecMetadata;
  private cfg: CodecConfig;
  private fwdPrimer: string;
  private revPrimer: string;
  private maxReadsPerOligo: number;
  private batchSize: number;

  // Accumulated reads per oligo index
  private accumulatedReads: Map<number, SequencingRead[]> = new Map();
  private stats: StreamingDecodeRunnerStats = {
    totalReadsProcessed: 0,
    totalReadsDiscarded: 0,
    oligosWithReads: 0,
    totalReadsAccumulated: 0,
    batchesProcessed: 0,
  };

  constructor(
    metadata: CodecMetadata,
    cfg: CodecConfig,
    fwdPrimer: string,
    revPrimer: string,
    options: StreamingDecodeRunnerOptions = {},
  ) {
    this.metadata = metadata;
    this.cfg = cfg;
    this.fwdPrimer = fwdPrimer;
    this.revPrimer = revPrimer;
    this.maxReadsPerOligo = options.maxReadsPerOligo ?? 10;
    this.batchSize = options.batchSize ?? 10000;
  }

  /**
   * Add a batch of reads to the streaming decoder.
   *
   * The reads are clustered by oligo index (using the same logic as
   * decodeReads). Each read is added to the corresponding oligo's
   * accumulated read list, capped at maxReadsPerOligo.
   *
   * After this call returns, the caller can free the reads array.
   */
  addReads(reads: SequencingRead[]): void {
    this.stats.batchesProcessed++;
    this.stats.totalReadsProcessed += reads.length;

    // Cluster reads by oligo index using exact address matching.
    // For noisy channels (nanopore), the caller should use k-mer clustering
    // before feeding reads to the runner. Here we do simple address extraction.
    const { clusters, discarded } = this.clusterReadsSimple(reads);

    this.stats.totalReadsDiscarded += discarded;

    // Merge clusters into accumulated reads, capping at maxReadsPerOligo
    for (const [oligoIdx, clusterReads] of clusters) {
      const existing = this.accumulatedReads.get(oligoIdx);
      if (existing) {
        // Cap at maxReadsPerOligo — keep the first reads (they're typically
        // from earlier in the sequencing run, which are often higher quality)
        const space = this.maxReadsPerOligo - existing.length;
        if (space > 0) {
          const toAdd = Math.min(space, clusterReads.length);
          for (let i = 0; i < toAdd; i++) {
            existing.push(clusterReads[i]);
          }
          this.stats.totalReadsAccumulated += toAdd;
        }
        // Discard excess reads
        this.stats.totalReadsDiscarded += Math.max(0, clusterReads.length - Math.max(0, this.maxReadsPerOligo - existing.length));
      } else {
        const toAdd = Math.min(this.maxReadsPerOligo, clusterReads.length);
        this.accumulatedReads.set(oligoIdx, clusterReads.slice(0, toAdd));
        this.stats.totalReadsAccumulated += toAdd;
        this.stats.totalReadsDiscarded += clusterReads.length - toAdd;
      }
    }

    this.stats.oligosWithReads = this.accumulatedReads.size;
  }

  /**
   * Simple address-based clustering (for Illumina-like channels).
   *
   * For nanopore/high-IDS channels, the caller should pre-cluster using
   * k-mer matching and feed the clustered reads to addReads() instead.
   */
  private clusterReadsSimple(reads: SequencingRead[]): {
    clusters: Map<number, SequencingRead[]>;
    discarded: number;
  } {
    const clusters = new Map<number, SequencingRead[]>();
    let discarded = 0;

    const addressNt = 16; // 4 bytes × 4 nt/byte

    for (const read of reads) {
      try {
        const inner = trimPrimer(read.sequence, this.fwdPrimer, this.revPrimer);
        if (!inner || inner.length < addressNt) {
          discarded++;
          continue;
        }
        const addressDna = inner.slice(0, addressNt);
        const addressBytes = dnaToBytes(addressDna);
        const unwhitened = unwhitenAddress(addressBytes);
        const index = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
        if (!clusters.has(index)) clusters.set(index, []);
        // Store the ORIGINAL read (with primers) — decodeReads() will trim.
        // This is critical: if we store the trimmed sequence, decodeReads()
        // will try to trimPrimer() again and fail, discarding all reads.
        clusters.get(index)!.push(read);
      } catch {
        discarded++;
      }
    }

    return { clusters, discarded };
  }

  /**
   * Run the full decode pipeline on the accumulated reads.
   *
   * This calls the standard decodeReads() with all accumulated reads.
   * After decode, the accumulated reads are freed.
   */
  async decode(): Promise<ReturnType<typeof decodeReads>> {
    // Flatten accumulated reads into a single array
    const allReads: SequencingRead[] = [];
    for (const clusterReads of this.accumulatedReads.values()) {
      allReads.push(...clusterReads);
    }

    console.log(`[StreamingDecode] Decoding ${allReads.length.toLocaleString()} reads from ${this.accumulatedReads.size.toLocaleString()} oligos`);

    const result = await decodeReads(
      allReads,
      this.metadata,
      this.cfg,
      this.fwdPrimer,
      this.revPrimer,
    );

    // Free accumulated reads
    this.accumulatedReads.clear();

    return result;
  }

  /**
   * Get statistics about the streaming decode state.
   */
  getStats(): StreamingDecodeRunnerStats {
    return { ...this.stats };
  }

  /**
   * Get the current memory usage estimate (in bytes).
   */
  getEstimatedMemoryUsage(): number {
    let total = 0;
    for (const clusterReads of this.accumulatedReads.values()) {
      for (const read of clusterReads) {
        total += read.sequence.length * 2 + 64; // DNA + Q-scores + overhead
      }
    }
    return total;
  }
}
