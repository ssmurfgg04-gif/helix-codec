/**
 * LSM-Tree Journal for Incremental Archive Compaction (LAB-DB pattern).
 *
 * Instead of rewriting the full archive on every mutation, we append
 * journal entries to an L0 buffer. When L0 fills (default: 64KB),
 * we compact it into L1 (sorted run). Lower levels compact during
 * scheduled low-traffic windows.
 *
 * Each 64KB packed block is a B+ tree leaf.
 * Journal entry: ~16 bytes per mutation (block_id + bit_offset + new_base + timestamp).
 *
 * Synthesis queue: Batches journal deltas into oligo synthesis orders.
 * The codec never touches a pipette.
 *
 * Reference:
 *   - LAB-DB: Log-Structured Append-Only B-Tree with incremental partial compaction
 *   - LevelDB/RocksDB: LSM architecture patterns
 */

// ---------------------------------------------------------------------------
// Interfaces
// ---------------------------------------------------------------------------

/** A single mutation entry in the journal. */
export interface JournalEntry {
  /** Block ID (which oligo block is mutated). */
  blockId: number;
  /** Bit offset within the block (0-based). */
  bitOffset: number;
  /** New base value at this position (0=A, 1=C, 2=G, 3=T). */
  newBase: number;
  /** Timestamp of the mutation (ms since epoch). */
  timestamp: number;
}

/** Options for the LSM journal. */
export interface LsmJournalOptions {
  /** Maximum size of the L0 buffer in bytes before compaction trigger. Default: 65536 (64KB). */
  l0Size?: number;
  /**
   * Compaction trigger ratio — compact L0→L1 when L0 usage exceeds this fraction.
   * Default: 0.75 (compact when L0 is 75% full).
   */
  compactionTrigger?: number;
  /** Maximum number of L1 sorted runs before triggering L1→L2 compaction. Default: 4. */
  maxL1Runs?: number;
}

/** Result of an L1→L2 compaction. */
export interface CompactionResult {
  /** Total entries that were merged across all L1 runs + existing L2. */
  entriesMerged: number;
  /** Entries removed due to tombstones or deduplication. */
  entriesDeleted: number;
  /** Number of L1 runs that were compacted. */
  runsCompacted: number;
}

/** A batch of mutations ready for synthesis (DNA oligo production). */
export interface SynthesisBatch {
  /** Unique batch ID. */
  batchId: string;
  /** Mutations in this batch. */
  entries: JournalEntry[];
  /** Total byte size of the serialized batch. */
  sizeBytes: number;
  /** Timestamp when the batch was created. */
  createdAt: number;
}

/** Statistics about the LSM journal state. */
export interface LsmJournalStats {
  /** Number of entries currently in L0 buffer. */
  l0Count: number;
  /** Current L0 buffer usage in bytes. */
  l0Bytes: number;
  /** Maximum L0 buffer size in bytes. */
  l0MaxBytes: number;
  /** L0 fill ratio (0..1). */
  l0FillRatio: number;
  /** Number of L1 sorted runs. */
  l1RunCount: number;
  /** Total entries across all L1 runs. */
  l1TotalEntries: number;
  /** Number of L2 compactions performed. */
  l2CompactionCount: number;
  /** Total entries ever appended. */
  totalAppended: number;
  /** Number of compactions performed (L0→L1). */
  compactionCount: number;
  /** Number of synthesis batches created. */
  synthesisBatchCount: number;
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/** Default L0 buffer size: 64KB. */
const DEFAULT_L0_SIZE = 64 * 1024;

/** Default compaction trigger: compact when L0 is 75% full. */
const DEFAULT_COMPACTION_TRIGGER = 0.75;

/** Default max L1 sorted runs before L1→L2 compaction. */
const DEFAULT_MAX_L1_RUNS = 4;

/** Sentinel value for tombstone entries (delete markers). Valid bases are 0-3. */
export const TOMBSTONE = 0xFF;

/** Size of a serialized journal entry in bytes: 4 + 4 + 1 + 8 = 17 bytes (padded to 16 with packing). */
const ENTRY_SERIALIZED_SIZE = 16;

// ---------------------------------------------------------------------------
// Internal: Sorted Run (L1)
// ---------------------------------------------------------------------------

/**
 * A sorted run of journal entries at L1.
 * Entries are sorted by (blockId, bitOffset) for efficient point queries.
 */
class SortedRun {
  entries: JournalEntry[] = [];

  /** Insert an entry, maintaining sort order by (blockId, bitOffset). */
  insert(entry: JournalEntry): void {
    // Binary search for insertion point.
    let lo = 0;
    let hi = this.entries.length;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      const e = this.entries[mid];
      if (e.blockId < entry.blockId || (e.blockId === entry.blockId && e.bitOffset < entry.bitOffset)) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    // If entry at same (blockId, bitOffset) exists, update it (last-write-wins).
    if (lo < this.entries.length) {
      const existing = this.entries[lo];
      if (existing.blockId === entry.blockId && existing.bitOffset === entry.bitOffset) {
        this.entries[lo] = entry;
        return;
      }
    }
    this.entries.splice(lo, 0, entry);
  }

  /** Point query: find the entry with the given (blockId, bitOffset). */
  get(blockId: number, bitOffset: number): JournalEntry | undefined {
    let lo = 0;
    let hi = this.entries.length;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      const e = this.entries[mid];
      const cmp = e.blockId - blockId || e.bitOffset - bitOffset;
      if (cmp < 0) {
        lo = mid + 1;
      } else if (cmp > 0) {
        hi = mid;
      } else {
        return e;
      }
    }
    return undefined;
  }
}

// ---------------------------------------------------------------------------
// LsmJournal
// ---------------------------------------------------------------------------

/**
 * LSM-Tree Journal for incremental compaction of .hlx archives.
 *
 * Architecture:
 *   L0 (memtable): In-memory buffer, append-only. Compacted to L1 when full.
 *   L1 (sorted runs): One or more sorted runs. Merged to L2 during low-traffic.
 *   L2 (compacted): Single sorted run. Fully compacted.
 *
 * Point query path: L0 → L1 (newest to oldest) → L2.
 * Last-write-wins: the most recent entry for a given (blockId, bitOffset) is returned.
 */
export class LsmJournal {
  // L0 buffer (append-only, unsorted)
  private l0: JournalEntry[] = [];
  private l0Bytes: number = 0;

  // L1 sorted runs
  private l1Runs: SortedRun[] = [];

  // L2 fully-compacted run
  private l2: SortedRun | null = null;

  // Configuration
  private readonly l0MaxSize: number;
  private readonly compactionTrigger: number;
  private readonly maxL1Runs: number;

  // Stats
  private totalAppended: number = 0;
  private compactionCount: number = 0;
  private l2CompactionCount: number = 0;
  private synthesisBatchCount: number = 0;

  constructor(options?: LsmJournalOptions) {
    this.l0MaxSize = options?.l0Size ?? DEFAULT_L0_SIZE;
    this.compactionTrigger = options?.compactionTrigger ?? DEFAULT_COMPACTION_TRIGGER;
    this.maxL1Runs = options?.maxL1Runs ?? DEFAULT_MAX_L1_RUNS;
  }

  // -------------------------------------------------------------------------
  // Append
  // -------------------------------------------------------------------------

  /**
   * Append a mutation entry to the L0 buffer.
   * If L0 exceeds the compaction trigger ratio, automatically compact to L1.
   *
   * @param entry The journal entry to append.
   */
  append(entry: JournalEntry): void {
    this.l0.push(entry);
    this.l0Bytes += ENTRY_SERIALIZED_SIZE;
    this.totalAppended++;

    // Auto-compact if L0 exceeds trigger.
    if (this.l0Bytes >= this.l0MaxSize * this.compactionTrigger) {
      this.flush();
    }
  }

  // -------------------------------------------------------------------------
  // Delete (tombstone)
  // -------------------------------------------------------------------------

  /**
   * Mark a key as deleted by appending a tombstone entry.
   * Tombstones are resolved during L1→L2 compaction.
   *
   * @param blockId    Block ID to delete.
   * @param bitOffset  Bit offset within the block.
   */
  delete(blockId: number, bitOffset: number): void {
    this.append({
      blockId,
      bitOffset,
      newBase: TOMBSTONE,
      timestamp: Date.now(),
    });
  }

  // -------------------------------------------------------------------------
  // Flush (L0 → L1 compaction)
  // -------------------------------------------------------------------------

  /**
   * Compact L0 buffer into a new L1 sorted run.
   * Returns the serialized journal (for persistence).
   *
   * @returns Serialized journal bytes (L0 entries in binary format).
   */
  flush(): Uint8Array {
    if (this.l0.length === 0) {
      return new Uint8Array(0);
    }

    // Sort L0 entries by (blockId, bitOffset) for the sorted run.
    const sorted = [...this.l0].sort((a, b) => {
      return a.blockId - b.blockId || a.bitOffset - b.bitOffset;
    });

    // Create new L1 sorted run.
    const run = new SortedRun();
    for (const entry of sorted) {
      run.insert(entry);
    }
    this.l1Runs.push(run);

    // Serialize L0 entries to binary format.
    const serialized = this.serializeEntries(this.l0);

    // Clear L0.
    this.l0 = [];
    this.l0Bytes = 0;
    this.compactionCount++;

    return serialized;
  }

  // -------------------------------------------------------------------------
  // Point Query
  // -------------------------------------------------------------------------

  /**
   * Point query: find the value at (blockId, bitOffset).
   * Searches L0 first (most recent), then L1 runs (newest first), then L2.
   * Returns undefined if no entry exists for this key.
   *
   * @param blockId    Block ID to query.
   * @param bitOffset  Bit offset within the block.
   * @returns The base value (0-3) or undefined if not found.
   */
  get(blockId: number, bitOffset: number): number | undefined {
    // Search L0 (most recent first — iterate in reverse).
    for (let i = this.l0.length - 1; i >= 0; i--) {
      const e = this.l0[i];
      if (e.blockId === blockId && e.bitOffset === bitOffset) {
        // Tombstone means deleted — stop searching lower levels.
        return e.newBase === TOMBSTONE ? undefined : e.newBase;
      }
    }

    // Search L1 runs (newest first).
    for (let i = this.l1Runs.length - 1; i >= 0; i--) {
      const entry = this.l1Runs[i].get(blockId, bitOffset);
      if (entry !== undefined) {
        return entry.newBase === TOMBSTONE ? undefined : entry.newBase;
      }
    }

    // Search L2.
    if (this.l2 !== null) {
      const entry = this.l2.get(blockId, bitOffset);
      if (entry !== undefined) {
        return entry.newBase === TOMBSTONE ? undefined : entry.newBase;
      }
    }

    return undefined;
  }

  // -------------------------------------------------------------------------
  // Synthesis Queue
  // -------------------------------------------------------------------------

  /**
   * Batch pending journal deltas into synthesis orders for DNA oligo production.
   * Each batch represents a set of mutations that should be synthesized together.
   *
   * @returns Array of synthesis batches ready for production.
   */
  getSynthesisQueue(): SynthesisBatch[] {
    const batches: SynthesisBatch[] = [];

    // Collect all pending mutations: L0 buffer + L1 runs not yet compacted to L2.
    const pending: JournalEntry[] = [...this.l0];
    for (const run of this.l1Runs) {
      pending.push(...run.entries);
    }

    if (pending.length === 0) {
      return batches;
    }

    // Create batches grouped by blockId (each block maps to one or more oligos).
    const byBlock = new Map<number, JournalEntry[]>();
    for (const entry of pending) {
      let list = byBlock.get(entry.blockId);
      if (!list) {
        list = [];
        byBlock.set(entry.blockId, list);
      }
      list.push(entry);
    }

    const blockIds = Array.from(byBlock.keys());
    for (const blockId of blockIds) {
      const entries = byBlock.get(blockId)!;
      const batch: SynthesisBatch = {
        batchId: `synth_${blockId}_${Date.now()}`,
        entries,
        sizeBytes: entries.length * ENTRY_SERIALIZED_SIZE,
        createdAt: Date.now(),
      };
      batches.push(batch);
      this.synthesisBatchCount++;
    }

    return batches;
  }

  // -------------------------------------------------------------------------
  // L1 → L2 Compaction
  // -------------------------------------------------------------------------

  /**
   * Trigger L1→L2 compaction during low-traffic window.
   * Merges all L1 sorted runs into a single L2 run:
   *   1. Collect all entries from L1 runs + existing L2
   *   2. Deduplicate by (blockId, bitOffset) — keep only the latest (highest timestamp)
   *   3. Remove tombstoned entries (delete markers) entirely
   *   4. Sort result by (blockId, bitOffset)
   *   5. Replace L1 runs with the merged L2 run
   *
   * @returns Compaction statistics.
   */
  compact(): CompactionResult {
    const runsCompacted = this.l1Runs.length;

    if (this.l1Runs.length === 0 && this.l2 === null) {
      return { entriesMerged: 0, entriesDeleted: 0, runsCompacted: 0 };
    }

    // Step 1: Collect all entries from L1 runs + existing L2.
    const allEntries: JournalEntry[] = [];
    if (this.l2 !== null) {
      allEntries.push(...this.l2.entries);
    }
    for (const run of this.l1Runs) {
      allEntries.push(...run.entries);
    }

    const entriesMerged = allEntries.length;

    // Step 2: Deduplicate by (blockId, bitOffset) — keep only the latest timestamp.
    // We use a map keyed by "blockId:bitOffset" to pick the newest entry per key.
    const latestByKey = new Map<string, JournalEntry>();
    for (const entry of allEntries) {
      const key = `${entry.blockId}:${entry.bitOffset}`;
      const existing = latestByKey.get(key);
      if (!existing || entry.timestamp > existing.timestamp) {
        latestByKey.set(key, entry);
      }
    }

    // Step 3: Remove tombstoned entries (delete markers) entirely.
    let entriesDeleted = 0;
    const surviving: JournalEntry[] = [];
    for (const entry of latestByKey.values()) {
      if (entry.newBase === TOMBSTONE) {
        entriesDeleted++;
      } else {
        surviving.push(entry);
      }
    }

    // Step 4: Sort surviving entries by (blockId, bitOffset).
    surviving.sort((a, b) => a.blockId - b.blockId || a.bitOffset - b.bitOffset);

    // Step 5: Build new L2 sorted run and replace L1 runs.
    const merged = new SortedRun();
    for (const entry of surviving) {
      merged.insert(entry);
    }

    this.l2 = merged;
    this.l1Runs = [];
    this.l2CompactionCount++;

    return { entriesMerged, entriesDeleted, runsCompacted };
  }

  // -------------------------------------------------------------------------
  // Stats
  // -------------------------------------------------------------------------

  /**
   * Get current journal statistics.
   *
   * @returns Snapshot of journal state.
   */
  getStats(): LsmJournalStats {
    const l1TotalEntries = this.l1Runs.reduce((sum, run) => sum + run.entries.length, 0);

    return {
      l0Count: this.l0.length,
      l0Bytes: this.l0Bytes,
      l0MaxBytes: this.l0MaxSize,
      l0FillRatio: this.l0Bytes / this.l0MaxSize,
      l1RunCount: this.l1Runs.length,
      l1TotalEntries,
      l2CompactionCount: this.l2CompactionCount,
      totalAppended: this.totalAppended,
      compactionCount: this.compactionCount,
      synthesisBatchCount: this.synthesisBatchCount,
    };
  }

  // -------------------------------------------------------------------------
  // Serialization
  // -------------------------------------------------------------------------

  /**
   * Serialize journal entries to a compact binary format.
   *
   * Layout per entry (16 bytes):
   *   [0..3]  blockId    (uint32 LE)
   *   [4..7]  bitOffset  (uint32 LE)
   *   [8]     newBase    (uint8, 0-3)
   *   [9..15] timestamp  (uint56 LE, truncated from float64 — ms since epoch fits in 6 bytes until year 10889)
   *
   * @param entries Entries to serialize.
   * @returns Compact binary representation.
   */
  private serializeEntries(entries: JournalEntry[]): Uint8Array {
    const buf = new Uint8Array(entries.length * ENTRY_SERIALIZED_SIZE);
    const view = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);

    for (let i = 0; i < entries.length; i++) {
      const e = entries[i];
      const off = i * ENTRY_SERIALIZED_SIZE;
      view.setUint32(off, e.blockId, true);        // LE
      view.setUint32(off + 4, e.bitOffset, true);  // LE
      view.setUint8(off + 8, e.newBase);
      // Store timestamp as float64 for simplicity (8 bytes, but we only have 7 remaining).
      // We use 7 bytes for timestamp (uint56), which covers until year ~10889.
      const ts = e.timestamp;
      // Byte 9-15: uint56 LE
      const low = ts & 0xFFFFFF;
      const high = (ts / 0x1000000) & 0xFFFFFFFF;
      view.setUint24(off + 9, low, true);
      view.setUint32(off + 12, high, true);
    }

    return buf;
  }

  /**
   * Deserialize journal entries from the compact binary format.
   *
   * @param buf Binary data produced by serializeEntries.
   * @returns Deserialized journal entries.
   */
  static deserializeEntries(buf: Uint8Array): JournalEntry[] {
    const entries: JournalEntry[] = [];
    const view = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
    const count = Math.floor(buf.length / ENTRY_SERIALIZED_SIZE);

    for (let i = 0; i < count; i++) {
      const off = i * ENTRY_SERIALIZED_SIZE;
      const blockId = view.getUint32(off, true);
      const bitOffset = view.getUint32(off + 4, true);
      const newBase = view.getUint8(off + 8);
      const low = view.getUint24(off + 9, true);
      const high = view.getUint32(off + 12, true);
      const timestamp = low + high * 0x1000000;

      entries.push({ blockId, bitOffset, newBase, timestamp });
    }

    return entries;
  }
}

// ---------------------------------------------------------------------------
// DataView extensions for uint24
// ---------------------------------------------------------------------------

/** Read an unsigned 24-bit integer (3 bytes, little-endian). */
declare global {
  interface DataView {
    getUint24(byteOffset: number, littleEndian?: boolean): number;
    setUint24(byteOffset: number, value: number, littleEndian?: boolean): void;
  }
}

// Polyfill getUint24 / setUint24 on DataView prototype.
DataView.prototype.getUint24 = function (byteOffset: number, littleEndian?: boolean): number {
  const b0 = this.getUint8(byteOffset);
  const b1 = this.getUint8(byteOffset + 1);
  const b2 = this.getUint8(byteOffset + 2);
  return littleEndian ? b0 | (b1 << 8) | (b2 << 16) : (b0 << 16) | (b1 << 8) | b2;
};

DataView.prototype.setUint24 = function (byteOffset: number, value: number, littleEndian?: boolean): void {
  if (littleEndian) {
    this.setUint8(byteOffset, value & 0xFF);
    this.setUint8(byteOffset + 1, (value >>> 8) & 0xFF);
    this.setUint8(byteOffset + 2, (value >>> 16) & 0xFF);
  } else {
    this.setUint8(byteOffset, (value >>> 16) & 0xFF);
    this.setUint8(byteOffset + 1, (value >>> 8) & 0xFF);
    this.setUint8(byteOffset + 2, value & 0xFF);
  }
};
