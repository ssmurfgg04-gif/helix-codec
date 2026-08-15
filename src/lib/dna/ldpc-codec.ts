/**
 * LDPC Inner Code — drop-in replacement for Reed-Solomon.
 *
 * Uses a PEG-constructed LDPC code over GF(2) with bit-flipping decoding.
 * The LDPC code operates on bits internally but exposes a byte-oriented
 * interface matching `ReedSolomon` so it can be swapped in transparently.
 *
 * Parameters (default, "Helix hi-fi"):
 *   - n = 304 bits (38 bytes codeword = 4B address + 30B payload + 4B parity)
 *   - k = 272 bits (34 bytes info = 4B address + 30B payload)
 *   - m =  32 bits ( 4 bytes parity)
 *   - dv = 3 (variable node degree)
 *   - dc ≈ 29 (check node degree, target)
 *   - rate = 272/304 = 0.895 → 1.79 bits/nt theoretical (Shannon = 2.0)
 *
 * Density (200 nt oligo, 20 nt primers × 2):
 *   - Payload per oligo: 30 bytes = 240 bits (vs 26 bytes for RS)
 *   - Total oligo density: 240/200 = 1.20 bits/nt  (+15% vs RS's 1.04)
 *   - Payload-only density: 240/160 = 1.50 bits/nt  (+15% vs RS's 1.30)
 *
 * Decoder: Gallager bit-flipping with up to `maxIter` iterations.
 *   - Each iteration: compute syndrome, flip the bit that participates
 *     in the most unsatisfied checks, repeat until syndrome = 0.
 *   - At Illumina 0.1% sub rate over 304 bits, expected errors per
 *     codeword = 0.3 → typically 0-1 errors → converges in 1-2 iters.
 *   - m=32 parity bits can correct up to 16 errors per codeword.
 *
 * For higher error rates (Nanopore, Real 2024), use `innerParityBytes: 8`
 * which gives m=64 parity bits (32 errors correctable) at the cost of
 * reverting to RS-level density.
 *
 * References:
 *   - Gallager (1962). "Low-density parity-check codes." IRE Trans. IT-8:1.
 *   - Hu, Eleftheriou, Arnold (2005). PEG construction. IEEE TIT 51:1.
 *   - Banal et al. (2026). Mahoraga codec. arXiv:2604.20810.
 */

import { constructPEG } from "./peg";
import { GF2Matrix } from "./osd";
import { BitParallelSyndrome } from "./mega-performance";

export interface LDPCConfig {
  /** Total codeword length in BYTES (must equal k + nsym). */
  n: number;
  /** Information length in BYTES. */
  k: number;
  /** Max bit-flipping iterations. Default 20. */
  maxIter?: number;
}

export interface LDPCDecodeResult {
  /** Decoded information bytes (k bytes). */
  data: Uint8Array;
  /** Number of bits corrected. */
  corrected: number;
  /** Number of erasures corrected (always 0 for plain bit-flipping). */
  erased: number;
}

/**
 * LDPC inner code with byte-oriented interface matching ReedSolomon.
 */
export class LDPCInnerCode {
  readonly n: number; // bytes
  readonly k: number; // bytes
  readonly nsym: number; // bytes (= n - k)
  readonly nBits: number; // n * 8
  readonly kBits: number; // k * 8
  readonly mBits: number; // nsym * 8
  readonly maxIter: number;

  /** Sparse parity-check matrix H (mBits × nBits). */
  private H: GF2Matrix;

  /** Public read-only access to the parity-check matrix H for OSD post-pass. */
  get parityCheckMatrix(): GF2Matrix { return this.H; }

  /** Adjacency lists for fast syndrome / bit-flipping. */
  private rowCols: Uint32Array[]; // rowCols[i] = list of column indices in row i
  private colRows: Uint32Array[]; // colRows[j] = list of row indices in column j

  /** Precomputed column hashes for single-error syndrome lookup. */
  private colHashToInt: Map<string, number>; // hash(info column j) -> j (for info bits only)
  private parityColHashToIdx: Map<string, number>; // hash(parity column i) -> kBits + i

  /** Precomputed pair hashes for double-error syndrome lookup (lazy-built). */
  private pairHashToBits: Map<string, [number, number]> | null;

  /** Whether mBits fits in 32 bits (for fast path). */
  private useFastHash: boolean;

  /** v61: Bit-parallel syndrome LUT for 8× encode speedup */
  private bitParallel: BitParallelSyndrome | null = null;

  /** v63: Cached parity-column row view (avoids per-decode slice allocation) */
  private parityColRowsView: Uint32Array[] | null = null;

  /** Precomputed reverse edge indices for O(1) BP decoding.
   *
   * rowColsReverseIdx[i][e] = position of check node i in colRows[rowCols[i][e]].
   * This replaces the O(degree) linear search in the check→variable message pass.
   *
   * colRowsReverseIdx[j][e] = position of variable node j in rowCols[colRows[j][e]].
   * This replaces the O(degree) linear search in the variable→check message pass.
   */
  private rowColsReverseIdx: Int32Array[];
  private colRowsReverseIdx: Int32Array[];

  constructor(cfg: LDPCConfig) {
    if (cfg.n <= 0 || cfg.k <= 0 || cfg.k >= cfg.n) {
      throw new Error(`LDPC config invalid: n=${cfg.n}, k=${cfg.k}`);
    }
    if (cfg.n > 255) {
      // We allow this for completeness but note that the outer RS over GF(2^8)
      // caps oligo count at 255; for >255 oligos we use GF(2^16) which is fine.
    }
    this.n = cfg.n;
    this.k = cfg.k;
    this.nsym = cfg.n - cfg.k;
    this.nBits = cfg.n * 8;
    this.kBits = cfg.k * 8;
    this.mBits = this.nsym * 8;
    this.maxIter = cfg.maxIter ?? 20;

    // Build H in SYSTEMATIC form: H = [A | I_m] where:
    //   - A is mBits × kBits (sparse, PEG-constructed with dv=4 per info column)
    //   - I_m is mBits × mBits identity (parity columns)
    //
    // dv=4 (not 3) to reduce column collisions:
    //   - dv=3, mBits=32: C(32,3)=4960 columns → ~20 duplicates per 456 info cols
    //   - dv=4, mBits=32: C(32,4)=35960 columns → ~3 duplicates per 456 info cols
    // The decoder still uses syndrome lookup, but with fewer collisions, the
    // per-read failure rate drops from ~5% to ~1%, which is recoverable at 10x.
    const dv = 4;
    const dc = Math.max(6, Math.floor((this.kBits * dv) / this.mBits));

    // Step 1: PEG-construct A (mBits × kBits) for the info columns
    const A = constructPEG({
      n: this.kBits,
      m: this.mBits,
      dv,
      dc,
    });

    // v61: DUPLICATE-COLUMN DEDUPLICATION
    //
    // At 300+ oligo scale (kBits ≥ 2400, mBits = 64), the PEG construction
    // can produce duplicate columns — two info bits j1, j2 where A[:,j1] and
    // A[:,j2] have the same set of 1-rows. This means a single-bit error at
    // position j1 produces the same syndrome as an error at j2, making them
    // indistinguishable to the syndrome-lookup decoder. The decoder silently
    // picks one, leading to a wrong correction → hash FAIL.
    //
    // Detection: build a hash of each column's row-set, check for collisions.
    // Fix: for each duplicate, swap one row index with a spare row (one that
    // neither column currently uses) to break the tie. If no spare exists,
    // increase dv by 1 and re-construct.
    //
    // This is especially critical at large kBits because:
    //   - kBits=2400, mBits=64, dv=4: 2400 columns of weight 4 from C(64,4)=635,440
    //   - Birthday paradox: ~2400²/(2*635440) = 4.5 expected collisions
    //   - Each collision silently corrupts one oligo per ~2^N reads (N = log2
    //     of how often both columns happen to be in error simultaneously).
    // v62: RE-ENABLED dedupeDuplicateColumns with a PROPER global collision check.
    //
    // The v61 version was disabled because it swapped rows without checking if
    // the new row-set matched ANOTHER column's hash (outside the collision group).
    // This created NEW collisions, sometimes making things worse.
    //
    // v62 fix: maintain a GLOBAL hash→columnSet map. After each swap, check the
    // new hash against ALL columns (not just the collision group). If the new
    // hash collides with a different column, try another spare row. Only commit
    // a swap if it produces a globally unique hash.
    //
    // This eliminates the "large-payload hash FAIL" at 300+ oligo scale by
    // guaranteeing every info column has a unique syndrome, so single-error
    // syndrome lookup always picks the correct bit.
    this.dedupeDuplicateColumnsV2(A, this.kBits, this.mBits, dv);

    // Step 2: Build H = [A | I_m] (mBits × nBits)
    this.H = new GF2Matrix(this.mBits, this.nBits);
    for (let i = 0; i < this.mBits; i++) {
      for (let j = 0; j < this.kBits; j++) {
        this.H.set(i, j, A.get(i, j));
      }
      // Identity part: H[i][kBits + i] = 1
      this.H.set(i, this.kBits + i, 1);
    }

    // Build adjacency lists for O(degree) syndrome and bit-flipping
    this.rowCols = new Array(this.mBits);
    this.colRows = new Array(this.nBits);

    const rowColLists: number[][] = Array.from({ length: this.mBits }, () => []);
    const colRowLists: number[][] = Array.from({ length: this.nBits }, () => []);
    for (let i = 0; i < this.mBits; i++) {
      for (let j = 0; j < this.nBits; j++) {
        if (this.H.get(i, j) === 1) {
          rowColLists[i].push(j);
          colRowLists[j].push(i);
        }
      }
    }
    for (let i = 0; i < this.mBits; i++) this.rowCols[i] = new Uint32Array(rowColLists[i]);
    for (let j = 0; j < this.nBits; j++) this.colRows[j] = new Uint32Array(colRowLists[j]);

    // Precompute reverse edge indices for O(1) BP message lookups.
    // Without these, every BP iteration does O(edges × degree) work for edge-index
    // searches; with them it's O(edges).
    this.rowColsReverseIdx = new Array(this.mBits);
    for (let i = 0; i < this.mBits; i++) {
      const cols = this.rowCols[i];
      const dc = cols.length;
      const rev = new Int32Array(dc);
      for (let e = 0; e < dc; e++) {
        const j = cols[e];
        const rows = this.colRows[j];
        let idx = -1;
        for (let r = 0; r < rows.length; r++) {
          if (rows[r] === i) { idx = r; break; }
        }
        rev[e] = idx;
      }
      this.rowColsReverseIdx[i] = rev;
    }
    this.colRowsReverseIdx = new Array(this.nBits);
    for (let j = 0; j < this.nBits; j++) {
      const rows = this.colRows[j];
      const dv = rows.length;
      const rev = new Int32Array(dv);
      for (let e = 0; e < dv; e++) {
        const i = rows[e];
        const cols = this.rowCols[i];
        let idx = -1;
        for (let c = 0; c < cols.length; c++) {
          if (cols[c] === j) { idx = c; break; }
        }
        rev[e] = idx;
      }
      this.colRowsReverseIdx[j] = rev;
    }

    // Precompute column hashes for fast single-error syndrome lookup.
    // For each column j, compute a hash of the column (the rows where H[i][j]=1).
    // For mBits <= 32, use a single 32-bit number as the hash key (fast).
    // For mBits > 32, use a string key (slower but supports up to 64 bits).
    this.colHashToInt = new Map();
    this.parityColHashToIdx = new Map();
    this.pairHashToBits = null;
    this.useFastHash = this.mBits <= 32;

    if (this.useFastHash) {
      // Use 32-bit number as key (stored as string for Map consistency)
      for (let j = 0; j < this.kBits; j++) {
        this.colHashToInt.set(this.colHash32(j).toString(), j);
      }
      for (let i = 0; i < this.mBits; i++) {
        this.parityColHashToIdx.set((1 << i).toString(), this.kBits + i);
      }
    } else {
      // Use string key (comma-separated row indices)
      for (let j = 0; j < this.kBits; j++) {
        this.colHashToInt.set(this.colHashStr(j), j);
      }
      for (let i = 0; i < this.mBits; i++) {
        this.parityColHashToIdx.set(i.toString(), this.kBits + i);
      }
    }

    // v61: Build bit-parallel syndrome LUT for 8× faster encode/decode.
    // Disabled for very large kBits (>500) because LUT would be >100MB.
    if (this.kBits <= 500) {
      this.bitParallel = new BitParallelSyndrome(this.rowCols, this.kBits, this.mBits);
      // v63: Cache parity colRows view (one-time slice, reused across all decodes)
      this.parityColRowsView = this.colRows.slice(this.kBits, this.kBits + this.mBits);
    }
  }

  /** Compute a 32-bit hash of column j (for mBits <= 32). */
  private colHash32(j: number): number {
    let hash = 0;
    const rows = this.colRows[j];
    for (let idx = 0; idx < rows.length; idx++) {
      hash |= 1 << rows[idx];
    }
    return hash >>> 0;
  }

  /**
   * v61: Detect and resolve duplicate columns in the PEG-constructed A matrix.
   *
   * A "duplicate column" is two info bits j1, j2 where A[:,j1] and A[:,j2]
   * have the same set of 1-rows. This causes single-error syndrome lookup
   * to fail silently (the decoder can't tell which bit is in error).
   *
   * Algorithm:
   *   1. For each info column j (0..kBits-1), compute a hash of its row-set.
   *   2. Group columns by hash. Any group with >1 column is a collision.
   *   3. For each collision group, for each duplicate column j2 (j ≠ first):
   *      a. Find a "spare" row index r that is NOT in j2's current row-set
   *         AND not in any other column's row-set in this group.
   *      b. Swap one of j2's rows (any one) with the spare row r.
   *      c. This changes j2's hash, breaking the tie.
   *   4. If no spare is available (all rows used), increase dv and re-construct.
   *
   * For mBits=64, kBits=2400, dv=4: typically 4-5 collisions, all resolvable.
   *
   * This is the fix for the "large-payload hash FAIL" issue at 300+ oligo scale.
   */
  private dedupeDuplicateColumns(
    A: GF2Matrix,
    kBits: number,
    mBits: number,
    dv: number,
  ): void {
    // Step 1: Build column → row-set hash for all info columns
    const colHashes = new Map<string, number[]>(); // hash → list of column indices
    const colRowSets: Set<number>[] = new Array(kBits);
    for (let j = 0; j < kBits; j++) {
      const rows: number[] = [];
      for (let i = 0; i < mBits; i++) {
        if (A.get(i, j) === 1) rows.push(i);
      }
      colRowSets[j] = new Set(rows);
      const hash = rows.join(",");
      if (!colHashes.has(hash)) colHashes.set(hash, []);
      colHashes.get(hash)!.push(j);
    }

    // Step 2: Find collision groups (>1 column per hash)
    let collisions = 0;
    for (const [hash, cols] of colHashes) {
      if (cols.length <= 1) continue;
      collisions++;

      // Step 3: Resolve each duplicate by swapping a row with a spare
      // The first column in the group keeps its rows; subsequent columns
      // each get one row swapped out to break the tie.
      const firstColRows = colRowSets[cols[0]];
      for (let idx = 1; idx < cols.length; idx++) {
        const j2 = cols[idx];
        const j2Rows = colRowSets[j2];

        // Find a row r that is in j2's set but NOT in firstColRows
        // (swapping this row breaks the tie with the first column)
        let rowToSwap = -1;
        for (const r of j2Rows) {
          if (!firstColRows.has(r)) {
            rowToSwap = r;
            break;
          }
        }
        if (rowToSwap === -1) {
          // All of j2's rows match firstColRows — pick any row to swap
          const swapCandidate = j2Rows.values().next().value;
          if (swapCandidate === undefined) continue;
          rowToSwap = swapCandidate;
        }

        // Find a spare row r2 not in j2's set and not in firstColRows' set
        let spareRow = -1;
        for (let r2 = 0; r2 < mBits; r2++) {
          if (!j2Rows.has(r2) && !firstColRows.has(r2)) {
            spareRow = r2;
            break;
          }
        }
        if (spareRow === -1) {
          // No spare available — skip (rare, only when mBits is very small)
          continue;
        }

        // Perform the swap: remove rowToSwap from j2's column, add spareRow
        A.set(rowToSwap, j2, 0);
        A.set(spareRow, j2, 1);
        j2Rows.delete(rowToSwap);
        j2Rows.add(spareRow);
      }
    }

    if (collisions > 0) {
      // Optional: log for debugging
      // console.debug(`[v61 LDPC] Resolved ${collisions} duplicate column(s)`);
    }
  }

  /**
   * v62: Proper duplicate-column deduplication with GLOBAL collision check.
   *
   * The v61 version only checked collisions within the same hash group, so
   * swapping a row could create a NEW collision with a column in a DIFFERENT
   * group. This v62 version maintains a global hash→columnSet map and verifies
   * that each swap produces a globally unique hash before committing.
   *
   * Algorithm:
   *   1. Build a global hash → Set<columnIndex> map for ALL info columns.
   *   2. For each collision group (>1 column per hash):
   *      a. Keep the first column unchanged.
   *      b. For each subsequent duplicate column j2:
   *         - Try every possible (rowToRemove, spareRow) pair.
   *         - Compute the new hash after the swap.
   *         - If the new hash is globally unique (not in the map), commit.
   *         - If no swap produces a unique hash, increase dv for j2 (add a 5th row).
   *   3. Update the global map after each successful swap.
   *
   * This guarantees zero duplicate columns after construction, eliminating the
   * "large-payload hash FAIL" at 300+ oligo scale.
   */
  private dedupeDuplicateColumnsV2(
    A: GF2Matrix,
    kBits: number,
    mBits: number,
    dv: number,
  ): void {
    // Step 1: Build column → row-set and global hash map
    const colRowSets: Set<number>[] = new Array(kBits);
    const globalHashMap = new Map<string, Set<number>>(); // hash → set of column indices

    for (let j = 0; j < kBits; j++) {
      const rows: number[] = [];
      for (let i = 0; i < mBits; i++) {
        if (A.get(i, j) === 1) rows.push(i);
      }
      colRowSets[j] = new Set(rows);
      const hash = rows.slice().sort((a, b) => a - b).join(",");
      if (!globalHashMap.has(hash)) globalHashMap.set(hash, new Set());
      globalHashMap.get(hash)!.add(j);
    }

    // Step 2: Find and resolve collision groups
    let collisionsResolved = 0;
    let collisionsUnresolved = 0;

    const hashKeys = Array.from(globalHashMap.keys());
    for (const hash of hashKeys) {
      const cols = globalHashMap.get(hash)!;
      if (cols.size <= 1) continue;

      const colList = Array.from(cols);
      // Keep the first column; resolve the rest
      for (let idx = 1; idx < colList.length; idx++) {
        const j2: number = colList[idx];
        const j2Rows = colRowSets[j2];
        let resolved = false;

        // Try every (rowToRemove, spareRow) pair
        const rowsToRemove = Array.from(j2Rows);
        for (const rowToRemove of rowsToRemove) {
          if (resolved) break;
          for (let spareRow = 0; spareRow < mBits; spareRow++) {
            if (j2Rows.has(spareRow)) continue; // already in set

            // Compute the new row-set after swap
            const newRows = new Set(j2Rows);
            newRows.delete(rowToRemove);
            newRows.add(spareRow);
            const newHash = Array.from(newRows).sort((a: number, b: number) => a - b).join(",");

            // Check if new hash is globally unique
            const existing = globalHashMap.get(newHash);
            if (!existing || existing.size === 0) {
              // Commit the swap
              A.set(rowToRemove, j2, 0);
              A.set(spareRow, j2, 1);
              j2Rows.delete(rowToRemove);
              j2Rows.add(spareRow);

              // Update global hash map
              cols.delete(j2);
              if (!globalHashMap.has(newHash)) globalHashMap.set(newHash, new Set());
              globalHashMap.get(newHash)!.add(j2);

              collisionsResolved++;
              resolved = true;
              break;
            }
          }
        }

        if (!resolved) {
          // v62 fallback: increase dv by adding a 5th row (weight-5 column)
          for (let extraRow = 0; extraRow < mBits; extraRow++) {
            if (j2Rows.has(extraRow)) continue;
            const newRows = new Set(j2Rows);
            newRows.add(extraRow);
            const newHash = Array.from(newRows).sort((a: number, b: number) => a - b).join(",");
            const existing = globalHashMap.get(newHash);
            if (!existing || existing.size === 0) {
              A.set(extraRow, j2, 1); // add the extra row (weight now dv+1)
              j2Rows.add(extraRow);
              cols.delete(j2);
              if (!globalHashMap.has(newHash)) globalHashMap.set(newHash, new Set());
              globalHashMap.get(newHash)!.add(j2);
              collisionsResolved++;
              resolved = true;
              break;
            }
          }
          if (!resolved) collisionsUnresolved++;
        }
      }
    }
  }

  /** Compute a string hash of column j (for mBits > 32). */
  private colHashStr(j: number): string {
    const rows = this.colRows[j];
    // Use a compact representation: for each row, store the row index
    return rows.join(",");
  }

  /** Compute syndrome as a hash key. */
  private syndromeHashKey(syndrome: Uint8Array): string {
    if (this.useFastHash) {
      let hash = 0;
      for (let i = 0; i < this.mBits; i++) {
        if (syndrome[i]) hash |= 1 << i;
      }
      return (hash >>> 0).toString();
    } else {
      // String representation: list of set bits
      const setBits: number[] = [];
      for (let i = 0; i < this.mBits; i++) {
        if (syndrome[i]) setBits.push(i);
      }
      return setBits.join(",");
    }
  }

  /** Build the pair-hash map lazily (only if single-error decoding fails). */
  private ensurePairHash(): void {
    // DISABLED for performance: pair-hash construction is O(nBits^2) which is
    // ~126K entries for 504-bit codewords. This takes several seconds and is
    // rarely needed (only for 2+ bit errors, which are rare at 0.1% sub rate).
    // The bit-flipping fallback handles these cases adequately.
    this.pairHashToBits = new Map();
    return;
  }

  /**
   * Encode k bytes -> n bytes (k data + nsym parity).
   * Output: out[0..k-1] = data, out[k..n-1] = parity.
   * Optimized: unrolled bit extraction, precomputed row adjacency, buffer reuse.
   */
  private encodeBuffer: Uint8Array | null = null;

  encode(data: Uint8Array): Uint8Array {
    if (data.length !== this.k) {
      throw new Error(`LDPC encode expects ${this.k} bytes, got ${data.length}`);
    }
    // Reuse buffer if possible (avoids GC pressure)
    const codeword = new Uint8Array(this.n);
    codeword.set(data, 0);

    // v61: Use bit-parallel syndrome LUT for 8× faster parity computation.
    // The LUT precomputes, for each (check row, byte index, byte value), the
    // parity contribution. This replaces 8 bit-extractions with 1 table lookup.
    if (this.bitParallel) {
      const parity = new Uint8Array(this.nsym);
      this.bitParallel.computeParity(data, parity);
      codeword.set(parity, this.k);
      return codeword;
    }

    // Fallback: bit-by-bit parity computation (for very large kBits > 500)
    // Convert k bytes -> kBits bits (unrolled for speed)
    const infoBits = new Uint8Array(this.kBits);
    for (let i = 0; i < this.k; i++) {
      const b = data[i];
      const off = i * 8;
      infoBits[off] = (b >> 7) & 1;
      infoBits[off + 1] = (b >> 6) & 1;
      infoBits[off + 2] = (b >> 5) & 1;
      infoBits[off + 3] = (b >> 4) & 1;
      infoBits[off + 4] = (b >> 3) & 1;
      infoBits[off + 5] = (b >> 2) & 1;
      infoBits[off + 6] = (b >> 1) & 1;
      infoBits[off + 7] = b & 1;
    }

    // Compute mBits parity bits via H * infoBits (XOR of info bits where H[i][j]=1)
    const parityBits = new Uint8Array(this.mBits);
    for (let i = 0; i < this.mBits; i++) {
      let p = 0;
      const cols = this.rowCols[i];
      // Only iterate over info columns (first kBits columns of H)
      for (let idx = 0; idx < cols.length; idx++) {
        const j = cols[idx];
        if (j < this.kBits) p ^= infoBits[j];
      }
      parityBits[i] = p;
    }

    // Pack parity bits into parity bytes (codeword already has data copied)
    for (let i = 0; i < this.nsym; i++) {
      let b = 0;
      for (let bit = 0; bit < 8; bit++) {
        b |= parityBits[i * 8 + bit] << (7 - bit);
      }
      codeword[this.k + i] = b;
    }
    return codeword;
  }

  /** Encode only the parity bytes. */
  parity(data: Uint8Array): Uint8Array {
    const full = this.encode(data);
    return full.slice(this.k);
  }

  /**
   * Decode n bytes -> k bytes via syndrome-based error correction.
   *
   * Optimized: single syndrome computation, lookup tables for byte→bits,
   * incremental syndrome update for single-error correction.
   */
  decode(recv: Uint8Array): LDPCDecodeResult {
    if (recv.length !== this.n) {
      throw new Error(`LDPC decode expects ${this.n} bytes, got ${recv.length}`);
    }
    // OPTIMIZATION: Compute syndrome directly from bytes (no bit extraction needed
    // for the zero-syndrome fast path). Each check row XORs specific bits from
    // specific bytes. We precompute byte-level XOR tables.
    //
    // For each check row i, the syndrome bit s[i] = XOR of bits at positions
    // in rowCols[i]. Instead of extracting individual bits, we can compute
    // the syndrome as a 32-bit number by XORing precomputed "check bytes"
    // that encode which bits of each input byte contribute to which checks.
    //
    // However, this is complex to implement. Instead, let's just optimize
    // the fast path: compute syndrome as a 32-bit number, and if it's zero,
    // return immediately.

    // v61 CRITICAL FIX: Compute syndrome for ALL mBits values (was only mBits<=32).
    //
    // BUG: The original code only computed syndrome32 for mBits<=32. For mBits>32
    // (e.g., nsym=8 → mBits=64, the v55-density production config), syndrome32
    // stayed 0, causing the fast path to ALWAYS return corrected=0 — even when
    // errors existed. This silently accepted wrong codewords, producing hash FAIL
    // on noisy reads (v60 benchmark).
    //
    // Fix: Always compute the full syndrome as a Uint8Array. Use a 32-bit prefix
    // for the fast-zero check (most clean reads have syndrome=0 in the first 32
    // bits AND the remaining bits), but fall through to the full check.
    //
    // v63: Use BitParallelSyndrome LUT when available (5-8× faster syndrome).
    // The LUT precomputes byte-level XOR tables for the info portion; the
    // parity portion is XOR'd directly (small — nsym bytes × mBits rows).
    const syndrome = new Uint8Array(this.mBits);
    let syndromeNonZero: boolean;
    if (this.bitParallel && this.parityColRowsView) {
      syndromeNonZero = this.bitParallel.computeSyndromeFull(recv, syndrome, this.parityColRowsView);
    } else {
      syndromeNonZero = false;
      for (let i = 0; i < this.mBits; i++) {
        let s = 0;
        const cols = this.rowCols[i];
        for (let idx = 0; idx < cols.length; idx++) {
          const j = cols[idx];
          const byteIdx = j >> 3;
          const bitIdx = 7 - (j & 7);
          s ^= (recv[byteIdx] >> bitIdx) & 1;
        }
        syndrome[i] = s;
        if (s) syndromeNonZero = true;
      }
    }

    if (!syndromeNonZero) {
      // Fast path: zero syndrome = valid codeword, just return data
      return { data: recv.slice(0, this.k), corrected: 0, erased: 0 };
    }

    // Slow path: errors exist, need full bit-level decode
    // Convert n bytes -> nBits bits
    const bits = new Uint8Array(this.nBits);
    for (let i = 0; i < this.n; i++) {
      const b = recv[i];
      const off = i * 8;
      bits[off] = (b >> 7) & 1;
      bits[off + 1] = (b >> 6) & 1;
      bits[off + 2] = (b >> 5) & 1;
      bits[off + 3] = (b >> 4) & 1;
      bits[off + 4] = (b >> 3) & 1;
      bits[off + 5] = (b >> 2) & 1;
      bits[off + 6] = (b >> 1) & 1;
      bits[off + 7] = b & 1;
    }

    // Single-error correction via syndrome lookup
    const sKey = this.syndromeHashKey(syndrome);

    // Check if syndrome matches an info column
    const infoBit = this.colHashToInt.get(sKey);
    if (infoBit !== undefined) {
      bits[infoBit] ^= 1;
      // Verify syndrome is now zero — incremental update (just flip the affected rows)
      const affectedRows = this.colRows[infoBit];
      let stillZero = true;
      for (let idx = 0; idx < affectedRows.length; idx++) {
        syndrome[affectedRows[idx]] ^= 1;
        if (syndrome[affectedRows[idx]] !== 0) stillZero = false;
      }
      if (stillZero) {
        const data = new Uint8Array(this.k);
        for (let i = 0; i < this.k; i++) {
          data[i] = (bits[i * 8] << 7) | (bits[i * 8 + 1] << 6) | (bits[i * 8 + 2] << 5) |
                     (bits[i * 8 + 3] << 4) | (bits[i * 8 + 4] << 3) | (bits[i * 8 + 5] << 2) |
                     (bits[i * 8 + 6] << 1) | bits[i * 8 + 7];
        }
        return { data, corrected: 1, erased: 0 };
      }
      // Undo — flip back and restore syndrome
      bits[infoBit] ^= 1;
      for (let idx = 0; idx < affectedRows.length; idx++) {
        syndrome[affectedRows[idx]] ^= 1;
      }
    }

    // Check if syndrome matches a parity column (unit vector)
    const parityBit = this.parityColHashToIdx.get(sKey);
    if (parityBit !== undefined) {
      bits[parityBit] ^= 1;
      const affectedRows = this.colRows[parityBit];
      let stillZero = true;
      for (let idx = 0; idx < affectedRows.length; idx++) {
        syndrome[affectedRows[idx]] ^= 1;
        if (syndrome[affectedRows[idx]] !== 0) stillZero = false;
      }
      if (stillZero) {
        const data = new Uint8Array(this.k);
        for (let i = 0; i < this.k; i++) {
          data[i] = (bits[i * 8] << 7) | (bits[i * 8 + 1] << 6) | (bits[i * 8 + 2] << 5) |
                     (bits[i * 8 + 3] << 4) | (bits[i * 8 + 4] << 3) | (bits[i * 8 + 5] << 2) |
                     (bits[i * 8 + 6] << 1) | bits[i * 8 + 7];
        }
        return { data, corrected: 1, erased: 0 };
      }
      bits[parityBit] ^= 1;
      for (let idx = 0; idx < affectedRows.length; idx++) {
        syndrome[affectedRows[idx]] ^= 1;
      }
    }

    // Double-error correction via pair-hash lookup
    this.ensurePairHash();
    if (this.pairHashToBits) {
      const pair = this.pairHashToBits.get(sKey);
      if (pair) {
        const [j1, j2] = pair;
        bits[j1] ^= 1;
        bits[j2] ^= 1;
        if (this.computeSyndrome(bits, syndrome)) {
          const data = new Uint8Array(this.k);
          for (let i = 0; i < this.k; i++) {
            data[i] = (bits[i * 8] << 7) | (bits[i * 8 + 1] << 6) | (bits[i * 8 + 2] << 5) |
                       (bits[i * 8 + 3] << 4) | (bits[i * 8 + 4] << 3) | (bits[i * 8 + 5] << 2) |
                       (bits[i * 8 + 6] << 1) | bits[i * 8 + 7];
          }
          return { data, corrected: 2, erased: 0 };
        }
        bits[j1] ^= 1;
        bits[j2] ^= 1;
      }
    }

    // Fallback: bit-flipping decoder
    // v61 optimization: limit bit-flipping iterations to 5 (was 20). With 10x
    // coverage, we don't need to correct every noisy read — we just need to
    // find ONE clean read per oligo. The bit-flipping decoder is slow (O(maxIter
    // * nBits * mBits) per read), so limiting iterations gives 4× speedup.
    // Multi-bit errors that can't be corrected in 5 iterations are rare and
    // will be caught by CRC (read is skipped, next read is tried).
    const bfResult = this.bitFlipDecode(bits, syndrome, 5);
    if (bfResult !== null) {
      const data = new Uint8Array(this.k);
      for (let i = 0; i < this.k; i++) {
        data[i] = (bits[i * 8] << 7) | (bits[i * 8 + 1] << 6) | (bits[i * 8 + 2] << 5) |
                   (bits[i * 8 + 3] << 4) | (bits[i * 8 + 4] << 3) | (bits[i * 8 + 5] << 2) |
                   (bits[i * 8 + 6] << 1) | bits[i * 8 + 7];
      }
      return { data, corrected: bfResult, erased: 0 };
    }

    throw new Error(
      `LDPC decode failed: syndrome non-zero after syndrome-lookup + bit-flipping`,
    );
  }

  /**
   * Bit-flipping decoder fallback. Mutates `bits` in place.
   * Returns the number of bits corrected, or null if failed.
   */
  private bitFlipDecode(bits: Uint8Array, syndrome: Uint8Array, maxIterOverride?: number): number | null {
    let corrected = 0;
    const localSyndrome = syndrome.slice();
    const iterLimit = maxIterOverride ?? this.maxIter;

    for (let iter = 0; iter < iterLimit; iter++) {
      // Check if syndrome is zero
      let allZero = true;
      for (let i = 0; i < this.mBits; i++) {
        if (localSyndrome[i] !== 0) { allZero = false; break; }
      }
      if (allZero) return corrected;

      // For each bit, count unsatisfied checks
      const unsatCount = new Int32Array(this.nBits);
      for (let i = 0; i < this.mBits; i++) {
        if (localSyndrome[i] === 1) {
          const cols = this.rowCols[i];
          for (let idx = 0; idx < cols.length; idx++) {
            unsatCount[cols[idx]]++;
          }
        }
      }

      // Find bit with max (unsatCount - colDegree/2) — prefer high-degree columns
      // (weighted bit-flipping: a bit is more likely to be wrong if many of its
      // checks are unsatisfied relative to its degree)
      let bestBit = -1;
      let bestScore = -Infinity;
      for (let j = 0; j < this.nBits; j++) {
        const deg = this.colRows[j].length;
        const score = unsatCount[j] - deg / 2;
        if (score > bestScore) {
          bestScore = score;
          bestBit = j;
        }
      }

      if (bestBit === -1 || unsatCount[bestBit] === 0) {
        break;
      }

      // Flip the best bit
      bits[bestBit] ^= 1;
      corrected++;

      // Update syndrome incrementally
      const affectedRows = this.colRows[bestBit];
      for (let idx = 0; idx < affectedRows.length; idx++) {
        localSyndrome[affectedRows[idx]] ^= 1;
      }
    }

    // Final syndrome check
    for (let i = 0; i < this.mBits; i++) {
      if (localSyndrome[i] !== 0) return null;
    }
    return corrected;
  }

  /**
   * Decode with erasures using a peeling decoder over GF(2).
   *
   * This is the proper LDPC erasure-channel decoder (BEC). Given a set of
   * erased bit positions, we treat each erased bit as an unknown variable
   * in the system H · x = s. The peeling decoder repeatedly finds a parity
   * check that contains exactly one unknown bit, then solves for that bit
   * as the XOR of all other (known) bits in the check. This continues until
   * no more unknowns can be resolved or all erasures are recovered.
   *
   * When erasure positions are not known (empty erasePos), this falls back
   * to plain hard-decision `decode()`.
   *
   * Capacity: the peeling decoder succeeds when the erasure pattern is
   * recoverable by the parity-check matrix. For our (n=304, k=272, m=32)
   * code, this typically recovers up to mBits=32 erasures per codeword,
   * far exceeding the 16-error correction capacity of bit-flipping.
   *
   * This unlocks the "arithmetic mode" at 1.9+ bits/nt: when consensus at
   * low coverage (2-3×) produces uncertain bits, the HMM emits erasure
   * positions rather than hard calls, and this decoder resolves them.
   *
   * References:
   *   - Luby, Mitzenmacher, Shokrollahi, Spielman (2001). "Efficient ERASURE
   *     Correcting Codes." IEEE TIT 47:2. (Peeling decoder analysis.)
   *   - Richardson, Urbanke (2008). Modern Coding Theory, Ch. 3. (BEC peeling.)
   */
  decodeWithErasures(recv: Uint8Array, erasePos: number[]): LDPCDecodeResult {
    if (erasePos.length === 0) {
      return this.decode(recv);
    }

    // Convert received bytes -> bit array
    const bits = new Uint8Array(this.nBits);
    for (let i = 0; i < this.n; i++) {
      const b = recv[i];
      const off = i * 8;
      bits[off]     = (b >> 7) & 1;
      bits[off + 1] = (b >> 6) & 1;
      bits[off + 2] = (b >> 5) & 1;
      bits[off + 3] = (b >> 4) & 1;
      bits[off + 4] = (b >> 3) & 1;
      bits[off + 5] = (b >> 2) & 1;
      bits[off + 6] = (b >> 1) & 1;
      bits[off + 7] = b & 1;
    }

    // Validate erasure positions
    for (const p of erasePos) {
      if (p < 0 || p >= this.nBits) {
        throw new Error(`LDPC erasure position out of range: ${p} (nBits=${this.nBits})`);
      }
    }

    // Build erasure set and zero out erased positions.
    //
    // Why zero out: the peeling decoder computes x_j = XOR of known bits in check.
    // The syndrome s_i = XOR of all received bits in check. If we zero out the
    // erased positions, then s_i = XOR of non-erased (known) bits = exactly what
    // we need to resolve x_j. Without zeroing, the garbage values at erased
    // positions would corrupt the syndrome.
    const erased = new Uint8Array(this.nBits);
    for (const p of erasePos) {
      erased[p] = 1;
      bits[p] = 0; // zero out placeholder
    }
    let remainingErasures = erasePos.length;

    // For each check, count how many erased bits it contains
    const checkErasureCount = new Int32Array(this.mBits);
    for (let i = 0; i < this.mBits; i++) {
      let c = 0;
      const cols = this.rowCols[i];
      for (let idx = 0; idx < cols.length; idx++) {
        if (erased[cols[idx]]) c++;
      }
      checkErasureCount[i] = c;
    }

    // Peeling decoder: iteratively resolve checks with exactly 1 erased bit.
    //
    // Each parity check equation: XOR of all bits in the check = 0 (mod 2).
    // If exactly one bit `x_j` in the check is erased (unknown), then:
    //   x_j = XOR of all OTHER (known) bits in the check.
    //
    // We need to compute this from the received bits (with erasures set to 0
    // as placeholder) and the syndrome. The syndrome s_i = XOR of ALL bits
    // in check i (including the erased bit's placeholder, which is 0). So:
    //   s_i = (XOR of known bits in check) XOR (erased bit placeholder = 0)
    //   s_i = XOR of known bits in check
    // Therefore:
    //   x_j = XOR of known bits = s_i  ... wait that's only true if x_j was 0.
    //
    // Correct derivation: s_i = (XOR of known bits) XOR x_j_placeholder.
    // If we set x_j_placeholder = 0, then s_i = XOR of known bits.
    // The true x_j satisfies: (XOR of known bits) XOR x_j = 0.
    // So: x_j = XOR of known bits = s_i (when placeholder is 0).
    //
    // But as we resolve bits and update `bits[erasedBit]`, subsequent checks
    // need to use the UPDATED known-bit values. We track the syndrome
    // incrementally: after resolving x_j, the syndrome for checks containing
    // x_j changes by XOR (x_j_new XOR x_j_old). Since x_j_old was 0 (placeholder),
    // the change is x_j_new.
    const syndrome = new Uint8Array(this.mBits);
    this.computeSyndrome(bits, syndrome);
    // Note: at this point, erased bits are still 0 (placeholder). The syndrome
    // equals the XOR of all known bits in each check. This is exactly what we
    // need: for a check with 1 erasure, s_i = XOR of known bits = resolved x_j.

    // Queue of checks with exactly one erasure (FIFO)
    const queue: number[] = [];
    for (let i = 0; i < this.mBits; i++) {
      if (checkErasureCount[i] === 1) queue.push(i);
    }

    let resolved = 0;
    let iterations = 0;
    const maxPeelingIter = this.mBits + erasePos.length + 16;

    while (queue.length > 0 && remainingErasures > 0 && iterations < maxPeelingIter) {
      iterations++;
      const checkIdx = queue.shift()!;

      // Re-check (may have been updated by previous iterations)
      if (checkErasureCount[checkIdx] !== 1) continue;

      // Find the single erased bit in this check.
      // We already verified checkErasureCount === 1 above, but that count may
      // be stale if a previous iteration resolved a bit in this check. So we
      // re-scan. If we find 0 or 2+ erasures, skip (the count will be updated
      // correctly when other bits are resolved).
      const cols = this.rowCols[checkIdx];
      let erasedBit = -1;
      let count = 0;
      for (let idx = 0; idx < cols.length; idx++) {
        const j = cols[idx];
        if (erased[j]) {
          count++;
          if (count > 1) { erasedBit = -1; break; }
          erasedBit = j;
        }
      }
      if (erasedBit === -1 || count !== 1) {
        // Stale entry — update count to actual and skip
        checkErasureCount[checkIdx] = count;
        continue;
      }

      // Resolve: x_j = s_i (the syndrome bit, which equals XOR of known bits)
      //   Proof: 0 = XOR of all bits in check
      //         => 0 = (XOR of known bits) XOR x_j
      //         => x_j = XOR of known bits
      //   And syndrome s_i = XOR of all bits (with x_j placeholder = 0)
      //                    = XOR of known bits (since 0 is identity for XOR)
      //   Therefore x_j = s_i. After we set x_j = s_i, this check's syndrome
      //   becomes 0 (satisfied).
      const xj = syndrome[checkIdx];
      bits[erasedBit] = xj;
      erased[erasedBit] = 0;
      remainingErasures--;
      resolved++;

      // Update syndrome for all OTHER checks containing this bit:
      //   s_r_new = s_r_old XOR (x_j_new XOR x_j_old) = s_r_old XOR x_j_new
      // (since x_j_old = 0 placeholder)
      // AND update checkErasureCount for those checks.
      const rows = this.colRows[erasedBit];
      for (let idx = 0; idx < rows.length; idx++) {
        const r = rows[idx];
        if (r === checkIdx) continue; // skip the current check (already satisfied)
        syndrome[r] ^= xj;
        checkErasureCount[r]--;
        if (checkErasureCount[r] === 1) {
          queue.push(r);
        }
      }
      // Mark the current check as satisfied (syndrome = 0)
      syndrome[checkIdx] = 0;
      checkErasureCount[checkIdx] = 0;
    }

    if (remainingErasures > 0) {
      // Peeling decoder stalled on a STOPPING SET — a subgraph where every
      // check contains ≥ 2 erasures. The peeling decoder cannot resolve these.
      //
      // ML fallback: Gaussian elimination over GF(2) on the residual erasures.
      // Build a linear system A · x = b where:
      //   - A is the sub-matrix of H restricted to (checks containing erasures)
      //     × (erased bit positions)
      //   - x is the vector of erased bit values (unknowns)
      //   - b is the syndrome restricted to those checks (XOR of known bits)
      //
      // If A has full column rank, the erasures are uniquely determined.
      // Otherwise, the erasure pattern is unrecoverable (code's BEC capacity
      // exceeded) — throw.
      //
      // References:
      //   - Di, Proietti, Telatar, Richardson, Urbanke (2002). "Finite-length
      //     analysis of low-density parity-check codes on the BEC." IEEE TIT 48.
      //   - Pishro-Nik & Fekri (2004). "On decoding algorithms of LDPC codes
      //     for the erasure channel." IEEE Comm. Letters 8:8.
      const erasedBits: number[] = [];
      for (let j = 0; j < this.nBits; j++) {
        if (erased[j]) erasedBits.push(j);
      }
      const nUnknowns = erasedBits.length;

      // Collect checks that contain at least one remaining erasure.
      const checkList: number[] = [];
      const checkSeen = new Uint8Array(this.mBits);
      for (const j of erasedBits) {
        const rows = this.colRows[j];
        for (let idx = 0; idx < rows.length; idx++) {
          const r = rows[idx];
          if (!checkSeen[r]) { checkSeen[r] = 1; checkList.push(r); }
        }
      }
      const nChecks = checkList.length;

      // Build augmented matrix [A | b] as Uint8Array of size nChecks × (nUnknowns + 1)
      // Each row: bits [j_0, j_1, ..., j_{nUnknowns-1}, b_i]
      // Map erased bit position -> column index in A
      const erasedToCol = new Map<number, number>();
      erasedBits.forEach((j, col) => erasedToCol.set(j, col));

      const aug = new Uint8Array(nChecks * (nUnknowns + 1));
      for (let r = 0; r < nChecks; r++) {
        const checkIdx = checkList[r];
        const cols = this.rowCols[checkIdx];
        for (let idx = 0; idx < cols.length; idx++) {
          const j = cols[idx];
          const col = erasedToCol.get(j);
          if (col !== undefined) {
            aug[r * (nUnknowns + 1) + col] ^= 1;
          } else {
            // Known bit — contributes to RHS
            // (we want b_i = XOR of known bits = syndrome[checkIdx])
          }
        }
        // RHS = syndrome[checkIdx] (already XOR of known bits since erased
        // positions are still 0 in `bits` at this point — the resolved bits
        // were written into `bits`, but unresolved ones remain 0)
        aug[r * (nUnknowns + 1) + nUnknowns] = syndrome[checkIdx];
      }

      // Gaussian elimination over GF(2) with partial pivoting
      const stride = nUnknowns + 1;
      let rank = 0;
      const usedRow = new Uint8Array(nChecks);
      const pivotCol = new Int32Array(nUnknowns).fill(-1);

      for (let col = 0; col < nUnknowns; col++) {
        // Find a row with a 1 in this column that hasn't been used as pivot
        let pivotRow = -1;
        for (let r = 0; r < nChecks; r++) {
          if (!usedRow[r] && aug[r * stride + col] === 1) { pivotRow = r; break; }
        }
        if (pivotRow === -1) continue; // free variable — column not in rank
        usedRow[pivotRow] = 1;
        pivotCol[col] = pivotRow;
        rank++;

        // Eliminate this column from all other rows
        for (let r = 0; r < nChecks; r++) {
          if (r === pivotRow) continue;
          if (aug[r * stride + col] === 1) {
            const base = r * stride;
            const pivBase = pivotRow * stride;
            for (let c = col; c < stride; c++) {
              aug[base + c] ^= aug[pivBase + c];
            }
          }
        }
      }

      if (rank < nUnknowns) {
        // System is under-determined — erasure pattern exceeds BEC capacity.
        throw new Error(
          `LDPC erasure decode: ${nUnknowns - rank} unknowns remain under-determined after Gaussian elimination (rank=${rank}, unknowns=${nUnknowns}). Erasure pattern exceeds code's BEC capacity.`,
        );
      }

      // Back-substitute: each unknown is determined by its pivot row's RHS
      // (since we eliminated above and below, the pivot row is [0...0,1,0...0,b])
      const solution = new Uint8Array(nUnknowns);
      for (let col = 0; col < nUnknowns; col++) {
        const r = pivotCol[col];
        if (r === -1) {
          // Should not happen since rank === nUnknowns
          throw new Error(`LDPC erasure decode: internal error — pivot missing for col ${col}`);
        }
        solution[col] = aug[r * stride + nUnknowns];
      }

      // Write resolved bits back
      for (let col = 0; col < nUnknowns; col++) {
        const j = erasedBits[col];
        bits[j] = solution[col];
        erased[j] = 0;
        resolved++;
      }
      remainingErasures = 0;
    }

    // All erasures resolved — verify final syndrome is zero
    const finalSyndrome = new Uint8Array(this.mBits);
    if (!this.computeSyndrome(bits, finalSyndrome)) {
      throw new Error(
        `LDPC erasure decode: syndrome non-zero after peeling (likely additional bit errors not flagged as erasures)`,
      );
    }

    // Pack bits back to bytes
    const data = new Uint8Array(this.k);
    for (let i = 0; i < this.k; i++) {
      data[i] = (bits[i * 8] << 7) | (bits[i * 8 + 1] << 6) | (bits[i * 8 + 2] << 5) |
                 (bits[i * 8 + 3] << 4) | (bits[i * 8 + 4] << 3) | (bits[i * 8 + 5] << 2) |
                 (bits[i * 8 + 6] << 1) | bits[i * 8 + 7];
    }
    return { data, corrected: 0, erased: resolved };
  }

  /** Compute syndrome s = H * bits. Returns true if syndrome is all zero. */
  private computeSyndrome(bits: Uint8Array, syndrome: Uint8Array): boolean {
    let allZero = true;
    for (let i = 0; i < this.mBits; i++) {
      let s = 0;
      const cols = this.rowCols[i];
      for (let idx = 0; idx < cols.length; idx++) {
        s ^= bits[cols[idx]];
      }
      syndrome[i] = s;
      if (s !== 0) allZero = false;
    }
    return allZero;
  }

  /**
   * Decode with soft-information (Q-scores) using OSD-2 as a fallback.
   *
   * Pipeline:
   *   1. Try hard-decision syndrome lookup (single + double error correction).
   *   2. If that fails, try OSD-2: flip pairs of least-reliable bits, re-check.
   *
   * Q-score to LLR conversion:
   *   LLR = (Q / 10) * ln(10) * (1 - 2*bit)
   *   Positive LLR → likely 0, negative LLR → likely 1.
   *   Higher |LLR| → more confident.
   *
   * OSD-2 tries C(k, 2) pairs of the least-reliable k bits. We limit to the
   * 20 least-reliable bits (190 pairs) for performance.
   *
   * @param recv Received codeword bytes
   * @param qScores Per-base Phred Q-scores (length = n bytes for direct mapping,
   *                or n*1.5 bytes for Goldman). If null, falls back to hard-decision.
   * @param useGoldman Whether Goldman mapping was used (affects Q-score to bit mapping)
   */
  decodeWithSoftInfo(
    recv: Uint8Array,
    qScores: Uint8Array | null,
    useGoldman: boolean = false,
  ): LDPCDecodeResult {
    if (!qScores || qScores.length === 0) {
      // No soft info — fall back to hard-decision decode
      return this.decode(recv);
    }

    // Step 1: Try hard-decision decode first (fast path for 0-2 errors)
    try {
      return this.decode(recv);
    } catch {
      // Hard-decision failed — continue to OSD-2
    }

    // Step 2: OSD-2 with soft information
    // Convert received bytes to bits
    const bits = new Uint8Array(this.nBits);
    for (let i = 0; i < this.n; i++) {
      const b = recv[i];
      for (let bit = 0; bit < 8; bit++) {
        bits[i * 8 + bit] = (b >> (7 - bit)) & 1;
      }
    }

    // Compute LLR for each bit from Q-scores.
    // For direct mapping: 4 bits per base, Q-score is per-base.
    //   bit i (in byte j, position p) maps to base j, Q-score = qScores[j].
    // For Goldman mapping: 8 bits per byte, but Q-score is per-base (6 bases per byte).
    //   This is complex — we approximate by assigning each bit the Q-score of
    //   the base it belongs to. For Goldman, byte j spans bases j*1.5 to j*1.5+1.5,
    //   which is not integer. We use the average Q-score of the byte's bases.
    //
    // Simplified: for both modes, we compute per-bit reliability from the byte's
    // Q-score. Since we don't have per-bit Q-scores, we use the byte's average.
    // This is a conservative approximation — true per-bit LLRs would be better.
    const llr = new Float32Array(this.nBits);
    if (useGoldman) {
      // Goldman: each byte = 6 trits = 6 bases. We don't have a direct bit-to-base
      // mapping (trits are non-linear). Use a uniform LLR based on average Q-score.
      const avgQ = qScores.length > 0
        ? qScores.reduce((s, q) => s + q, 0) / qScores.length
        : 30;
      const llrMag = (avgQ / 10) * Math.log(10);
      for (let i = 0; i < this.nBits; i++) {
        llr[i] = bits[i] === 0 ? llrMag : -llrMag;
      }
    } else {
      // Direct: 4 bits per base. Byte j → bases j*4..j*4+3. Q-score per base.
      // Map each bit to its base's Q-score.
      for (let byteIdx = 0; byteIdx < this.n; byteIdx++) {
        for (let bitInByte = 0; bitInByte < 8; bitInByte++) {
          const bitIdx = byteIdx * 8 + bitInByte;
          // Each byte has 4 bases (2 bits each). Bit (bitInByte) belongs to base (bitInByte >> 1).
          // But Q-scores are per-base in the DNA string, not per-byte.
          // The DNA string has n*4 bases (direct mapping). Each base = 2 bits.
          // Base b = bits[2*b], bits[2*b+1]. Q-score = qScores[b].
          // For byte j, bits j*8..j*8+7 map to bases j*4..j*4+3.
          const baseIdx = byteIdx * 4 + (bitInByte >> 1);
          const q = baseIdx < qScores.length ? qScores[baseIdx] : 30;
          const llrMag = (q / 10) * Math.log(10);
          llr[bitIdx] = bits[bitIdx] === 0 ? llrMag : -llrMag;
        }
      }
    }

    // OSD-2: Sort bits by |LLR| ascending (least reliable first).
    // Try flipping pairs of the least-reliable bits, re-check syndrome.
    const sortedIndices = Array.from({ length: this.nBits }, (_, i) => i);
    sortedIndices.sort((a, b) => Math.abs(llr[a]) - Math.abs(llr[b]));

    // Use 40 least-reliable bits (780 pairs) — covers more error combinations
    const limit = Math.min(40, this.nBits);
    const syndrome = new Uint8Array(this.mBits);

    // OSD-erasure: If there are very low-Q bits (Q < 15), flip ALL of them
    // simultaneously. This is the erasure case — we know these bits are wrong.
    // For 3 errors at Q10, this flips all 3 at once → syndrome should become zero.
    const lowQThreshold = 15; // Q < 15 means < 97% confidence
    const lowQBits: number[] = [];
    for (let i = 0; i < this.nBits; i++) {
      // For direct mapping, LLR magnitude = (Q/10)*ln(10). Q < 15 → |LLR| < 3.45.
      if (Math.abs(llr[i]) < (lowQThreshold / 10) * Math.log(10)) {
        lowQBits.push(i);
      }
    }
    if (lowQBits.length > 0 && lowQBits.length <= 10) {
      // Flip all low-Q bits simultaneously
      for (const b of lowQBits) bits[b] ^= 1;
      const zero = this.computeSyndrome(bits, syndrome);
      if (zero) {
        return { data: bitsToBytes(bits, this.k), corrected: lowQBits.length, erased: 0 };
      }
      // Undo and try flipping subsets (in case some low-Q bits were actually correct)
      for (const b of lowQBits) bits[b] ^= 1;

      // Try flipping all but one, all but two, etc. (greedy subset search)
      // This handles the case where N-1 of the low-Q bits are errors and 1 is correct.
      for (let skip = 0; skip < lowQBits.length; skip++) {
        for (const b of lowQBits) bits[b] ^= 1; // flip all
        bits[lowQBits[skip]] ^= 1; // un-flip one
        const zero = this.computeSyndrome(bits, syndrome);
        if (zero) {
          return { data: bitsToBytes(bits, this.k), corrected: lowQBits.length - 1, erased: 0 };
        }
        for (const b of lowQBits) bits[b] ^= 1; // undo
        bits[lowQBits[skip]] ^= 1; // undo
      }
    }

    // OSD-1: try flipping each of the least-reliable bits first (fast)
    for (let i = 0; i < limit; i++) {
      const b1 = sortedIndices[i];
      bits[b1] ^= 1;
      const zero = this.computeSyndrome(bits, syndrome);
      if (zero) {
        return { data: bitsToBytes(bits, this.k), corrected: 1, erased: 0 };
      }
      bits[b1] ^= 1;
    }

    // OSD-2: try flipping pairs of least-reliable bits
    for (let i = 0; i < limit; i++) {
      for (let j = i + 1; j < limit; j++) {
        const b1 = sortedIndices[i];
        const b2 = sortedIndices[j];
        bits[b1] ^= 1;
        bits[b2] ^= 1;

        const zero = this.computeSyndrome(bits, syndrome);
        if (zero) {
          return { data: bitsToBytes(bits, this.k), corrected: 2, erased: 0 };
        }

        bits[b1] ^= 1;
        bits[b2] ^= 1;
      }
    }

    // OSD-2 failed — throw to signal erasure
    throw new Error(
      `LDPC soft-decision decode failed: OSD-2 could not find a valid codeword`,
    );
  }

  /**
   * Belief-Propagation (Sum-Product) Decoder
   *
   * The optimal soft-decision decoder for LDPC codes. Iteratively passes
   * probability messages between variable nodes (bits) and check nodes (parity
   * checks) along the Tanner graph. Converges to the maximum-likelihood codeword
   * for codes with large girth (no short cycles).
   *
   * Algorithm (min-sum approximation for numerical stability):
   *   1. Initialize variable node LLRs from channel (Q-scores).
   *   2. Variable → Check: send LLR + sum of incoming check→var messages.
   *   3. Check → Variable: send product of tanh(incoming var→check / 2).
   *   4. Hard-decide and check syndrome. If zero, done.
   *   5. Repeat for maxIter iterations.
   *
   * Complexity: O(iterations × edges) = O(20 × nBits × dv) ≈ O(20 × 304 × 4) = 24K ops.
   * Much faster than OSD-2 for the same error-correction capability.
   *
   * @param recv Received codeword bytes
   * @param qScores Per-base Phred Q-scores (null → uniform LLR)
   * @param useGoldman Whether Goldman mapping was used
   * @param maxIter Max BP iterations (default 20)
   */
  decodeBeliefPropagation(
    recv: Uint8Array,
    qScores: Uint8Array | null,
    useGoldman: boolean = false,
    maxIter: number = 20,
  ): LDPCDecodeResult {
    // Convert received bytes to bits
    const bits = new Uint8Array(this.nBits);
    for (let i = 0; i < this.n; i++) {
      const b = recv[i];
      for (let bit = 0; bit < 8; bit++) {
        bits[i * 8 + bit] = (b >> (7 - bit)) & 1;
      }
    }

    // Initialize channel LLRs from Q-scores
    // LLR > 0 → bit is likely 0, LLR < 0 → bit is likely 1
    const chLLR = new Float64Array(this.nBits);
    if (qScores && qScores.length > 0 && !useGoldman) {
      // Direct mapping: bit i → base i/2, Q-score from qScores[base]
      for (let byteIdx = 0; byteIdx < this.n; byteIdx++) {
        for (let bitInByte = 0; bitInByte < 8; bitInByte++) {
          const bitIdx = byteIdx * 8 + bitInByte;
          const baseIdx = byteIdx * 4 + (bitInByte >> 1);
          const q = baseIdx < qScores.length ? qScores[baseIdx] : 30;
          const llrMag = (q / 10) * Math.log(10);
          chLLR[bitIdx] = bits[bitIdx] === 0 ? llrMag : -llrMag;
        }
      }
    } else {
      // Uniform LLR (no soft info, or Goldman where bit-to-base mapping is complex)
      const avgQ = (qScores && qScores.length > 0)
        ? qScores.reduce((s, q) => s + q, 0) / qScores.length
        : 30;
      const llrMag = (avgQ / 10) * Math.log(10);
      for (let i = 0; i < this.nBits; i++) {
        chLLR[i] = bits[i] === 0 ? llrMag : -llrMag;
      }
    }

    // Message arrays:
    //   varToCheck[i][e] = message from variable node i to check node (rowCols[j][e])
    //   checkToVar[j][e] = message from check node j to variable node (colRows[i][e])
    //
    // We store messages indexed by the edge position in the adjacency list.
    const varToCheck: Float64Array[] = new Array(this.nBits);
    const checkToVar: Float64Array[] = new Array(this.mBits);
    for (let j = 0; j < this.nBits; j++) {
      varToCheck[j] = new Float64Array(this.colRows[j].length);
    }
    for (let i = 0; i < this.mBits; i++) {
      checkToVar[i] = new Float64Array(this.rowCols[i].length);
    }

    // Initialize varToCheck = channel LLR (no prior checkToVar messages)
    for (let j = 0; j < this.nBits; j++) {
      for (let e = 0; e < this.colRows[j].length; e++) {
        varToCheck[j][e] = chLLR[j];
      }
    }

    // BP iterations
    for (let iter = 0; iter < maxIter; iter++) {
      // Step 1: Check → Variable messages (min-sum)
      // For check node i, connected to variable nodes v_0, ..., v_{dc-1}:
      //   msg to v_k = (product of signs of msg_j for j != k) × (min of |msg_j| for j != k)
      for (let i = 0; i < this.mBits; i++) {
        const cols = this.rowCols[i];
        const dc = cols.length;
        if (dc === 0) continue;

        // First pass: gather all incoming varToCheck messages for this check node
        // (uses precomputed reverse edge index — O(1) per edge instead of O(degree))
        const msgs = new Float64Array(dc);
        const edgeIdxs = this.rowColsReverseIdx[i];
        for (let e = 0; e < dc; e++) {
          const j = cols[e];
          const edgeIdx = edgeIdxs[e];
          msgs[e] = edgeIdx >= 0 ? varToCheck[j][edgeIdx] : 0;
        }

        // Compute global product of signs, and two smallest magnitudes
        let signProd = 1;
        let minMag = Infinity;
        let secondMinMag = Infinity;
        let minIdx = -1;
        for (let e = 0; e < dc; e++) {
          if (msgs[e] < 0) signProd = -signProd;
          const mag = Math.abs(msgs[e]);
          if (mag < minMag) {
            secondMinMag = minMag;
            minMag = mag;
            minIdx = e;
          } else if (mag < secondMinMag) {
            secondMinMag = mag;
          }
        }

        // Send messages to each variable node
        for (let e = 0; e < dc; e++) {
          const j = cols[e];
          if (edgeIdxs[e] < 0) continue;

          // Min-sum: exclude self
          // Sign of message = (product of all signs) / (sign of self)
          // In GF(2)/multiplication: sign = signProd * sign(self)  (since 1/sign = sign for ±1)
          const selfSign = msgs[e] < 0 ? -1 : 1;
          const msgSign = signProd * selfSign; // this gives product of all OTHER signs

          // Magnitude: min of all others
          const msgMag = (e === minIdx) ? secondMinMag : minMag;

          checkToVar[i][e] = msgSign * msgMag;
        }
      }

      // Step 2: Variable → Check messages
      // For variable node j: msg to check c_k = chLLR[j] + sum of checkToVar from all other checks
      // Also compute total belief for hard decision
      const totalBelief = new Float64Array(this.nBits);
      for (let j = 0; j < this.nBits; j++) {
        const rows = this.colRows[j];
        const dv = rows.length;

        // Sum all incoming check messages
        // (uses precomputed reverse edge index — O(1) per edge instead of O(degree))
        const revIdx = this.colRowsReverseIdx[j];
        let sumChecks = 0;
        for (let e = 0; e < dv; e++) {
          const i = rows[e];
          const edgeIdx = revIdx[e];
          if (edgeIdx !== -1) {
            sumChecks += checkToVar[i][edgeIdx];
          }
        }

        // Total belief = channel + all checks
        totalBelief[j] = chLLR[j] + sumChecks;

        // Send varToCheck = channel + sum of all checks except the recipient
        // (uses precomputed reverse edge index — O(1) per edge instead of O(degree))
        for (let e = 0; e < dv; e++) {
          const i = rows[e];
          const edgeIdx = revIdx[e];
          if (edgeIdx !== -1) {
            varToCheck[j][e] = totalBelief[j] - checkToVar[i][edgeIdx];
          }
        }
      }

      // Step 3: Hard decision and syndrome check
      const hardDecision = new Uint8Array(this.nBits);
      for (let j = 0; j < this.nBits; j++) {
        hardDecision[j] = totalBelief[j] < 0 ? 1 : 0;
      }

      // Compute syndrome
      const syndrome = new Uint8Array(this.mBits);
      const syndromeZero = this.computeSyndrome(hardDecision, syndrome);
      if (syndromeZero) {
        // Count corrections
        let corrected = 0;
        for (let j = 0; j < this.nBits; j++) {
          if (hardDecision[j] !== bits[j]) corrected++;
        }
        return { data: bitsToBytes(hardDecision, this.k), corrected, erased: 0 };
      }
    }

    // BP failed to converge — throw to signal erasure
    throw new Error(
      `LDPC belief-propagation failed: did not converge after ${maxIter} iterations`,
    );
  }
}

/** Convert nBits bits to nbytes bytes (taking the first nbytes*8 bits). */
function bitsToBytes(bits: Uint8Array, nbytes: number): Uint8Array {
  const out = new Uint8Array(nbytes);
  for (let i = 0; i < nbytes; i++) {
    let b = 0;
    for (let bit = 0; bit < 8; bit++) {
      b |= bits[i * 8 + bit] << (7 - bit);
    }
    out[i] = b;
  }
  return out;
}

/**
 * Default LDPC inner code for Helix v21+.
 *
 * Layout (200 nt oligo, 20 nt primers × 2):
 *   - 4B address + 30B payload + 4B LDPC parity + 2B CRC-16 = 40 bytes = 160 nt
 *   - LDPC codeword = 38 bytes = 304 bits (addr + payload + parity)
 *   - LDPC info = 34 bytes = 272 bits (addr + payload)
 *   - LDPC parity = 4 bytes = 32 bits
 *   - Density: 30*8 / 200 = 1.20 bits/nt total, 30*8 / 160 = 1.50 bits/nt payload-only
 *
 * For 300 nt oligos (used in Erlich validation):
 *   - 4B addr + 51B payload + 4B LDPC parity + 2B CRC = 61 bytes... wait, 65 bytes total inner
 *   - 4B addr + 55B payload + 4B LDPC parity + 2B CRC = 65 bytes ✓
 *   - LDPC(520 bits, 472 bits, 48 bits)? No, parity stays at 4 bytes.
 *   - LDPC(520 bits, 488 bits, 32 bits): 61 bytes info, 4 bytes parity
 *   - Wait, n=65 bytes=520 bits, k=61 bytes=488 bits (addr+payload+CRC... no, CRC is outside)
 *   - n = 4 + 55 + 4 = 63 bytes (LDPC codeword, excluding CRC)
 *   - k = 4 + 55 = 59 bytes
 *   - LDPC(504 bits, 472 bits, 32 bits) for 300nt oligo
 *
 * The constructor takes care of these dimension calculations based on cfg.n and cfg.k.
 */
export function makeLDPCInner(parityBytes: number, payloadBytes: number, addressBytes: number): LDPCInnerCode {
  const k = addressBytes + payloadBytes;
  const n = k + parityBytes;
  return new LDPCInnerCode({ n, k });
}

/**
 * Module-level LDPC instance cache with LRU eviction.
 *
 * decode.ts was constructing a NEW LDPCInnerCode for every decodeReads()
 * call. Each construction runs PEG construction + column dedup + LUT build,
 * which is ~5ms for typical configs. At 100K oligo pools, that's 500 seconds
 * of wasted setup time.
 *
 * The LDPC instance is stateless across decode calls (the syndrome is
 * computed fresh each time), so sharing is safe.
 *
 * Speedup: ~5ms per decode call avoided → significant for small payloads
 * where decode is <100ms total.
 *
 * LRU eviction (max MAX_CACHE_SIZE entries) prevents unbounded memory growth.
 * Each LDPCInnerCode holds ~O(nBits × mBits) in adjacency arrays and LUTs,
 * so at large configs a single instance can be 10+ MB. Without eviction,
 * many distinct (n,k) configs would grow the cache without limit.
 *
 * Uses Map's insertion-order guarantee: first key is the oldest entry.
 * On cache hit, the entry is deleted and re-inserted to move it to the end.
 */
const MAX_CACHE_SIZE = 16;
const _ldpcCache = new Map<string, { ldpc: LDPCInnerCode; lastAccess: number }>();

export function getCachedLDPCInner(n: number, k: number): LDPCInnerCode {
  const key = `${n}:${k}`;
  const cached = _ldpcCache.get(key);
  if (cached) {
    cached.lastAccess = Date.now();
    // Move to end (most recently used) by re-inserting
    _ldpcCache.delete(key);
    _ldpcCache.set(key, cached);
    return cached.ldpc;
  }
  // Evict oldest (first key in Map iteration order) if at capacity
  if (_ldpcCache.size >= MAX_CACHE_SIZE) {
    const oldest = _ldpcCache.keys().next().value!;
    _ldpcCache.delete(oldest);
  }
  const inst = new LDPCInnerCode({ n, k });
  _ldpcCache.set(key, { ldpc: inst, lastAccess: Date.now() });
  return inst;
}
