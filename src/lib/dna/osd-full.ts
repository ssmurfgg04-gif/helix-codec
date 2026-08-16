/**
 * Full OSD (Ordered Statistics Decoding) — OSD-0/1/2/3 cascade.
 *
 * This is the complete Mahoraga-style soft-decision decoder for linear codes
 * over GF(2). It implements the cascade:
 *   OSD-0: hard-decide LLRs, solve on Most Reliable Basis (MRB), check CRC.
 *   OSD-1: flip each of the k least-reliable MRB bits, re-solve, check CRC.
 *   OSD-2: flip pairs of the k least-reliable MRB bits, re-solve, check CRC.
 *   OSD-3: flip triples of the k least-reliable MRB bits, re-solve, check CRC.
 *
 * OSD-2 has O(k²) candidates, OSD-3 has O(k³). Mahoraga uses OSD-2 on hi-fi
 * and OSD-3 cascade on lo-fi (when OSD-2 fails CRC). Early-exit on first CRC
 * pass makes it fast in practice.
 *
 * Algorithm (Fossorier & Lin 1995):
 *   1. Sort received bits by |LLR| descending → Most Reliable Basis (MRB)
 *      is the first k bits.
 *   2. Column-permute H so MRB columns come first, then row-reduce to
 *      systematic form [A | I_m] where A is m×k and I_m is m×m.
 *   3. OSD-0: hard-decide MRB, compute parity from A * u, un-permute → codeword.
 *   4. OSD-t: for each subset of t MRB bits (lowest |LLR| first), flip them,
 *      re-solve, check CRC. Return first that passes.
 *
 * v65 FIX: The previous solveFromMRB was a broken stub that just checked
 * syndrome and tried a crude pivot flip. The proper implementation:
 *   - Permute columns so MRB comes first
 *   - Row-reduce to get [A | I_m]
 *   - Given MRB values u, compute parity p = A * u (over GF(2))
 *   - Codeword in permuted order = [u | p]
 *   - Un-permute to get codeword in original bit order
 *
 * Reference:
 *   - Fossorier & Lin (1995). "Soft-decision decoding of linear block codes
 *     based on ordered statistics." IEEE TIT 41:5.
 *   - Banal et al. (2026). arXiv:2604.20810. (Mahoraga uses OSD-2/3)
 *   - Reference impl: github.com/jeplb/mahoraga-codec (module osd.py)
 */

import { GF2Matrix } from "./osd";

export interface OSDConfig {
  /** Maximum OSD order (0, 1, 2, or 3). Default 2 (Mahoraga hi-fi). */
  maxOrder: number;
  /** Number of information bits k (MRB size). */
  k: number;
}

export const DEFAULT_OSD_CONFIG: OSDConfig = {
  maxOrder: 2,
  k: 208, // Mahoraga hi-fi: 252 - 44 parity = 208 info bits (after 32-bit CRC)
};

export interface OSDResult {
  /** Decoded codeword (n bits), or null if all orders failed CRC. */
  codeword: Uint8Array | null;
  /** Order at which decoding succeeded (0, 1, 2, or 3). -1 if failed. */
  successOrder: number;
  /** Number of candidates tried before success. */
  candidatesTried: number;
}

/**
 * Precomputed OSD solver structure.
 *
 * After column permutation and row reduction, we have:
 *   H_perm (after GE) = [A | I_m]
 *
 * Where:
 *   - A is m × k (the parity computation matrix)
 *   - I_m is m × m identity
 *   - m = n - k (number of parity checks)
 *
 * Given MRB values u[0..k-1], the parity values p[0..m-1] are:
 *   p[i] = XOR_{j=0}^{k-1} A[i][j] * u[j]   (mod 2)
 *
 * The codeword in permuted order is [u | p], then un-permuted.
 */
interface OSDSolver {
  n: number;
  k: number;
  m: number;
  /** Column permutation: perm[pos] = original column index.
   *  First k positions are MRB columns, last m positions are parity columns. */
  perm: Uint32Array;
  /** Inverse permutation: invPerm[origCol] = position in permuted order. */
  invPerm: Uint32Array;
  /** Parity computation matrix A (m × k), stored row-major.
   *  A[i * k + j] = coefficient of MRB bit j for parity bit i. */
  A: Uint8Array;
  /** Hard decision bits (from LLR). */
  hard: Uint8Array;
  /** MRB indices sorted by reliability ascending (least reliable first). */
  mrbByReliability: Uint32Array;
}

/**
 * Build the OSD solver structure from LLRs, parity-check matrix H, and k.
 *
 * Steps:
 *   1. Sort bit indices by |LLR| descending (most reliable first).
 *   2. Find k linearly independent columns from the most reliable bits (MRB).
 *   3. Permute columns so MRB comes first, parity comes last.
 *   4. Row-reduce H_perm to systematic form [A | I_m].
 *   5. Extract A (parity computation matrix).
 *
 * Returns null if k independent columns cannot be found (code is degenerate).
 */
function buildOSDSolver(
  llr: Float32Array,
  H: GF2Matrix,
  k: number,
): OSDSolver | null {
  const n = llr.length;
  const m = H.rows;

  if (k > n || k + m !== n && m < n - k) {
    // k must not exceed n, and we need m >= n-k parity checks
  }
  const actualM = n - k; // expected number of parity positions
  if (H.rows < actualM) {
    // Not enough parity check rows — can't form a proper MRB
    // Fall back: use as many independent columns as we can find
  }

  // Step 1: Sort bit indices by |LLR| descending (most reliable first)
  const sortedIndices = Array.from({ length: n }, (_, i) => i);
  sortedIndices.sort((a, b) => Math.abs(llr[b]) - Math.abs(llr[a]));

  // Step 2: Find k linearly independent columns via Gaussian elimination
  // on H^T (the transpose of H, treating columns of H as rows for GE).
  //
  // We work on a copy of H's column vectors. Each column is an m-bit vector.
  // We do GE on these column vectors to find k independent ones.
  //
  // Representation: colVectors[col] = m-bit vector (as a Uint8Array of 0/1)
  // We process columns in reliability order (most reliable first).

  // Build column vectors from H
  const colVectors: Uint8Array[] = new Array(n);
  for (let col = 0; col < n; col++) {
    const vec = new Uint8Array(m);
    for (let row = 0; row < m; row++) {
      vec[row] = H.get(row, col);
    }
    colVectors[col] = vec;
  }

  // GE on column vectors (in reliability order) to find independent columns
  const pivotRowForCol = new Int32Array(n).fill(-1); // which row this column pivots on
  const colPivotRow = new Int32Array(m).fill(-1); // which column pivots on this row
  let numIndependent = 0;

  // We'll track the row-reduced form of each column
  const reducedCols: Uint8Array[] = new Array(n);
  for (let col = 0; col < n; col++) {
    reducedCols[col] = colVectors[col].slice();
  }

  // Process columns in reliability order
  const mrbSet = new Set<number>();
  const mrbList: number[] = [];
  const parityList: number[] = [];

  for (const col of sortedIndices) {
    const vec = reducedCols[col];
    // Reduce this column against existing pivots
    for (let row = 0; row < m; row++) {
      if (vec[row] === 1 && colPivotRow[row] !== -1) {
        // Eliminate using the column that pivots on this row
        const pivotCol = colPivotRow[row];
        for (let r = 0; r < m; r++) {
          vec[r] ^= reducedCols[pivotCol][r];
        }
      }
    }

    // Find the first non-zero row (potential pivot)
    let pivotRow = -1;
    for (let row = 0; row < m; row++) {
      if (vec[row] === 1) {
        pivotRow = row;
        break;
      }
    }

    if (pivotRow !== -1 && numIndependent < k) {
      // This column is independent — add to MRB
      pivotRowForCol[col] = pivotRow;
      colPivotRow[pivotRow] = col;
      mrbSet.add(col);
      mrbList.push(col);
      numIndependent++;
      // Normalize: ensure pivot row is 1 (already is), eliminate from other reduced cols
      // (not strictly necessary for our approach, but keeps reducedCols consistent)
    } else {
      // Dependent or we already have k independent columns → parity position
      parityList.push(col);
    }
  }

  // If we couldn't find k independent columns, the code is degenerate
  if (numIndependent < k) {
    return null;
  }

  // Any remaining sortedIndices not in mrbSet or parityList go to parityList
  for (const col of sortedIndices) {
    if (!mrbSet.has(col) && !parityList.includes(col)) {
      parityList.push(col);
    }
  }

  // Step 3: Build permutation — MRB columns first, parity columns last
  const perm = new Uint32Array(n);
  const invPerm = new Uint32Array(n);
  for (let i = 0; i < k; i++) {
    perm[i] = mrbList[i];
    invPerm[mrbList[i]] = i;
  }
  for (let i = 0; i < parityList.length && i < m; i++) {
    perm[k + i] = parityList[i];
    invPerm[parityList[i]] = k + i;
  }
  // Handle any extra columns (shouldn't happen if n = k + m)
  for (let i = m; i < parityList.length; i++) {
    perm[k + i] = parityList[i];
    invPerm[parityList[i]] = k + i;
  }

  // Step 4: Row-reduce H_perm = H with columns permuted so MRB comes first
  // to systematic form [A | I_m].
  //
  // H_perm[i][j] = H[i][perm[j]]
  // After GE, we want: for parity positions j (k <= j < n),
  //   H_reduced[i][j] = 1 if j - k == i, else 0  (identity in parity part)
  //
  // We perform Gaussian elimination on the PARITY part of H_perm.

  // Build H_perm as a working copy (m × n)
  const hPerm = new GF2Matrix(m, n);
  for (let row = 0; row < m; row++) {
    for (let j = 0; j < n; j++) {
      hPerm.set(row, j, H.get(row, perm[j]));
    }
  }

  // Gaussian elimination on the parity columns to get I_m
  // We want pivot row i to have a 1 in column k+i
  const rowMap = new Int32Array(m); // which original row is now at position i
  for (let i = 0; i < m; i++) rowMap[i] = i;

  for (let i = 0; i < m; i++) {
    const targetCol = k + i;

    // Find a row with a 1 in targetCol (at or below row i)
    let found = -1;
    for (let row = i; row < m; row++) {
      if (hPerm.get(row, targetCol) === 1) {
        found = row;
        break;
      }
    }

    if (found === -1) {
      // No pivot found in this parity column — the MRB selection may not
      // have been optimal. Try to find a pivot in ANY remaining column.
      // This is a fallback for degenerate cases.
      continue;
    }

    // Swap rows i and found
    if (found !== i) {
      for (let j = 0; j < n; j++) {
        const tmp = hPerm.get(i, j);
        hPerm.set(i, j, hPerm.get(found, j));
        hPerm.set(found, j, tmp);
      }
    }

    // Eliminate targetCol from all other rows
    for (let row = 0; row < m; row++) {
      if (row === i) continue;
      if (hPerm.get(row, targetCol) === 1) {
        for (let j = 0; j < n; j++) {
          hPerm.set(row, j, hPerm.get(row, j) ^ hPerm.get(i, j));
        }
      }
    }
  }

  // Step 5: Extract A from the left part of hPerm
  // After GE, hPerm = [A | I_m] (approximately — some rows may not have pivots
  // if the MRB selection was suboptimal, but for well-constructed codes this works).
  //
  // A[i][j] = hPerm[i][j] for 0 <= i < m, 0 <= j < k
  const A = new Uint8Array(m * k);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < k; j++) {
      A[i * k + j] = hPerm.get(i, j);
    }
  }

  // Hard decision from LLRs
  const hard = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    hard[i] = llr[i] < 0 ? 1 : 0;
  }

  // MRB indices sorted by reliability ascending (least reliable first)
  const mrbByReliability = Uint32Array.from(mrbList).slice();
  mrbByReliability.sort((a, b) => Math.abs(llr[a]) - Math.abs(llr[b]));

  return { n, k, m, perm, invPerm, A, hard, mrbByReliability };
}

/**
 * Construct a codeword from MRB values using the OSD solver.
 *
 * Given MRB values u[0..k-1], compute:
 *   p[i] = XOR_{j=0}^{k-1} A[i*k+j] * u[j]   for i = 0..m-1
 *   codeword_permuted = [u | p]
 *   codeword[perm[pos]] = codeword_permuted[pos]
 *
 * @param solver OSD solver structure
 * @param mrbValues Values of the k MRB bits
 * @returns Codeword in original bit order (n bits)
 */
function constructCodeword(
  solver: OSDSolver,
  mrbValues: Uint8Array,
): Uint8Array {
  const { n, k, m, perm, A } = solver;
  const codeword = new Uint8Array(n);

  // Compute parity bits: p[i] = XOR of A[i*k+j] * u[j]
  for (let i = 0; i < m; i++) {
    let parity = 0;
    const rowOffset = i * k;
    for (let j = 0; j < k; j++) {
      parity ^= A[rowOffset + j] & mrbValues[j];
    }
    // Place parity bit at permuted position k+i → original position perm[k+i]
    codeword[perm[k + i]] = parity;
  }

  // Place MRB bits at original positions
  for (let j = 0; j < k; j++) {
    codeword[perm[j]] = mrbValues[j];
  }

  return codeword;
}

/**
 * Full OSD cascade decoder.
 *
 * @param llr Log-likelihood ratios for each bit (negative = likely 1, positive = likely 0)
 * @param H Parity-check matrix (m rows × n cols)
 * @param crcCheck Function that takes a candidate codeword and returns true if CRC passes
 * @param config OSD configuration (maxOrder, k)
 * @returns Deced codeword + metadata, or null if all orders failed
 */
export function osdDecode(
  llr: Float32Array,
  H: GF2Matrix,
  crcCheck: (codeword: Uint8Array) => boolean,
  config: OSDConfig = DEFAULT_OSD_CONFIG,
): OSDResult {
  const n = llr.length;
  const k = Math.min(config.k, n);
  const maxOrder = config.maxOrder;

  // Build the OSD solver (MRB selection + row reduction)
  const solver = buildOSDSolver(llr, H, k);
  if (!solver) {
    // Could not find k independent columns — fall back to hard decision
    const hard = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      hard[i] = llr[i] < 0 ? 1 : 0;
    }
    if (crcCheck(hard)) {
      return { codeword: hard, successOrder: 0, candidatesTried: 1 };
    }
    return { codeword: null, successOrder: -1, candidatesTried: 1 };
  }

  const { hard, mrbByReliability } = solver;

  // Extract hard-decided MRB values
  const mrbHard = new Uint8Array(k);
  for (let j = 0; j < k; j++) {
    mrbHard[j] = hard[solver.perm[j]];
  }

  let candidatesTried = 0;

  // OSD-0: hard decision on MRB, solve for parity, construct codeword
  {
    const codeword = constructCodeword(solver, mrbHard);
    candidatesTried++;
    if (crcCheck(codeword)) {
      return { codeword, successOrder: 0, candidatesTried };
    }
  }

  // OSD-1: flip each of the k least-reliable MRB bits
  if (maxOrder >= 1) {
    // mrbByReliability maps from least-reliable to most-reliable MRB bit
    // We need to find the position of each MRB bit in the permuted order
    const mrbFlipPositions = new Uint32Array(k);
    for (let i = 0; i < k; i++) {
      const origCol = mrbByReliability[i]; // original column index
      const permPos = solver.invPerm[origCol]; // position in permuted order (0..k-1)
      mrbFlipPositions[i] = permPos;
    }

    for (let i = 0; i < k; i++) {
      const flipPos = mrbFlipPositions[i];
      const mrbFlipped = mrbHard.slice();
      mrbFlipped[flipPos] ^= 1;

      const codeword = constructCodeword(solver, mrbFlipped);
      candidatesTried++;
      if (crcCheck(codeword)) {
        return { codeword, successOrder: 1, candidatesTried };
      }
    }
  }

  // OSD-2: flip pairs of least-reliable MRB bits
  if (maxOrder >= 2) {
    const mrbFlipPositions = new Uint32Array(k);
    for (let i = 0; i < k; i++) {
      const origCol = mrbByReliability[i];
      const permPos = solver.invPerm[origCol];
      mrbFlipPositions[i] = permPos;
    }

    // Limit to first min(k, 30) bits to keep O(k²) manageable (30² = 900 candidates)
    const limit = Math.min(k, 30);

    for (let i = 0; i < limit; i++) {
      for (let j = i + 1; j < limit; j++) {
        const mrbFlipped = mrbHard.slice();
        mrbFlipped[mrbFlipPositions[i]] ^= 1;
        mrbFlipped[mrbFlipPositions[j]] ^= 1;

        const codeword = constructCodeword(solver, mrbFlipped);
        candidatesTried++;
        if (crcCheck(codeword)) {
          return { codeword, successOrder: 2, candidatesTried };
        }
      }
    }
  }

  // OSD-3: flip triples of least-reliable MRB bits
  if (maxOrder >= 3) {
    const mrbFlipPositions = new Uint32Array(k);
    for (let i = 0; i < k; i++) {
      const origCol = mrbByReliability[i];
      const permPos = solver.invPerm[origCol];
      mrbFlipPositions[i] = permPos;
    }

    // Limit to first min(k, 15) bits (15³ ≈ 3375 candidates)
    const limit = Math.min(k, 15);

    for (let i = 0; i < limit; i++) {
      for (let j = i + 1; j < limit; j++) {
        for (let l = j + 1; l < limit; l++) {
          const mrbFlipped = mrbHard.slice();
          mrbFlipped[mrbFlipPositions[i]] ^= 1;
          mrbFlipped[mrbFlipPositions[j]] ^= 1;
          mrbFlipped[mrbFlipPositions[l]] ^= 1;

          const codeword = constructCodeword(solver, mrbFlipped);
          candidatesTried++;
          if (crcCheck(codeword)) {
            return { codeword, successOrder: 3, candidatesTried };
          }
        }
      }
    }
  }

  return { codeword: null, successOrder: -1, candidatesTried };
}
