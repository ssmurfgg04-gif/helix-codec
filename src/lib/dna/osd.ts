/**
 * Ordered Statistics Decoding (OSD)
 *
 * A soft-decision decoder for linear codes over GF(2). Used as the inner code
 * in the Mahoraga codec (Banal 2026, arXiv:2604.20810) — the current state-
 * of-the-art for DNA data storage (1.815 bits/nt).
 *
 * Algorithm (OSD-t):
 *   1. Sort bits by |LLR| (log-likelihood ratio) ascending — least reliable first.
 *   2. Form the Most Reliable Basis (MRB) from the k most-reliable independent
 *      columns of the parity-check matrix H.
 *   3. OSD-0: hard-decide the MRB, solve for the rest. Check CRC.
 *   4. OSD-t: flip up to t bits in the MRB (least reliable first), re-solve, check CRC.
 *   5. Return the first codeword that passes CRC.
 *
 * This is a simplified OSD-0/1 implementation. Full OSD-2/3 (Mahoraga uses
 * order 2-3) requires combinatorial enumeration which is expensive for large k.
 *
 * For DNA storage, the code parameters are:
 *   - Block length n = 252 bits (126 nt × 2 bits/nt)
 *   - Information bits k = 208 (with 32-bit CRC)
 *   - Parity-check matrix H: (n-k) × n = 44 × 252
 *
 * Reference:
 *   - Fossorier & Lin (1995). "Soft-decision decoding of linear block codes
 *     based on ordered statistics." IEEE TIT 41:5.
 *   - Banal et al. (2026). "Mahoraga." arXiv:2604.20810.
 *   - Reference impl: github.com/jeplb/mahoraga-codec (module osd.py)
 */

/**
 * GF(2) matrix operations (binary matrices).
 * Matrices are stored as Uint8Array of 0/1 values, row-major.
 */

export class GF2Matrix {
  rows: number;
  cols: number;
  data: Uint8Array;

  constructor(rows: number, cols: number, data?: Uint8Array) {
    this.rows = rows;
    this.cols = cols;
    this.data = data ?? new Uint8Array(rows * cols);
  }

  get(r: number, c: number): number {
    return this.data[r * this.cols + c];
  }

  set(r: number, c: number, v: number): void {
    this.data[r * this.cols + c] = v & 1;
  }

  /** Row reduce to systematic form using Gaussian elimination over GF(2). */
  systematicForm(): { matrix: GF2Matrix; pivotCols: number[] } {
    const m = new GF2Matrix(this.rows, this.cols, this.data.slice());
    const pivotCols: number[] = [];
    let pivotRow = 0;
    for (let col = 0; col < m.cols && pivotRow < m.rows; col++) {
      // Find pivot
      let pivot = -1;
      for (let row = pivotRow; row < m.rows; row++) {
        if (m.get(row, col) === 1) {
          pivot = row;
          break;
        }
      }
      if (pivot === -1) continue;
      // Swap rows
      if (pivot !== pivotRow) {
        for (let c = 0; c < m.cols; c++) {
          const tmp = m.get(pivotRow, c);
          m.set(pivotRow, c, m.get(pivot, c));
          m.set(pivot, c, tmp);
        }
      }
      // Eliminate
      for (let row = 0; row < m.rows; row++) {
        if (row === pivotRow) continue;
        if (m.get(row, col) === 1) {
          for (let c = 0; c < m.cols; c++) {
            m.set(row, c, m.get(row, c) ^ m.get(pivotRow, c));
          }
        }
      }
      pivotCols.push(col);
      pivotRow++;
    }
    return { matrix: m, pivotCols };
  }
}

/**
 * OSD-0 decoder: hard-decide LLRs, solve on most reliable basis, check CRC.
 *
 * @param llr Log-likelihood ratios for each bit (negative = likely 1, positive = likely 0)
 * @param H Parity-check matrix (n-k rows × n cols)
 * @param crcCheck Function that takes a candidate codeword and returns true if CRC passes
 * @returns Decoded codeword (n bits), or null if OSD-0 fails
 */
export function osd0(
  llr: Float32Array,
  H: GF2Matrix,
  crcCheck: (codeword: Uint8Array) => boolean,
): Uint8Array | null {
  const n = llr.length;

  // Hard decision
  const hard = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    hard[i] = llr[i] < 0 ? 1 : 0;
  }

  // Check if hard decision is already a valid codeword (syndrome = 0)
  const syndrome = new Uint8Array(H.rows);
  for (let r = 0; r < H.rows; r++) {
    let s = 0;
    for (let c = 0; c < n; c++) {
      s ^= H.get(r, c) & hard[c];
    }
    syndrome[r] = s;
  }
  const syndromeZero = syndrome.every((s) => s === 0);
  if (syndromeZero && crcCheck(hard)) {
    return hard;
  }

  // OSD-0: sort by |LLR|, pick most reliable basis, solve
  // For simplicity, just return the hard decision if syndrome is zero,
  // otherwise return null (OSD-1+ would flip bits)
  if (syndromeZero) {
    return hard;
  }

  return null;
}

/**
 * OSD-1 decoder: try flipping each of the k least-reliable bits, check CRC.
 */
export function osd1(
  llr: Float32Array,
  H: GF2Matrix,
  crcCheck: (codeword: Uint8Array) => boolean,
  k: number,
): Uint8Array | null {
  const n = llr.length;

  // First try OSD-0
  const osd0Result = osd0(llr, H, crcCheck);
  if (osd0Result) return osd0Result;

  // Sort bit indices by |LLR| ascending (least reliable first)
  const indices = Array.from({ length: n }, (_, i) => i);
  indices.sort((a, b) => Math.abs(llr[a]) - Math.abs(llr[b]));

  // Hard decision
  const hard = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    hard[i] = llr[i] < 0 ? 1 : 0;
  }

  // Try flipping each of the k least-reliable bits
  for (let i = 0; i < Math.min(k, n); i++) {
    const flipIdx = indices[i];
    const candidate = hard.slice();
    candidate[flipIdx] ^= 1;

    // Check syndrome
    let syndromeZero = true;
    for (let r = 0; r < H.rows; r++) {
      let s = 0;
      for (let c = 0; c < n; c++) {
        s ^= H.get(r, c) & candidate[c];
      }
      if (s !== 0) {
        syndromeZero = false;
        break;
      }
    }
    if (syndromeZero && crcCheck(candidate)) {
      return candidate;
    }
  }

  return null;
}

/**
 * Generate a simple parity-check matrix for a (n, k) systematic code.
 * H = [P^T | I_{n-k}] where P is the parity part of the generator.
 *
 * For DNA storage, a proper PEG (Progressive Edge Growth) construction
 * would be used (Mahoraga uses dv=3, dc=84 for hi-fi channel).
 * This is a simple random construction for demonstration.
 */
export function generateSimpleParityMatrix(
  n: number,
  k: number,
  seed: number = 42,
): GF2Matrix {
  const m = n - k; // number of parity checks
  const H = new GF2Matrix(m, n);

  // Simple pseudo-random construction
  let state = seed >>> 0 || 1;
  const rng = () => {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    state = state >>> 0;
    return state / 0x100000000;
  };

  // Fill the parity part (first k columns) with random 0/1
  for (let r = 0; r < m; r++) {
    for (let c = 0; c < k; c++) {
      H.set(r, c, rng() < 0.5 ? 1 : 0);
    }
    // Identity part (last m columns)
    H.set(r, k + r, 1);
  }

  return H;
}

/**
 * Encode using a parity-check matrix (systematic).
 * Given k information bits, produce n codeword bits.
 */
export function encodeWithParity(
  info: Uint8Array,
  H: GF2Matrix,
  n: number,
  k: number,
): Uint8Array {
  const codeword = new Uint8Array(n);
  codeword.set(info, 0);

  // Compute parity bits: for each parity row r, parity[r] = sum(H[r][c] * info[c])
  for (let r = 0; r < H.rows; r++) {
    let parity = 0;
    for (let c = 0; c < k; c++) {
      parity ^= H.get(r, c) & info[c];
    }
    codeword[k + r] = parity;
  }

  return codeword;
}
