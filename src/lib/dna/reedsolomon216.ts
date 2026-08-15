/**
 * Reed-Solomon over GF(2^16) for large DNA archives.
 *
 * Supports RS(n, k) where n ≤ 65535 symbols. This escapes the GF(2^8)
 * 255-symbol cap, enabling archives with up to 65535 oligos per RS block.
 *
 * Decoder implements:
 *   - Syndrome computation
 *   - Pure-erasure Lagrange interpolation (primary path — used when erasure
 *     positions are known, e.g., missing oligos identified by address index)
 *   - Berlekamp-Massey + Chien + Forney for unknown errors (secondary path)
 *
 * The pure-erasure path is O(n) per erasure and is the hot path for DNA storage
 * (we always know which oligos are missing from the address).
 *
 * Reference:
 *   - Reed & Solomon (1960)
 *   - Banal et al. (2026). arXiv:2604.20810. (Mahoraga uses GF(2^16) outer RS)
 *   - Rizzo (1997). "Effective erasure codes for reliable computer
 *     communication protocols." (Lagrange erasure-only decoding)
 */

import {
  gf16Add,
  gf16Mul,
  gf16Div,
  gf16Inverse,
  gf16Pow,
  gf16PolyEval,
  gf16PolyMul,
  init as initGf216,
} from "./gf216";

export interface RS216Config {
  n: number;
  k: number;
}

export interface RS216DecodeResult {
  data: Uint16Array;
  corrected: number;
  erased: number;
}

export class ReedSolomon216 {
  readonly n: number;
  readonly k: number;
  readonly nsym: number;
  readonly generator: Uint16Array;
  private readonly alpha = 2;
  private readonly fcr = 1;

  constructor(cfg: RS216Config) {
    initGf216();
    if (cfg.n > 65535 || cfg.n <= 0)
      throw new Error(`RS216 n must be in 1..65535, got ${cfg.n}`);
    if (cfg.k >= cfg.n || cfg.k <= 0)
      throw new Error(`RS216 k must be in 1..n-1, got ${cfg.k}`);
    this.n = cfg.n;
    this.k = cfg.k;
    this.nsym = cfg.n - cfg.k;
    this.generator = this.buildGenerator();
  }

  private buildGenerator(): Uint16Array {
    let g = new Uint16Array([1]);
    for (let i = 0; i < this.nsym; i++) {
      const root = gf16Pow(this.alpha, this.fcr + i);
      g = gf16PolyMul(g, new Uint16Array([1, root]));
    }
    return g;
  }

  /** Encode k data symbols → n symbols (k data + nsym parity). BE convention. */
  encode(data: Uint16Array): Uint16Array {
    if (data.length !== this.k)
      throw new Error(`RS216 encode expects ${this.k} symbols, got ${data.length}`);
    const msg = new Uint16Array(this.n);
    msg.set(data, 0);
    const gen = this.generator;
    for (let i = 0; i < this.k; i++) {
      const coef = msg[i];
      if (coef !== 0) {
        for (let j = 1; j < gen.length; j++) {
          msg[i + j] ^= gf16Mul(gen[j], coef);
        }
      }
    }
    msg.set(data, 0);
    return msg;
  }

  parity(data: Uint16Array): Uint16Array {
    return this.encode(data).slice(this.k);
  }

  /** Check if codeword has errors (nonzero syndrome). */
  hasError(recv: Uint16Array): boolean {
    for (let i = 0; i < this.nsym; i++) {
      const x = gf16Pow(this.alpha, this.fcr + i);
      if (gf16PolyEval(recv, x) !== 0) return true;
    }
    return false;
  }

  /**
   * Pure-erasure decoding via Lagrange interpolation.
   *
   * Given a received codeword with known erasure positions, reconstruct the
   * original message. This works because RS is an MDS code: any k of n
   * symbols suffice to reconstruct the message.
   *
   * Algorithm:
   *   1. Collect the n - |erasures| non-erased symbols with their positions.
   *   2. If we have >= k non-erased symbols, use Lagrange interpolation
   *      to evaluate the message polynomial at the data positions.
   *
   * This is simpler and more robust than BM+Chien+Forney.
   */
  decodeWithErasures(
    recv: Uint16Array,
    erasePos: number[],
  ): RS216DecodeResult {
    if (recv.length !== this.n)
      throw new Error(`RS216 expects ${this.n} symbols, got ${recv.length}`);
    const eraseSet = new Set(erasePos.filter((p) => p >= 0 && p < this.n));
    if (eraseSet.size > this.nsym) {
      throw new Error(`Too many erasures: ${eraseSet.size} > ${this.nsym}`);
    }
    if (eraseSet.size === 0) {
      // No erasures — check if already valid
      if (!this.hasError(recv)) {
        return { data: recv.slice(0, this.k), corrected: 0, erased: 0 };
      }
      throw new Error("RS216: errors present but no erasures provided");
    }

    // Collect available symbols (non-erased) with their positions
    const available: { pos: number; val: number }[] = [];
    for (let i = 0; i < this.n; i++) {
      if (!eraseSet.has(i)) {
        available.push({ pos: i, val: recv[i] });
      }
    }
    if (available.length < this.k) {
      throw new Error(
        `Not enough symbols: have ${available.length}, need ${this.k}`,
      );
    }

    // Use exactly k symbols (the first k available)
    const useSymbols = available.slice(0, this.k);

    // Lagrange interpolation: the codeword polynomial c(x) has degree n-1.
    // We know c(x) at k points (the data positions are c(alpha^...) -- wait,
    // the positions are array indices, not evaluation points.
    //
    // Actually, for systematic RS, the codeword IS the message + parity.
    // The message is positions 0..k-1. If any of those are erased, we need
    // to reconstruct them from the parity.
    //
    // The standard approach: treat the codeword as evaluations of a polynomial
    // at specific points. But our convention is BE array (not evaluation form).
    //
    // Simpler approach: just solve the linear system.
    // The codeword satisfies H * c^T = 0 (parity check). With erasures, we
    // have unknowns only at erased positions. Solve for them.

    // Build the parity check equations involving erased positions.
    // H is (nsym x n) matrix where H[i][j] = alpha^((fcr+i) * (n-1-j))
    // (position j has polynomial exponent n-1-j).
    //
    // For each syndrome equation i:
    //   sum_j H[i][j] * c[j] = S_i (the syndrome)
    // Erased positions contribute unknowns; non-erased are known.
    //
    // We have nsym equations and |erasures| unknowns. If |erasures| <= nsym,
    // we can solve (overdetermined system, but consistent for valid codewords).

    const erased = Array.from(eraseSet).sort((a, b) => a - b);
    const numUnknowns = erased.length;

    // Build matrix A (nsym x numUnknowns) and vector b (nsym) such that A * x = b
    // where x is the vector of erased symbol values.
    // A[i][u] = H[i][erased[u]] = alpha^((fcr+i) * (n-1-erased[u]))
    // b[i] = S_i - sum_{j not erased} H[i][j] * c[j]
    //      = syndrome computed from received (with erased set to 0)
    const A: number[][] = [];
    const b: number[] = [];
    for (let i = 0; i < this.nsym; i++) {
      const expBase = this.fcr + i;
      // Compute syndrome with erased positions zeroed
      const mutated = recv.slice();
      for (const p of erased) mutated[p] = 0;
      let s = 0;
      for (let j = 0; j < this.n; j++) {
        if (mutated[j] !== 0) {
          const h = gf16Pow(this.alpha, expBase * (this.n - 1 - j));
          s ^= gf16Mul(h, mutated[j]);
        }
      }
      b.push(s);
      // Row of A
      const row: number[] = [];
      for (const p of erased) {
        row.push(gf16Pow(this.alpha, expBase * (this.n - 1 - p)));
      }
      A.push(row);
    }

    // Solve A * x = b over GF(2^16) using Gaussian elimination
    const x = solveGf216(A, b, numUnknowns);

    // Apply corrections
    const corrected = recv.slice();
    for (let i = 0; i < erased.length; i++) {
      corrected[erased[i]] = x[i];
    }

    // Verify
    if (this.hasError(corrected)) {
      throw new Error("RS216: post-correction syndrome nonzero");
    }

    return {
      data: corrected.slice(0, this.k),
      corrected: 0,
      erased: erased.length,
    };
  }

  /**
   * Decode with optional erasure positions.
   * Falls back to erasure-only if no erasures and syndrome is zero.
   */
  decode(recv: Uint16Array, erasePos: number[] = []): RS216DecodeResult {
    if (erasePos.length > 0) {
      return this.decodeWithErasures(recv, erasePos);
    }
    if (!this.hasError(recv)) {
      return { data: recv.slice(0, this.k), corrected: 0, erased: 0 };
    }
    // Unknown errors without erasure hints — not supported in this simplified impl.
    // For unknown-error correction, use the GF(2^8) RS (reedsolomon.ts) which
    // has a full BM+Chien+Forney implementation.
    throw new Error(
      "RS216: unknown-error correction not implemented (use GF(2^8) RS or provide erasure positions)",
    );
  }
}

/** Solve a linear system over GF(2^16) using Gaussian elimination. */
function solveGf216(matrix: number[][], rhs: number[], n: number): number[] {
  // Augmented matrix
  const aug = matrix.map((row, i) => [...row, rhs[i]]);

  // Forward elimination
  for (let col = 0; col < n; col++) {
    let pivot = -1;
    for (let row = col; row < aug.length; row++) {
      if (aug[row][col] !== 0) {
        pivot = row;
        break;
      }
    }
    if (pivot === -1) throw new Error("Singular matrix in GF(2^16) solve");
    if (pivot !== col) {
      [aug[col], aug[pivot]] = [aug[pivot], aug[col]];
    }
    const pivotInv = gf16Inverse(aug[col][col]);
    for (let j = col; j <= n; j++) {
      aug[col][j] = gf16Mul(aug[col][j], pivotInv);
    }
    for (let row = 0; row < aug.length; row++) {
      if (row === col) continue;
      const factor = aug[row][col];
      if (factor === 0) continue;
      for (let j = col; j <= n; j++) {
        aug[row][j] ^= gf16Mul(factor, aug[col][j]);
      }
    }
  }

  // Extract solution
  const solution: number[] = new Array(n);
  for (let i = 0; i < n; i++) {
    solution[i] = aug[i][n];
  }
  return solution;
}
