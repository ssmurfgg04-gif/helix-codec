/**
 * PEG (Progressive Edge Growth) LDPC Construction
 *
 * Constructs a sparse parity-check matrix H for LDPC codes with good girth
 * (no short cycles). Used by the Mahoraga codec (Banal 2026) for its inner
 * LDPC code at n=252 bits.
 *
 * Algorithm:
 *   For each variable node (column), for each check node (row) to connect:
 *     1. BFS from the variable node to find check nodes at increasing distance.
 *     2. Pick the check node with the lowest degree among those at maximum
 *        distance (to spread edges and avoid short cycles).
 *     3. Add the edge.
 *
 * This is a simplified PEG that doesn't do the full BFS but instead picks
 * the lowest-degree check node (a common heuristic that works well in practice).
 *
 * Reference:
 *   - Hu, Eleftheriou, Arnold (2005). "Regular and irregular progressive
 *     edge-growth Tanner graphs." IEEE TIT 51:1.
 *   - Banal et al. (2026). arXiv:2604.20810.
 */

import { GF2Matrix } from "./osd";

export interface PEGConfig {
  /** Block length n (number of variable nodes / columns). */
  n: number;
  /** Number of check nodes (rows). */
  m: number;
  /** Variable node degree dv (edges per column). */
  dv: number;
  /** Check node degree dc (edges per row). Target, not strict. */
  dc: number;
}

/**
 * Construct an LDPC parity-check matrix using simplified PEG with
 * column-uniqueness enforcement.
 *
 * Two improvements over the naive "pick lowest-degree check" heuristic:
 *
 * 1. GLOBAL minimum-degree selection: scan ALL check nodes to find the one
 *    with the lowest degree (not the first one below dc). This distributes
 *    edges evenly.
 *
 * 2. Column-uniqueness: track the set of checks already used for each prior
 *    column (as a sorted string key). When building a new column, if the
 *    chosen set of dv checks would duplicate an existing column's pattern,
 *    swap one check for the next-best alternative. This eliminates stopping
 *    sets caused by duplicate columns — the dominant failure mode of the
 *    peeling + Gaussian-elimination LDPC erasure decoder.
 *
 * This is critical for the LDPC erasure decoder (`decodeWithErasures`) to
 * achieve its theoretical BEC capacity. Without column uniqueness, two
 * identical columns share all their checks, making both bits unrecoverable
 * when simultaneously erased (a rank-1 stopping set of size 2).
 *
 * Reference:
 *   - Hu, Eleftheriou, Arnold (2005). "Regular and irregular PEG Tanner
 *     graphs." IEEE TIT 51:1. (Original PEG with girth maximization.)
 *   - Di, Proietti, Telatar, Richardson, Urbanke (2002). "Finite-length
 *     analysis on the BEC." IEEE TIT 48. (Stopping sets.)
 *
 * @returns H as a GF2Matrix (m rows × n columns)
 */
export function constructPEG(cfg: PEGConfig): GF2Matrix {
  const { n, m, dv, dc } = cfg;
  const H = new GF2Matrix(m, n);

  // Track degree of each check node (row)
  const checkDegree = new Uint32Array(m);

  // Track existing column signatures to prevent duplicates
  const existingSignatures = new Set<string>();

  // For each variable node (column), add dv edges
  for (let col = 0; col < n; col++) {
    const usedChecks = new Set<number>(); // don't add duplicate edges within this column

    // Try up to (dv + 4) rounds to find a unique column signature
    let attempts = 0;
    const maxAttempts = dv * 4;

    while (attempts < maxAttempts) {
      attempts++;
      usedChecks.clear();

      // Build this column's edges using greedy lowest-degree selection
      // with a small randomization on ties (to escape duplicate patterns)
      const candidateChecks: number[] = [];
      for (let edge = 0; edge < dv; edge++) {
        // Find the check node with GLOBAL minimum degree
        let bestCheck = -1;
        let bestDegree = Infinity;
        const tiedChecks: number[] = [];

        for (let row = 0; row < m; row++) {
          if (usedChecks.has(row)) continue;
          if (checkDegree[row] >= dc + 2) continue; // don't over-fill
          if (checkDegree[row] < bestDegree) {
            bestDegree = checkDegree[row];
            bestCheck = row;
            tiedChecks.length = 0;
            tiedChecks.push(row);
          } else if (checkDegree[row] === bestDegree) {
            tiedChecks.push(row);
          }
        }

        if (bestCheck === -1) break;

        // If multiple ties, pick pseudo-randomly based on col + edge + attempts
        // (deterministic but varies across attempts to escape duplicates)
        const chosen = tiedChecks.length > 1
          ? tiedChecks[(col * 7 + edge * 13 + attempts * 3) % tiedChecks.length]
          : bestCheck;

        usedChecks.add(chosen);
        candidateChecks.push(chosen);
      }

      if (candidateChecks.length !== dv) continue; // retry

      // Check signature uniqueness
      const sig = Array.from(usedChecks).sort((a, b) => a - b).join(",");
      if (!existingSignatures.has(sig)) {
        // Accept this column
        for (const row of candidateChecks) {
          H.set(row, col, 1);
          checkDegree[row]++;
        }
        existingSignatures.add(sig);
        break;
      }
      // Otherwise retry with different randomization
    }

    // If we exhausted attempts, just place edges without uniqueness check
    // (better than leaving column empty)
    let colHasEdge = false;
    for (let r = 0; r < m; r++) {
      if (H.get(r, col) === 1) { colHasEdge = true; break; }
    }
    if (attempts >= maxAttempts && !colHasEdge) {
      for (let edge = 0; edge < dv; edge++) {
        let bestCheck = -1;
        let bestDegree = Infinity;
        for (let row = 0; row < m; row++) {
          if (usedChecks.has(row)) continue;
          if (checkDegree[row] < bestDegree) {
            bestDegree = checkDegree[row];
            bestCheck = row;
          }
        }
        if (bestCheck === -1) break;
        H.set(bestCheck, col, 1);
        checkDegree[bestCheck]++;
        usedChecks.add(bestCheck);
      }
    }
  }

  return H;
}

/**
 * Generate the Mahoraga hi-fi LDPC matrix.
 *
 * Parameters (from arXiv:2604.20810):
 *   n = 252 bits (126 nt × 2 bits/nt)
 *   dv = 3 (variable node degree)
 *   dc = 84 (check node degree)
 *   k = 243 information bits (after 32-bit CRC → 208-bit user payload)
 *   m = 9 base parity checks + redundant rows for BP enrichment
 *
 * Returns H with (n-k) = 9 rows for the base code.
 * (Mahoraga adds extra redundant rows; we omit those for simplicity.)
 */
export function mahoragaHiFiMatrix(): GF2Matrix {
  const n = 252;
  const k = 243;
  const m = n - k; // 9
  const dv = 3;
  const dc = 84;
  return constructPEG({ n, m, dv, dc });
}

/**
 * Generate the Mahoraga lo-fi LDPC matrix.
 *
 * Parameters (from arXiv:2604.20810):
 *   n = 252 bits
 *   dv = 3
 *   dc = 21
 *   k = 216 (after 32-bit CRC → 176-bit user payload)
 *   m = 36 base parity checks
 */
export function mahoragaLoFiMatrix(): GF2Matrix {
  const n = 252;
  const k = 216;
  const m = n - k; // 36
  const dv = 3;
  const dc = 21;
  return constructPEG({ n, m, dv, dc });
}

/**
 * Compute the rate of an LDPC code: k/n where k = n - rank(H).
 * For a systematic code, rank(H) = m (number of independent rows).
 */
export function ldpcRate(H: GF2Matrix): number {
  const { matrix, pivotCols } = H.systematicForm();
  const rank = pivotCols.length;
  return (H.cols - rank) / H.cols;
}
