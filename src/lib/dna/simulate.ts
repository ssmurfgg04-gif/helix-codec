/**
 * Mutation / sequencing error simulation for DNA storage.
 *
 * Models three primary error modes observed in real DNA synthesis and sequencing:
 *   1. Substitutions: a base is replaced by another (A->C, etc.).
 *   2. Insertions: an extra base is added at a position.
 *   3. Deletions: a base is removed.
 *
 * Each error type has an independent per-position rate. The simulator also
 * models sequencing COVERAGE: each oligo is "read" multiple times (coverage
 * depth), with independent errors per read. The decoder then clusters reads
 * by index and takes a consensus to recover the original.
 *
 * Realistic defaults (per Illumina DNA-storage runs, Chandak 2018):
 *   - Substitution: 1e-3 per position
 *   - Insertion:     5e-4 per position
 *   - Deletion:      1e-3 per position
 *   - Coverage:      20x
 *
 * For Nanopore (ONT): indels dominate, ~5-10% total. We provide a preset.
 *
 * Reference:
 *   - Chandak et al. (2018), "Improved read/write cost tradeoff in a DNA-based
 *     storage system using rotating fountain codes", ISMB.
 *   - Organick et al. (2018), "Random access in large-scale DNA data storage",
 *     Nature Biotechnology 36:242-248.
 */

import { Oligo } from "./types";

export interface MutationConfig {
  /** Per-position substitution rate (0..1). */
  substitutionRate: number;
  /** Per-position insertion rate (0..1). */
  insertionRate: number;
  /** Per-position deletion rate (0..1). */
  deletionRate: number;
  /** Sequencing coverage depth (each oligo read N times). */
  coverage: number;
  /** Fraction of oligos that are completely lost (synthesis failure). */
  dropoutRate: number;
  /** Random seed for reproducibility (0 = non-deterministic). */
  seed: number;
}

export const PRESET_ILLUMINA: MutationConfig = {
  // Chandak 2018 (ISMB) — Illumina DNA-storage run, slightly conservative.
  substitutionRate: 0.001,
  insertionRate: 0.0005,
  deletionRate: 0.001,
  coverage: 20,
  dropoutRate: 0.0,
  seed: 0,
};

/**
 * PRESET_REAL_2024 — measured real-world rates from Preuss et al. (2026, Nature Sci Rep).
 * Deletions dominate by ~5x over substitutions. This is the most realistic preset.
 * Per-nt rates: del 0.082, sub 0.025, ins 0.016.
 */
export const PRESET_REAL_2024: MutationConfig = {
  substitutionRate: 0.025,
  insertionRate: 0.016,
  deletionRate: 0.082,
  coverage: 25,
  dropoutRate: 0.02,
  seed: 0,
};

/**
 * PRESET_BANAL_2026 — soft-information-aware simulator matching Banal et al. (2026).
 * Same underlying error rates as REAL_2024 but with higher coverage to model
 * the profile-HMM + ordered-statistics decoding pipeline.
 */
export const PRESET_BANAL_2026: MutationConfig = {
  substitutionRate: 0.025,
  insertionRate: 0.016,
  deletionRate: 0.082,
  coverage: 30,
  dropoutRate: 0.01,
  seed: 0,
};

export const PRESET_NANOPORE: MutationConfig = {
  // ONT R10.4.1 — indel-heavy, ~9% total raw error
  substitutionRate: 0.02,
  insertionRate: 0.03,
  deletionRate: 0.04,
  coverage: 15,
  dropoutRate: 0.05,
  seed: 0,
};

export const PRESET_PACBIO: MutationConfig = {
  // PacBio HiFi (CCS-corrected) — insertion-dominant
  substitutionRate: 0.005,
  insertionRate: 0.05,
  deletionRate: 0.03,
  coverage: 10,
  dropoutRate: 0.02,
  seed: 0,
};

export const PRESET_CLEAN: MutationConfig = {
  substitutionRate: 0,
  insertionRate: 0,
  deletionRate: 0,
  coverage: 1,
  dropoutRate: 0,
  seed: 0,
};

/** A single noisy read of an oligo. */
export interface SequencingRead {
  /** Index of the source oligo (extracted from address). */
  oligoIndex: number;
  /** Noisy DNA sequence (with substitutions, insertions, deletions). */
  sequence: string;
  /** Per-base Phred quality scores (Q10 = 90% accuracy, Q20 = 99%, Q30 = 99.9%). */
  quality?: Uint8Array;
  /** Number of substitutions applied. */
  substitutions: number;
  /** Number of insertions applied. */
  insertions: number;
  /** Number of deletions applied. */
  deletions: number;
}

/** Result of simulating sequencing on an encoded file. */
export interface SimulationResult {
  /** All reads across all oligos (length = oligoCount * coverage). */
  reads: SequencingRead[];
  /** Indices of oligos that were dropped (synthesis failure). */
  droppedOligos: number[];
  /** Per-oligo error stats (averaged across reads). */
  perOligoStats: {
    index: number;
    avgSubstitutions: number;
    avgInsertions: number;
    avgDeletions: number;
    readCount: number;
  }[];
  /** Aggregate stats. */
  totalReads: number;
  avgReadLength: number;
  totalErrors: number;
  simulationTimeMs: number;
}

// Xorshift32 PRNG (deterministic if seed != 0)
class Rng {
  private state: number;
  constructor(seed: number) {
    this.state = (seed >>> 0) || 1;
  }
  next(): number {
    this.state ^= this.state << 13;
    this.state ^= this.state >>> 17;
    this.state ^= this.state << 5;
    this.state = this.state >>> 0;
    return this.state / 0x100000000;
  }
  nextInt(maxExclusive: number): number {
    return Math.floor(this.next() * maxExclusive);
  }
}

const BASES = "ACGT";

/**
 * Apply a single substitution at a random position.
 * @deprecated O(n²) due to string slicing. Use simulateRead() which uses array-based mutation.
 */
function applySubstitution(dna: string, rng: Rng): { result: string; pos: number } {
  if (dna.length === 0) return { result: dna, pos: -1 };
  const pos = rng.nextInt(dna.length);
  const original = dna[pos];
  let newBase: string;
  do {
    newBase = BASES[rng.nextInt(4)];
  } while (newBase === original);
  return {
    result: dna.slice(0, pos) + newBase + dna.slice(pos + 1),
    pos,
  };
}

/**
 * Apply a single insertion at a random position.
 * @deprecated O(n²) due to string slicing. Use simulateRead() which uses array-based mutation.
 */
function applyInsertion(dna: string, rng: Rng): { result: string; pos: number } {
  if (dna.length === 0) {
    const b = BASES[rng.nextInt(4)];
    return { result: b, pos: 0 };
  }
  const pos = rng.nextInt(dna.length + 1); // can insert at end
  const newBase = BASES[rng.nextInt(4)];
  return {
    result: dna.slice(0, pos) + newBase + dna.slice(pos),
    pos,
  };
}

/**
 * Apply a single deletion at a random position.
 * @deprecated O(n²) due to string slicing. Use simulateRead() which uses array-based mutation.
 */
function applyDeletion(dna: string, rng: Rng): { result: string; pos: number } {
  if (dna.length === 0) return { result: dna, pos: -1 };
  const pos = rng.nextInt(dna.length);
  return {
    result: dna.slice(0, pos) + dna.slice(pos + 1),
    pos,
  };
}

/**
 * Simulate sequencing errors on a single oligo.
 * Returns one noisy read WITH per-base Phred quality scores.
 *
 * Q-score model:
 *   - Correct base: Q30-Q40 (99.9%-99.99% confidence)
 *   - Substituted base: Q5-Q15 (low confidence — should be treated as erasure)
 *   - Inserted base: Q2-Q8 (very low confidence — strong erasure signal)
 *   - Deleted base: (no Q-score — base is absent)
 *
 * The Q-scores allow the decoder to use soft-information: low-Q positions are
 * marked as erasures for RS decoding, doubling correction capacity.
 * Reference: Banal et al. 2026 (arXiv:2604.20810) — soft-info passthrough.
 */
function simulateRead(
  oligo: Oligo,
  cfg: MutationConfig,
  rng: Rng,
): SequencingRead {
  const original = oligo.sequence;
  // Pre-allocate arrays for performance (avoid O(n²) string concatenation)
  const resultParts: string[] = [];
  const qualities: number[] = [];
  let subs = 0;
  let ins = 0;
  let dels = 0;

  // First pass: decide deletions (mark which original positions survive)
  const survived = new Array<boolean>(original.length).fill(true);
  if (cfg.deletionRate > 0) {
    for (let i = 0; i < original.length; i++) {
      if (rng.next() < cfg.deletionRate) {
        survived[i] = false;
        dels++;
      }
    }
  }

  // Second pass: walk through original, emit surviving bases (with possible
  // substitution and insertion), tracking Q-scores.
  for (let i = 0; i < original.length; i++) {
    if (!survived[i]) continue;

    const origBase = original[i];
    let emitBase = origBase;
    let qScore: number;

    // Substitution?
    if (rng.next() < cfg.substitutionRate) {
      let newBase: string;
      do {
        newBase = BASES[rng.nextInt(4)];
      } while (newBase === origBase);
      emitBase = newBase;
      // Low Q-score for substituted base: Q5-Q15
      qScore = 5 + rng.nextInt(11); // 5..15
      subs++;
    } else {
      // Correct base: Q30-Q40
      qScore = 30 + rng.nextInt(11); // 30..40
    }

    resultParts.push(emitBase);
    qualities.push(qScore);

    // Insertion after this position?
    if (rng.next() < cfg.insertionRate) {
      const insBase = BASES[rng.nextInt(4)];
      resultParts.push(insBase);
      // Very low Q-score for inserted base: Q2-Q8
      qualities.push(2 + rng.nextInt(7)); // 2..8
      ins++;
    }
  }

  // Possible insertion at the very end
  if (rng.next() < cfg.insertionRate) {
    resultParts.push(BASES[rng.nextInt(4)]);
    qualities.push(2 + rng.nextInt(7));
    ins++;
  }

  return {
    oligoIndex: oligo.index,
    sequence: resultParts.join(""),
    quality: new Uint8Array(qualities),
    substitutions: subs,
    insertions: ins,
    deletions: dels,
  };
}

/**
 * Simulate full sequencing of an encoded file.
 */
export function simulate(
  oligos: Oligo[],
  cfg: MutationConfig,
): SimulationResult {
  const t0 = Date.now();
  const rng = new Rng(cfg.seed || Math.floor(Math.random() * 0xffffffff));
  const reads: SequencingRead[] = [];
  const droppedOligos: number[] = [];
  const perOligoStats: SimulationResult["perOligoStats"] = [];

  for (const oligo of oligos) {
    // Synthesis dropout
    if (rng.next() < cfg.dropoutRate) {
      droppedOligos.push(oligo.index);
      perOligoStats.push({
        index: oligo.index,
        avgSubstitutions: 0,
        avgInsertions: 0,
        avgDeletions: 0,
        readCount: 0,
      });
      continue;
    }

    let totalSubs = 0;
    let totalIns = 0;
    let totalDels = 0;
    for (let r = 0; r < cfg.coverage; r++) {
      const read = simulateRead(oligo, cfg, rng);
      reads.push(read);
      totalSubs += read.substitutions;
      totalIns += read.insertions;
      totalDels += read.deletions;
    }
    perOligoStats.push({
      index: oligo.index,
      avgSubstitutions: totalSubs / cfg.coverage,
      avgInsertions: totalIns / cfg.coverage,
      avgDeletions: totalDels / cfg.coverage,
      readCount: cfg.coverage,
    });
  }

  const totalReads = reads.length;
  const avgReadLength =
    totalReads === 0 ? 0 : reads.reduce((s, r) => s + r.sequence.length, 0) / totalReads;
  const totalErrors = reads.reduce((s, r) => s + r.substitutions + r.insertions + r.deletions, 0);

  return {
    reads,
    droppedOligos,
    perOligoStats,
    totalReads,
    avgReadLength,
    totalErrors,
    simulationTimeMs: Date.now() - t0,
  };
}

/**
 * Compare two DNA strings position-by-position, returning a diff for visualization.
 * Uses simple alignment (no gap handling — for visualization only).
 * @deprecated Broken for indels — a single insertion at position 0 marks every
 * subsequent base as a substitution. Use Needleman-Wunsch alignment instead.
 */
export interface PositionDiff {
  position: number;
  original: string;
  mutated: string;
  type: "sub" | "ins" | "del";
}

export function diffSequences(original: string, mutated: string): PositionDiff[] {
  const diffs: PositionDiff[] = [];
  const maxLen = Math.max(original.length, mutated.length);
  for (let i = 0; i < maxLen; i++) {
    const o = original[i] ?? "";
    const m = mutated[i] ?? "";
    if (o === m) continue;
    if (o === "") diffs.push({ position: i, original: o, mutated: m, type: "ins" });
    else if (m === "") diffs.push({ position: i, original: o, mutated: m, type: "del" });
    else diffs.push({ position: i, original: o, mutated: m, type: "sub" });
  }
  return diffs;
}
