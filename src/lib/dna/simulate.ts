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
import {
  simulateWetlab,
  WetlabConfig,
  WetlabResult,
  ILLUMINA_WETLAB,
  NANOPORE_WETLAB,
  PACBIO_WETLAB,
} from "./dt4dds-simulate";
// v69: napi-rs native simulation hot path (FIRST PRIORITY for bulk read simulation)
import { getNativeAddon } from "./native/helix-napi";

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
  /**
   * Simulator backend:
   *   - "basic"  — Simple uniform per-position model (this file)
   *   - "dt4dds" — Parametric wet-lab pipeline (dt4dds-simulate.ts)
   * Default: "dt4dds" (since v3.0)
   */
  simulator?: "basic" | "dt4dds";
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
  simulator: "basic", // Clean preset uses basic for fast, deterministic test paths
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

// ---------------------------------------------------------------------------
// Bridge: MutationConfig → WetlabConfig (dt4dds)
// ---------------------------------------------------------------------------

/**
 * Convert a simple MutationConfig to a WetlabConfig for the dt4dds pipeline.
 *
 * Maps the flat error rates to the full parametric model:
 *   - Synthesis: uses substitution/insertion/deletion rates with position-dependent scaling
 *   - PCR: defaults (15 cycles, 0.85 duplication, 0.0001 sub rate, 0.1 GC bias)
 *   - Aging: no aging by default (days = 0)
 *   - Sequencing: uses coverage/dropout + platform-appropriate error rates
 *
 * The channel is inferred from the error profile:
 *   - deletionRate > 0.02 → nanopore
 *   - insertionRate > 0.02 → pacbio
 *   - otherwise → illumina
 */
export function mutationConfigToWetlabConfig(cfg: MutationConfig): WetlabConfig {
  // Infer platform from error profile
  let platform: "illumina" | "nanopore" | "pacbio";
  let base: WetlabConfig;
  if (cfg.deletionRate > 0.02 || cfg.insertionRate > 0.03) {
    platform = "nanopore";
    base = NANOPORE_WETLAB;
  } else if (cfg.insertionRate > 0.02) {
    platform = "pacbio";
    base = PACBIO_WETLAB;
  } else {
    platform = "illumina";
    base = ILLUMINA_WETLAB;
  }

  // Override synthesis rates from MutationConfig
  base = structuredClone(base);
  base.synthesis.substitutionRate = cfg.substitutionRate;
  base.synthesis.insertionRate = cfg.insertionRate;
  base.synthesis.deletionRate = cfg.deletionRate;

  // Override sequencing rates + coverage + dropout from MutationConfig
  base.sequencing.platform = platform;
  base.sequencing.substitutionRate = cfg.substitutionRate;
  base.sequencing.insertionRate = cfg.insertionRate;
  base.sequencing.deletionRate = cfg.deletionRate;
  base.sequencing.coverage = cfg.coverage;
  base.sequencing.dropoutRate = cfg.dropoutRate;

  return base;
}

/**
 * Convert a WetlabResult from dt4dds into a SimulationResult for downstream consumers.
 */
export function wetlabResultToSimulationResult(
  result: WetlabResult,
  simulationTimeMs: number,
): SimulationResult {
  const { reads, stats } = result;
  const droppedOligos: number[] = [];

  // Identify dropped oligos (those in the original set with no reads)
  const readIndices = new Set(reads.map((r) => r.oligoIndex));
  // We don't have the original oligo list here, so we derive dropped from stats
  // The dt4dds pipeline already handles dropout; we leave droppedOligos empty
  // since the reads already exclude dropped oligos.

  // Compute per-oligo stats
  const byIndex = new Map<number, SequencingRead[]>();
  for (const read of reads) {
    const arr = byIndex.get(read.oligoIndex) ?? [];
    arr.push(read);
    byIndex.set(read.oligoIndex, arr);
  }

  const perOligoStats: SimulationResult["perOligoStats"] = [];
  const entries = Array.from(byIndex.entries());
  for (const [index, oligoReads] of entries) {
    const avgSubs = oligoReads.reduce((s, r) => s + r.substitutions, 0) / oligoReads.length;
    const avgIns = oligoReads.reduce((s, r) => s + r.insertions, 0) / oligoReads.length;
    const avgDels = oligoReads.reduce((s, r) => s + r.deletions, 0) / oligoReads.length;
    perOligoStats.push({
      index,
      avgSubstitutions: avgSubs,
      avgInsertions: avgIns,
      avgDeletions: avgDels,
      readCount: oligoReads.length,
    });
  }

  const totalReads = reads.length;
  const avgReadLength =
    totalReads === 0 ? 0 : reads.reduce((s, r) => s + r.sequence.length, 0) / totalReads;
  const totalErrors = reads.reduce(
    (s, r) => s + r.substitutions + r.insertions + r.deletions,
    0,
  );

  return {
    reads,
    droppedOligos,
    perOligoStats,
    totalReads,
    avgReadLength,
    totalErrors,
    simulationTimeMs,
  };
}

/**
 * Simulate full sequencing of an encoded file using the dt4dds parametric pipeline.
 */
export function simulateDt4dds(
  oligos: Oligo[],
  cfg: MutationConfig,
): SimulationResult {
  const t0 = Date.now();
  const wetlabCfg = mutationConfigToWetlabConfig(cfg);
  const result = simulateWetlab(oligos, wetlabCfg, cfg.seed);
  return wetlabResultToSimulationResult(result, Date.now() - t0);
}

/**
 * Simulate full sequencing of an encoded file.
 *
 * Delegates to the dt4dds parametric pipeline by default (since v3.0).
 * Set `cfg.simulator = "basic"` to use the simple uniform model.
 */
export function simulate(
  oligos: Oligo[],
  cfg: MutationConfig,
): SimulationResult {
  const useDt4dds = cfg.simulator !== "basic"; // default is dt4dds
  if (useDt4dds) {
    return simulateDt4dds(oligos, cfg);
  }
  const t0 = Date.now();
  const rng = new Rng(cfg.seed || Math.floor(Math.random() * 0xffffffff));
  const reads: SequencingRead[] = [];
  const droppedOligos: number[] = [];
  const perOligoStats: SimulationResult["perOligoStats"] = [];

  // v69: napi-rs native simulation FIRST PRIORITY (bulk path).
  // For each oligo, call the native simulate_oligo_reads which returns a flat
  // array of [coverage_u32, r0_len_u32, r0_bytes, r1_len_u32, r1_bytes, ...]
  // We parse it into SequencingRead objects on the JS side.
  const nativeAddon = getNativeAddon();
  if (nativeAddon && oligos.length > 0) {
    for (const oligo of oligos) {
      // Synthesis dropout (use the same RNG as the JS path)
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

      try {
        // Strip primers (simulate only the inner oligo region)
        const seq = oligo.sequence;
        const nativeResult = nativeAddon.simulateOligoReads(seq, {
          substitutionRate: cfg.substitutionRate,
          insertionRate: cfg.insertionRate,
          deletionRate: cfg.deletionRate,
          coverage: cfg.coverage,
          dropoutRate: 0, // already handled above
          seed: cfg.seed || 0,
          positionDependent: false, // basic simulator: uniform rates
          fivePrimeMult: 1.0,
          threePrimeMult: 1.0,
        });

        // Parse flat output: [coverage_u32, r0_len_u32, r0_bytes, ...]
        const view = new DataView(nativeResult.buffer, nativeResult.byteOffset, nativeResult.byteLength);
        let pos = 0;
        if (pos + 4 > nativeResult.byteLength) continue;
        const actualCoverage = view.getUint32(pos, true); // little-endian
        pos += 4;

        let totalSubs = 0;
        let totalIns = 0;
        let totalDels = 0;
        for (let r = 0; r < actualCoverage; r++) {
          if (pos + 4 > nativeResult.byteLength) break;
          const rlen = view.getUint32(pos, true);
          pos += 4;
          if (pos + rlen > nativeResult.byteLength) break;
          // Decode bytes back to string
          let readSeq = '';
          for (let i = 0; i < rlen; i++) {
            readSeq += String.fromCharCode(nativeResult[pos + i]);
          }
          pos += rlen;
          // Estimate error counts from length delta (best-effort)
          const origLen = seq.length;
          const newLen = readSeq.length;
          const ins = Math.max(0, newLen - origLen);
          const dels = Math.max(0, origLen - newLen);
          let subs = 0;
          const minLen = Math.min(origLen, newLen);
          for (let i = 0; i < minLen; i++) {
            if (readSeq[i] !== seq[i]) subs++;
          }
          totalSubs += subs;
          totalIns += ins;
          totalDels += dels;
          reads.push({
            oligoIndex: oligo.index,
            readIndex: r,
            sequence: readSeq,
            quality: new Uint8Array(rlen).fill(20),
            substitutions: subs,
            insertions: ins,
            deletions: dels,
          } as any);
        }
        perOligoStats.push({
          index: oligo.index,
          avgSubstitutions: totalSubs / Math.max(1, actualCoverage),
          avgInsertions: totalIns / Math.max(1, actualCoverage),
          avgDeletions: totalDels / Math.max(1, actualCoverage),
          readCount: actualCoverage,
        });
      } catch {
        // Native sim failed for this oligo — fall back to JS simulateRead
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

  // JS fallback (original path)
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
} // end basic simulator path

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
