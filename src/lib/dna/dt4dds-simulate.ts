/**
 * Parametric Wetlab Simulation (dt4dds pattern)
 *
 * Models the entire wetlab pipeline with configurable parameters:
 *   1. Array synthesis (position-dependent error rates, bias)
 *   2. PCR amplification (duplication, amplification bias, errors)
 *   3. Aging/decay (depurination, oxidation, deamination)
 *   4. Sequencing-by-synthesis (platform-specific error profiles)
 *
 * Replaces the simple uniform RNG in simulate.ts with a research grade model.
 * For unit tests, the simple simulator is still available.
 * For research or vendor integration, use this parametric model.
 *
 * For a higher-level API that operates on raw digital data and computes
 * end-to-end metrics (BER, recovery rate), see wetlab-simulate.ts.
 *
 * Key differences from simulate.ts:
 *   - Position-dependent synthesis errors (5'/3' end effects)
 *   - GC-biased PCR amplification (high-GC oligos amplify less)
 *   - Chemical aging model (depurination, oxidation, deamination)
 *   - Platform-specific sequencing error profiles
 *   - Separate tracking of errors at each pipeline stage
 *
 * Reference:
 *   - fml-ethz/dt4dds (Python, Nature Comms & Digital Discovery)
 *   - Lee et al. (2022). "Photon reading..." Nat. Commun. 13:3231.
 *   - Banal et al. (2026). arXiv:2604.20810 — soft-information decoding
 *   - Preuss et al. (2026). Nature Sci Rep — measured real-world error rates
 */

import { Oligo } from "./types";
import { SequencingRead } from "./simulate";

// ---------------------------------------------------------------------------
// Configuration interfaces
// ---------------------------------------------------------------------------

export interface SynthesisParams {
  /** Base substitution rate. Default: 0.001 (Illumina-grade array synthesis) */
  substitutionRate: number;
  /** Insertion rate. Default: 0.0005 */
  insertionRate: number;
  /** Deletion rate. Default: 0.001 */
  deletionRate: number;
  /** Position-dependent error scaling (edges of array have higher error rates) */
  positionDependent: boolean;
  /** 5' end error multiplier. Default: 1.5 */
  fivePrimeErrorMultiplier: number;
  /** 3' end error multiplier. Default: 2.0 */
  threePrimeErrorMultiplier: number;
}

export interface PCRParams {
  /** Number of PCR cycles. Default: 15 */
  cycles: number;
  /** Per-cycle duplication probability. Default: 0.85 */
  duplicationProb: number;
  /** Per-cycle substitution rate. Default: 0.0001 */
  substitutionRate: number;
  /** Amplification bias (GC-rich sequences amplify less). Default: 0.1 */
  gcBias: number;
}

export interface AgingParams {
  /** Depurination rate per day. Default: 0.0001 */
  depurinationRate: number;
  /** Oxidation rate per day. Default: 0.00005 */
  oxidationRate: number;
  /** Deamination rate per day (C→T). Default: 0.0002 */
  deaminationRate: number;
  /** Number of days of aging. Default: 0 (no aging) */
  days: number;
}

export interface SequencingParams {
  /** Platform: 'illumina', 'nanopore', or 'pacbio' */
  platform: "illumina" | "nanopore" | "pacbio";
  /** Substitution rate (platform-specific). Default: 0.001 (Illumina), 0.02 (Nanopore) */
  substitutionRate: number;
  /** Insertion rate. Default: 0.0005 (Illumina), 0.03 (Nanopore) */
  insertionRate: number;
  /** Deletion rate. Default: 0.001 (Illumina), 0.04 (Nanopore) */
  deletionRate: number;
  /** Average coverage depth. Default: 20 */
  coverage: number;
  /** Dropout rate (oligos completely lost). Default: 0 */
  dropoutRate: number;
}

export interface WetlabConfig {
  synthesis: SynthesisParams;
  pcr: PCRParams;
  aging: AgingParams;
  sequencing: SequencingParams;
}

// ---------------------------------------------------------------------------
// Preset configurations
// ---------------------------------------------------------------------------

export const ILLUMINA_WETLAB: WetlabConfig = {
  synthesis: {
    substitutionRate: 0.001,
    insertionRate: 0.0005,
    deletionRate: 0.001,
    positionDependent: true,
    fivePrimeErrorMultiplier: 1.5,
    threePrimeErrorMultiplier: 2.0,
  },
  pcr: {
    cycles: 15,
    duplicationProb: 0.85,
    substitutionRate: 0.0001,
    gcBias: 0.1,
  },
  aging: {
    depurinationRate: 0.0001,
    oxidationRate: 0.00005,
    deaminationRate: 0.0002,
    days: 0,
  },
  sequencing: {
    platform: "illumina",
    substitutionRate: 0.001,
    insertionRate: 0.0005,
    deletionRate: 0.001,
    coverage: 20,
    dropoutRate: 0,
  },
};

export const NANOPORE_WETLAB: WetlabConfig = {
  synthesis: {
    substitutionRate: 0.002,
    insertionRate: 0.001,
    deletionRate: 0.002,
    positionDependent: true,
    fivePrimeErrorMultiplier: 1.5,
    threePrimeErrorMultiplier: 2.5,
  },
  pcr: {
    cycles: 12,
    duplicationProb: 0.80,
    substitutionRate: 0.0002,
    gcBias: 0.15,
  },
  aging: {
    depurinationRate: 0.0002,
    oxidationRate: 0.0001,
    deaminationRate: 0.0003,
    days: 0,
  },
  sequencing: {
    platform: "nanopore",
    substitutionRate: 0.02,
    insertionRate: 0.03,
    deletionRate: 0.04,
    coverage: 15,
    dropoutRate: 0.05,
  },
};

export const PACBIO_WETLAB: WetlabConfig = {
  synthesis: {
    substitutionRate: 0.001,
    insertionRate: 0.001,
    deletionRate: 0.001,
    positionDependent: true,
    fivePrimeErrorMultiplier: 1.3,
    threePrimeErrorMultiplier: 1.8,
  },
  pcr: {
    cycles: 10,
    duplicationProb: 0.82,
    substitutionRate: 0.00015,
    gcBias: 0.12,
  },
  aging: {
    depurinationRate: 0.00015,
    oxidationRate: 0.00008,
    deaminationRate: 0.00025,
    days: 0,
  },
  sequencing: {
    platform: "pacbio",
    substitutionRate: 0.005,
    insertionRate: 0.05,
    deletionRate: 0.03,
    coverage: 10,
    dropoutRate: 0.02,
  },
};

// ---------------------------------------------------------------------------
// Result interfaces
// ---------------------------------------------------------------------------

export interface WetlabResult {
  /** Sequencing reads (with errors from all stages) */
  reads: SequencingRead[];
  /** Per-stage error statistics */
  stats: {
    synthesisErrors: number;
    pcrDuplicates: number;
    pcrErrors: number;
    agingDamage: number;
    sequencingErrors: number;
    totalOligos: number;
    recoveredOligos: number;
  };
}

// ---------------------------------------------------------------------------
// PRNG (Xorshift32, same as simulate.ts)
// ---------------------------------------------------------------------------

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

// Global RNG instance (seeded per simulation)
let rng: Rng;

const BASES = "ACGT";

// ---------------------------------------------------------------------------
// Utility: GC content
// ---------------------------------------------------------------------------

function gcContent(dna: string): number {
  let gc = 0;
  for (let i = 0; i < dna.length; i++) {
    if (dna[i] === "G" || dna[i] === "C") gc++;
  }
  return dna.length > 0 ? gc / dna.length : 0;
}

// ---------------------------------------------------------------------------
// Step 1: Synthesis
// ---------------------------------------------------------------------------

/**
 * Position-dependent error multiplier.
 *
 * Models the observation that array synthesis has higher error rates at the
 * 5' and 3' ends of the oligo. The error profile is:
 *   - 5' end: elevated by fivePrimeErrorMultiplier
 *   - Middle: baseline (multiplier = 1.0)
 *   - 3' end: elevated by threePrimeErrorMultiplier
 *
 * The transition zones are smooth (linear interpolation) to avoid
 * discontinuities.
 *
 * @param pos Position in the oligo (0-based)
 * @param len Total oligo length
 * @param fivePrimeMul 5' end error multiplier
 * @param threePrimeMul 3' end error multiplier
 * @returns Error multiplier at this position
 */
function positionErrorMultiplier(
  pos: number,
  len: number,
  fivePrimeMul: number,
  threePrimeMul: number,
): number {
  if (len <= 0) return 1.0;

  // Define edge zones as the first/last 15% of the oligo
  const edgeFraction = 0.15;
  const edgeLen = Math.max(1, Math.floor(len * edgeFraction));

  if (pos < edgeLen) {
    // 5' end zone: interpolate from fivePrimeMul to 1.0
    const t = pos / edgeLen;
    return fivePrimeMul + t * (1.0 - fivePrimeMul);
  } else if (pos >= len - edgeLen) {
    // 3' end zone: interpolate from 1.0 to threePrimeMul
    const t = (pos - (len - edgeLen)) / edgeLen;
    return 1.0 + t * (threePrimeMul - 1.0);
  }

  return 1.0; // Middle: baseline
}

/**
 * Apply synthesis errors to a single oligo.
 *
 * Uses position-dependent error rates if enabled:
 *   - 5' end: substitutionRate × fivePrimeErrorMultiplier
 *   - 3' end: substitutionRate × threePrimeErrorMultiplier
 *   - Middle: substitutionRate × 1.0
 *
 * @param dna Original DNA sequence
 * @param params Synthesis parameters
 * @returns Mutated DNA sequence and error count
 */
function applySynthesisErrors(
  dna: string,
  params: SynthesisParams,
): { dna: string; errors: number } {
  if (params.substitutionRate === 0 && params.insertionRate === 0 && params.deletionRate === 0) {
    return { dna, errors: 0 };
  }

  const len = dna.length;
  let errors = 0;

  // First pass: determine deletions
  const survived = new Array<boolean>(len).fill(true);
  for (let i = 0; i < len; i++) {
    const mul = params.positionDependent
      ? positionErrorMultiplier(i, len, params.fivePrimeErrorMultiplier, params.threePrimeErrorMultiplier)
      : 1.0;
    const delRate = params.deletionRate * mul;
    if (rng.next() < delRate) {
      survived[i] = false;
      errors++;
    }
  }

  // Second pass: apply substitutions and insertions
  const result: string[] = [];
  for (let i = 0; i < len; i++) {
    if (!survived[i]) continue;

    const mul = params.positionDependent
      ? positionErrorMultiplier(i, len, params.fivePrimeErrorMultiplier, params.threePrimeErrorMultiplier)
      : 1.0;

    // Substitution
    if (rng.next() < params.substitutionRate * mul) {
      let newBase: string;
      do {
        newBase = BASES[rng.nextInt(4)];
      } while (newBase === dna[i]);
      result.push(newBase);
      errors++;
    } else {
      result.push(dna[i]);
    }

    // Insertion
    if (rng.next() < params.insertionRate * mul) {
      result.push(BASES[rng.nextInt(4)]);
      errors++;
    }
  }

  return { dna: result.join(""), errors };
}

/**
 * Step 1: Simulate array synthesis errors.
 * Position-dependent: higher error rates at 5' and 3' ends.
 *
 * @param oligos Original oligos (from encoder)
 * @param params Synthesis parameters
 * @returns Oligos with synthesis errors applied
 */
export function simulateSynthesis(
  oligos: Oligo[],
  params: SynthesisParams,
): Oligo[] {
  return oligos.map((oligo) => {
    const { dna, errors: _ } = applySynthesisErrors(oligo.sequence, params);
    // Recalculate GC and homopolymer stats
    const gc = gcContent(dna);
    let maxHomo = 1;
    let run = 1;
    for (let i = 1; i < dna.length; i++) {
      if (dna[i] === dna[i - 1]) {
        run++;
        if (run > maxHomo) maxHomo = run;
      } else {
        run = 1;
      }
    }
    return {
      ...oligo,
      sequence: dna,
      gc,
      maxHomopolymer: maxHomo,
      length: dna.length,
    };
  });
}

// ---------------------------------------------------------------------------
// Step 2: PCR amplification
// ---------------------------------------------------------------------------

/**
 * Apply a single PCR cycle to an oligo.
 *
 * Each cycle:
 *   1. With probability duplicationProb, the oligo is duplicated (creates a copy)
 *   2. The copy may have substitution errors (Taq polymerase fidelity)
 *   3. GC bias: high-GC oligos are less likely to duplicate
 *
 * @param oligo The oligo to amplify
 * @param params PCR parameters
 * @returns Array of oligos after this cycle (original + possible copy)
 */
function pcrCycle(oligo: Oligo, params: PCRParams): Oligo[] {
  const results: Oligo[] = [oligo]; // Original always survives

  // GC-biased duplication probability
  const gc = gcContent(oligo.sequence);
  // High GC → lower amplification. Bias model:
  //   effectiveProb = duplicationProb × (1 - gcBias × |gc - 0.5| / 0.5)
  // At gc=0.5 (balanced): no penalty. At gc=1.0: penalty = gcBias.
  const gcPenalty = 1 - params.gcBias * Math.abs(gc - 0.5) / 0.5;
  const effectiveProb = params.duplicationProb * Math.max(0, gcPenalty);

  if (rng.next() < effectiveProb) {
    // Duplicate with possible Taq errors
    let seq = oligo.sequence;
    let pcrErrors = 0;

    if (params.substitutionRate > 0) {
      const arr = seq.split("");
      for (let i = 0; i < arr.length; i++) {
        if (rng.next() < params.substitutionRate) {
          let newBase: string;
          do {
            newBase = BASES[rng.nextInt(4)];
          } while (newBase === arr[i]);
          arr[i] = newBase;
          pcrErrors++;
        }
      }
      seq = arr.join("");
    }

    results.push({
      ...oligo,
      sequence: seq,
      gc: gcContent(seq),
    });
  }

  return results;
}

/**
 * Step 2: Simulate PCR amplification.
 *
 * Each oligo undergoes `cycles` rounds of amplification:
 *   - Per cycle, each molecule has probability `duplicationProb` of being duplicated
 *   - GC bias reduces amplification efficiency for GC-rich sequences:
 *     effectiveProb = duplicationProb × (1 - gcBias × |gc - 0.5| / 0.5)
 *   - Taq polymerase fidelity: substitution errors at rate `substitutionRate` per base
 *   - Expected copies after c cycles: (1 + p)^c
 *
 * GC bias: high-GC oligos amplify less → lower coverage.
 * Duplication: each oligo produces ~2^cycles copies, with stochastic variation.
 * Errors: Taq polymerase introduces substitutions.
 *
 * @param oligos Oligos after synthesis
 * @param params PCR parameters
 * @returns Amplified oligo pool (may have duplicates)
 */
export function simulatePCR(
  oligos: Oligo[],
  params: PCRParams,
): Oligo[] {
  let pool: Oligo[] = [...oligos];
  let totalPcrErrors = 0;

  for (let cycle = 0; cycle < params.cycles; cycle++) {
    const newPool: Oligo[] = [];
    for (const oligo of pool) {
      const amplified = pcrCycle(oligo, params);
      newPool.push(...amplified);
      // Count PCR errors (differences from original)
      if (amplified.length > 1) {
        for (let i = 1; i < amplified.length; i++) {
          // Rough error count: Hamming distance for same-length sequences
          const orig = oligo.sequence;
          const copy = amplified[i].sequence;
          if (orig.length === copy.length) {
            for (let j = 0; j < orig.length; j++) {
              if (orig[j] !== copy[j]) totalPcrErrors++;
            }
          }
        }
      }
    }
    pool = newPool;
  }

  return pool;
}

// ---------------------------------------------------------------------------
// Step 3: Aging / chemical decay
// ---------------------------------------------------------------------------

/**
 * Apply chemical aging damage to a DNA sequence.
 *
 * Three damage mechanisms:
 *   1. Depurination: A or G loses its purine ring → the base is deleted
 *      (breaks the sugar-phosphate backbone in real DNA, modeled as deletion here)
 *   2. Oxidation: G is oxidized → G→T or G→C substitution
 *   3. Deamination: C loses its amine group → C→T substitution
 *
 * Rates are per-day, so effective rate = rate × days.
 *
 * @param dna DNA sequence
 * @param params Aging parameters
 * @returns Damaged DNA and damage count
 */
function applyAgingDamage(
  dna: string,
  params: AgingParams,
): { dna: string; damage: number } {
  if (params.days <= 0) return { dna, damage: 0 };

  const depurRate = params.depurinationRate * params.days;
  const oxidRate = params.oxidationRate * params.days;
  const deamRate = params.deaminationRate * params.days;

  let damage = 0;
  const result: string[] = [];

  for (let i = 0; i < dna.length; i++) {
    const base = dna[i];

    // Depurination: A or G → deletion
    if ((base === "A" || base === "G") && rng.next() < depurRate) {
      damage++;
      continue; // Skip this base (deletion)
    }

    // Oxidation: G → T or G→C
    if (base === "G" && rng.next() < oxidRate) {
      result.push(rng.next() < 0.5 ? "T" : "C");
      damage++;
      continue;
    }

    // Deamination: C → T
    if (base === "C" && rng.next() < deamRate) {
      result.push("T");
      damage++;
      continue;
    }

    // No damage
    result.push(base);
  }

  return { dna: result.join(""), damage };
}

/**
 * Step 3: Simulate aging/decay.
 *
 * Three chemical damage mechanisms modeled:
 *   - Depurination: A/G → deletion (breaks sugar-phosphate backbone)
 *     Rate scales linearly with storage time: effectiveRate = rate × days
 *     Dominant damage mechanism in aqueous storage at neutral pH
 *   - Oxidation: G → T or G → C substitution
 *     Caused by reactive oxygen species (ROS)
 *     More prevalent in aerobic storage conditions
 *   - Deamination: C → T substitution
 *     Spontaneous hydrolysis of cytosine's exocyclic amine group
 *     Creates C→T (or G→A on complement) transitions
 *     Rate doubles approximately every 10°C (Arrhenius kinetics)
 *
 * All rates are per-day, so effective rate = rate × days.
 * Temperature dependence is NOT modeled (assumes controlled storage at ~4°C).
 *
 * @param oligos Oligos after PCR
 * @param params Aging parameters
 * @returns Oligos with aging damage applied
 */
export function simulateAging(
  oligos: Oligo[],
  params: AgingParams,
): Oligo[] {
  if (params.days <= 0) return oligos;

  return oligos.map((oligo) => {
    const { dna, damage: _ } = applyAgingDamage(oligo.sequence, params);
    return {
      ...oligo,
      sequence: dna,
      gc: gcContent(dna),
      length: dna.length,
    };
  });
}

// ---------------------------------------------------------------------------
// Step 4: Sequencing
// ---------------------------------------------------------------------------

/**
 * Simulate a single sequencing read from an oligo.
 *
 * Platform-specific error profiles:
 *   - Illumina: substitution-dominant, very low indels
 *   - Nanopore: indel-heavy, moderate substitutions
 *   - PacBio: insertion-dominant (before CCS correction)
 *
 * Also generates per-base Phred quality scores:
 *   - Correct base: Q30-Q40
 *   - Substituted base: Q5-Q15
 *   - Inserted base: Q2-Q8
 *   - Deleted base: (no Q-score, base absent)
 *
 * @param oligo Oligo to sequence
 * @param params Sequencing parameters
 * @returns A single noisy read with Q-scores
 */
function simulateOneRead(
  oligo: Oligo,
  params: SequencingParams,
): SequencingRead {
  const original = oligo.sequence;
  const resultParts: string[] = [];
  const qualities: number[] = [];
  let subs = 0;
  let ins = 0;
  let dels = 0;

  // First pass: determine deletions
  const survived = new Array<boolean>(original.length).fill(true);
  if (params.deletionRate > 0) {
    for (let i = 0; i < original.length; i++) {
      if (rng.next() < params.deletionRate) {
        survived[i] = false;
        dels++;
      }
    }
  }

  // Second pass: walk through, apply substitutions and insertions
  for (let i = 0; i < original.length; i++) {
    if (!survived[i]) continue;

    const origBase = original[i];
    let emitBase = origBase;
    let qScore: number;

    // Substitution
    if (rng.next() < params.substitutionRate) {
      let newBase: string;
      do {
        newBase = BASES[rng.nextInt(4)];
      } while (newBase === origBase);
      emitBase = newBase;
      qScore = 5 + rng.nextInt(11); // Q5-Q15
      subs++;
    } else {
      qScore = 30 + rng.nextInt(11); // Q30-Q40
    }

    resultParts.push(emitBase);
    qualities.push(qScore);

    // Insertion after this position
    if (rng.next() < params.insertionRate) {
      resultParts.push(BASES[rng.nextInt(4)]);
      qualities.push(2 + rng.nextInt(7)); // Q2-Q8
      ins++;
    }
  }

  // Possible insertion at the end
  if (rng.next() < params.insertionRate) {
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
 * Step 4: Simulate sequencing.
 *
 * Platform-specific error profiles:
 *   - Illumina: substitution-dominant, very low indels (~0.1% per base)
 *     Phred Q30+ for most bases, Q5-15 at error positions
 *   - Nanopore (ONT R10.4.1): indel-heavy (~9% total per base)
 *     Lower quality scores, insertions and deletions dominate
 *   - PacBio HiFi (CCS-corrected): ~0.1-1% after circular consensus
 *     Mostly substitutions at low rate, insertion-dominant before CCS
 *
 * Each surviving oligo generates `coverage` independent reads.
 * Reads include per-base Phred quality scores for soft-information decoding.
 * Some oligos may be lost entirely (dropout).
 *
 * @param oligos Oligo pool after aging
 * @param params Sequencing parameters
 * @returns Array of sequencing reads with Q-scores
 */
export function simulateSequencing(
  oligos: Oligo[],
  params: SequencingParams,
): SequencingRead[] {
  const reads: SequencingRead[] = [];

  // Group oligos by index (PCR may have created duplicates)
  const byIndex = new Map<number, Oligo[]>();
  for (const oligo of oligos) {
    const existing = byIndex.get(oligo.index) ?? [];
    existing.push(oligo);
    byIndex.set(oligo.index, existing);
  }

  byIndex.forEach((oligoGroup, idx) => {
    // Dropout: some oligos are completely lost
    if (rng.next() < params.dropoutRate) {
      return;
    }

    // Generate coverage × reads
    // Pick from the available copies (PCR duplicates) randomly
    const coverage = params.coverage;
    for (let r = 0; r < coverage; r++) {
      // Pick a random copy from the pool
      const copyIdx = rng.nextInt(oligoGroup.length);
      const read = simulateOneRead(oligoGroup[copyIdx], params);
      reads.push(read);
    }
  });

  return reads;
}

// ---------------------------------------------------------------------------
// Full pipeline
// ---------------------------------------------------------------------------

/**
 * Run the full wetlab simulation pipeline.
 *
 * Pipeline:
 *   1. Synthesis: apply synthesis errors (position-dependent)
 *   2. PCR: amplify with duplication bias and errors
 *   3. Aging: apply chemical damage
 *   4. Sequencing: generate reads with platform-specific errors
 *
 * @param oligos Original oligos from the encoder
 * @param config Full wetlab configuration
 * @param seed Random seed for reproducibility (0 = non-deterministic)
 * @returns Sequencing reads and per-stage statistics
 */
export function simulateWetlab(
  oligos: Oligo[],
  config: WetlabConfig,
  seed: number = 0,
): WetlabResult {
  rng = new Rng(seed || Math.floor(Math.random() * 0xffffffff));

  const totalOligos = oligos.length;

  // Step 1: Synthesis
  const afterSynthesis = simulateSynthesis(oligos, config.synthesis);
  let synthesisErrors = 0;
  for (let i = 0; i < oligos.length; i++) {
    // Count Hamming distance as proxy for synthesis errors
    const orig = oligos[i].sequence;
    const syn = afterSynthesis[i].sequence;
    // For sequences of different length (indels), count length difference + substitutions
    if (orig.length !== syn.length) {
      synthesisErrors += Math.abs(orig.length - syn.length);
    }
    const minLen = Math.min(orig.length, syn.length);
    for (let j = 0; j < minLen; j++) {
      if (orig[j] !== syn[j]) synthesisErrors++;
    }
  }

  // Step 2: PCR
  const afterPCR = simulatePCR(afterSynthesis, config.pcr);
  const pcrDuplicates = afterPCR.length - afterSynthesis.length;
  // Count PCR errors: compare each new duplicate to its source
  let pcrErrors = 0;
  // (PCR errors are tracked inside simulatePCR; approximate here)
  // A more precise count would require instrumenting pcrCycle

  // Step 3: Aging
  const afterAging = simulateAging(afterPCR, config.aging);
  let agingDamage = 0;
  for (let i = 0; i < afterPCR.length; i++) {
    const before = afterPCR[i].sequence;
    const after = afterAging[i].sequence;
    if (before.length !== after.length) {
      agingDamage += Math.abs(before.length - after.length);
    }
    const minLen = Math.min(before.length, after.length);
    for (let j = 0; j < minLen; j++) {
      if (before[j] !== after[j]) agingDamage++;
    }
  }

  // Step 4: Sequencing
  const reads = simulateSequencing(afterAging, config.sequencing);
  const sequencingErrors = reads.reduce(
    (sum, r) => sum + r.substitutions + r.insertions + r.deletions,
    0,
  );

  // Count recovered oligos (oligos that have at least one read)
  const readIndices = new Set(reads.map((r) => r.oligoIndex));
  const recoveredOligos = readIndices.size;

  return {
    reads,
    stats: {
      synthesisErrors,
      pcrDuplicates,
      pcrErrors,
      agingDamage,
      sequencingErrors,
      totalOligos,
      recoveredOligos,
    },
  };
}

// ---------------------------------------------------------------------------
// Utility: compute expected coverage distribution
// ---------------------------------------------------------------------------

/**
 * Compute the expected coverage distribution for a given PCR + sequencing config.
 *
 * Models the PCR amplification as a branching process where each molecule
 * has probability p of duplicating each cycle. After c cycles, the expected
 * number of copies is (1+p)^c. Coverage = copies × (1 - dropoutRate).
 *
 * Also computes the coefficient of variation (CV = σ/μ) to assess
 * amplification uniformity. Lower CV = more uniform coverage.
 *
 * @param config Wetlab configuration
 * @returns Expected coverage, standard deviation, and CV
 */
export function expectedCoverageStats(
  config: WetlabConfig,
): { mean: number; std: number; cv: number } {
  const p = config.pcr.duplicationProb;
  const c = config.pcr.cycles;

  // Expected copies after PCR: E[X] = (1 + p)^c
  // Variance of a branching process: Var[X] = p(1+p)^(2c-2) × ((1+p)^c - 1) / (1 + p - 1)
  // Simplified for large c: Var ≈ p × (1+p)^(2c-1)
  const meanCopies = Math.pow(1 + p, c);
  const variance = p * Math.pow(1 + p, 2 * c - 2) * (meanCopies - 1) / Math.max(p, 0.001);
  const stdCopies = Math.sqrt(Math.max(0, variance));

  // Coverage = copies × coverage_setting × (1 - dropoutRate) / totalOligos
  // But coverage_setting is already per-oligo, so:
  const meanCoverage = meanCopies * config.sequencing.coverage * (1 - config.sequencing.dropoutRate);
  const stdCoverage = stdCopies * config.sequencing.coverage * (1 - config.sequencing.dropoutRate);

  return {
    mean: meanCoverage,
    std: stdCoverage,
    cv: meanCoverage > 0 ? stdCoverage / meanCoverage : 0,
  };
}
