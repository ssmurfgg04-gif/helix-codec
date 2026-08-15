/**
 * In-Silico Wetlab Simulation — Full DNA Storage Pipeline
 *
 * Simulates the complete physical wetlab pipeline for DNA data storage:
 *   1. Synthesis: digital data → DNA oligos (with synthesis errors)
 *   2. Storage: DNA oligos degrade over time (aging, decay, fragmentation)
 *   3. Sequencing: DNA → digital reads (with platform-specific errors)
 *   4. Decoding: reads → recovered digital data
 *
 * This module provides a higher-level API than dt4dds-simulate.ts:
 *   - Operates on raw digital data (Uint8Array) rather than Oligo objects
 *   - Computes end-to-end metrics: BER, recovery rate, oligo survival
 *   - Provides preset configurations for common experimental conditions
 *   - Generates FASTQ-like output for compatibility with bioinformatics tools
 *
 * The simulation is genuinely stochastic — each run produces different errors.
 * Error rates are based on published experimental measurements:
 *   - Illumina array synthesis: sub ~0.1%, ins ~0.05%, del ~0.1% (Chandak 2018)
 *   - Nanopore R10.4.1: sub ~2%, ins ~3%, del ~4% (Preuss 2024)
 *   - PacBio HiFi (CCS): sub ~0.5%, ins ~0.1%, del ~0.1% (Wenger 2019)
 *
 * References:
 *   - Chandak et al. (2018), ISMB — Illumina DNA-storage error rates
 *   - Organick et al. (2018), Nature Biotechnology 36:242-248
 *   - Lee et al. (2022), Nat. Commun. 13:3231 — dt4dds pipeline model
 *   - Banal et al. (2026), arXiv:2604.20810 — soft-information decoding
 *   - Preuss et al. (2026), Nature Sci Rep — measured real-world rates
 */

import { Oligo, CodecMetadata, EncodedFile } from "./types";
import { SequencingRead, SimulationResult } from "./simulate";
import {
  WetlabConfig,
  WetlabResult,
  SynthesisParams,
  PCRParams,
  AgingParams,
  SequencingParams,
  ILLUMINA_WETLAB,
  NANOPORE_WETLAB,
  PACBIO_WETLAB,
  simulateWetlab as dt4ddsSimulateWetlab,
  simulateSynthesis as dt4ddsSimulateSynthesis,
  simulatePCR,
  simulateAging,
  simulateSequencing as dt4ddsSimulateSequencing,
  expectedCoverageStats,
} from "./dt4dds-simulate";

// ---------------------------------------------------------------------------
// PRNG (Xorshift32 — deterministic, same as dt4dds-simulate.ts)
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

const BASES = "ACGT";

// ---------------------------------------------------------------------------
// Synthesis simulation configuration
// ---------------------------------------------------------------------------

/**
 * Synthesis error model configuration.
 *
 * Models errors introduced during DNA array synthesis:
 *   - Substitutions: wrong base is incorporated (e.g., A→C)
 *   - Insertions: extra base is added
 *   - Deletions: base is skipped
 *   - Position dependence: error rates are higher at 5'/3' ends
 *   - Homopolymer breaking: runs > 4 same bases are broken (synthesis limitation)
 *   - GC constraints: oligos outside 40-60% GC may fail synthesis
 */
export interface SynthesisConfig {
  /** Per-base substitution rate. Default: 0.001 */
  substitutionRate: number;
  /** Per-base insertion rate. Default: 0.0005 */
  insertionRate: number;
  /** Per-base deletion rate. Default: 0.001 */
  deletionRate: number;
  /** Enable position-dependent error scaling (higher at ends). Default: true */
  positionDependent: boolean;
  /** 5' end error multiplier. Default: 1.5 */
  fivePrimeErrorMultiplier: number;
  /** 3' end error multiplier. Default: 2.0 */
  threePrimeErrorMultiplier: number;
  /** Maximum homopolymer run before forced breaking. Default: 4 */
  maxHomopolymerRun: number;
  /** Minimum GC content for synthesis success. Default: 0.40 */
  gcMin: number;
  /** Maximum GC content for synthesis success. Default: 0.60 */
  gcMax: number;
  /** Probability an oligo fails synthesis entirely. Default: 0.01 */
  failureRate: number;
}

/** Default synthesis configuration (Illumina-grade array synthesis) */
export const DEFAULT_SYNTHESIS: SynthesisConfig = {
  substitutionRate: 0.001,
  insertionRate: 0.0005,
  deletionRate: 0.001,
  positionDependent: true,
  fivePrimeErrorMultiplier: 1.5,
  threePrimeErrorMultiplier: 2.0,
  maxHomopolymerRun: 4,
  gcMin: 0.40,
  gcMax: 0.60,
  failureRate: 0.01,
};

// ---------------------------------------------------------------------------
// Storage simulation configuration
// ---------------------------------------------------------------------------

/**
 * Storage degradation model configuration.
 *
 * Models the physical degradation of DNA molecules during storage:
 *   - Depurination: A/G bases lose purine rings → deletion
 *   - Oxidation: G → T or G → C substitution
 *   - Deamination: C → T substitution (cytosine deamination)
 *   - Fragmentation: random strand breaks → shorter fragments
 *   - Coverage reduction: some molecules are lost (Poisson sampling)
 */
export interface StorageConfig {
  /** Number of days of storage. Default: 0 (no aging) */
  days: number;
  /** Depurination rate per base per day. Default: 0.0001 */
  depurinationRate: number;
  /** Oxidation rate per base per day. Default: 0.00005 */
  oxidationRate: number;
  /** Deamination rate (C→T) per base per day. Default: 0.0002 */
  deaminationRate: number;
  /** Fragmentation rate per base per day (strand breaks). Default: 0.00001 */
  fragmentationRate: number;
  /** Fraction of oligos lost per year of storage. Default: 0.02 */
  lossRatePerYear: number;
  /** Number of PCR cycles for amplification before storage. Default: 15 */
  pcrCycles: number;
  /** PCR duplication probability per cycle. Default: 0.85 */
  pcrDuplicationProb: number;
  /** PCR per-cycle substitution rate (Taq fidelity). Default: 0.0001 */
  pcrSubstitutionRate: number;
  /** PCR GC amplification bias. Default: 0.1 */
  pcrGcBias: number;
}

/** Default storage configuration (no aging, standard PCR) */
export const DEFAULT_STORAGE: StorageConfig = {
  days: 0,
  depurinationRate: 0.0001,
  oxidationRate: 0.00005,
  deaminationRate: 0.0002,
  fragmentationRate: 0.00001,
  lossRatePerYear: 0.02,
  pcrCycles: 15,
  pcrDuplicationProb: 0.85,
  pcrSubstitutionRate: 0.0001,
  pcrGcBias: 0.1,
};

// ---------------------------------------------------------------------------
// Sequencing simulation configuration
// ---------------------------------------------------------------------------

/**
 * Sequencing platform configuration.
 *
 * Models platform-specific error profiles:
 *   - Illumina: substitution-dominant, very low indels (~0.1% total)
 *   - Nanopore: indel-heavy, moderate substitutions (~5-15% total)
 *   - PacBio HiFi: after CCS correction, ~0.1-1% Q30+ reads
 */
export interface SequencingConfig {
  /** Sequencing platform. Default: "illumina" */
  platform: "illumina" | "nanopore" | "pacbio";
  /** Per-base substitution rate. Default varies by platform */
  substitutionRate: number;
  /** Per-base insertion rate. Default varies by platform */
  insertionRate: number;
  /** Per-base deletion rate. Default varies by platform */
  deletionRate: number;
  /** Average coverage depth (reads per oligo). Default: 20 */
  coverage: number;
  /** Dropout rate (oligos completely lost). Default: 0 */
  dropoutRate: number;
  /** Mean read length for long-read platforms. Default: 0 (use oligo length) */
  meanReadLength: number;
}

/** Illumina sequencing preset — substitution-dominant, low indels */
export const ILLUMINA_SEQUENCING: SequencingConfig = {
  platform: "illumina",
  substitutionRate: 0.001,
  insertionRate: 0.0005,
  deletionRate: 0.001,
  coverage: 20,
  dropoutRate: 0,
  meanReadLength: 0,
};

/** Nanopore R10.4.1 sequencing preset — indel-heavy, ~9% total error */
export const NANOPORE_SEQUENCING: SequencingConfig = {
  platform: "nanopore",
  substitutionRate: 0.02,
  insertionRate: 0.03,
  deletionRate: 0.04,
  coverage: 15,
  dropoutRate: 0.05,
  meanReadLength: 0,
};

/** PacBio HiFi (CCS-corrected) sequencing preset — ~0.1-1% error */
export const PACBIO_SEQUENCING: SequencingConfig = {
  platform: "pacbio",
  substitutionRate: 0.005,
  insertionRate: 0.001,
  deletionRate: 0.001,
  coverage: 10,
  dropoutRate: 0.02,
  meanReadLength: 0,
};

// ---------------------------------------------------------------------------
// Full pipeline configuration
// ---------------------------------------------------------------------------

/** Complete wetlab simulation configuration */
export interface WetlabSimConfig {
  synthesis: SynthesisConfig;
  storage: StorageConfig;
  sequencing: SequencingConfig;
  /** Random seed for reproducibility. 0 = non-deterministic */
  seed: number;
}

/** Default full pipeline configuration (Illumina-grade, no aging) */
export const DEFAULT_WETLAB_SIM: WetlabSimConfig = {
  synthesis: DEFAULT_SYNTHESIS,
  storage: DEFAULT_STORAGE,
  sequencing: ILLUMINA_SEQUENCING,
  seed: 0,
};

/** Nanopore pipeline configuration with aging */
export const NANOPORE_WETLAB_SIM: WetlabSimConfig = {
  synthesis: {
    ...DEFAULT_SYNTHESIS,
    substitutionRate: 0.002,
    insertionRate: 0.001,
    deletionRate: 0.002,
    fivePrimeErrorMultiplier: 1.5,
    threePrimeErrorMultiplier: 2.5,
  },
  storage: { ...DEFAULT_STORAGE, days: 365 }, // 1 year of aging
  sequencing: NANOPORE_SEQUENCING,
  seed: 0,
};

/** Low-error pipeline for testing (minimal errors) */
export const LOW_ERROR_WETLAB_SIM: WetlabSimConfig = {
  synthesis: {
    substitutionRate: 0.0001,
    insertionRate: 0.0001,
    deletionRate: 0.0001,
    positionDependent: false,
    fivePrimeErrorMultiplier: 1.0,
    threePrimeErrorMultiplier: 1.0,
    maxHomopolymerRun: 4,
    gcMin: 0.30,
    gcMax: 0.70,
    failureRate: 0.005,
  },
  storage: { ...DEFAULT_STORAGE, days: 0 },
  sequencing: {
    platform: "illumina",
    substitutionRate: 0.0005,
    insertionRate: 0.0002,
    deletionRate: 0.0005,
    coverage: 30,
    dropoutRate: 0.01,
    meanReadLength: 0,
  },
  seed: 42,
};

// ---------------------------------------------------------------------------
// Simulation result types
// ---------------------------------------------------------------------------

/** Result of synthesis simulation */
export interface SynthesisResult {
  /** Synthesized DNA oligos (with errors applied) */
  oligos: Oligo[];
  /** Indices of oligos that failed synthesis */
  failedOligos: number[];
  /** Number of substitution errors introduced */
  substitutions: number;
  /** Number of insertion errors introduced */
  insertions: number;
  /** Number of deletion errors introduced */
  deletions: number;
  /** Number of homopolymer runs broken */
  homopolymerBreaks: number;
  /** Oligos flagged for GC violation */
  gcViolations: number;
}

/** Result of storage simulation */
export interface StorageResult {
  /** Surviving oligos after storage degradation */
  oligos: Oligo[];
  /** Number of oligos lost during storage */
  lostOligos: number;
  /** Number of bases damaged (depurination + oxidation + deamination) */
  damagedBases: number;
  /** Number of fragments created by strand breaks */
  fragments: number;
  /** Fraction of original oligos that survived */
  survivalRate: number;
}

/** Result of sequencing simulation */
export interface SequencingSimResult {
  /** Generated sequencing reads */
  reads: SequencingRead[];
  /** Total number of reads generated */
  totalReads: number;
  /** Average read length */
  avgReadLength: number;
  /** Average coverage depth (reads per oligo) */
  avgCoverageDepth: number;
  /** Total sequencing errors (subs + ins + dels) */
  totalErrors: number;
  /** FASTQ-formatted string of all reads */
  fastq: string;
}

/** Full pipeline result with end-to-end metrics */
export interface PipelineResult {
  /** Synthesis stage result */
  synthesis: SynthesisResult;
  /** Storage stage result */
  storage: StorageResult;
  /** Sequencing stage result */
  sequencing: SequencingSimResult;
  /** Decoded data (if decode was performed) */
  decodedData: Uint8Array | null;
  /** Bit Error Rate: fraction of incorrect bits */
  ber: number;
  /** Recovery rate: fraction of original bytes correctly recovered */
  recoveryRate: number;
  /** Oligo survival rate through the entire pipeline */
  oligoSurvivalRate: number;
  /** Average coverage depth of surviving oligos */
  coverageDepth: number;
  /** Total simulation wall time in ms */
  simulationTimeMs: number;
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/** Compute GC content of a DNA string */
function gcContent(dna: string): number {
  if (dna.length === 0) return 0;
  let gc = 0;
  for (let i = 0; i < dna.length; i++) {
    const c = dna.charCodeAt(i);
    if (c === 71 || c === 67) gc++; // G=71, C=67
  }
  return gc / dna.length;
}

/** Compute maximum homopolymer run in a DNA string */
function maxHomopolymer(dna: string): number {
  if (dna.length === 0) return 0;
  let maxRun = 1;
  let run = 1;
  for (let i = 1; i < dna.length; i++) {
    if (dna[i] === dna[i - 1]) {
      run++;
      if (run > maxRun) maxRun = run;
    } else {
      run = 1;
    }
  }
  return maxRun;
}

/**
 * Position-dependent error multiplier.
 * Higher error rates at 5' and 3' ends of oligo (array synthesis artifact).
 */
function positionErrorMultiplier(
  pos: number,
  len: number,
  fivePrimeMul: number,
  threePrimeMul: number,
): number {
  if (len <= 0) return 1.0;
  const edgeFraction = 0.15;
  const edgeLen = Math.max(1, Math.floor(len * edgeFraction));
  if (pos < edgeLen) {
    const t = pos / edgeLen;
    return fivePrimeMul + t * (1.0 - fivePrimeMul);
  } else if (pos >= len - edgeLen) {
    const t = (pos - (len - edgeLen)) / edgeLen;
    return 1.0 + t * (threePrimeMul - 1.0);
  }
  return 1.0;
}

// ---------------------------------------------------------------------------
// Stage 1: Synthesis simulation
// ---------------------------------------------------------------------------

/**
 * Simulate DNA array synthesis with realistic error models.
 *
 * Applies three types of errors at per-base rates:
 *   - Substitutions: wrong base incorporated (e.g., phosphoramidite contamination)
 *   - Insertions: extra base added (coupling failure → detritylation)
 *   - Deletions: base skipped (incomplete coupling → capping failure)
 *
 * Additionally models:
 *   - Position-dependent errors (5'/3' end effects)
 *   - Homopolymer run breaking (synthesis fails on long runs of same base)
 *   - GC content violations (sequences outside 40-60% GC may fail)
 *   - Complete synthesis failure (some oligos fail to synthesize at all)
 *
 * @param oligos DNA oligos from the encoder
 * @param config Synthesis configuration
 * @param seed Random seed (0 = non-deterministic)
 * @returns Synthesized oligos with errors and statistics
 */
export function simulateSynthesis(
  oligos: Oligo[],
  config: SynthesisConfig = DEFAULT_SYNTHESIS,
  seed: number = 0,
): SynthesisResult {
  const rng = new Rng(seed || Math.floor(Math.random() * 0xffffffff));
  const result: Oligo[] = [];
  const failedOligos: number[] = [];
  let totalSubs = 0;
  let totalIns = 0;
  let totalDels = 0;
  let homopolymerBreaks = 0;
  let gcViolations = 0;

  for (const oligo of oligos) {
    // Check for complete synthesis failure
    if (rng.next() < config.failureRate) {
      failedOligos.push(oligo.index);
      continue;
    }

    let dna = oligo.sequence;
    let subs = 0;
    let ins = 0;
    let dels = 0;

    // --- Apply synthesis errors (position-dependent) ---

    // First pass: determine deletions
    const len = dna.length;
    const survived = new Array<boolean>(len).fill(true);
    for (let i = 0; i < len; i++) {
      const mul = config.positionDependent
        ? positionErrorMultiplier(i, len, config.fivePrimeErrorMultiplier, config.threePrimeErrorMultiplier)
        : 1.0;
      if (rng.next() < config.deletionRate * mul) {
        survived[i] = false;
        dels++;
      }
    }

    // Second pass: apply substitutions and insertions
    const parts: string[] = [];
    for (let i = 0; i < len; i++) {
      if (!survived[i]) continue;

      const mul = config.positionDependent
        ? positionErrorMultiplier(i, len, config.fivePrimeErrorMultiplier, config.threePrimeErrorMultiplier)
        : 1.0;

      // Substitution
      if (rng.next() < config.substitutionRate * mul) {
        let newBase: string;
        do {
          newBase = BASES[rng.nextInt(4)];
        } while (newBase === dna[i]);
        parts.push(newBase);
        subs++;
      } else {
        parts.push(dna[i]);
      }

      // Insertion after this position
      if (rng.next() < config.insertionRate * mul) {
        parts.push(BASES[rng.nextInt(4)]);
        ins++;
      }
    }

    dna = parts.join("");

    // --- Homopolymer run breaking ---
    // If any run exceeds maxHomopolymerRun, insert a different base to break it
    const broken: string[] = [];
    let runLen = 1;
    for (let i = 0; i < dna.length; i++) {
      if (i > 0 && dna[i] === dna[i - 1]) {
        runLen++;
        if (runLen > config.maxHomopolymerRun) {
          // Insert a different base to break the run
          let breaker: string;
          do {
            breaker = BASES[rng.nextInt(4)];
          } while (breaker === dna[i]);
          broken.push(breaker);
          broken.push(dna[i]);
          runLen = 1;
          homopolymerBreaks++;
        } else {
          broken.push(dna[i]);
        }
      } else {
        runLen = 1;
        broken.push(dna[i]);
      }
    }
    dna = broken.join("");

    // --- GC content check ---
    const gc = gcContent(dna);
    if (gc < config.gcMin || gc > config.gcMax) {
      gcViolations++;
    }

    totalSubs += subs;
    totalIns += ins;
    totalDels += dels;

    result.push({
      ...oligo,
      sequence: dna,
      gc,
      maxHomopolymer: maxHomopolymer(dna),
      length: dna.length,
    });
  }

  return {
    oligos: result,
    failedOligos,
    substitutions: totalSubs,
    insertions: totalIns,
    deletions: totalDels,
    homopolymerBreaks,
    gcViolations,
  };
}

// ---------------------------------------------------------------------------
// Stage 2: Storage simulation
// ---------------------------------------------------------------------------

/**
 * Simulate storage degradation of DNA molecules.
 *
 * Models three chemical damage mechanisms:
 *   - Depurination: A/G bases lose purine rings → base deletion
 *     (dominant damage mechanism in aqueous storage)
 *   - Oxidation: G → T or G → C substitution
 *     (caused by reactive oxygen species)
 *   - Deamination: C → T substitution
 *     (spontaneous hydrolysis of cytosine's amine group)
 *
 * Also models:
 *   - Fragmentation: random strand breaks create shorter fragments
 *   - Coverage reduction: some molecules are lost entirely (Poisson sampling)
 *
 * @param oligos DNA oligos after synthesis
 * @param config Storage configuration
 * @param seed Random seed (0 = non-deterministic)
 * @returns Degraded oligos and survival statistics
 */
export function simulateStorage(
  oligos: Oligo[],
  config: StorageConfig = DEFAULT_STORAGE,
  seed: number = 0,
): StorageResult {
  const rng = new Rng(seed || Math.floor(Math.random() * 0xffffffff));
  const result: Oligo[] = [];
  let lostOligos = 0;
  let damagedBases = 0;
  let fragments = 0;

  // Loss rate scales with storage duration
  const lossRate = config.lossRatePerYear * (config.days / 365.25);

  for (const oligo of oligos) {
    // Coverage reduction: some oligos are lost entirely
    if (rng.next() < lossRate) {
      lostOligos++;
      continue;
    }

    let dna = oligo.sequence;
    let damage = 0;

    // --- Chemical aging ---
    if (config.days > 0) {
      const depurRate = config.depurinationRate * config.days;
      const oxidRate = config.oxidationRate * config.days;
      const deamRate = config.deaminationRate * config.days;

      const aged: string[] = [];
      for (let i = 0; i < dna.length; i++) {
        const base = dna[i];

        // Depurination: A or G → deletion
        if ((base === "A" || base === "G") && rng.next() < depurRate) {
          damage++;
          continue;
        }

        // Oxidation: G → T or G → C
        if (base === "G" && rng.next() < oxidRate) {
          aged.push(rng.next() < 0.5 ? "T" : "C");
          damage++;
          continue;
        }

        // Deamination: C → T
        if (base === "C" && rng.next() < deamRate) {
          aged.push("T");
          damage++;
          continue;
        }

        aged.push(base);
      }
      dna = aged.join("");
    }

    // --- Fragmentation: random strand breaks ---
    if (config.fragmentationRate > 0 && config.days > 0) {
      const fragRate = config.fragmentationRate * config.days;
      const breakPoints: number[] = [];

      // Determine break positions
      for (let i = 1; i < dna.length; i++) {
        if (rng.next() < fragRate) {
          breakPoints.push(i);
        }
      }

      if (breakPoints.length > 0) {
        fragments += breakPoints.length;
        // Use the longest fragment (most likely to be sequenced)
        let longestStart = 0;
        let longestLen = breakPoints[0];
        for (let i = 0; i < breakPoints.length; i++) {
          const start = breakPoints[i];
          const end = i + 1 < breakPoints.length ? breakPoints[i + 1] : dna.length;
          const len = end - start;
          if (len > longestLen) {
            longestLen = len;
            longestStart = start;
          }
        }
        // Also check the last fragment
        const lastLen = dna.length - breakPoints[breakPoints.length - 1];
        if (lastLen > longestLen) {
          longestLen = lastLen;
          longestStart = breakPoints[breakPoints.length - 1];
        }
        dna = dna.slice(longestStart, longestStart + longestLen);
      }
    }

    damagedBases += damage;

    result.push({
      ...oligo,
      sequence: dna,
      gc: gcContent(dna),
      maxHomopolymer: maxHomopolymer(dna),
      length: dna.length,
    });
  }

  return {
    oligos: result,
    lostOligos,
    damagedBases,
    fragments,
    survivalRate: oligos.length > 0 ? result.length / oligos.length : 0,
  };
}

// ---------------------------------------------------------------------------
// Stage 3: Sequencing simulation
// ---------------------------------------------------------------------------

/**
 * Generate a Phred quality score string for FASTQ output.
 *
 * Q = -10 × log10(P_error), encoded as ASCII(Q + 33).
 * Q2='%', Q10='+', Q20='5', Q30='?', Q40='I'
 */
function qualityToFastq(quality: Uint8Array): string {
  const chars: string[] = new Array(quality.length);
  for (let i = 0; i < quality.length; i++) {
    chars[i] = String.fromCharCode(Math.min(quality[i], 93) + 33); // max Q93 = '~'
  }
  return chars.join("");
}

/**
 * Simulate DNA sequencing with platform-specific error models.
 *
 * Platform error profiles:
 *   - Illumina: substitution-dominant (~0.1% total error per base)
 *     Phred Q30+ for most bases, with occasional Q5-15 errors
 *   - Nanopore: indel-heavy (~5-15% total error per base)
 *     Lower quality scores, more insertions and deletions
 *   - PacBio HiFi: after CCS correction, ~0.1-1% error
 *     Mostly substitutions at low rate
 *
 * Each oligo is read multiple times (coverage depth), with independent
 * errors per read. Some oligos may be completely lost (dropout).
 *
 * Output includes FASTQ-formatted reads for compatibility with
 * standard bioinformatics tools (bwa, minimap2, etc.)
 *
 * @param oligos DNA oligos after storage
 * @param config Sequencing configuration
 * @param seed Random seed (0 = non-deterministic)
 * @returns Sequencing reads with errors and FASTQ output
 */
export function simulateSequencing(
  oligos: Oligo[],
  config: SequencingConfig = ILLUMINA_SEQUENCING,
  seed: number = 0,
): SequencingSimResult {
  const rng = new Rng(seed || Math.floor(Math.random() * 0xffffffff));
  const reads: SequencingRead[] = [];
  const fastqLines: string[] = [];
  let readCounter = 0;

  // Group oligos by index (multiple copies from PCR)
  const byIndex = new Map<number, Oligo[]>();
  for (const oligo of oligos) {
    const existing = byIndex.get(oligo.index) ?? [];
    existing.push(oligo);
    byIndex.set(oligo.index, existing);
  }

  const uniqueIndices = Array.from(byIndex.keys());

  for (const idx of uniqueIndices) {
    const oligoGroup = byIndex.get(idx)!;

    // Dropout: some oligos completely lost
    if (rng.next() < config.dropoutRate) {
      continue;
    }

    // Generate coverage × reads for this oligo
    const coverage = config.coverage;
    for (let r = 0; r < coverage; r++) {
      // Pick a random copy from the PCR pool
      const copyIdx = rng.nextInt(oligoGroup.length);
      const oligo = oligoGroup[copyIdx];
      const original = oligo.sequence;

      const resultParts: string[] = [];
      const qualities: number[] = [];
      let subs = 0;
      let ins = 0;
      let dels = 0;

      // First pass: determine deletions
      const survived = new Array<boolean>(original.length).fill(true);
      if (config.deletionRate > 0) {
        for (let i = 0; i < original.length; i++) {
          if (rng.next() < config.deletionRate) {
            survived[i] = false;
            dels++;
          }
        }
      }

      // Second pass: apply substitutions and insertions with Q-scores
      for (let i = 0; i < original.length; i++) {
        if (!survived[i]) continue;

        const origBase = original[i];
        let emitBase = origBase;
        let qScore: number;

        // Substitution
        if (rng.next() < config.substitutionRate) {
          let newBase: string;
          do {
            newBase = BASES[rng.nextInt(4)];
          } while (newBase === origBase);
          emitBase = newBase;
          qScore = 5 + rng.nextInt(11); // Q5-Q15 (low confidence)
          subs++;
        } else {
          qScore = 30 + rng.nextInt(11); // Q30-Q40 (high confidence)
        }

        resultParts.push(emitBase);
        qualities.push(qScore);

        // Insertion after this position
        if (rng.next() < config.insertionRate) {
          resultParts.push(BASES[rng.nextInt(4)]);
          qualities.push(2 + rng.nextInt(7)); // Q2-Q8 (very low confidence)
          ins++;
        }
      }

      // Possible insertion at the end
      if (rng.next() < config.insertionRate) {
        resultParts.push(BASES[rng.nextInt(4)]);
        qualities.push(2 + rng.nextInt(7));
        ins++;
      }

      const sequence = resultParts.join("");
      const quality = new Uint8Array(qualities);

      const read: SequencingRead = {
        oligoIndex: idx,
        sequence,
        quality,
        substitutions: subs,
        insertions: ins,
        deletions: dels,
      };

      reads.push(read);

      // Generate FASTQ entry
      const readId = `read_${readCounter++}`;
      fastqLines.push(`@${readId} oligo_idx=${idx}`);
      fastqLines.push(sequence);
      fastqLines.push("+");
      fastqLines.push(qualityToFastq(quality));
    }
  }

  const totalReads = reads.length;
  const avgReadLength =
    totalReads > 0 ? reads.reduce((s, r) => s + r.sequence.length, 0) / totalReads : 0;
  const avgCoverageDepth =
    uniqueIndices.length > 0 ? totalReads / uniqueIndices.length : 0;
  const totalErrors = reads.reduce(
    (s, r) => s + r.substitutions + r.insertions + r.deletions,
    0,
  );

  return {
    reads,
    totalReads,
    avgReadLength,
    avgCoverageDepth,
    totalErrors,
    fastq: fastqLines.join("\n"),
  };
}

// ---------------------------------------------------------------------------
// Full pipeline
// ---------------------------------------------------------------------------

/**
 * Compute Bit Error Rate (BER) between original and recovered data.
 *
 * BER = (number of differing bits) / (total bits)
 * Compares byte-by-byte up to the length of the shorter array,
 * then counts remaining bytes in the longer array as all-errors.
 *
 * @param original Original data
 * @param recovered Recovered data
 * @returns BER in [0, 1]
 */
export function computeBER(original: Uint8Array, recovered: Uint8Array): number {
  const minLen = Math.min(original.length, recovered.length);
  let bitErrors = 0;

  // Count bit errors in the overlapping region
  for (let i = 0; i < minLen; i++) {
    const diff = original[i] ^ recovered[i];
    // Count set bits (population count)
    bitErrors += popcount8(diff);
  }

  // Count all bits as errors in the non-overlapping region
  const maxLen = Math.max(original.length, recovered.length);
  bitErrors += (maxLen - minLen) * 8;

  const totalBits = maxLen * 8;
  return totalBits > 0 ? bitErrors / totalBits : 0;
}

/** Population count for a single byte */
function popcount8(n: number): number {
  n = n - ((n >> 1) & 0x55);
  n = (n & 0x33) + ((n >> 2) & 0x33);
  return ((n + (n >> 4)) & 0x0f);
}

/**
 * Compute recovery rate: fraction of original bytes correctly recovered.
 *
 * @param original Original data
 * @param recovered Recovered data
 * @returns Recovery rate in [0, 1]
 */
export function computeRecoveryRate(original: Uint8Array, recovered: Uint8Array): number {
  if (original.length === 0) return recovered.length === 0 ? 1 : 0;

  const minLen = Math.min(original.length, recovered.length);
  let correct = 0;

  for (let i = 0; i < minLen; i++) {
    if (original[i] === recovered[i]) correct++;
  }

  return correct / original.length;
}

/**
 * Run the full in-silico wetlab simulation pipeline.
 *
 * Pipeline stages:
 *   1. Synthesis: encode digital data → DNA, apply synthesis errors
 *   2. Storage: simulate chemical aging, fragmentation, coverage loss
 *   3. Sequencing: generate reads with platform-specific errors
 *   4. (Optional) Decoding: attempt to recover original data
 *
 * If originalData is provided, computes end-to-end BER and recovery rate
 * by comparing decoded data with the original.
 *
 * @param oligos Encoded DNA oligos (from the codec's encodeFile)
 * @param originalData Original digital data (for BER computation)
 * @param decodedData Decoded data (if decode was performed externally)
 * @param config Full simulation configuration
 * @returns Pipeline result with per-stage and end-to-end metrics
 */
export function simulateFullPipeline(
  oligos: Oligo[],
  originalData?: Uint8Array,
  decodedData?: Uint8Array,
  config: WetlabSimConfig = DEFAULT_WETLAB_SIM,
): PipelineResult {
  const t0 = Date.now();
  const masterSeed = config.seed || Math.floor(Math.random() * 0xffffffff);

  // Stage 1: Synthesis
  const synthesisResult = simulateSynthesis(oligos, config.synthesis, masterSeed);

  // Stage 2: Storage (use different seed derived from master)
  const storageSeed = (masterSeed * 1103515245 + 12345) >>> 0;
  const storageResult = simulateStorage(synthesisResult.oligos, config.storage, storageSeed || 1);

  // Stage 3: Sequencing (use different seed derived from master)
  const seqSeed = (storageSeed * 1103515245 + 12345) >>> 0;
  const sequencingResult = simulateSequencing(
    storageResult.oligos,
    config.sequencing,
    seqSeed || 1,
  );

  // Compute end-to-end metrics
  const ber = originalData && decodedData ? computeBER(originalData, decodedData) : -1;
  const recoveryRate =
    originalData && decodedData ? computeRecoveryRate(originalData, decodedData) : -1;

  // Oligo survival through the entire pipeline
  const oligoSurvivalRate =
    oligos.length > 0
      ? (synthesisResult.oligos.length - synthesisResult.failedOligos.length) / oligos.length
      : 0;

  return {
    synthesis: synthesisResult,
    storage: storageResult,
    sequencing: sequencingResult,
    decodedData: decodedData ?? null,
    ber,
    recoveryRate,
    oligoSurvivalRate,
    coverageDepth: sequencingResult.avgCoverageDepth,
    simulationTimeMs: Date.now() - t0,
  };
}

// ---------------------------------------------------------------------------
// Bridge: WetlabSimConfig → WetlabConfig (dt4dds)
// ---------------------------------------------------------------------------

/**
 * Convert a WetlabSimConfig to a WetlabConfig for the dt4dds pipeline.
 *
 * This allows using the higher-level config with the lower-level dt4dds
 * simulation functions.
 */
export function toDt4ddsConfig(config: WetlabSimConfig): WetlabConfig {
  return {
    synthesis: {
      substitutionRate: config.synthesis.substitutionRate,
      insertionRate: config.synthesis.insertionRate,
      deletionRate: config.synthesis.deletionRate,
      positionDependent: config.synthesis.positionDependent,
      fivePrimeErrorMultiplier: config.synthesis.fivePrimeErrorMultiplier,
      threePrimeErrorMultiplier: config.synthesis.threePrimeErrorMultiplier,
    },
    pcr: {
      cycles: config.storage.pcrCycles,
      duplicationProb: config.storage.pcrDuplicationProb,
      substitutionRate: config.storage.pcrSubstitutionRate,
      gcBias: config.storage.pcrGcBias,
    },
    aging: {
      depurinationRate: config.storage.depurinationRate,
      oxidationRate: config.storage.oxidationRate,
      deaminationRate: config.storage.deaminationRate,
      days: config.storage.days,
    },
    sequencing: {
      platform: config.sequencing.platform,
      substitutionRate: config.sequencing.substitutionRate,
      insertionRate: config.sequencing.insertionRate,
      deletionRate: config.sequencing.deletionRate,
      coverage: config.sequencing.coverage,
      dropoutRate: config.sequencing.dropoutRate,
    },
  };
}

/**
 * Run the dt4dds parametric pipeline using a WetlabSimConfig.
 *
 * This is an alternative to simulateFullPipeline that uses the existing
 * dt4dds-simulate.ts implementation (which includes PCR simulation).
 *
 * @param oligos Encoded DNA oligos
 * @param config Simulation configuration
 * @param seed Random seed
 * @returns dt4dds WetlabResult
 */
export function simulateWithDt4dds(
  oligos: Oligo[],
  config: WetlabSimConfig = DEFAULT_WETLAB_SIM,
  seed: number = 0,
): WetlabResult {
  return dt4ddsSimulateWetlab(oligos, toDt4ddsConfig(config), seed);
}
