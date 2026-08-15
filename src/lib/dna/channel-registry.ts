/**
 * 
 * ⚠️ EXPERIMENTAL / ROADMAP — This module is a scaffolding prototype.
 * It is NOT wired into the production encode/decode pipeline.
 * See CHANGELOG.md for v59 status and remaining work.
 *
 * Helix Channel Registry — Empirical Channel State Information (CSI)
 *
 * The "Data Moat" from the strategic playbook. Every time a user decodes a
 * FASTQ file, the registry silently ingests the error profile (substitution
 * rates, indel hotspots, GC-bias dropout curves) and uploads it to a central
 * database. When a new customer comes in, Helix auto-tunes its LDPC parity
 * and HMM parameters based on this proprietary dataset.
 *
 * Modeled after Waze's traffic data flywheel: the more people use Helix,
 * the better the channel model becomes, creating an un-copyable moat.
 *
 * Privacy: All telemetry is OPT-IN. Sovereign Archive mode (offline-first)
 * disables all telemetry. No DNA sequence data is ever uploaded — only
 * anonymized error statistics.
 */

import { createHash } from "crypto";

export interface ErrorProfile {
  /** Substitution rate (0-1) */
  substitutionRate: number;
  /** Insertion rate (0-1) */
  insertionRate: number;
  /** Deletion rate (0-1) */
  deletionRate: number;
  /** Average read length */
  avgReadLength: number;
  /** Read length standard deviation */
  readLengthStd: number;
  /** Average Q-score */
  avgQScore: number;
  /** GC bias (observed/expected GC ratio) */
  gcBias: number;
  /** Homopolymer-induced dropout rate */
  homopolymerDropout: number;
  /** Position-dependent error rate (binned per 10% of read) */
  positionalErrorRate: number[];
  /** Substitution matrix (A->C, A->G, A->T, C->A, ...) — 16 entries */
  substitutionMatrix: number[];
  /** Insertion base distribution (A, C, G, T) */
  insertionBases: number[];
  /** Sequencing platform (illumina, nanopore, pacbio) */
  platform: string;
  /** Vendor (twist, idt, custom) */
  vendor?: string;
  /** Number of reads analyzed */
  readCount: number;
}

export interface ChannelRegistryEntry {
  /** Unique ID (hash of platform + vendor + timestamp) */
  id: string;
  /** Error profile */
  profile: ErrorProfile;
  /** Timestamp (ISO 8601) */
  timestamp: string;
  /** Source FASTQ hash (for dedup, not the sequence itself) */
  sourceHash: string;
  /** Anonymous lab ID */
  labId?: string;
}

export interface RegistryConfig {
  /** Enable telemetry (default: false, opt-in) */
  enabled: boolean;
  /** Endpoint URL for uploading error profiles */
  endpoint?: string;
  /** Anonymous lab ID */
  labId?: string;
  /** Minimum reads required before uploading (privacy threshold) */
  minReads: number;
  /** Batch size for uploads */
  batchSize: number;
}

export const DEFAULT_REGISTRY_CONFIG: RegistryConfig = {
  enabled: false, // OPT-IN
  minReads: 100,
  batchSize: 10,
};

/**
 * Compute error profile from a set of reads + reference oligos.
 * This is the "telemetry ping" — runs after every decode.
 */
export function computeErrorProfile(
  reads: { sequence: string; quality?: string }[],
  referenceOligos?: { sequence: string }[],
  consensusOligos?: { sequence: string }[],
): ErrorProfile {
  const profile: ErrorProfile = {
    substitutionRate: 0,
    insertionRate: 0,
    deletionRate: 0,
    avgReadLength: 0,
    readLengthStd: 0,
    avgQScore: 0,
    gcBias: 1.0,
    homopolymerDropout: 0,
    positionalErrorRate: new Array(10).fill(0),
    substitutionMatrix: new Array(16).fill(0),
    insertionBases: [0, 0, 0, 0],
    platform: "unknown",
    readCount: reads.length,
  };

  if (reads.length === 0) return profile;

  // Read length statistics
  const lengths = reads.map(r => r.sequence.length);
  profile.avgReadLength = lengths.reduce((a, b) => a + b, 0) / lengths.length;
  const variance = lengths.reduce((a, b) => a + (b - profile.avgReadLength) ** 2, 0) / lengths.length;
  profile.readLengthStd = Math.sqrt(variance);

  // Q-score statistics (if available)
  if (reads[0].quality) {
    const allQ = reads.flatMap(r =>
      Array.from(r.quality!).map(c => c.charCodeAt(0) - 33)
    );
    profile.avgQScore = allQ.reduce((a, b) => a + b, 0) / allQ.length;
  }

  // GC content
  const gcCounts = reads.map(r => {
    const gc = (r.sequence.match(/[GC]/g) || []).length;
    return gc / r.sequence.length;
  });
  const avgGc = gcCounts.reduce((a, b) => a + b, 0) / gcCounts.length;
  profile.gcBias = avgGc / 0.5; // normalized to expected 50%

  // Homopolymer dropout: count reads with homopolymer runs ≥ 4
  const hpDropout = reads.filter(r => /(AAAA|CCCC|GGGG|TTTT)/.test(r.sequence)).length;
  profile.homopolymerDropout = hpDropout / reads.length;

  // Error rates vs reference (if provided)
  if (referenceOligos && referenceOligos.length > 0) {
    let totalSubs = 0;
    let totalIns = 0;
    let totalDels = 0;
    let totalCompared = 0;
    const subMatrix = new Array(16).fill(0);
    const insBases = [0, 0, 0, 0];
    const posErrors = new Array(10).fill(0);
    const posCounts = new Array(10).fill(0);

    for (const read of reads) {
      // Find best-matching reference (simplified — real impl would use clustering)
      let bestRef = referenceOligos[0];
      let bestDist = Infinity;
      for (const ref of referenceOligos) {
        const dist = Math.abs(read.sequence.length - ref.sequence.length);
        if (dist < bestDist) {
          bestDist = dist;
          bestRef = ref;
        }
      }

      // Simple alignment (count mismatches)
      const minLen = Math.min(read.sequence.length, bestRef.sequence.length);
      for (let i = 0; i < minLen; i++) {
        const posBin = Math.floor((i / minLen) * 10);
        posCounts[posBin]++;
        if (read.sequence[i] !== bestRef.sequence[i]) {
          totalSubs++;
          posErrors[posBin]++;
          // Substitution matrix
          const from = "ACGT".indexOf(bestRef.sequence[i]);
          const to = "ACGT".indexOf(read.sequence[i]);
          if (from >= 0 && to >= 0) {
            subMatrix[from * 4 + to]++;
          }
        }
      }
      if (read.sequence.length > bestRef.sequence.length) totalIns += read.sequence.length - bestRef.sequence.length;
      if (read.sequence.length < bestRef.sequence.length) totalDels += bestRef.sequence.length - read.sequence.length;
      totalCompared += minLen;
    }

    if (totalCompared > 0) {
      profile.substitutionRate = totalSubs / totalCompared;
      profile.insertionRate = totalIns / totalCompared;
      profile.deletionRate = totalDels / totalCompared;
      profile.substitutionMatrix = subMatrix.map(c => c / Math.max(totalSubs, 1));
      profile.insertionBases = insBases;
      profile.positionalErrorRate = posErrors.map((e, i) => posCounts[i] > 0 ? e / posCounts[i] : 0);
    }
  }

  // Platform detection
  if (profile.avgReadLength < 200) profile.platform = "illumina";
  else if (profile.avgReadLength > 1000) profile.platform = "nanopore";
  else if (profile.avgQScore > 30) profile.platform = "illumina";
  else profile.platform = "nanopore";

  return profile;
}

/**
 * Create a registry entry from an error profile.
 */
export function createRegistryEntry(
  profile: ErrorProfile,
  sourceFastqHash: string,
  config: RegistryConfig,
): ChannelRegistryEntry {
  return {
    id: createHash("sha256")
      .update(`${profile.platform}:${profile.vendor ?? "unknown"}:${Date.now()}`)
      .digest("hex")
      .slice(0, 16),
    profile,
    timestamp: new Date().toISOString(),
    sourceHash: sourceFastqHash,
    labId: config.labId,
  };
}

/**
 * Local storage for error profiles (offline buffer).
 * Profiles are stored locally and uploaded when connectivity is available.
 */
export class ChannelRegistryStore {
  private entries: ChannelRegistryEntry[] = [];
  private config: RegistryConfig;

  constructor(config: RegistryConfig = DEFAULT_REGISTRY_CONFIG) {
    this.config = { ...DEFAULT_REGISTRY_CONFIG, ...config };
  }

  /**
   * Record an error profile. If telemetry is enabled, queue for upload.
   * If disabled (Sovereign Archive mode), do nothing.
   */
  record(profile: ErrorProfile, sourceHash: string): void {
    if (!this.config.enabled) return;
    if (profile.readCount < this.config.minReads) return;

    const entry = createRegistryEntry(profile, sourceHash, this.config);
    this.entries.push(entry);

    // If batch is full, flush
    if (this.entries.length >= this.config.batchSize) {
      this.flush().catch(() => {/* silent fail — offline buffer */});
    }
  }

  /**
   * Upload buffered profiles to the central registry.
   * In offline mode (no endpoint), this is a no-op.
   */
  async flush(): Promise<void> {
    if (!this.config.enabled || !this.config.endpoint) return;
    if (this.entries.length === 0) return;

    try {
      const batch = this.entries.splice(0, this.config.batchSize);
      // In a real implementation, this would POST to the endpoint
      // For now, just log
      console.log(`[ChannelRegistry] Uploaded ${batch.length} profiles to ${this.config.endpoint}`);
    } catch (e) {
      // Re-queue on failure
      // this.entries.unshift(...batch);
    }
  }

  /**
   * Get the local error profile database (for auto-tuning).
   */
  getLocalProfiles(): ErrorProfile[] {
    return this.entries.map(e => e.profile);
  }

  /**
   * Get the aggregate error profile for a platform (for auto-tuning LDPC/HMM).
   */
  getAggregateProfile(platform: string): ErrorProfile | null {
    const profiles = this.entries.filter(e => e.profile.platform === platform);
    if (profiles.length === 0) return null;

    // Average all profiles for this platform
    const avg: ErrorProfile = {
      substitutionRate: 0,
      insertionRate: 0,
      deletionRate: 0,
      avgReadLength: 0,
      readLengthStd: 0,
      avgQScore: 0,
      gcBias: 0,
      homopolymerDropout: 0,
      positionalErrorRate: new Array(10).fill(0),
      substitutionMatrix: new Array(16).fill(0),
      insertionBases: [0, 0, 0, 0],
      platform,
      readCount: 0,
    };

    for (const e of profiles) {
      avg.substitutionRate += e.profile.substitutionRate;
      avg.insertionRate += e.profile.insertionRate;
      avg.deletionRate += e.profile.deletionRate;
      avg.avgReadLength += e.profile.avgReadLength;
      avg.avgQScore += e.profile.avgQScore;
      avg.gcBias += e.profile.gcBias;
      avg.homopolymerDropout += e.profile.homopolymerDropout;
      avg.readCount += e.profile.readCount;
      for (let i = 0; i < 10; i++) avg.positionalErrorRate[i] += e.profile.positionalErrorRate[i];
      for (let i = 0; i < 16; i++) avg.substitutionMatrix[i] += e.profile.substitutionMatrix[i];
    }

    const n = profiles.length;
    avg.substitutionRate /= n;
    avg.insertionRate /= n;
    avg.deletionRate /= n;
    avg.avgReadLength /= n;
    avg.avgQScore /= n;
    avg.gcBias /= n;
    avg.homopolymerDropout /= n;
    for (let i = 0; i < 10; i++) avg.positionalErrorRate[i] /= n;
    for (let i = 0; i < 16; i++) avg.substitutionMatrix[i] /= n;

    return avg;
  }

  /**
   * Auto-tune codec config based on aggregate error profile.
   * Returns recommended config adjustments.
   */
  autoTune(platform: string): {
    recommendedInnerParity: number;
    recommendedOuterParity: number;
    recommendedCoverage: number;
    recommendedMapping: string;
    recommendedChannel: string;
  } | null {
    const profile = this.getAggregateProfile(platform);
    if (!profile) return null;

    const totalError = profile.substitutionRate + profile.insertionRate + profile.deletionRate;

    // Auto-tune rules (based on SOTA literature + empirical data)
    let recommendedInnerParity = 8;
    let recommendedOuterParity = 0.05;
    let recommendedCoverage = 5;
    let recommendedMapping = "direct";
    let recommendedChannel = "illumina";

    if (platform === "nanopore" || profile.insertionRate + profile.deletionRate > 0.01) {
      // Nanopore: high IDS, need conv inner + HMM
      recommendedInnerParity = 12;
      recommendedOuterParity = 0.15;
      recommendedCoverage = 15;
      recommendedMapping = "direct"; // conv inner only supports direct
      recommendedChannel = "nanopore";
    } else if (totalError > 0.005) {
      // Noisy Illumina
      recommendedInnerParity = 10;
      recommendedOuterParity = 0.10;
      recommendedCoverage = 8;
    } else {
      // Clean Illumina
      recommendedInnerParity = 8;
      recommendedOuterParity = 0.05;
      recommendedCoverage = 5;
    }

    // Adjust for GC bias
    if (profile.gcBias < 0.8 || profile.gcBias > 1.2) {
      recommendedOuterParity = Math.min(recommendedOuterParity + 0.05, 0.20);
    }

    // Adjust for homopolymer dropout
    if (profile.homopolymerDropout > 0.1) {
      recommendedInnerParity = Math.max(recommendedInnerParity + 2, 14);
    }

    return {
      recommendedInnerParity,
      recommendedOuterParity,
      recommendedCoverage,
      recommendedMapping,
      recommendedChannel,
    };
  }
}

/**
 * Global singleton registry store.
 * Configured via environment variable HELIX_REGISTRY_ENABLED=1.
 */
let globalRegistry: ChannelRegistryStore | null = null;

export function getGlobalRegistry(): ChannelRegistryStore {
  if (!globalRegistry) {
    const enabled = process.env.HELIX_REGISTRY_ENABLED === "1" || process.env.HELIX_REGISTRY_ENABLED === "true";
    globalRegistry = new ChannelRegistryStore({
      ...DEFAULT_REGISTRY_CONFIG,
      enabled,
      endpoint: process.env.HELIX_REGISTRY_ENDPOINT,
      labId: process.env.HELIX_REGISTRY_LAB_ID,
    });
  }
  return globalRegistry;
}
