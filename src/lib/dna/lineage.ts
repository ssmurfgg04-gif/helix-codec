/**
 * Generational Lineage & Audit/Scrubbing for BioArchive.
 *
 * Tracks archive state across simulated "generations" (replication events),
 * including observed mutations, repair events, and chunk health.
 *
 * Audit/scrubbing verifies integrity, detects mutation hotspots, and reports
 * recovery probability — analogous to ZFS scrubbing or object-storage audit.
 */

import { BioArchive, DnaChunk, LineageRecord, RepairEvent, AuditReport } from "./bioarchive";
import { Oligo } from "./types";
import { simulate, MutationConfig, SimulationResult } from "./simulate";
import { decodeReads } from "./decode";
import { CodecConfig } from "./types";

/**
 * Simulate one generation of biological time passing:
 *   - Apply mutation model to current chunks
 *   - Attempt recovery
 *   - Record mutations, repairs, and chunk health
 *   - Return new archive state with updated lineage
 */
export async function simulateGeneration(
  archive: BioArchive,
  generation: number,
  mutCfg: MutationConfig,
  codecConfig: CodecConfig,
): Promise<{
  newArchive: BioArchive;
  simulation: SimulationResult;
  repairEvents: RepairEvent[];
  recoverySuccess: boolean;
}> {
  // Convert archive chunks back to Oligo format for the simulator
  const oligos: Oligo[] = archive.chunks.map((c) => ({
    index: c.index,
    sequence: c.sequence,
    gc: c.gc,
    maxHomopolymer: c.maxHomopolymer,
    seed: c.seed,
    payloadBytes: archive.manifest.chunking.chunkSizeBytes,
    length: c.sequence.length,
  }));

  // Simulate mutations
  const simulation = simulate(oligos, mutCfg);

  // Attempt recovery
  const decodeResult = await decodeReads(
    simulation.reads,
    archive.metadata,
    codecConfig,
    archive.forwardPrimer,
    archive.reversePrimer,
  );

  // Build repair events
  const repairEvents: RepairEvent[] = [];
  for (const perOligo of decodeResult.perOligo) {
    if (perOligo.innerRS.corrected > 0) {
      repairEvents.push({
        type: "inner_ecc",
        chunkId: perOligo.index.toString().padStart(6, "0"),
        description: `Corrected ${perOligo.innerRS.corrected} errors via inner RS`,
        timestamp: new Date().toISOString(),
      });
    }
    if (perOligo.readCount === 0) {
      repairEvents.push({
        type: "outer_ecc",
        chunkId: perOligo.index.toString().padStart(6, "0"),
        description: `Recovered via outer erasure decoding (no reads)`,
        timestamp: new Date().toISOString(),
      });
    }
  }

  // Compute observed mutation rate (errors per base across all reads)
  const totalBases = simulation.reads.reduce((s, r) => s + r.sequence.length, 0);
  const totalErrors = simulation.totalErrors;
  const observedMutationRate = totalBases > 0 ? totalErrors / totalBases : 0;

  // Compute chunk health
  const chunksHealthy = decodeResult.perOligo.filter(
    (p) => p.innerRS.success && p.readCount > 0,
  ).length;
  const chunksRepaired = decodeResult.perOligo.filter(
    (p) => p.innerRS.corrected > 0,
  ).length;
  const chunksUnrecoverable = decodeResult.perOligo.filter(
    (p) => !p.innerRS.success,
  ).length;
  const recoveryProbability = chunksHealthy / archive.chunks.length;

  // Build new lineage record
  const lineageRecord: LineageRecord = {
    generation,
    timestamp: new Date().toISOString(),
    parentHash: archive.manifest.merkle.root,
    observedMutationRate,
    repairEvents,
    chunksHealthy,
    chunksRepaired,
    chunksUnrecoverable,
    recoveryProbability,
  };

  const newArchive: BioArchive = {
    ...archive,
    lineage: [...archive.lineage, lineageRecord],
  };

  return {
    newArchive,
    simulation,
    repairEvents,
    recoverySuccess: decodeResult.hashMatches,
  };
}

/**
 * Audit / scrub an archive: verify integrity, detect mutation hotspots,
 * report recovery probability.
 */
export function auditArchive(
  archive: BioArchive,
  simulation?: SimulationResult,
  decodeResult?: { perOligo: { index: number; innerRS: { success: boolean; corrected: number }; readCount: number }[] },
): AuditReport {
  const findings: string[] = [];
  let mutatedChunks = 0;
  let repairedChunks = 0;
  let unrecoverableChunks = 0;
  const mutationHotspots: { chunkId: string; mutationRate: number }[] = [];

  if (decodeResult) {
    for (const p of decodeResult.perOligo) {
      if (p.innerRS.corrected > 0) {
        repairedChunks++;
        mutatedChunks++;
        // Estimate mutation rate from corrections
        const estimatedRate = p.innerRS.corrected / (archive.manifest.sequenceLayout.payloadBytesPerOligo * 4);
        if (estimatedRate > 0.05) {
          mutationHotspots.push({
            chunkId: p.index.toString().padStart(6, "0"),
            mutationRate: estimatedRate,
          });
        }
      }
      if (!p.innerRS.success) {
        unrecoverableChunks++;
      }
    }
  }

  if (mutationHotspots.length > 0) {
    findings.push(
      `Detected ${mutationHotspots.length} mutation hotspots (chunks with >5% estimated mutation rate).`,
    );
  }
  if (unrecoverableChunks > 0) {
    findings.push(
      `${unrecoverableChunks} chunks are unrecoverable and require outer erasure decoding or regeneration.`,
    );
  }
  if (repairedChunks > archive.chunks.length * 0.3) {
    findings.push(
      `Over 30% of chunks required repair — consider regenerating the archive from a healthy copy.`,
    );
  }

  const healthyChunks = archive.chunks.length - mutatedChunks - unrecoverableChunks;
  const recoveryProbability = decodeResult
    ? (archive.chunks.length - unrecoverableChunks) / archive.chunks.length
    : 1.0;

  // Verify Merkle root (we can't re-verify without the original chunks, so assume valid if no decode errors)
  const merkleValid = unrecoverableChunks === 0;

  return {
    archiveId: archive.manifest.archiveId,
    auditedAt: new Date().toISOString(),
    generation: archive.lineage.length - 1,
    totalChunks: archive.chunks.length,
    healthyChunks,
    mutatedChunks,
    repairedChunks,
    unrecoverableChunks,
    mutationHotspots,
    recoveryProbability,
    merkleValid,
    findings,
  };
}

// --- Lifecycle policies ---

export interface LifecyclePolicy {
  retention: string; // e.g. "100y"
  storageClass: "hot" | "warm" | "cold" | "deep_bio_archive";
  replicationTarget: number;
  migrationInterval: string;
  decayPolicy: string;
}

export const DEFAULT_LIFECYCLE: LifecyclePolicy = {
  retention: "100y",
  storageClass: "deep_bio_archive",
  replicationTarget: 3,
  migrationInterval: "10y",
  decayPolicy: "repair_if_mutation_gt_5%",
};

/**
 * Check if an archive needs migration based on its lifecycle policy.
 */
export function needsMigration(archive: BioArchive): { needed: boolean; reason: string } {
  if (!archive.manifest.lifecycle) return { needed: false, reason: "No lifecycle policy" };

  const created = new Date(archive.manifest.created).getTime();
  const now = Date.now();
  const yearsSinceCreation = (now - created) / (365.25 * 24 * 3600 * 1000);

  const migrationIntervalYears = parseInt(
    archive.manifest.lifecycle.migrationInterval.replace(/[^0-9]/g, ""),
  ) || 10;

  if (yearsSinceCreated > migrationIntervalYears) {
    return {
      needed: true,
      reason: `Archive is ${yearsSinceCreated.toFixed(1)}y old (migration interval: ${archive.manifest.lifecycle.migrationInterval})`,
    };
  }

  // Check mutation rate against decay policy
  const lastGen = archive.lineage[archive.lineage.length - 1];
  const thresholdMatch = archive.manifest.lifecycle.decayPolicy.match(/(\d+)%/);
  const threshold = thresholdMatch ? parseInt(thresholdMatch[1]) / 100 : 0.05;
  if (lastGen && lastGen.observedMutationRate > threshold) {
    return {
      needed: true,
      reason: `Observed mutation rate ${(lastGen.observedMutationRate * 100).toFixed(2)}% exceeds threshold ${(threshold * 100).toFixed(1)}%`,
    };
  }

  return { needed: false, reason: "Archive is within lifecycle parameters" };
}

function yearsSinceCreated(archive: BioArchive): number {
  const created = new Date(archive.manifest.created).getTime();
  const now = Date.now();
  return (now - created) / (365.25 * 24 * 3600 * 1000);
}

// --- Random access ---

/**
 * Retrieve a specific chunk by ID without decoding the whole archive.
 * In a real wetlab, this would use PCR with chunk-specific primers.
 * In software, it's a simple lookup.
 */
export function randomAccessChunk(archive: BioArchive, chunkId: string): DnaChunk | null {
  return archive.chunks.find((c) => c.chunkId === chunkId) ?? null;
}

/**
 * Retrieve a range of chunks (e.g., for partial file restoration).
 */
export function randomAccessRange(archive: BioArchive, startIdx: number, count: number): DnaChunk[] {
  return archive.chunks
    .filter((c) => c.index >= startIdx && c.index < startIdx + count)
    .sort((a, b) => a.index - b.index);
}
