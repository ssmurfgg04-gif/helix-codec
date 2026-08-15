/**
 * 
 * ⚠️ EXPERIMENTAL / ROADMAP — This module is a scaffolding prototype.
 * It is NOT wired into the production encode/decode pipeline.
 * See CHANGELOG.md for v59 status and remaining work.
 *
 * Helix Bio-DFS — Biological Distributed File System
 *
 * Evolves Helix from a "file encoder" into a Biological Distributed File System.
 * Think Hadoop/IPFS but data shards live in physical DNA pools (micro-vials)
 * instead of hard drives.
 *
 * A 1 Petabyte archive at 1.76 bits/nt requires ~5 grams of DNA across thousands
 * of micro-tubes. Bio-DFS manages this as a distributed file system:
 *   - Splits large files into "pools" (each pool = one test tube)
 *   - Each pool has a manifest (which oligos it contains)
 *   - Pools are organized in a directory hierarchy (pool groups)
 *   - Physical labels (QR codes) map to pool IDs
 *   - Metadata (file system tree) is stored separately (small DNA archive or digital)
 *
 * This is the foundation for the "African Sovereign Archive" — solar-powered,
 * off-grid biological biobanks that store national archives in DNA.
 */

import { createHash } from "crypto";

export interface PoolManifest {
  /** Unique pool ID (deterministic from archive ID + pool index) */
  poolId: string;
  /** Archive ID (parent) */
  archiveId: string;
  /** Pool index within the archive */
  poolIndex: number;
  /** Total pools in the archive */
  totalPools: number;
  /** Oligo count in this pool */
  oligoCount: number;
  /** Data capacity per oligo (bytes) */
  payloadBytesPerOligo: number;
  /** Total data in this pool (bytes) */
  poolDataSize: number;
  /** Oligo length (nt) */
  oligoLength: number;
  /** Primer sequences (forward + reverse) */
  forwardPrimer: string;
  reversePrimer: string;
  /** LDPC/RS parameters */
  innerN: number;
  innerK: number;
  outerN: number;
  outerK: number;
  /** Pool hash (SHA-256 of all oligo sequences, for integrity) */
  poolHash: string;
  /** Physical label (QR code data) */
  physicalLabel: string;
  /** Storage conditions */
  storageConditions: {
    temperature: number; // °C
    humidity: number; // %
    desiccant: boolean;
  };
}

export interface ArchiveManifest {
  /** Unique archive ID */
  archiveId: string;
  /** Original file name */
  fileName: string;
  /** Original file size (bytes) */
  fileSize: number;
  /** File hash (SHA-256) */
  fileHash: string;
  /** Number of pools */
  poolCount: number;
  /** Total oligos across all pools */
  totalOligos: number;
  /** Oligo length (nt) */
  oligoLength: number;
  /** Net density (bits/nt) */
  netDensity: number;
  /** Physical DNA weight estimate (grams) */
  dnaWeightGrams: number;
  /** Pool manifests (one per pool) */
  pools: PoolManifest[];
  /** Creation timestamp */
  createdAt: string;
  /** Expiry (for sovereign archives — "never" by default) */
  expiresAt: string;
  /** Encryption status */
  encrypted: boolean;
  /** Sovereign mode flag (offline, no telemetry) */
  sovereign: boolean;
}

export interface BioDFSConfig {
  /** Maximum oligos per pool (physical limit of synthesis) */
  maxOligosPerPool: number;
  /** Oligo length (nt) */
  oligoLength: number;
  /** Primer length (nt) */
  primerLength: number;
  /** Payload bytes per oligo (computed from layout) */
  payloadBytesPerOligo: number;
  /** Outer RS parity fraction */
  outerParityFraction: number;
  /** Physical storage conditions */
  storageConditions: {
    temperature: number;
    humidity: number;
    desiccant: boolean;
  };
  /** Sovereign mode (offline, no telemetry, no cloud) */
  sovereign: boolean;
  /** Pool naming convention */
  poolNamingConvention: "sequential" | "hash-based";
}

export const DEFAULT_BIO_DFS_CONFIG: BioDFSConfig = {
  maxOligosPerPool: 100000, // ~100K oligos per pool (typical synthesis limit)
  oligoLength: 500,
  primerLength: 20,
  payloadBytesPerOligo: 96, // after v53 -3 fix
  outerParityFraction: 0.05,
  storageConditions: {
    temperature: 25, // room temperature
    humidity: 30,
    desiccant: true,
  },
  sovereign: true, // default to sovereign mode (African context)
  poolNamingConvention: "hash-based",
};

/**
 * Compute the number of pools needed for a given file size.
 *
 * At 1.76 bits/nt, 1 PB requires ~9 Trillion oligos.
 * At 100K oligos/pool, that's 90 Million pools.
 * At 5 grams DNA total, each pool has ~55 ng DNA.
 */
export function computePoolCount(
  fileSize: number,
  config: BioDFSConfig = DEFAULT_BIO_DFS_CONFIG,
): {
  poolCount: number;
  oligosPerPool: number;
  totalOligos: number;
  dnaWeightGrams: number;
  netDensity: number;
} {
  const payloadBits = fileSize * 8;
  const bitsPerOligo = config.payloadBytesPerOligo * 8;
  const outerEfficiency = 1 / (1 + config.outerParityFraction);
  const infoBitsPerOligo = bitsPerOligo * outerEfficiency;
  const totalOligos = Math.ceil(payloadBits / infoBitsPerOligo);
  const poolCount = Math.ceil(totalOligos / config.maxOligosPerPool);
  const oligosPerPool = Math.ceil(totalOligos / poolCount);

  // DNA weight: ~1 gram dsDNA ≈ 215 PB theoretical, ~57 EB practical
  // Use practical: 1 gram = 57 EB = 57 * 10^18 bytes
  const dnaWeightGrams = fileSize / (57e18);

  // Net density
  const totalNt = totalOligos * config.oligoLength;
  const netDensity = payloadBits / totalNt;

  return {
    poolCount,
    oligosPerPool,
    totalOligos,
    dnaWeightGrams,
    netDensity,
  };
}

/**
 * Generate a pool ID from archive ID + pool index.
 * Uses hash-based naming for sovereignty (no sequential correlation).
 */
export function generatePoolId(archiveId: string, poolIndex: number): string {
  const hash = createHash("sha256")
    .update(`${archiveId}:${poolIndex}`)
    .digest("hex")
    .slice(0, 12);
  return `pool_${hash}`;
}

/**
 * Generate a physical label (QR code data) for a pool.
 * Format: helix://<archiveId>/<poolId>/<poolIndex>/<totalPools>
 */
export function generatePhysicalLabel(
  archiveId: string,
  poolId: string,
  poolIndex: number,
  totalPools: number,
): string {
  return `helix://${archiveId}/${poolId}/${poolIndex + 1}/${totalPools}`;
}

/**
 * Create the full archive manifest for a large file.
 * This is the "directory" of the Bio-DFS — maps file → pools → oligos.
 */
export function createArchiveManifest(
  fileName: string,
  fileSize: number,
  fileHash: string,
  config: BioDFSConfig = DEFAULT_BIO_DFS_CONFIG,
): ArchiveManifest {
  const { poolCount, oligosPerPool, totalOligos, dnaWeightGrams, netDensity } =
    computePoolCount(fileSize, config);

  const archiveId = createHash("sha256")
    .update(`${fileName}:${fileSize}:${fileHash}:${Date.now()}`)
    .digest("hex")
    .slice(0, 16);

  const pools: PoolManifest[] = [];
  for (let i = 0; i < poolCount; i++) {
    const poolId = generatePoolId(archiveId, i);
    const isLast = i === poolCount - 1;
    const oligosThisPool = isLast
      ? totalOligos - (poolCount - 1) * oligosPerPool
      : oligosPerPool;
    const poolDataSize = oligosThisPool * config.payloadBytesPerOligo;

    pools.push({
      poolId,
      archiveId,
      poolIndex: i,
      totalPools: poolCount,
      oligoCount: oligosThisPool,
      payloadBytesPerOligo: config.payloadBytesPerOligo,
      poolDataSize,
      oligoLength: config.oligoLength,
      forwardPrimer: "AGTACTGCTGTGTATGTACG", // default primers
      reversePrimer: "GATATATGCGATACATATC", // (would be customized per archive)
      innerN: 106,
      innerK: 98,
      outerN: Math.ceil(oligosThisPool * (1 + config.outerParityFraction)),
      outerK: oligosThisPool,
      poolHash: "", // computed after oligos are generated
      physicalLabel: generatePhysicalLabel(archiveId, poolId, i, poolCount),
      storageConditions: config.storageConditions,
    });
  }

  return {
    archiveId,
    fileName,
    fileSize,
    fileHash,
    poolCount,
    totalOligos,
    oligoLength: config.oligoLength,
    netDensity,
    dnaWeightGrams,
    pools,
    createdAt: new Date().toISOString(),
    expiresAt: "never",
    encrypted: false,
    sovereign: config.sovereign,
  };
}

/**
 * Serialize an archive manifest to JSON (for digital storage or small DNA archive).
 */
export function serializeManifest(manifest: ArchiveManifest): string {
  return JSON.stringify(manifest, null, 2);
}

/**
 * Deserialize an archive manifest from JSON.
 */
export function deserializeManifest(json: string): ArchiveManifest {
  return JSON.parse(json) as ArchiveManifest;
}

/**
 * Compute the physical storage requirements for an archive.
 * Useful for planning lab space, shipping containers, etc.
 */
export function computeStorageRequirements(manifest: ArchiveManifest): {
  poolCount: number;
  tubeCount: number;
  rackCount: number;
  shippingContainerCount: number;
  estimatedWeightKg: number;
  estimatedVolumeLiters: number;
} {
  const poolCount = manifest.poolCount;
  const tubeCount = poolCount; // 1 tube per pool
  const rackCount = Math.ceil(tubeCount / 96); // 96-well plates
  const shippingContainerCount = Math.ceil(rackCount / 100); // 100 racks per container
  const estimatedWeightKg = manifest.dnaWeightGrams * 1000 + tubeCount * 0.002; // DNA + tubes
  const estimatedVolumeLiters = tubeCount * 0.0001 + rackCount * 0.001; // tubes + racks

  return {
    poolCount,
    tubeCount,
    rackCount,
    shippingContainerCount,
    estimatedWeightKg,
    estimatedVolumeLiters,
  };
}

/**
 * Format file size for human-readable display.
 */
export function formatFileSize(bytes: number): string {
  const units = ["B", "KB", "MB", "GB", "TB", "PB", "EB"];
  let size = bytes;
  let unitIndex = 0;
  while (size >= 1024 && unitIndex < units.length - 1) {
    size /= 1024;
    unitIndex++;
  }
  return `${size.toFixed(2)} ${units[unitIndex]}`;
}

/**
 * Format DNA weight for human-readable display.
 */
export function formatDnaWeight(grams: number): string {
  if (grams < 1e-6) return `${(grams * 1e9).toFixed(2)} ng`;
  if (grams < 1e-3) return `${(grams * 1e6).toFixed(2)} µg`;
  if (grams < 1) return `${(grams * 1e3).toFixed(2)} mg`;
  return `${grams.toFixed(2)} g`;
}

/**
 * Generate a summary report for an archive (for lab documentation).
 */
export function generateArchiveReport(manifest: ArchiveManifest): string {
  const storage = computeStorageRequirements(manifest);
  return `
╔══════════════════════════════════════════════════════════════╗
║  HELIX BIO-DFS ARCHIVE REPORT                                ║
╠══════════════════════════════════════════════════════════════╣
  Archive ID:     ${manifest.archiveId}
  File Name:      ${manifest.fileName}
  File Size:      ${formatFileSize(manifest.fileSize)}
  File Hash:      ${manifest.fileHash.slice(0, 32)}...
  Created:        ${manifest.createdAt}
  Expires:        ${manifest.expiresAt}
  Encrypted:      ${manifest.encrypted ? "Yes" : "No"}
  Sovereign:      ${manifest.sovereign ? "Yes (offline, no telemetry)" : "No"}

  ─── Physical Layout ───
  Pool Count:     ${manifest.poolCount}
  Total Oligos:   ${manifest.totalOligos.toLocaleString()}
  Oligo Length:   ${manifest.oligoLength} nt
  Net Density:    ${manifest.netDensity.toFixed(3)} bits/nt
  DNA Weight:     ${formatDnaWeight(manifest.dnaWeightGrams)}

  ─── Storage Requirements ───
  Tubes:          ${storage.tubeCount.toLocaleString()}
  Racks (96-well): ${storage.rackCount.toLocaleString()}
  Ship Containers: ${storage.shippingContainerCount.toLocaleString()}
  Est. Weight:    ${storage.estimatedWeightKg.toFixed(2)} kg
  Est. Volume:    ${storage.estimatedVolumeLiters.toFixed(2)} L

  ─── Pool Distribution ───
  Pool 1:         ${manifest.pools[0]?.poolId} (${manifest.pools[0]?.oligoCount.toLocaleString()} oligos)
  Pool 2:         ${manifest.pools[1]?.poolId ?? "N/A"}
  ...
  Pool ${manifest.poolCount}: ${manifest.pools[manifest.poolCount - 1]?.poolId}

╚══════════════════════════════════════════════════════════════╝
`;
}
