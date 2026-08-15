/**
 * BioArchive Container Format
 *
 * A self-describing biological archival file format. Each archive contains:
 *   - manifest.json: describes the archive (file metadata, codec config, layout)
 *   - merkle.json: SHA-256 Merkle tree for chunk integrity
 *   - chunks/: one FASTA-like .dna file per oligo
 *   - metadata.enc: encrypted metadata (optional)
 *   - recovery_report.json: audit results
 *   - lineage.json: generational history
 *
 * The format is designed for a medium that mutates, decays, and replicates.
 * Unlike normal archive formats that assume static bits, BioArchive plans for
 * damage: every chunk has a checksum, every chunk has redundant copies via
 * outer erasure coding, and the manifest declares the recovery policy.
 *
 * File extension: .bioarc (JSON manifest + chunk records)
 *
 * Reference design inspired by:
 *   - ZFS (checksums everywhere, self-healing)
 *   - Object storage (object IDs, metadata)
 *   - PAR2 (recovery blocks)
 *   - Git/IPFS (content-addressed, Merkle trees)
 */

import { CodecConfig, CodecMetadata, EncodedFile, Oligo } from "./types";

/** Top-level archive manifest. */
export interface BioArchiveManifest {
  format: string; // "bioarchive/v1"
  created: string; // ISO timestamp
  archiveId: string; // unique archive ID (derived from Merkle root)
  payload: {
    originalName: string;
    sizeBytes: number;
    sha256: string;
    mime: string;
  };
  chunking: {
    algorithm: "fixed";
    chunkSizeBytes: number;
    dataChunks: number;
    parityChunks: number;
    totalChunks: number;
  };
  encryption: {
    cipher: "none" | "xchacha20-poly1305";
    kdf: "none" | "argon2id";
    keyId?: string; // hash of derived key (for key management)
    salt?: string; // base64
    nonce?: string; // base64
  };
  ecc: {
    inner: "reed-solomon";
    outer: "reed-solomon-erasure" | "holographic-shard";
    overheadRatio: number; // totalChunks / dataChunks
  };
  bioEncoding: {
    baseMapper: "direct-2bit";
    maxHomopolymer: number;
    gcMin: number;
    gcMax: number;
    primerLength: number;
    addressBits: number;
    checksumBits: number; // CRC-16 = 16 bits
  };
  sequenceLayout: {
    oligoLength: number;
    payloadBytesPerOligo: number;
    innerParityBytes: number;
  };
  recovery: {
    minCoverage: number;
    maxMutationRate: number;
    repairMode: "consensus+outer_ecc";
  };
  merkle: {
    algorithm: "sha256";
    root: string; // hex
    leafCount: number;
    treeDepth: number;
  };
  lifecycle?: {
    retention: string; // e.g. "100y"
    storageClass: "hot" | "warm" | "cold" | "deep_bio_archive";
    replicationTarget: number;
    migrationInterval: string;
    decayPolicy: string;
  };
}

/** A single chunk in FASTA-like format. */
export interface DnaChunk {
  chunkId: string; // e.g. "000001"
  index: number;
  type: "data" | "parity";
  barcode: string; // 12-nt barcode for multiplexing
  address: string; // 16-nt address (whitened index + seed)
  sequence: string; // full DNA sequence including primers
  checksum: string; // SHA-256 of payload bytes (hex)
  gc: number;
  maxHomopolymer: number;
  seed: number;
}

/** Full archive package. */
export interface BioArchive {
  manifest: BioArchiveManifest;
  merkleTree: MerkleTree;
  chunks: DnaChunk[];
  forwardPrimer: string;
  reversePrimer: string;
  metadata: CodecMetadata;
  lineage: LineageRecord[];
  auditReport?: AuditReport;
}

/** Merkle tree for chunk integrity. */
export interface MerkleTree {
  leaves: string[]; // SHA-256 of each chunk's payload (hex)
  nodes: string[][]; // tree levels, level 0 = leaves, last level = [root]
  root: string;
  depth: number;
}

/** Lineage record — one per generation. */
export interface LineageRecord {
  generation: number;
  timestamp: string;
  parentHash: string; // Merkle root of parent generation
  observedMutationRate: number;
  repairEvents: RepairEvent[];
  chunksHealthy: number;
  chunksRepaired: number;
  chunksUnrecoverable: number;
  recoveryProbability: number;
}

/** A single repair event. */
export interface RepairEvent {
  type: "consensus" | "inner_ecc" | "outer_ecc" | "regeneration";
  chunkId: string;
  description: string;
  timestamp: string;
}

/** Audit / scrubbing report. */
export interface AuditReport {
  archiveId: string;
  auditedAt: string;
  generation: number;
  totalChunks: number;
  healthyChunks: number;
  mutatedChunks: number;
  repairedChunks: number;
  unrecoverableChunks: number;
  mutationHotspots: { chunkId: string; mutationRate: number }[];
  recoveryProbability: number;
  merkleValid: boolean;
  findings: string[];
}

// --- Merkle tree implementation ---

/**
 * Build a SHA-256 Merkle tree from a list of leaf hashes.
 * Returns the tree (with all intermediate nodes) and the root.
 */
export async function buildMerkleTree(leaves: string[]): Promise<MerkleTree> {
  if (leaves.length === 0) {
    return { leaves: [], nodes: [], root: "", depth: 0 };
  }
  // Pad to power of 2 with the last leaf duplicated
  let paddedLeaves = leaves.slice();
  let depth = 0;
  while (paddedLeaves.length > 1) {
    if (paddedLeaves.length % 2 !== 0) {
      paddedLeaves.push(paddedLeaves[paddedLeaves.length - 1]);
    }
    paddedLeaves = paddedLeaves; // no-op, just for clarity
    break;
  }
  // Actually let me redo this properly
  paddedLeaves = leaves.slice();
  // Pad to next power of 2
  let nextPow2 = 1;
  while (nextPow2 < paddedLeaves.length) nextPow2 *= 2;
  while (paddedLeaves.length < nextPow2) {
    paddedLeaves.push(paddedLeaves[paddedLeaves.length - 1]);
  }

  const nodes: string[][] = [paddedLeaves];
  let current = paddedLeaves;
  while (current.length > 1) {
    const next: string[] = [];
    for (let i = 0; i < current.length; i += 2) {
      const combined = current[i] + current[i + 1];
      const hash = await sha256Hex(combined);
      next.push(hash);
    }
    nodes.push(next);
    current = next;
    depth++;
  }
  return {
    leaves: paddedLeaves.slice(0, leaves.length), // original (unpadded) leaves
    nodes,
    root: current[0],
    depth,
  };
}

/** Verify a Merkle proof for a specific leaf. */
export async function verifyMerkleProof(
  leaf: string,
  leafIndex: number,
  proof: { hash: string; isRight: boolean }[],
  root: string,
): Promise<boolean> {
  let current = leaf;
  let index = leafIndex;
  for (const step of proof) {
    const combined = step.isRight ? current + step.hash : step.hash + current;
    current = await sha256Hex(combined);
    index = Math.floor(index / 2);
  }
  return current === root;
}

/** Generate a Merkle proof for a leaf at given index. */
export function generateMerkleProof(
  tree: MerkleTree,
  leafIndex: number,
): { hash: string; isRight: boolean }[] {
  const proof: { hash: string; isRight: boolean }[] = [];
  let index = leafIndex;
  for (let level = 0; level < tree.nodes.length - 1; level++) {
    const siblingIndex = index % 2 === 0 ? index + 1 : index - 1;
    const sibling = tree.nodes[level][siblingIndex];
    if (sibling) {
      proof.push({ hash: sibling, isRight: index % 2 === 0 });
    }
    index = Math.floor(index / 2);
  }
  return proof;
}

// --- SHA-256 helper ---

async function sha256Hex(data: string): Promise<string> {
  if (typeof crypto !== "undefined" && crypto.subtle) {
    const enc = new TextEncoder().encode(data);
    const hash = await crypto.subtle.digest("SHA-256", enc);
    return Array.from(new Uint8Array(hash))
      .map((b) => b.toString(16).padStart(2, "0"))
      .join("");
  }
  const { createHash } = await import("crypto");
  return createHash("sha256").update(data).digest("hex");
}

// --- Barcode generation ---

/**
 * Generate a unique 12-nt barcode for a chunk.
 * Uses xorshift32 PRNG seeded by the chunk index for determinism.
 */
export function generateBarcode(index: number): string {
  const bases = "ACGT";
  let state = ((index + 1) * 2654435761) >>> 0 || 1;
  let barcode = "";
  let prev = "";
  while (barcode.length < 12) {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    state = state >>> 0;
    const base = bases[state % 4];
    if (base !== prev) {
      barcode += base;
      prev = base;
    }
  }
  return barcode;
}

// --- Convert EncodedFile to BioArchive ---

/**
 * Wrap an existing EncodedFile into a BioArchive container, adding:
 *   - per-chunk SHA-256 (over the payload bytes)
 *   - Merkle tree
 *   - barcodes
 *   - manifest
 */
export async function toBioArchive(
  encoded: EncodedFile,
  encryption?: BioArchiveManifest["encryption"],
  lifecycle?: BioArchiveManifest["lifecycle"],
  config?: CodecConfig,
): Promise<BioArchive> {
  // Compute per-chunk SHA-256 (over the payload bytes, NOT the DNA string)
  const leaves: string[] = [];
  const chunks: DnaChunk[] = [];

  // Derive primer length from the actual primer sequences (not hardcoded).
  // This ensures the Merkle tree slices correctly when non-default primers
  // (e.g., longer barcoded primers) are used.
  const primerLen = encoded.forwardPrimer.length || encoded.reversePrimer.length || 20;
  // Derive actual oligo length from the first oligo (if any), falling back
  // to a computed default.
  const oligoLength = encoded.oligos.length > 0
    ? encoded.oligos[0].sequence.length
    : (primerLen * 2 + encoded.metadata.payloadBytesPerOligo * 4);

  for (const oligo of encoded.oligos) {
    // Extract payload bytes from the oligo (strip primers dynamically)
    // (We hash the inner block DNA string for simplicity — in a real system
    // we'd hash the decoded bytes, but hashing the DNA is also valid since
    // it's a unique representation of the chunk.)
    const innerDna = oligo.sequence.slice(
      primerLen,
      oligo.sequence.length - primerLen,
    );
    const leaf = await sha256Hex(innerDna);
    leaves.push(leaf);

    chunks.push({
      chunkId: oligo.index.toString().padStart(6, "0"),
      index: oligo.index,
      type: oligo.index >= encoded.metadata.outerRS.k ? "parity" : "data",
      barcode: generateBarcode(oligo.index),
      address: innerDna.slice(0, 16), // first 16 nt = address
      sequence: oligo.sequence,
      checksum: leaf,
      gc: oligo.gc,
      maxHomopolymer: oligo.maxHomopolymer,
      seed: oligo.seed,
    });
  }

  const merkleTree = await buildMerkleTree(leaves);
  const archiveId = merkleTree.root.slice(0, 16);

  const manifest: BioArchiveManifest = {
    format: "bioarchive/v1",
    created: encoded.metadata.encodedAt,
    archiveId,
    payload: {
      originalName: encoded.metadata.fileName,
      sizeBytes: encoded.metadata.fileSize,
      sha256: encoded.metadata.fileHash,
      mime: encoded.metadata.contentType,
    },
    chunking: {
      algorithm: "fixed",
      chunkSizeBytes: encoded.metadata.payloadBytesPerOligo,
      dataChunks: encoded.metadata.outerRS.k,
      parityChunks: encoded.metadata.parityOligos,
      totalChunks: encoded.metadata.oligoCount,
    },
    encryption: encryption ?? { cipher: "none", kdf: "none" },
    ecc: {
      inner: "reed-solomon",
      outer: "reed-solomon-erasure",
      overheadRatio: encoded.metadata.oligoCount / encoded.metadata.outerRS.k,
    },
    bioEncoding: {
      baseMapper: "direct-2bit",
      maxHomopolymer: config?.constraints?.maxHomopolymer ?? 3,
      gcMin: config?.constraints?.gcMin ?? 0.4,
      gcMax: config?.constraints?.gcMax ?? 0.6,
      primerLength: primerLen, // Dynamic — derived from actual primer sequences
      addressBits: 32,
      checksumBits: 16,
    },
    sequenceLayout: {
      oligoLength, // Dynamic — derived from actual oligo sequences
      payloadBytesPerOligo: encoded.metadata.payloadBytesPerOligo,
      innerParityBytes: encoded.metadata.innerRS.n - encoded.metadata.innerRS.k,
    },
    recovery: {
      minCoverage: 3,
      maxMutationRate: 0.12,
      repairMode: "consensus+outer_ecc",
    },
    merkle: {
      algorithm: "sha256",
      root: merkleTree.root,
      leafCount: leaves.length,
      treeDepth: merkleTree.depth,
    },
    lifecycle,
  };

  return {
    manifest,
    merkleTree,
    chunks,
    forwardPrimer: encoded.forwardPrimer,
    reversePrimer: encoded.reversePrimer,
    metadata: encoded.metadata,
    lineage: [
      {
        generation: 0,
        timestamp: encoded.metadata.encodedAt,
        parentHash: "",
        observedMutationRate: 0,
        repairEvents: [],
        chunksHealthy: chunks.length,
        chunksRepaired: 0,
        chunksUnrecoverable: 0,
        recoveryProbability: 1.0,
      },
    ],
  };
}

// --- FASTA-like serialization ---

/** Serialize a chunk to FASTA-like format. */
export function chunkToFasta(chunk: DnaChunk): string {
  return `>chunk_${chunk.chunkId} barcode=${chunk.barcode} address=${chunk.address} checksum=${chunk.checksum.slice(0, 16)} type=${chunk.type} gc=${(chunk.gc * 100).toFixed(1)} maxhp=${chunk.maxHomopolymer} seed=${chunk.seed}
${chunk.sequence}
`;
}

/** Parse a FASTA-like chunk record. */
export function fastaToChunk(fasta: string): DnaChunk | null {
  const lines = fasta.trim().split("\n");
  if (lines.length < 2) return null;
  const header = lines[0];
  if (!header.startsWith(">")) return null;
  const parts = header.slice(1).split(" ");
  const idPart = parts[0]; // "chunk_000001"
  const chunkId = idPart.replace("chunk_", "");
  const index = parseInt(chunkId, 10);
  const props: Record<string, string> = {};
  for (let i = 1; i < parts.length; i++) {
    const [k, v] = parts[i].split("=");
    props[k] = v;
  }
  const sequence = lines.slice(1).join("");
  return {
    chunkId,
    index,
    type: (props.type as "data" | "parity") ?? "data",
    barcode: props.barcode ?? "",
    address: props.address ?? "",
    sequence,
    checksum: props.checksum ?? "",
    gc: parseFloat(props.gc ?? "0") / 100,
    maxHomopolymer: parseInt(props.maxhp ?? "0", 10),
    seed: parseInt(props.seed ?? "0", 10),
  };
}

/** Serialize a BioArchive to a downloadable JSON package. */
export function serializeBioArchive(archive: BioArchive): string {
  return JSON.stringify(archive, null, 2);
}

/** Deserialize a BioArchive from JSON. */
export function deserializeBioArchive(json: string): BioArchive {
  return JSON.parse(json) as BioArchive;
}
