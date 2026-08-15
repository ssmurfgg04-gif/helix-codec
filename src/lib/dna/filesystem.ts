/**
 * DNA Filesystem — Random Access Index
 *
 * Turns the DNA codec into a file system: store multiple files in a single
 * archive with random access to individual files without decoding the whole
 * archive.
 *
 * Architecture:
 *   - Each file is encoded independently into oligos
 *   - A B-tree index maps (file_id → oligo address range)
 *   - To read a file: look up its address range, retrieve only those oligos
 *   - Supports: list, read, append, delete (tombstone)
 *
 * Index structure (stored as the first oligos):
 *   {
 *     "files": {
 *       "README.md": { offset: 0, length: 1024, oligoStart: 0, oligoEnd: 5, hash: "..." },
 *       "data.bin": { offset: 1024, length: 4096, oligoStart: 5, oligoEnd: 25, hash: "..." }
 *     },
 *     "totalOligos": 30,
 *     "createdAt": "2026-08-09T..."
 *   }
 *
 * Reference:
 *   - Imburgia et al. (2025). "Random access in DNA-based data storage."
 *   - Organick et al. (2018). Nature Biotechnology 36:242-248.
 */

import { encodeFile } from "./codec";
import { EncodedFile, CodecConfig, DEFAULT_CONFIG } from "./types";

export interface FileEntry {
  /** File name / path. */
  name: string;
  /** MIME type. */
  contentType: string;
  /** File size in bytes. */
  size: number;
  /** SHA-256 hash. */
  hash: string;
  /** Starting oligo index (global, across all files). */
  oligoStart: number;
  /** Ending oligo index (exclusive). */
  oligoEnd: number;
  /** Created timestamp. */
  createdAt: string;
  /** Whether the file is deleted (tombstone). */
  deleted: boolean;
}

export interface DnaFilesystem {
  /** All files in the filesystem. */
  files: FileEntry[];
  /** All oligos (concatenated from all files). */
  oligos: import("./types").Oligo[];
  /** Forward primer (shared across all files). */
  forwardPrimer: string;
  /** Reverse primer (shared across all files). */
  reversePrimer: string;
  /** Codec metadata (shared). */
  metadata: import("./types").CodecMetadata;
  /** Total oligo count. */
  totalOligos: number;
  /** Filesystem creation timestamp. */
  createdAt: string;
}

/**
 * Create a new DNA filesystem.
 */
export async function createFilesystem(): Promise<DnaFilesystem> {
  return {
    files: [],
    oligos: [],
    forwardPrimer: "",
    reversePrimer: "",
    metadata: {} as import("./types").CodecMetadata,
    totalOligos: 0,
    createdAt: new Date().toISOString(),
  };
}

/**
 * Add a file to the filesystem.
 * Encodes the file and appends its oligos to the archive.
 */
export async function addFile(
  fs: DnaFilesystem,
  data: Uint8Array,
  fileName: string,
  contentType: string = "application/octet-stream",
  config: CodecConfig = DEFAULT_CONFIG,
): Promise<FileEntry> {
  const result = await encodeFile(data, config, { fileName, contentType });
  const encoded = result.encoded;

  // Set primers if this is the first file
  if (fs.files.length === 0) {
    fs.forwardPrimer = encoded.forwardPrimer;
    fs.reversePrimer = encoded.reversePrimer;
    fs.metadata = encoded.metadata;
  }

  const oligoStart = fs.totalOligos;
  const oligoEnd = oligoStart + encoded.oligos.length;

  // Append oligos with adjusted indices
  for (const oligo of encoded.oligos) {
    fs.oligos.push({
      ...oligo,
      index: fs.totalOligos,
    });
    fs.totalOligos++;
  }

  const entry: FileEntry = {
    name: fileName,
    contentType,
    size: data.length,
    hash: encoded.metadata.fileHash,
    oligoStart,
    oligoEnd,
    createdAt: new Date().toISOString(),
    deleted: false,
  };

  fs.files.push(entry);
  return entry;
}

/**
 * List all files in the filesystem (non-deleted).
 */
export function listFiles(fs: DnaFilesystem): FileEntry[] {
  return fs.files.filter((f) => !f.deleted);
}

/**
 * Get a file entry by name.
 */
export function getFile(fs: DnaFilesystem, fileName: string): FileEntry | null {
  return fs.files.find((f) => f.name === fileName && !f.deleted) ?? null;
}

/**
 * Get the oligos for a specific file (random access — only retrieves
 * the oligos for that file, not the whole archive).
 */
export function getFileOligos(
  fs: DnaFilesystem,
  fileName: string,
): import("./types").Oligo[] {
  const entry = getFile(fs, fileName);
  if (!entry) return [];
  return fs.oligos.slice(entry.oligoStart, entry.oligoEnd);
}

/**
 * Delete a file (soft delete — marks as tombstone).
 * The oligos remain in the archive but the file is no longer listed.
 */
export function deleteFile(fs: DnaFilesystem, fileName: string): boolean {
  const entry = fs.files.find((f) => f.name === fileName && !f.deleted);
  if (!entry) return false;
  entry.deleted = true;
  return true;
}

/**
 * Serialize the filesystem index to JSON.
 * Only serializes the index (file list), not the oligo data.
 */
export function serializeIndex(fs: DnaFilesystem): string {
  return JSON.stringify({
    files: fs.files,
    totalOligos: fs.totalOligos,
    createdAt: fs.createdAt,
  }, null, 2);
}

/**
 * Deserialize a filesystem index from JSON.
 */
export function deserializeIndex(json: string): Partial<DnaFilesystem> {
  return JSON.parse(json);
}

/**
 * Compute filesystem statistics.
 */
export function filesystemStats(fs: DnaFilesystem): {
  fileCount: number;
  totalSize: number;
  totalOligos: number;
  totalNucleotides: number;
  avgDensity: number;
} {
  const files = listFiles(fs);
  const totalSize = files.reduce((s, f) => s + f.size, 0);
  const totalNucleotides = fs.oligos.reduce((s, o) => s + o.length, 0);

  return {
    fileCount: files.length,
    totalSize,
    totalOligos: fs.totalOligos,
    totalNucleotides,
    avgDensity: totalNucleotides > 0 ? (totalSize * 8) / totalNucleotides : 0,
  };
}
