/**
 * Read Until / Adaptive Sampling Simulator
 *
 * Simulates Oxford Nanopore's "Read Until" API where the sequencer ejects
 * a DNA strand if the first ~400 bases don't match the target file header.
 *
 * This enables PCR-free, instant random access to specific files in a DNA
 * pool without sequencing the entire archive.
 *
 * Pipeline:
 *   1. Simulate a Nanopore run: strands enter the pore in random order
 *   2. For each strand, read the first 400 bases (simulating real-time)
 *   3. Align the partial read to known file headers using DTW
 *   4. If match: keep sequencing (this strand contains target data)
 *   5. If no match: "eject" the strand (reverse voltage)
 *   6. Measure: time saved, reagent cost saved, coverage achieved
 *
 * Reference:
 *   - Loose, Malla, Stout (2016). "Real-time selective sequencing with
 *     ReadUntil." Nature Methods 13:751-754.
 *   - Oxford Nanopore Read Until API documentation
 */

import { dtwDistance } from "./squiggle";

export interface FileHeader {
  fileName: string;
  /** Header DNA sequence (first 400+ bases of the file's oligos). */
  headerSequence: string;
  /** Expected squiggle pattern for the header. */
  headerSquiggle: Float32Array;
}

export interface ReadUntilConfig {
  /** Number of bases to read before making eject decision. */
  decisionBases: number;
  /** DTW threshold for header matching. */
  matchThreshold: number;
  /** Target file to search for. */
  targetFile: string;
  /** Total strands in the pool. */
  totalStrands: number;
  /** Sequencing speed (bases per second). */
  sequencingSpeed: number;
  /** Cost per base in USD (for reagent cost calculation). */
  costPerBase: number;
}

export const DEFAULT_READUNTIL_CONFIG: ReadUntilConfig = {
  decisionBases: 400,
  matchThreshold: 50,
  targetFile: "",
  totalStrands: 100000,
  sequencingSpeed: 450, // 450 bases/sec per pore
  costPerBase: 0.00001, // ~$10 per 1M bases (Nanopore pricing)
};

export interface StrandResult {
  strandId: number;
  headerMatched: boolean;
  matchedFile: string | null;
  dtwDistance: number;
  decisionTimeMs: number;
  ejected: boolean;
  basesSequenced: number;
}

export interface ReadUntilResult {
  totalStrands: number;
  strandsKept: number;
  strandsEjected: number;
  totalBasesSequenced: number;
  basesSaved: number;
  timeSavedSeconds: number;
  costSavedUSD: number;
  coverageOfTarget: number;
  results: StrandResult[];
  summary: string;
}

/**
 * Simulate a Read Until run.
 *
 * Given a pool of strands (each belonging to a file), simulates the Nanopore
 * Read Until API: reads the first N bases of each strand, matches against
 * file headers, and ejects non-matching strands.
 *
 * @param strands Array of strand sequences (each from a specific file)
 * @param fileMap Map of strand index → file name
 * @param headers Array of file headers to match against
 * @param config Read Until configuration
 */
export function simulateReadUntil(
  strands: string[],
  fileMap: Map<number, string>,
  headers: FileHeader[],
  config: ReadUntilConfig = DEFAULT_READUNTIL_CONFIG,
): ReadUntilResult {
  const results: StrandResult[] = [];
  let strandsKept = 0;
  let strandsEjected = 0;
  let totalBasesSequenced = 0;
  let totalBasesThatWouldHaveBeenSequenced = 0;

  for (let i = 0; i < strands.length; i++) {
    const strand = strands[i];
    const fileName = fileMap.get(i) ?? "unknown";

    // Read the first `decisionBases` bases
    const partialRead = strand.slice(0, config.decisionBases);
    const decisionTimeMs = (config.decisionBases / config.sequencingSpeed) * 1000;

    // Match against all file headers using DTW
    let bestMatch: string | null = null;
    let bestDistance = Infinity;

    for (const header of headers) {
      // Generate squiggle for the partial read
      const partialSquiggle = generateSquiggleSimple(partialRead);
      const dist = dtwDistance(partialSquiggle, header.headerSquiggle.slice(0, Math.min(header.headerSquiggle.length, partialSquiggle.length)), 20);

      if (dist < bestDistance) {
        bestDistance = dist;
        bestMatch = dist < config.matchThreshold ? header.fileName : null;
      }
    }

    const isTarget = fileName === config.targetFile;
    const headerMatched = bestMatch === config.targetFile;
    const shouldEject = !headerMatched;

    if (shouldEject) {
      strandsEjected++;
      totalBasesSequenced += config.decisionBases; // only read decision bases
      totalBasesThatWouldHaveBeenSequenced += strand.length;
    } else {
      strandsKept++;
      totalBasesSequenced += strand.length;
      totalBasesThatWouldHaveBeenSequenced += strand.length;
    }

    results.push({
      strandId: i,
      headerMatched,
      matchedFile: bestMatch,
      dtwDistance: bestDistance,
      decisionTimeMs,
      ejected: shouldEject,
      basesSequenced: shouldEject ? config.decisionBases : strand.length,
    });
  }

  const basesSaved = totalBasesThatWouldHaveBeenSequenced - totalBasesSequenced;
  const timeSavedSeconds = basesSaved / config.sequencingSpeed;
  const costSavedUSD = basesSaved * config.costPerBase;

  // Calculate coverage of target file
  const targetStrands = Array.from(fileMap.entries()).filter(([, name]) => name === config.targetFile).length;
  const targetStrandsKept = results.filter(r => !r.ejected && fileMap.get(r.strandId) === config.targetFile).length;
  const coverageOfTarget = targetStrands > 0 ? targetStrandsKept / targetStrands : 0;

  const ejectRate = (strandsEjected / strands.length * 100).toFixed(1);
  const summary = `Read Until: ejected ${strandsEjected}/${strands.length} strands (${ejectRate}%), saved ${basesSaved.toLocaleString()} bases, ${timeSavedSeconds.toFixed(1)}s, $${costSavedUSD.toFixed(2)}`;

  return {
    totalStrands: strands.length,
    strandsKept,
    strandsEjected,
    totalBasesSequenced,
    basesSaved,
    timeSavedSeconds,
    costSavedUSD,
    coverageOfTarget,
    results: results.slice(0, 100), // return first 100 for inspection
    summary,
  };
}

/**
 * Simple squiggle generation (for Read Until simulation).
 */
function generateSquiggleSimple(sequence: string): Float32Array {
  const signals = new Float32Array(sequence.length);
  for (let i = 0; i < sequence.length; i++) {
    const base = sequence[i];
    const mean = base === "A" ? 55 : base === "C" ? 65 : base === "G" ? 75 : 45;
    signals[i] = mean + (Math.random() - 0.5) * 4;
  }
  return signals;
}

/**
 * Generate file headers for a DNA archive.
 * Each file's header is the first 400+ bases of its first oligo.
 */
export function generateFileHeaders(
  files: { name: string; oligos: { sequence: string }[] }[],
): FileHeader[] {
  return files.map(file => {
    const firstOligo = file.oligos[0];
    const headerSequence = firstOligo ? firstOligo.sequence.slice(0, 400) : "";
    return {
      fileName: file.name,
      headerSequence,
      headerSquiggle: generateSquiggleSimple(headerSequence),
    };
  });
}

/**
 * Generate a simulated DNA pool with multiple files.
 */
export function generateDnaPool(
  files: { name: string; numOligos: number; oligoLength: number }[],
): { strands: string[]; fileMap: Map<number, string> } {
  const strands: string[] = [];
  const fileMap = new Map<number, string>();
  const bases = "ACGT";

  for (const file of files) {
    for (let i = 0; i < file.numOligos; i++) {
      let seq = "";
      // Add file-specific header prefix (first 20 bases = file ID)
      const fileId = file.name.charCodeAt(0) + i;
      for (let b = 0; b < file.oligoLength; b++) {
        seq += bases[(fileId + b * 7) % 4];
      }
      fileMap.set(strands.length, file.name);
      strands.push(seq);
    }
  }

  // Shuffle strands (random order in the pool)
  for (let i = strands.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [strands[i], strands[j]] = [strands[j], strands[i]];
    // Swap file map entries too
    const tmp = fileMap.get(i);
    fileMap.set(i, fileMap.get(j));
    fileMap.set(j, tmp);
  }

  return { strands, fileMap };
}
