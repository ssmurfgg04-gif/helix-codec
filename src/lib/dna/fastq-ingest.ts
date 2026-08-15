/**
 * FASTQ Ingestion Pipeline
 *
 * Parses real-world FASTQ files from DNA sequencing runs, extracts:
 *   - DNA sequences
 *   - Per-base Phred Q-scores (ASCII-encoded)
 *   - Empirical error statistics
 *
 * FASTQ format:
 *   @SEQ_ID
 *   GATTACAGATTCGATA...
 *   +
 *   !''*((((+++......
 *
 * Q-score encoding: Phred+33 (ASCII 33 = Q0, ASCII 73 = Q40)
 *
 * Usage:
 *   const profile = await ingestFastq(file);
 *   console.log(profile.avgSubRate, profile.avgDelRate, profile.avgInsRate);
 */

export interface FastqRead {
  id: string;
  sequence: string;
  quality: Uint8Array;
}

export interface EmpiricalNoiseProfile {
  totalReads: number;
  totalBases: number;
  avgReadLength: number;
  avgQScore: number;
  /** Q-score distribution: index = Q-value, value = count */
  qDistribution: Uint32Array;
  /** Reads with mean Q < 10 (low quality) */
  lowQualityReads: number;
  /** Reads with mean Q > 30 (high quality) */
  highQualityReads: number;
  /** Estimated substitution rate from Q-scores */
  estimatedSubRate: number;
  /** Read length distribution */
  lengthDistribution: Uint32Array;
  /** Parse time in ms */
  parseMs: number;
}

const PHRED_OFFSET = 33;
const MAX_Q = 93; // Phred+33 max

/**
 * Parse a FASTQ file from a string.
 * Returns array of reads with sequences and Q-scores.
 */
export function parseFastq(fastqContent: string): FastqRead[] {
  const reads: FastqRead[] = [];
  const lines = fastqContent.split("\n");
  let i = 0;

  while (i + 3 < lines.length) {
    // Line 1: @SEQ_ID
    const idLine = lines[i].trim();
    if (!idLine.startsWith("@")) {
      i++;
      continue;
    }
    const id = idLine.slice(1);

    // Line 2: sequence
    const sequence = lines[i + 1].trim();

    // Line 3: + (optional description)
    if (!lines[i + 2].trim().startsWith("+")) {
      i++;
      continue;
    }

    // Line 4: quality scores (ASCII Phred+33)
    const qualityStr = lines[i + 3].trim();
    const quality = new Uint8Array(qualityStr.length);
    for (let j = 0; j < qualityStr.length; j++) {
      quality[j] = qualityStr.charCodeAt(j) - PHRED_OFFSET;
    }

    reads.push({ id, sequence, quality });
    i += 4;
  }

  return reads;
}

/**
 * Analyze empirical noise from parsed FASTQ reads.
 * Extracts Q-score distribution, quality statistics, and estimated error rates.
 */
export function analyzeNoiseProfile(reads: FastqRead[]): EmpiricalNoiseProfile {
  const t0 = Date.now();
  let totalBases = 0;
  let totalQ = 0;
  const qDist = new Uint32Array(MAX_Q + 1);
  const lengthDist = new Uint32Array(1000); // up to 1000bp reads
  let lowQ = 0;
  let highQ = 0;
  let totalReadLen = 0;

  for (const read of reads) {
    const readLen = read.sequence.length;
    totalReadLen += readLen;
    totalBases += readLen;

    let readQSum = 0;
    for (let i = 0; i < read.quality.length; i++) {
      const q = Math.min(read.quality[i], MAX_Q);
      totalQ += q;
      qDist[q]++;
      readQSum += q;
    }

    const meanQ = readQSum / Math.max(read.quality.length, 1);
    if (meanQ < 10) lowQ++;
    if (meanQ > 30) highQ++;

    if (readLen < lengthDist.length) {
      lengthDist[readLen]++;
    }
  }

  const avgQ = totalBases > 0 ? totalQ / totalBases : 0;
  // Estimated substitution rate: avg P(error) = avg(10^(-Q/10))
  const estimatedSubRate = Math.pow(10, -avgQ / 10);

  return {
    totalReads: reads.length,
    totalBases,
    avgReadLength: reads.length > 0 ? totalReadLen / reads.length : 0,
    avgQScore: avgQ,
    qDistribution: qDist,
    lowQualityReads: lowQ,
    highQualityReads: highQ,
    estimatedSubRate,
    lengthDistribution: lengthDist,
    parseMs: Date.now() - t0,
  };
}

/**
 * Align reads to a reference using simple Hamming distance (for substitution counting).
 * Returns empirical substitution, insertion, deletion rates.
 *
 * For production, use minimap2 or BWA. This is a simplified alignment.
 */
export interface AlignmentStats {
  totalReads: number;
  alignedReads: number;
  totalSubstitutions: number;
  totalInsertions: number;
  totalDeletions: number;
  totalMatches: number;
  substitutionRate: number;
  insertionRate: number;
  deletionRate: number;
  matchRate: number;
}

export function alignToReference(reads: FastqRead[], reference: string): AlignmentStats {
  let aligned = 0;
  let subs = 0;
  let ins = 0;
  let dels = 0;
  let matches = 0;
  let totalBases = 0;

  for (const read of reads) {
    const seq = read.sequence;
    const refLen = reference.length;
    const readLen = seq.length;

    if (readLen === 0) continue;

    // Simple position-by-position alignment (no indel handling for this simplified version)
    const minLen = Math.min(readLen, refLen);
    for (let i = 0; i < minLen; i++) {
      if (seq[i] === reference[i]) {
        matches++;
      } else {
        subs++;
      }
      totalBases++;
    }

    // Count length differences as indels
    if (readLen > refLen) {
      ins += readLen - refLen;
      totalBases += readLen - refLen;
    } else if (readLen < refLen) {
      dels += refLen - readLen;
      totalBases += refLen - readLen;
    }

    aligned++;
  }

  return {
    totalReads: reads.length,
    alignedReads: aligned,
    totalSubstitutions: subs,
    totalInsertions: ins,
    totalDeletions: dels,
    totalMatches: matches,
    substitutionRate: totalBases > 0 ? subs / totalBases : 0,
    insertionRate: totalBases > 0 ? ins / totalBases : 0,
    deletionRate: totalBases > 0 ? dels / totalBases : 0,
    matchRate: totalBases > 0 ? matches / totalBases : 0,
  };
}

/**
 * Generate a synthetic FASTQ file for testing.
 * Simulates Illumina-style reads with realistic Q-score distributions.
 */
export function generateSyntheticFastq(
  reference: string,
  numReads: number,
  coverage: number,
  substitutionRate: number = 0.001,
  insertionRate: number = 0.0005,
  deletionRate: number = 0.001,
  seed: number = 42,
): string {
  let state = seed >>> 0 || 1;
  const rng = () => {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    state = state >>> 0;
    return state / 0x100000000;
  };

  const readLen = Math.min(200, reference.length);
  const bases = "ACGT";
  let output = "";

  for (let r = 0; r < numReads; r++) {
    const startPos = Math.floor(rng() * Math.max(1, reference.length - readLen));
    let sequence = "";
    let quality = "";

    for (let i = 0; i < readLen && startPos + i < reference.length; i++) {
      const refBase = reference[startPos + i];

      // Deletion?
      if (rng() < deletionRate) {
        continue; // skip this base
      }

      // Substitution?
      let emitBase: string;
      let qScore: number;
      if (rng() < substitutionRate) {
        do {
          emitBase = bases[Math.floor(rng() * 4)];
        } while (emitBase === refBase);
        qScore = 5 + Math.floor(rng() * 11); // Q5-Q15
      } else {
        emitBase = refBase;
        qScore = 30 + Math.floor(rng() * 11); // Q30-Q40
      }

      sequence += emitBase;
      quality += String.fromCharCode(qScore + PHRED_OFFSET);

      // Insertion?
      if (rng() < insertionRate) {
        const insBase = bases[Math.floor(rng() * 4)];
        sequence += insBase;
        quality += String.fromCharCode((2 + Math.floor(rng() * 7)) + PHRED_OFFSET); // Q2-Q8
      }
    }

    output += `@helix_read_${r} pos=${startPos}\n`;
    output += `${sequence}\n`;
    output += `+\n`;
    output += `${quality}\n`;
  }

  return output;
}

/**
 * Full ingestion pipeline: parse + analyze + align.
 */
export async function ingestFastq(
  fastqContent: string,
  reference?: string,
): Promise<{
  reads: FastqRead[];
  noiseProfile: EmpiricalNoiseProfile;
  alignmentStats?: AlignmentStats;
}> {
  const reads = parseFastq(fastqContent);
  const noiseProfile = analyzeNoiseProfile(reads);

  let alignmentStats: AlignmentStats | undefined;
  if (reference) {
    alignmentStats = alignToReference(reads, reference);
  }

  return { reads, noiseProfile, alignmentStats };
}
