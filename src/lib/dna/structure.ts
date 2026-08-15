/**
 * Secondary Structure Detection
 *
 * DNA sequences can fold into secondary structures (hairpins, stem-loops,
 * G-quadruplexes) that interfere with synthesis and sequencing. This module
 * detects problematic structures so the encoder can re-encode oligos that
 * would form them.
 *
 * Hairpin detection: a hairpin forms when a sequence contains a palindromic
 * region (stem) separated by a short loop. The stem is a reverse-complement
 * palindrome; the loop is typically 3-10 nt.
 *
 * Example hairpin:
 *   5'-ACGT....ACGT-3'  (stem = ACGT, loop = ....)
 *       ||||
 *   3'-TGCA....TGCA-5'
 *
 * The stem's two arms are reverse complements of each other.
 *
 * Reference:
 *   - SantaLucia & Hicks (2004). "The thermodynamics of DNA structural motifs."
 *     Annu Rev Biophys 33:415-440.
 *   - Zhou et al. (2023). "DNAfold: predicting DNA secondary structures."
 */

const COMPLEMENT: Record<string, string> = { A: "T", T: "A", C: "G", G: "C" };

/**
 * Get the reverse complement of a DNA string.
 */
export function reverseComplement(seq: string): string {
  let result = "";
  for (let i = seq.length - 1; i >= 0; i--) {
    result += COMPLEMENT[seq[i]] ?? "N";
  }
  return result;
}

export interface Hairpin {
  /** Start position of the stem (5' arm). */
  stemStart: number;
  /** Length of the stem. */
  stemLen: number;
  /** Length of the loop between the two stem arms. */
  loopLen: number;
  /** Total length of the hairpin structure. */
  totalLen: number;
  /** Estimated melting temperature (°C) — higher = more stable. */
  meltingTemp: number;
}

/**
 * Detect hairpins in a DNA sequence.
 *
 * A hairpin is defined as: a region where seq[i..i+stemLen] is the reverse
 * complement of seq[i+stemLen+loopLen..i+2*stemLen+loopLen].
 *
 * @param seq DNA sequence
 * @param minStemLen Minimum stem length to consider (default 4)
 * @param minLoopLen Minimum loop length (default 3)
 * @param maxLoopLen Maximum loop length (default 10)
 * @returns Array of detected hairpins
 */
export function detectHairpins(
  seq: string,
  minStemLen: number = 4,
  minLoopLen: number = 3,
  maxLoopLen: number = 10,
): Hairpin[] {
  const hairpins: Hairpin[] = [];
  const len = seq.length;

  for (let i = 0; i <= len - 2 * minStemLen - minLoopLen; i++) {
    for (let stemLen = minStemLen; stemLen <= Math.min(12, (len - i - minLoopLen) / 2); stemLen++) {
      for (let loopLen = minLoopLen; loopLen <= maxLoopLen; loopLen++) {
        const arm1Start = i;
        const arm1End = i + stemLen;
        const arm2Start = arm1End + loopLen;
        const arm2End = arm2Start + stemLen;

        if (arm2End > len) break;

        const arm1 = seq.slice(arm1Start, arm1End);
        const arm2 = seq.slice(arm2Start, arm2End);
        const arm2RevComp = reverseComplement(arm2);

        if (arm1 === arm2RevComp) {
          // Found a hairpin!
          const gc = (arm1.match(/[GC]/g)?.length ?? 0);
          const meltingTemp = estimateMeltingTemp(arm1, loopLen);
          hairpins.push({
            stemStart: i,
            stemLen,
            loopLen,
            totalLen: 2 * stemLen + loopLen,
            meltingTemp,
          });
        }
      }
    }
  }

  // Sort by melting temp (most stable first) and deduplicate overlapping
  hairpins.sort((a, b) => b.meltingTemp - a.meltingTemp);
  const result: Hairpin[] = [];
  for (const hp of hairpins) {
    const overlaps = result.some(
      (r) =>
        hp.stemStart < r.stemStart + r.totalLen &&
        hp.stemStart + hp.totalLen > r.stemStart,
    );
    if (!overlaps) result.push(hp);
  }

  return result;
}

/**
 * Estimate the melting temperature of a hairpin stem.
 * Uses the simple Wallace rule: Tm = 2*(A+T) + 4*(G+C).
 * More accurate: SantaLucia 1998 nearest-neighbor, but Wallace is fast.
 */
function estimateMeltingTemp(stem: string, loopLen: number): number {
  const at = (stem.match(/[AT]/g)?.length ?? 0);
  const gc = (stem.match(/[GC]/g)?.length ?? 0);
  const wallace = 2 * at + 4 * gc;
  // Longer loops reduce stability slightly
  const loopPenalty = Math.max(0, loopLen - 3) * 2;
  return wallace - loopPenalty;
}

/**
 * Check if a sequence has any problematic secondary structures.
 * Returns true if the sequence is "clean" (no hairpins with Tm > threshold).
 */
export function isStructureFree(
  seq: string,
  maxMeltingTemp: number = 30,
): boolean {
  const hairpins = detectHairpins(seq);
  return hairpins.every((hp) => hp.meltingTemp < maxMeltingTemp);
}

/**
 * Compute a "structure score" for a sequence — higher = more problematic.
 * Useful for ranking candidate encodings.
 */
export function structureScore(seq: string): number {
  const hairpins = detectHairpins(seq);
  return hairpins.reduce((sum, hp) => sum + hp.meltingTemp, 0);
}

/**
 * Detect G-quadruplex motifs (runs of 3+ Gs separated by short loops).
 * These are known to cause synthesis problems.
 *
 * Pattern: G3+ N1-7 G3+ N1-7 G3+ N1-7 G3+
 */
export function detectGQuadruplexes(seq: string): { start: number; length: number }[] {
  const results: { start: number; length: number }[] = [];
  const re = /G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}/g;
  let match;
  while ((match = re.exec(seq)) !== null) {
    results.push({ start: match.index, length: match[0].length });
  }
  return results;
}

/**
 * Full structural constraint check for an oligo.
 * Returns a list of all structural issues found.
 */
export interface StructureIssue {
  type: "hairpin" | "g_quadruplex" | "homopolymer" | "gc_bias";
  start: number;
  length: number;
  severity: number; // 0..1, higher = worse
  description: string;
}

export function checkStructureConstraints(
  seq: string,
  options: {
    maxHairpinTm?: number;
    maxHomopolymer?: number;
    gcMin?: number;
    gcMax?: number;
    windowSize?: number;
  } = {},
): StructureIssue[] {
  const {
    maxHairpinTm = 30,
    maxHomopolymer = 3,
    gcMin = 0.4,
    gcMax = 0.6,
    windowSize = 20,
  } = options;

  const issues: StructureIssue[] = [];

  // Hairpins
  const hairpins = detectHairpins(seq);
  for (const hp of hairpins) {
    if (hp.meltingTemp >= maxHairpinTm) {
      issues.push({
        type: "hairpin",
        start: hp.stemStart,
        length: hp.totalLen,
        severity: Math.min(1, hp.meltingTemp / 60),
        description: `Hairpin stem=${hp.stemLen}nt loop=${hp.loopLen}nt Tm=${hp.meltingTemp}°C`,
      });
    }
  }

  // G-quadruplexes
  const gquad = detectGQuadruplexes(seq);
  for (const gq of gquad) {
    issues.push({
      type: "g_quadruplex",
      start: gq.start,
      length: gq.length,
      severity: 0.8,
      description: `G-quadruplex motif (${gq.length}nt)`,
    });
  }

  // Homopolymers
  let runLen = 1;
  for (let i = 1; i < seq.length; i++) {
    if (seq[i] === seq[i - 1]) {
      runLen++;
    } else {
      if (runLen > maxHomopolymer) {
        issues.push({
          type: "homopolymer",
          start: i - runLen,
          length: runLen,
          severity: Math.min(1, (runLen - maxHomopolymer) / 4),
          description: `${seq[i - 1]} run of ${runLen}nt`,
        });
      }
      runLen = 1;
    }
  }
  if (runLen > maxHomopolymer) {
    issues.push({
      type: "homopolymer",
      start: seq.length - runLen,
      length: runLen,
      severity: Math.min(1, (runLen - maxHomopolymer) / 4),
      description: `${seq[seq.length - 1]} run of ${runLen}nt`,
    });
  }

  // GC bias in sliding windows
  for (let i = 0; i <= seq.length - windowSize; i += windowSize / 2) {
    const window = seq.slice(i, i + windowSize);
    const gc = (window.match(/[GC]/g)?.length ?? 0) / window.length;
    if (gc < gcMin || gc > gcMax) {
      issues.push({
        type: "gc_bias",
        start: i,
        length: windowSize,
        severity: Math.min(1, Math.abs(gc - 0.5) * 2),
        description: `GC ${(gc * 100).toFixed(0)}% in window [${i}, ${i + windowSize})`,
      });
    }
  }

  return issues;
}
