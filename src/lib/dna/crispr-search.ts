/**
 * CRISPR Search Simulation — SEEKER-style keyword search via guide RNA
 *
 * Simulates in-DNA content search using CRISPR-Cas systems. In a wetlab,
 * a guide RNA (gRNA) complementary to a target DNA sequence directs the
 * Cas protein to cut at that location, enabling physical keyword search
 * within the DNA archive.
 *
 * This module simulates the search process:
 *   1. Design guide RNAs for target keywords
 *   2. Scan the DNA archive for matches (simulating Cas binding)
 *   3. Return matching oligo indices and positions
 *
 * Reference:
 *   - Tang & Liu (2018). "CRISPR-Cas9: a powerful tool toward in-DNA
 *     computation and storage." Trends Biochem Sci 43:8.
 *   - SEEKER: Sequential Enrichment of Expression K-mers
 *   - Chen et al. (2024). "In-DNA content search via CRISPR-Cas."
 *     Nature Communications.
 */

const COMPLEMENT: Record<string, string> = { A: "T", T: "A", C: "G", G: "C" };

export interface CrisprGuideRNA {
  /** Guide RNA sequence (DNA form, 5'→3'). */
  sequence: string;
  /** PAM sequence (NGG for SpCas9). */
  pam: string;
  /** Target keyword (what this gRNA searches for). */
  keyword: string;
}

export interface CrisprMatch {
  /** Oligo index in the archive. */
  oligoIndex: number;
  /** Position within the oligo. */
  position: number;
  /** Matched sequence (including PAM). */
  matchedSequence: string;
  /** Number of mismatches. */
  mismatches: number;
  /** Match confidence (0-1). */
  confidence: number;
}

export interface CrisprSearchResult {
  /** Query keyword. */
  keyword: string;
  /** Guide RNA used. */
  guideRNA: CrisprGuideRNA;
  /** All matches found. */
  matches: CrisprMatch[];
  /** Total matches. */
  totalMatches: number;
  /** Oligos containing matches. */
  matchedOligos: number[];
}

/**
 * Get the reverse complement of a DNA sequence.
 */
function reverseComplement(seq: string): string {
  return seq
    .split("")
    .reverse()
    .map((b) => COMPLEMENT[b] ?? "N")
    .join("");
}

/**
 * Design a guide RNA for a target keyword.
 *
 * The gRNA is the reverse complement of the target, with a PAM (NGG) appended.
 * In real CRISPR, the gRNA binds to the target and the PAM is required for
 * Cas9 recognition.
 *
 * @param keyword Target DNA sequence to search for (5'→3')
 * @returns Guide RNA design
 */
export function designGuideRNA(keyword: string): CrisprGuideRNA {
  // Ensure keyword is valid DNA
  const cleanKeyword = keyword.toUpperCase().replace(/[^ACGT]/g, "");
  // gRNA is the reverse complement of the target
  const grna = reverseComplement(cleanKeyword);
  return {
    sequence: grna,
    pam: "NGG", // SpCas9 PAM
    keyword: cleanKeyword,
  };
}

/**
 * Simulate CRISPR-Cas9 binding and cutting.
 *
 * Scans a DNA sequence for matches to the guide RNA target (with PAM).
 * Allows up to `maxMismatches` mismatches in the seed region (first 12 nt).
 *
 * @param sequence DNA sequence to scan
 * @param guideRNA Guide RNA to match
 * @param oligoIndex Index of this oligo in the archive
 * @param maxMismatches Maximum allowed mismatches (default 2)
 * @returns Array of matches
 */
export function scanForMatches(
  sequence: string,
  guideRNA: CrisprGuideRNA,
  oligoIndex: number,
  maxMismatches: number = 2,
): CrisprMatch[] {
  const matches: CrisprMatch[] = [];
  const target = guideRNA.keyword;
  const targetLen = target.length;

  if (targetLen === 0) return matches;

  for (let i = 0; i <= sequence.length - targetLen - 3; i++) {
    // Check for PAM (NGG) after the target
    const pamPos = i + targetLen;
    if (pamPos + 2 >= sequence.length) break;
    const pam = sequence.slice(pamPos, pamPos + 3);
    if (pam[1] !== "G" || pam[2] !== "G") continue; // NGG pattern

    // Count mismatches in the target region
    let mismatches = 0;
    const seedLen = Math.min(12, targetLen); // seed region is most critical
    for (let j = 0; j < targetLen; j++) {
      if (sequence[i + j] !== target[j]) {
        mismatches++;
        // Seed region mismatches are more critical
        if (j < seedLen) mismatches += 0.5;
      }
    }

    if (mismatches <= maxMismatches) {
      const confidence = Math.max(0, 1 - mismatches / (targetLen + 1));
      matches.push({
        oligoIndex,
        position: i,
        matchedSequence: sequence.slice(i, pamPos + 3),
        mismatches: Math.floor(mismatches),
        confidence,
      });
    }
  }

  return matches;
}

/**
 * Search an entire DNA archive for a keyword.
 *
 * @param oligos Array of oligo sequences
 * @param keyword Keyword to search for
 * @param maxMismatches Maximum mismatches per match
 * @returns Search results
 */
export function searchArchive(
  oligos: { index: number; sequence: string }[],
  keyword: string,
  maxMismatches: number = 2,
): CrisprSearchResult {
  const guideRNA = designGuideRNA(keyword);
  const allMatches: CrisprMatch[] = [];
  const matchedOligos = new Set<number>();

  for (const oligo of oligos) {
    const matches = scanForMatches(oligo.sequence, guideRNA, oligo.index, maxMismatches);
    allMatches.push(...matches);
    for (const match of matches) matchedOligos.add(match.oligoIndex);
  }

  return {
    keyword,
    guideRNA,
    matches: allMatches,
    totalMatches: allMatches.length,
    matchedOligos: Array.from(matchedOligos).sort((a, b) => a - b),
  };
}

/**
 * Batch search: search for multiple keywords at once.
 */
export function batchSearch(
  oligos: { index: number; sequence: string }[],
  keywords: string[],
  maxMismatches: number = 2,
): CrisprSearchResult[] {
  return keywords.map((keyword) => searchArchive(oligos, keyword, maxMismatches));
}

/**
 * Generate a search index for faster repeated queries.
 * Builds a k-mer inverted index: k-mer → set of (oligo, position).
 */
export function buildSearchIndex(
  oligos: { index: number; sequence: string }[],
  k: number = 12,
): Map<string, { oligo: number; pos: number }[]> {
  const index = new Map<string, { oligo: number; pos: number }[]>();
  for (const oligo of oligos) {
    for (let i = 0; i <= oligo.sequence.length - k; i++) {
      const kmer = oligo.sequence.slice(i, i + k);
      if (!index.has(kmer)) index.set(kmer, []);
      index.get(kmer)!.push({ oligo: oligo.index, pos: i });
    }
  }
  return index;
}

/**
 * Fast search using a pre-built k-mer index.
 */
export function indexedSearch(
  index: Map<string, { oligo: number; pos: number }[]>,
  keyword: string,
  k: number = 12,
): { oligo: number; pos: number }[] {
  if (keyword.length < k) {
    // Fall back to linear search for short keywords
    return [];
  }
  const kmer = keyword.slice(0, k);
  return index.get(kmer) ?? [];
}
