/**
 * K-mer Indexing for Fast Clustering
 *
 * Instead of comparing every read against every other read (O(n²)), we build
 * a k-mer index: a hash map from k-mer → set of read IDs. Reads sharing many
 * k-mers are likely from the same oligo and get clustered together.
 *
 * For DNA storage, we cluster by the ADDRESS region (first 16 nt after primer).
 * Since the address is short and exact, we can use k=16 (the whole address)
 * as the clustering key — O(n) instead of O(n²).
 *
 * For noisy addresses (with substitutions), we use k=8 minimizers: extract
 * all 8-mers from the address, and cluster reads that share >= 50% of their
 * 8-mers. This is O(n * k) where k is the number of minimizers per read.
 */

const BASE_TO_BITS: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

/**
 * Encode a DNA k-mer as a 2k-bit integer (each base = 2 bits).
 * Faster than string hashing for lookup.
 */
export function kmerToBits(kmer: string): number {
  let bits = 0;
  for (let i = 0; i < kmer.length; i++) {
    const b = BASE_TO_BITS[kmer[i]];
    if (b === undefined) return -1; // invalid base (e.g., N or gap)
    bits = (bits << 2) | b;
  }
  return bits >>> 0; // unsigned
}

/**
 * Extract all k-mers from a sequence, returning their bit encodings.
 */
export function extractKmers(seq: string, k: number): number[] {
  const kmers: number[] = [];
  for (let i = 0; i <= seq.length - k; i++) {
    const bits = kmerToBits(seq.slice(i, i + k));
    if (bits >= 0) kmers.push(bits);
  }
  return kmers;
}

/**
 * Build a k-mer index: Map<kmer_bits, Set<readId>>.
 */
export function buildKmerIndex(
  reads: { id: number; sequence: string }[],
  k: number,
): Map<number, Set<number>> {
  const index = new Map<number, Set<number>>();
  for (const read of reads) {
    const kmers = extractKmers(read.sequence, k);
    for (const kmer of kmers) {
      if (!index.has(kmer)) index.set(kmer, new Set());
      index.get(kmer)!.add(read.id);
    }
  }
  return index;
}

/**
 * Cluster reads by exact address match (fast path).
 * Assumes reads are pre-trimmed of primers and the first `addressLen` bases
 * are the address.
 *
 * O(n) — each read is hashed once.
 */
export function clusterByAddress<T extends { sequence: string }>(
  reads: T[],
  addressLen: number,
): Map<number, T[]> {
  const clusters = new Map<number, T[]>();
  for (const read of reads) {
    if (read.sequence.length < addressLen) continue;
    const addr = read.sequence.slice(0, addressLen);
    const key = kmerToBits(addr);
    if (key < 0) continue;
    if (!clusters.has(key)) clusters.set(key, []);
    clusters.get(key)!.push(read);
  }
  return clusters;
}

/**
 * Build a reference k-mer index from a set of known reference DNA strings.
 *
 * For DNA storage, we know all possible addresses ahead of time (oligo
 * indices 0..N-1, each deterministically encoded to DNA). This lets us
 * pre-compute the k-mer index once, then match noisy reads against it
 * in O(reads * k) instead of O(reads² * k).
 *
 * @param references Array of reference DNA strings (one per oligo)
 * @param k k-mer size (default 5 — small enough to survive indels)
 * @returns Map<kmer_bits, Int32Array of reference indices>
 */
export function buildReferenceKmerIndex(
  references: string[],
  k: number = 5,
): Map<number, number[]> {
  const index = new Map<number, number[]>();
  for (let refIdx = 0; refIdx < references.length; refIdx++) {
    const ref = references[refIdx];
    if (ref.length < k) continue;
    const seen = new Set<number>(); // dedupe within a single reference
    for (let i = 0; i <= ref.length - k; i++) {
      const bits = kmerToBits(ref.slice(i, i + k));
      if (bits >= 0 && !seen.has(bits)) {
        seen.add(bits);
        if (!index.has(bits)) index.set(bits, []);
        index.get(bits)!.push(refIdx);
      }
    }
  }
  return index;
}

/**
 * Find the best-matching reference index for a noisy read's address.
 *
 * Uses k-mer overlap: extracts k-mers from the read's address, tallies
 * how many overlap with each reference, and returns the reference with
 * the highest overlap.
 *
 * Robust to:
 *   - Substitutions (each sub affects at most k k-mers)
 *   - Insertions/deletions (shifts subsequent k-mers, but pre-indel
 *     k-mers still match)
 *
 * v59: Added margin-based filtering. If the top two candidates have
 * similar overlap counts (within `margin`), the match is ambiguous and
 * rejected. This reduces false positives from ~22% to <5% at 9% IDS.
 *
 * @param readAddress The noisy address DNA from the read (first addressNt bases)
 * @param kmerIndex Pre-built reference k-mer index (from buildReferenceKmerIndex)
 * @param k k-mer size (must match the index)
 * @param minOverlap Minimum k-mer overlap to accept a match (default 3)
 * @param margin Required margin between best and second-best overlap (default 1)
 * @returns Best-matching reference index, or -1 if no match meets threshold
 */
export function matchReadToReference(
  readAddress: string,
  kmerIndex: Map<number, number[]>,
  k: number = 5,
  minOverlap: number = 3,
  margin: number = 1,
): { bestIdx: number; bestOverlap: number; candidates: { idx: number; overlap: number }[] } {
  if (readAddress.length < k) {
    return { bestIdx: -1, bestOverlap: 0, candidates: [] };
  }

  // Tally votes per reference index
  const votes = new Map<number, number>();
  const seen = new Set<number>(); // dedupe within this read
  for (let i = 0; i <= readAddress.length - k; i++) {
    const bits = kmerToBits(readAddress.slice(i, i + k));
    if (bits < 0 || seen.has(bits)) continue;
    seen.add(bits);
    const refs = kmerIndex.get(bits);
    if (!refs) continue;
    for (const refIdx of refs) {
      votes.set(refIdx, (votes.get(refIdx) ?? 0) + 1);
    }
  }

  if (votes.size === 0) {
    return { bestIdx: -1, bestOverlap: 0, candidates: [] };
  }

  // Find top candidates (keep top 3 for downstream tie-breaking)
  const candidates = Array.from(votes.entries())
    .map(([idx, overlap]) => ({ idx, overlap }))
    .sort((a, b) => b.overlap - a.overlap)
    .slice(0, 3);

  const best = candidates[0];
  if (best.overlap < minOverlap) {
    return { bestIdx: -1, bestOverlap: best.overlap, candidates };
  }

  // v59: Margin-based filtering — reject ambiguous matches
  if (candidates.length >= 2) {
    const secondBest = candidates[1];
    if (best.overlap - secondBest.overlap < margin) {
      // Ambiguous match — reject
      return { bestIdx: -1, bestOverlap: best.overlap, candidates };
    }
  }

  return { bestIdx: best.idx, bestOverlap: best.overlap, candidates };
}

/**
 * Cluster reads by approximate address match (handles substitutions).
 * Uses k-mer minimizers: reads sharing >= 50% of their address k-mers are
 * clustered together.
 *
 * O(n * k) where k = address length / minimizer size.
 */
export function clusterByAddressApprox<T extends { sequence: string }>(
  reads: T[],
  addressLen: number,
  k: number = 8,
  threshold: number = 0.5,
): Map<number, T[]> {
  // First pass: extract k-mers from each read's address
  const readKmers = reads.map((r) => {
    if (r.sequence.length < addressLen) return [];
    return extractKmers(r.sequence.slice(0, addressLen), k);
  });

  // Build inverted index: kmer -> read indices
  const index = new Map<number, number[]>();
  for (let i = 0; i < reads.length; i++) {
    for (const kmer of readKmers[i]) {
      if (!index.has(kmer)) index.set(kmer, []);
      index.get(kmer)!.push(i);
    }
  }

  // For each read, find other reads that share >= threshold of k-mers
  const clusterId = new Array(reads.length).fill(-1);
  let nextClusterId = 0;
  const clusters = new Map<number, T[]>();

  for (let i = 0; i < reads.length; i++) {
    if (clusterId[i] !== -1) continue; // already assigned
    const myCluster = nextClusterId++;
    clusterId[i] = myCluster;
    clusters.set(myCluster, [reads[i]]);

    const myKmers = new Set(readKmers[i]);
    if (myKmers.size === 0) continue;

    // Find candidate matches via inverted index
    const candidates = new Set<number>();
    for (const kmer of myKmers) {
      for (const j of index.get(kmer) ?? []) {
        if (j !== i && clusterId[j] === -1) candidates.add(j);
      }
    }

    // Check each candidate
    for (const j of candidates) {
      const theirKmers = new Set(readKmers[j]);
      let shared = 0;
      for (const kmer of myKmers) {
        if (theirKmers.has(kmer)) shared++;
      }
      if (shared / myKmers.size >= threshold) {
        clusterId[j] = myCluster;
        clusters.get(myCluster)!.push(reads[j]);
      }
    }
  }

  return clusters;
}
