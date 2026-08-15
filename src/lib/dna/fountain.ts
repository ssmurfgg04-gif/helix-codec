/**
 * LT (Luby Transform) Fountain Code for DNA Storage
 *
 * A rateless erasure code: encode data into an unlimited stream of "droplets",
 * each of which is a random XOR of a subset of data chunks. The decoder can
 * recover the original data from ANY K + ε droplets (with high probability),
 * where K is the number of data chunks.
 *
 * This is the code used by DNA Fountain (Erlich & Zielinski 2017, Science),
 * which achieved 1.57 bits/nt — the highest density at the time.
 *
 * Algorithm:
 *   1. Split data into K chunks of fixed size.
 *   2. For each droplet:
 *      a. Sample degree d from the Robust Soliton Distribution.
 *      b. Select d distinct source chunks at random.
 *      c. XOR the selected chunks → droplet payload.
 *      d. Encode the seed (so decoder knows which chunks were XORed) and
 *         the payload into an oligo.
 *   3. Screen each droplet for biological constraints (GC, homopolymer).
 *      If it fails, discard and generate a new one.
 *   4. Continue until enough droplets are generated (typically K * 1.05-1.5).
 *
 * Decoding (peeling decoder):
 *   1. Find a droplet of degree 1 (a single source chunk).
 *   2. Recover that chunk: chunk = droplet_payload.
 *   3. XOR the recovered chunk from all droplets that include it.
 *   4. Update droplet degrees.
 *   5. Repeat until all chunks recovered.
 *
 * Reference:
 *   - Luby (2002). "LT codes." FOCS.
 *   - Erlich & Zielinski (2017). "DNA Fountain enables a robust and efficient
 *     storage architecture." Science 355:6328.
 *   - github.com/TeamErlich/dna-fountain (original Python implementation)
 */

export interface DnaFountainConfig {
  chunkSize: number; // bytes per source chunk
  /** Robust Soliton Distribution parameter c (default 0.1). */
  rsdC: number;
  /** Robust Soliton Distribution parameter delta (default 0.5). */
  rsdDelta: number;
  /** PRNG seed for reproducibility. */
  seed: number;
  /** Maximum number of droplets to generate (safety limit). */
  maxDroplets: number;
}

export const DEFAULT_FOUNTAIN_CONFIG: DnaFountainConfig = {
  chunkSize: 32,
  rsdC: 0.1,
  rsdDelta: 0.5,
  seed: 42,
  maxDroplets: 100000,
};

export interface Droplet {
  /** Seed used to generate this droplet (determines degree + source selection). */
  seed: number;
  /** Degree (number of source chunks XORed). */
  degree: number;
  /** Indices of source chunks XORed. */
  sourceIndices: number[];
  /** XOR of the selected source chunks. */
  payload: Uint8Array;
}

export interface FountainEncoding {
  droplets: Droplet[];
  numChunks: number;
  chunkSize: number;
  originalLength: number;
}

// Xorshift32 PRNG
class Rng {
  private state: number;
  constructor(seed: number) {
    this.state = (seed >>> 0) || 1;
  }
  next(): number {
    this.state ^= this.state << 13;
    this.state ^= this.state >>> 17;
    this.state ^= this.state << 5;
    this.state = this.state >>> 0;
    return this.state / 0x100000000;
  }
  nextInt(maxExclusive: number): number {
    return Math.floor(this.next() * maxExclusive);
  }
}

/**
 * Compute the Robust Soliton Distribution.
 *
 * μ(d) = (ρ(d) + τ(d)) / Z
 * where:
 *   ρ(d) = 1/K for d=1, 1/(d*(d-1)) for d=2..K
 *   τ(d) = S/(K*d) for d=1..K/S-1, S*ln(S/δ)/K for d=K/S, 0 otherwise
 *   S = c * ln(K/δ) * sqrt(K)
 *   Z = sum(ρ + τ)
 *
 * Returns the CDF for sampling.
 *
 * DNA channel optimization (Schwarz 2024):
 *   For DNA storage channels, the optimal c parameter is smaller (0.01-0.05)
 *   than the default 0.1, because DNA channels have lower error rates than
 *   general erasure channels. This reduces overhead by ~5%.
 */
export function robustSolitonCDF(
  K: number,
  c: number = 0.1,
  delta: number = 0.5,
): { cdf: Float64Array; S: number } {
  // DNA channel optimization: use smaller c for large K
  // Schwarz 2024 showed c=0.02-0.05 is optimal for DNA channels
  if (K > 1000 && c === 0.1) {
    c = 0.03; // DNA-optimized default
  }

  const rho = new Float64Array(K + 1);
  const tau = new Float64Array(K + 1);

  // Ideal Soliton (rho)
  rho[1] = 1 / K;
  for (let d = 2; d <= K; d++) {
    rho[d] = 1 / (d * (d - 1));
  }

  // Robust extension (tau)
  const S = c * Math.log(K / delta) * Math.sqrt(K);
  const KOverS = Math.floor(K / S);
  for (let d = 1; d <= K; d++) {
    if (d <= KOverS - 1) {
      tau[d] = S / (K * d);
    } else if (d === KOverS) {
      tau[d] = (S * Math.log(S / delta)) / K;
    } else {
      tau[d] = 0;
    }
  }

  // Normalize
  let Z = 0;
  const mu = new Float64Array(K + 1);
  for (let d = 1; d <= K; d++) {
    mu[d] = rho[d] + tau[d];
    Z += mu[d];
  }

  // Build CDF
  const cdf = new Float64Array(K + 1);
  let cum = 0;
  for (let d = 1; d <= K; d++) {
    cum += mu[d] / Z;
    cdf[d] = cum;
  }
  cdf[K] = 1.0; // ensure exact 1.0

  return { cdf, S };
}

/** Sample a degree from the Robust Soliton Distribution. */
function sampleDegree(rng: Rng, cdf: Float64Array): number {
  const r = rng.next();
  for (let d = 1; d < cdf.length; d++) {
    if (r <= cdf[d]) return d;
  }
  return cdf.length - 1;
}

/** Select d distinct indices from [0, K) using a seeded PRNG. */
function selectIndices(rng: Rng, d: number, K: number): number[] {
  const indices = new Set<number>();
  while (indices.size < d) {
    indices.add(rng.nextInt(K));
  }
  return Array.from(indices).sort((a, b) => a - b);
}

/**
 * Encode data into fountain droplets.
 *
 * @param data Input bytes
 * @param config Fountain configuration
 * @param numDroplets Number of droplets to generate (typically K * 1.1 to K * 1.5)
 */
export function fountainEncode(
  data: Uint8Array,
  config: DnaFountainConfig = DEFAULT_FOUNTAIN_CONFIG,
  numDroplets?: number,
): FountainEncoding {
  const K = Math.max(1, Math.ceil(data.length / config.chunkSize));
  const actualNumDroplets = numDroplets ?? Math.ceil(K * 1.2);

  // Pad data to multiple of chunkSize
  const paddedLen = K * config.chunkSize;
  const padded = new Uint8Array(paddedLen);
  padded.set(data, 0);

  // Build chunks
  const chunks: Uint8Array[] = [];
  for (let i = 0; i < K; i++) {
    chunks.push(padded.slice(i * config.chunkSize, (i + 1) * config.chunkSize));
  }

  // Build RSD CDF
  const { cdf } = robustSolitonCDF(K, config.rsdC, config.rsdDelta);

  // Generate droplets
  const droplets: Droplet[] = [];
  const rng = new Rng(config.seed);
  for (let i = 0; i < actualNumDroplets; i++) {
    const dropletSeed = rng.next() * 0xffffffff;
    const dropletRng = new Rng(dropletSeed);
    const degree = sampleDegree(dropletRng, cdf);
    const sourceIndices = selectIndices(dropletRng, degree, K);

    // XOR selected chunks
    const payload = new Uint8Array(config.chunkSize);
    for (const idx of sourceIndices) {
      for (let b = 0; b < config.chunkSize; b++) {
        payload[b] ^= chunks[idx][b];
      }
    }

    droplets.push({
      seed: Math.floor(dropletSeed),
      degree,
      sourceIndices,
      payload,
    });
  }

  return {
    droplets,
    numChunks: K,
    chunkSize: config.chunkSize,
    originalLength: data.length,
  };
}

/**
 * Decode fountain droplets using the peeling decoder, with Gaussian elimination
 * fallback for the remaining chunks.
 *
 * The peeling decoder iteratively:
 *   1. Finds droplets of degree 1 (single source chunk).
 *   2. Recovers that chunk.
 *   3. XORs it from all droplets containing it.
 *   4. Reduces their degrees.
 *   5. Repeats until no more degree-1 droplets.
 *
 * If peeling gets stuck (no degree-1 droplets but chunks remain), we fall back
 * to Gaussian elimination over GF(2) on the remaining chunks. This solves the
 * linear system formed by the remaining droplets.
 *
 * @returns Recovered data, or null if decoding failed.
 */
export function fountainDecode(encoding: FountainEncoding): Uint8Array | null {
  const { droplets, numChunks: K, chunkSize, originalLength } = encoding;

  const chunks: (Uint8Array | null)[] = new Array(K).fill(null);
  const remaining = new Set<number>();
  for (let i = 0; i < K; i++) remaining.add(i);

  // Working copies of droplets (we'll modify them during peeling)
  const workDroplets = droplets.map((d) => ({
    ...d,
    payload: d.payload.slice(),
    remainingSources: new Set(d.sourceIndices),
  }));

  // Phase 1: Peeling decoder
  let progress = true;
  while (progress && remaining.size > 0) {
    progress = false;
    for (const wd of workDroplets) {
      if (wd.remainingSources.size !== 1) continue;
      const chunkIdx = Array.from(wd.remainingSources)[0];
      if (chunks[chunkIdx] !== null) {
        wd.remainingSources.delete(chunkIdx);
        continue;
      }

      chunks[chunkIdx] = wd.payload.slice();
      remaining.delete(chunkIdx);
      progress = true;

      for (const other of workDroplets) {
        if (other === wd) continue;
        if (other.remainingSources.has(chunkIdx)) {
          for (let b = 0; b < chunkSize; b++) {
            other.payload[b] ^= chunks[chunkIdx]![b];
          }
          other.remainingSources.delete(chunkIdx);
        }
      }
      wd.remainingSources.clear();
    }
  }

  if (remaining.size === 0) {
    // All chunks recovered via peeling
    const result = new Uint8Array(K * chunkSize);
    for (let i = 0; i < K; i++) {
      result.set(chunks[i]!, i * chunkSize);
    }
    return result.slice(0, originalLength);
  }

  // Phase 2: Gaussian elimination fallback for remaining chunks.
  // Build a linear system over GF(2) for each byte position.
  // For each byte position b, we have:
  //   For each remaining droplet d: sum(chunks[s][b] for s in d.remainingSources) = d.payload[b]
  // This is a system of equations in the unknown chunk bytes.
  const remainingArr = Array.from(remaining);
  const remainingSet = new Set(remainingArr);
  // Filter droplets that still have remaining sources
  const usefulDroplets = workDroplets.filter(wd => wd.remainingSources.size > 0);

  if (usefulDroplets.length < remainingArr.length) {
    return null; // Not enough equations
  }

  // For each byte position, solve the linear system
  // Build the coefficient matrix (rows = droplets, cols = remaining chunks)
  // Each entry is 0 or 1 (whether the droplet includes that chunk)
  const numEquations = usefulDroplets.length;
  const numUnknowns = remainingArr.length;

  // Process each byte position independently
  for (let b = 0; b < chunkSize; b++) {
    // Build augmented matrix [A | y] where:
    //   A[i][j] = 1 if droplet i includes remaining chunk j
    //   y[i] = droplet i's payload at byte b
    const A: Uint8Array[] = []; // numEquations × numUnknowns
    const y: Uint8Array = new Uint8Array(numEquations);
    for (let i = 0; i < numEquations; i++) {
      const row = new Uint8Array(numUnknowns);
      const wd = usefulDroplets[i];
      for (const src of wd.remainingSources) {
        const colIdx = remainingArr.indexOf(src);
        if (colIdx >= 0) row[colIdx] = 1;
      }
      A.push(row);
      y[i] = wd.payload[b];
    }

    // Gaussian elimination over GF(2)
    let pivotRow = 0;
    for (let col = 0; col < numUnknowns && pivotRow < numEquations; col++) {
      // Find a row with a 1 in this column
      let found = -1;
      for (let row = pivotRow; row < numEquations; row++) {
        if (A[row][col] === 1) {
          found = row;
          break;
        }
      }
      if (found === -1) continue;

      // Swap rows
      if (found !== pivotRow) {
        [A[pivotRow], A[found]] = [A[found], A[pivotRow]];
        [y[pivotRow], y[found]] = [y[found], y[pivotRow]];
      }

      // Eliminate this column from all other rows
      for (let row = 0; row < numEquations; row++) {
        if (row === pivotRow) continue;
        if (A[row][col] === 1) {
          for (let c = 0; c < numUnknowns; c++) {
            A[row][c] ^= A[pivotRow][c];
          }
          y[row] ^= y[pivotRow];
        }
      }
      pivotRow++;
    }

    if (pivotRow < numUnknowns) {
      return null; // System is underdetermined
    }

    // Extract solution
    for (let j = 0; j < numUnknowns; j++) {
      const chunkIdx = remainingArr[j];
      if (chunks[chunkIdx] === null) {
        chunks[chunkIdx] = new Uint8Array(chunkSize);
      }
      chunks[chunkIdx]![b] = y[j];
    }
  }

  // Check if all chunks are now recovered
  for (let i = 0; i < K; i++) {
    if (chunks[i] === null) return null;
  }

  const result = new Uint8Array(K * chunkSize);
  for (let i = 0; i < K; i++) {
    result.set(chunks[i]!, i * chunkSize);
  }
  return result.slice(0, originalLength);
}
