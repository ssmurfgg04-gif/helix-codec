/**
 * P1: Gungnir Hash-Based Single-Read Recovery
 * Based on HKU-BAL/Gungnir (Nature Communications 2026)
 * Proof-of-work: hash signature per fragment, guess until hash matches
 */

export interface GungnirOptions {
  maxGuesses?: number;
  hashBits?: number;
  maxSubstitutions?: number;
  maxIndels?: number;
  errorHints?: { likelyPositions?: number[]; homopolymerBias?: number; };
}

const DEFAULT_OPTIONS: Required<GungnirOptions> = {
  maxGuesses: 1000000, hashBits: 64, maxSubstitutions: 3, maxIndels: 2,
  errorHints: { likelyPositions: [], homopolymerBias: 0 },
};

// xxHash64 implementation
const P1 = 0x9E3779B185EBCA87n, P2 = 0xC2B2AE3D27D4EB4Fn;
const P3 = 0x165667B19E3779F9n, P4 = 0x85EBCA77C2B2AE63n, P5 = 0x27D4EB2F165667C5n;

function rotl64(x: bigint, n: number): bigint {
  return ((x << BigInt(n)) | (x >> BigInt(64 - n))) & 0xFFFFFFFFFFFFFFFFn;
}

function xxHash64(data: Uint8Array, seed: bigint = 0n): bigint {
  const len = data.length;
  let acc: bigint;
  if (len >= 32) {
    let acc1 = seed + P1 + P2, acc2 = seed + P2, acc3 = seed, acc4 = seed - P1;
    let i = 0;
    while (i + 32 <= len) {
      const read = (off: number) => {
        let v = 0n;
        for (let j = 0; j < 8; j++) v |= BigInt(data[off + j]) << BigInt(j * 8);
        return v;
      };
      const merge = (acc: bigint, val: bigint) => rotl64((acc ^ rotl64(val * P2, 31)) * P1, 27) * P4 + acc;
      acc1 = merge(acc1, read(i)); acc2 = merge(acc2, read(i + 8));
      acc3 = merge(acc3, read(i + 16)); acc4 = merge(acc4, read(i + 24));
      i += 32;
    }
    acc = rotl64(acc1, 1) + rotl64(acc2, 7) + rotl64(acc3, 12) + rotl64(acc4, 18);
  } else {
    acc = seed + P5;
  }
  acc += BigInt(len);
  let i = len >= 32 ? len - (len % 32) : 0;
  while (i + 8 <= len) {
    let val = 0n;
    for (let j = 0; j < 8; j++) val |= BigInt(data[i + j]) << BigInt(j * 8);
    acc = (rotl64(acc ^ rotl64(val * P2, 31), 27) * P1 + P4) & 0xFFFFFFFFFFFFFFFFn;
    i += 8;
  }
  while (i + 4 <= len) {
    let val = 0n;
    for (let j = 0; j < 4; j++) val |= BigInt(data[i + j]) << BigInt(j * 8);
    acc = (rotl64(acc ^ (val * P1), 23) * P2 + P3) & 0xFFFFFFFFFFFFFFFFn;
    i += 4;
  }
  while (i < len) {
    acc = (rotl64(acc ^ (BigInt(data[i]) * P5), 11) * P1) & 0xFFFFFFFFFFFFFFFFn;
    i++;
  }
  acc ^= acc >> 33n; acc = (acc * P2) & 0xFFFFFFFFFFFFFFFFn;
  acc ^= acc >> 29n; acc = (acc * P3) & 0xFFFFFFFFFFFFFFFFn;
  acc ^= acc >> 32n;
  return acc & 0xFFFFFFFFFFFFFFFFn;
}

const GUNGNIR_SEED = 0xBAD5060n;

export function computeFragmentHash(data: Uint8Array, hashBits: number = 64): bigint {
  const full = xxHash64(data, GUNGNIR_SEED);
  if (hashBits >= 64) return full;
  return full & ((1n << BigInt(hashBits)) - 1n);
}

export function gungnirEncode(fragments: Uint8Array[]): { fragments: Uint8Array[]; hashes: bigint[] } {
  return { fragments, hashes: fragments.map(f => computeFragmentHash(f)) };
}

export function* makeEducatedGuesses(corrupted: Uint8Array, errorHints?: GungnirOptions['errorHints']): Generator<Uint8Array> {
  yield corrupted; // 0 errors
  // 1 substitution
  for (let pos = 0; pos < corrupted.length; pos++) {
    for (let bit = 0; bit < 8; bit++) {
      const guess = new Uint8Array(corrupted);
      guess[pos] ^= (1 << bit);
      yield guess;
    }
  }
  // 1 deletion
  for (let pos = 0; pos < corrupted.length; pos++) {
    const guess = new Uint8Array(corrupted.length - 1);
    guess.set(corrupted.subarray(0, pos), 0);
    guess.set(corrupted.subarray(pos + 1), pos);
    yield guess;
  }
  // 1 insertion
  for (let pos = 0; pos <= corrupted.length; pos++) {
    for (let val = 0; val < 256; val += 64) {
      const guess = new Uint8Array(corrupted.length + 1);
      guess.set(corrupted.subarray(0, pos), 0);
      guess[pos] = val;
      guess.set(corrupted.subarray(pos), pos + 1);
      yield guess;
    }
  }
}

export function gungnirDecode(corrupted: Uint8Array, expectedHash: bigint, options?: GungnirOptions): Uint8Array | null {
  const opts = { ...DEFAULT_OPTIONS, ...options };
  const mask = opts.hashBits >= 64 ? 0xFFFFFFFFFFFFFFFFn : (1n << BigInt(opts.hashBits)) - 1n;
  let guesses = 0;
  for (const guess of makeEducatedGuesses(corrupted, opts.errorHints)) {
    if (++guesses > opts.maxGuesses) return null;
    if ((computeFragmentHash(guess, opts.hashBits) & mask) === (expectedHash & mask)) return guess;
  }
  return null;
}

export class GungnirDecoder {
  private options: Required<GungnirOptions>;
  public stats = { totalDecodes: 0, totalGuesses: 0, successes: 0, lastGuesses: 0 };
  constructor(options?: GungnirOptions) { this.options = { ...DEFAULT_OPTIONS, ...options }; }
  decode(corrupted: Uint8Array, expectedHash: bigint): Uint8Array | null {
    let guesses = 0; const mask = this.options.hashBits >= 64 ? 0xFFFFFFFFFFFFFFFFn : (1n << BigInt(this.options.hashBits)) - 1n;
    for (const guess of makeEducatedGuesses(corrupted, this.options.errorHints)) {
      if (++guesses > this.options.maxGuesses) break;
      if ((computeFragmentHash(guess, this.options.hashBits) & mask) === (expectedHash & mask)) {
        this.stats.totalDecodes++; this.stats.totalGuesses += guesses; this.stats.successes++; this.stats.lastGuesses = guesses;
        return guess;
      }
    }
    this.stats.totalDecodes++; this.stats.totalGuesses += guesses; this.stats.lastGuesses = guesses;
    return null;
  }
  encode(fragments: Uint8Array[]) { return gungnirEncode(fragments); }
  hash(data: Uint8Array) { return computeFragmentHash(data, this.options.hashBits); }
}
