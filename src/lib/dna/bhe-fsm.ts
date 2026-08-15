/**
 * P0: BHE Deterministic FSM Encoding
 * Based on Microsoft DNABoundedHomopolymerEncoding
 * Variable-base BigInt arithmetic: no seed-retry loop
 */

const BASES = ['A', 'C', 'G', 'T'] as const;
const BASE_MAP: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

export function validChoices(lastBase: number, runLength: number, k: number): number {
  return runLength >= k ? 3 : 4;
}

export function getForbiddenBase(lastBase: number, runLength: number, k: number): number {
  return runLength >= k ? lastBase : -1;
}

export function mapDigitToBase(digit: number, forbiddenBase: number): number {
  if (forbiddenBase < 0) return digit;
  let count = 0;
  for (let b = 0; b < 4; b++) {
    if (b === forbiddenBase) continue;
    if (count === digit) return b;
    count++;
  }
  return (forbiddenBase + 1) % 4;
}

export function baseToDigit(base: number, forbiddenBase: number): number {
  if (forbiddenBase < 0) return base;
  let count = 0;
  for (let b = 0; b < 4; b++) {
    if (b === forbiddenBase) continue;
    if (b === base) return count;
    count++;
  }
  return 0;
}

export function bheEncode(data: Uint8Array, maxHomopolymer: number): string {
  if (data.length === 0) return '';
  const marked = new Uint8Array(data.length + 2);
  marked[0] = 0xFF; marked[1] = data.length & 0xFF;
  marked.set(data, 2);
  let N = 0n;
  for (const b of marked) N = (N << 8n) | BigInt(b);
  if (N === 0n) return 'A';
  const k = maxHomopolymer;
  const digits: number[] = [];
  let lastBase = -1, runLength = 0;
  while (N > 0n) {
    const forbidden = getForbiddenBase(lastBase, runLength, k);
    const choices = BigInt(validChoices(lastBase, runLength, k));
    const digit = Number(N % choices);
    N = N / choices;
    const baseIdx = mapDigitToBase(digit, forbidden);
    if (baseIdx === lastBase) runLength++; else { lastBase = baseIdx; runLength = 1; }
    digits.push(baseIdx);
  }
  return digits.reverse().map(b => BASES[b]).join('');
}

export function bheDecode(dna: string, maxHomopolymer: number): Uint8Array {
  if (dna.length === 0) return new Uint8Array(0);
  const k = maxHomopolymer;
  let N = 0n;
  let lastBase = -1, runLength = 0;
  for (let i = dna.length - 1; i >= 0; i--) {
    const baseIdx = BASE_MAP[dna[i]] ?? 0;
    const forbidden = getForbiddenBase(lastBase, runLength, k);
    const choices = BigInt(validChoices(lastBase, runLength, k));
    const digit = BigInt(baseToDigit(baseIdx, forbidden));
    N = N * choices + digit;
    if (i > 0) {
      const prevIdx = BASE_MAP[dna[i - 1]] ?? 0;
      if (prevIdx === baseIdx) runLength++; else { lastBase = baseIdx; runLength = 1; }
    }
  }
  const bytes: number[] = [];
  while (N > 0n) { bytes.unshift(Number(N & 0xFFn)); N = N >> 8n; }
  if (bytes.length >= 2 && bytes[0] === 0xFF) bytes.shift();
  if (bytes.length > 0) bytes.shift();
  return new Uint8Array(bytes);
}

export function validateHomopolymer(dna: string, k: number): boolean {
  let run = 1;
  for (let i = 1; i < dna.length; i++) {
    if (dna[i] === dna[i - 1]) { run++; if (run > k) return false; } else run = 1;
  }
  return true;
}

export class BHECodebook {
  readonly k: number;
  readonly stateCount: number;
  readonly encodeTable: number[][];
  readonly decodeTable: number[][];

  constructor(k: number) {
    this.k = k;
    this.stateCount = 1 + 4 * k;
    this.encodeTable = [];
    this.decodeTable = [];
    for (let s = 0; s < this.stateCount; s++) {
      this.encodeTable[s] = [];
      this.decodeTable[s] = [];
    }
  }

  encode(data: Uint8Array): string { return bheEncode(data, this.k); }
  decode(dna: string): Uint8Array { return bheDecode(dna, this.k); }
}
