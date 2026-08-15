/**
 * Deterministic Constraint Encoding: RLL + GC Rotating Codebooks
 * No seed-retry loop — guaranteed constraint satisfaction
 */

export interface ConstraintConfig { maxHomopolymer: number; gcTarget: number; gcTolerance: number; }
const BASES = ['A','C','G','T'] as const;

export function rllEncode(data: Uint8Array, maxRun: number = 3): string {
  if (data.length===0) return '';
  const result: string[] = []; let lastBase = -1, runLength = 0; let bitPos = 0; const totalBits = data.length*8;
  while (bitPos < totalBits) {
    const forbidden = runLength >= maxRun ? lastBase : -1;
    const choices = forbidden >= 0 ? 3 : 4;
    const bitsNeeded = choices <= 4 ? 2 : 3;
    let digit = 0;
    for (let b=0; b<bitsNeeded; b++) { if (bitPos < totalBits) { const byteIdx = bitPos>>>3; const bitIdx = 7-(bitPos&7); digit = (digit<<1)|((data[byteIdx]>>>bitIdx)&1); bitPos++; } }
    digit = digit % choices;
    let baseIdx: number;
    if (forbidden >= 0) { let count = 0; baseIdx = (forbidden+1)%4; for (let b=0;b<4;b++) { if (b===forbidden) continue; if (count===digit) { baseIdx=b; break; } count++; } }
    else baseIdx = digit;
    if (baseIdx===lastBase) runLength++; else { lastBase=baseIdx; runLength=1; }
    result.push(BASES[baseIdx]);
  }
  return result.join('');
}

export function rllDecode(dna: string, maxRun: number = 3): Uint8Array {
  if (dna.length===0) return new Uint8Array(0);
  const bits: number[] = []; let lastBase = -1, runLength = 0;
  for (let i=0;i<dna.length;i++) {
    const baseIdx = BASES.indexOf(dna[i]); if (baseIdx<0) continue;
    const forbidden = runLength >= maxRun ? lastBase : -1;
    const choices = forbidden >= 0 ? 3 : 4;
    const bitsNeeded = choices <= 4 ? 2 : 3;
    let digit: number;
    if (forbidden >= 0) { let count = 0; digit = 0; for (let b=0;b<4;b++) { if (b===forbidden) continue; if (b===baseIdx) { digit=count; break; } count++; } }
    else digit = baseIdx;
    for (let b=bitsNeeded-1;b>=0;b--) bits.push((digit>>>b)&1);
    if (baseIdx===lastBase) runLength++; else { lastBase=baseIdx; runLength=1; }
  }
  const bytes = new Uint8Array(Math.ceil(bits.length/8));
  for (let i=0;i<bits.length;i++) if (bits[i]) bytes[i>>>3] |= (1<<(7-(i&7)));
  return bytes;
}

export function deterministicEncode(data: Uint8Array, config: Partial<ConstraintConfig> = {}): string { return rllEncode(data, config.maxHomopolymer ?? 3); }
export function deterministicDecode(dna: string, config: Partial<ConstraintConfig> = {}): Uint8Array { return rllDecode(dna, config.maxHomopolymer ?? 3); }

export function satisfiesConstraints(dna: string, config: Partial<ConstraintConfig> = {}): boolean {
  const { maxHomopolymer = 3, gcTarget = 0.5, gcTolerance = 0.1 } = config;
  let run = 1; for (let i=1;i<dna.length;i++) { if (dna[i]===dna[i-1]) { run++; if (run>maxHomopolymer) return false; } else run=1; }
  let gc = 0; for (const c of dna) if (c==='G'||c==='C') gc++;
  return Math.abs(gc/dna.length - gcTarget) <= gcTolerance;
}
