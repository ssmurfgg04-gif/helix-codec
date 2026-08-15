/**
 * P6: Yin-Yang Coding (YYC) for DNA Storage
 * Based on BGI-research/Chamaeleo YYC
 * Two binary streams → one DNA sequence via dual rules
 */

export interface YYCConfig { rotationIndex?: number; yinRule?: number[]; yangRule?: number[]; }
const DEFAULT_YIN = [0,1,2,3], DEFAULT_YANG = [3,2,1,0];
const BASES = ['A','C','G','T'] as const;
const BASE_MAP: Record<string,number> = {A:0,C:1,G:2,T:3};

function allPerms(): number[][] { const perms: number[][] = []; function hp(k: number, a: number[]) { if (k===1) { perms.push([...a]); return; } for (let i=0;i<k;i++) { hp(k-1,a); if (k%2===0) [a[i],a[k-1]]=[a[k-1],a[i]]; else [a[0],a[k-1]]=[a[k-1],a[0]]; } } hp(4,[0,1,2,3]); return perms; }
const ALL_PERMS = allPerms();

export function validateOrthogonality(yin: number[], yang: number[]): boolean { const pairs = new Set<string>(); for (let i=0;i<4;i++) { const key=`${yin[i]},${yang[i]}`; if (pairs.has(key)) return false; pairs.add(key); } return pairs.size===4; }

let _validPairs: Array<{yin:number[];yang:number[];index:number}>|null = null;
function getValidPairs() { if (_validPairs) return _validPairs; const pairs: Array<{yin:number[];yang:number[];index:number}> = []; let idx = 0; for (const yin of ALL_PERMS) for (const yang of ALL_PERMS) if (validateOrthogonality(yin,yang)) { pairs.push({yin,yang,index:idx}); idx++; } _validPairs = pairs; return pairs; }

export function generateRulePair(rotationIndex: number): {yin:number[];yang:number[]} { const pairs = getValidPairs(); const idx = ((rotationIndex % pairs.length) + pairs.length) % pairs.length; return { yin: pairs[idx].yin, yang: pairs[idx].yang }; }

function resolveConfig(config?: YYCConfig): {yin:number[];yang:number[]} { if (config?.yinRule && config?.yangRule) return {yin:config.yinRule,yang:config.yangRule}; if (config?.rotationIndex !== undefined) return generateRulePair(config.rotationIndex); return {yin:DEFAULT_YIN,yang:DEFAULT_YANG}; }

function getBit(data: Uint8Array, bitPos: number): number { const byteIdx = bitPos>>>3; const bitIdx = 7-(bitPos&7); if (byteIdx>=data.length) return 0; return (data[byteIdx]>>>bitIdx)&1; }
function setBit(data: Uint8Array, bitPos: number, value: number): void { const byteIdx = bitPos>>>3; const bitIdx = 7-(bitPos&7); if (byteIdx>=data.length) return; if (value) data[byteIdx]|=(1<<bitIdx); else data[byteIdx]&=~(1<<bitIdx); }

export function yycEncode(data: Uint8Array, config?: YYCConfig): string { if (data.length===0) return ''; const {yin,yang} = resolveConfig(config); const yinInv = new Array<number>(4), yangInv = new Array<number>(4); for (let i=0;i<4;i++) { yinInv[yin[i]]=i; yangInv[yang[i]]=i; } const totalBits = data.length*8; const result: string[] = []; for (let bitPos=0; bitPos+3<totalBits; bitPos+=4) { const yinVal = (getBit(data,bitPos)<<1)|getBit(data,bitPos+1); const yangVal = (getBit(data,bitPos+2)<<1)|getBit(data,bitPos+3); let foundBase = -1; for (let b=0;b<4;b++) if (yinInv[b]===yinVal && yangInv[b]===yangVal) { foundBase=b; break; } if (foundBase===-1) foundBase = yin[yinVal]; result.push(BASES[foundBase]); } return result.join(''); }

export function yycDecode(dna: string, config?: YYCConfig): Uint8Array { if (dna.length===0) return new Uint8Array(0); const {yin,yang} = resolveConfig(config); const totalBits = dna.length*4; const totalBytes = Math.ceil(totalBits/8); const result = new Uint8Array(totalBytes); for (let i=0;i<dna.length;i++) { const baseIdx = BASE_MAP[dna[i]]??0; const yinVal = yin[baseIdx]; const yangVal = yang[baseIdx]; const bitPos = i*4; setBit(result,bitPos,(yinVal>>>1)&1); setBit(result,bitPos+1,yinVal&1); setBit(result,bitPos+2,(yangVal>>>1)&1); setBit(result,bitPos+3,yangVal&1); } return result; }
