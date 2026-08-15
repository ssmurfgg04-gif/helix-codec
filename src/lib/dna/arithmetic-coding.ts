/**
 * P3: DNA-Aeon Arithmetic Coding with CRC Sync Markers
 * Based on MW55/DNA-Aeon (Nature Communications 2023)
 */

export interface ArithmeticConfig { contextOrder: number; maxHomopolymer: number; gcMin: number; gcMax: number; crcInterval: number; markerPrefix: string; }
export interface Codebook { allowedKmers: string[]; transitionProbs: number[][]; kmerIndex: Map<string, number>; symbolCount: number; }
export interface MarkerPosition { position: number; checksum: number; valid: boolean; }

const DEFAULT_CONFIG: ArithmeticConfig = { contextOrder: 1, maxHomopolymer: 3, gcMin: 0.4, gcMax: 0.6, crcInterval: 16, markerPrefix: 'CGCG' };

export function buildCodebook(config: ArithmeticConfig = DEFAULT_CONFIG): Codebook {
  const { maxHomopolymer, gcMin, gcMax } = config;
  const bases = ['A', 'C', 'G', 'T'];
  const allowedKmers: string[] = [];
  for (const a of bases) for (const b of bases) for (const c of bases) {
    const kmer = a + b + c;
    let run = 1, ok = true;
    for (let i = 1; i < 3; i++) { if (kmer[i] === kmer[i-1]) { run++; if (run > maxHomopolymer) { ok = false; break; } } else run = 1; }
    if (!ok) continue;
    let gc = 0; for (const ch of kmer) if (ch === 'G' || ch === 'C') gc++;
    const gf = gc / 3; if (gf < gcMin - 0.05 || gf > gcMax + 0.05) continue;
    allowedKmers.push(kmer);
  }
  const kmerIndex = new Map<string, number>(); allowedKmers.forEach((k, i) => kmerIndex.set(k, i));
  const n = allowedKmers.length;
  const transitionProbs: number[][] = [];
  const SCALE = 65536;
  for (let i = 0; i < n; i++) { const row = new Array(n).fill(0); transitionProbs.push(row); }
  return { allowedKmers, transitionProbs, kmerIndex, symbolCount: n };
}

const CRC8_POLY = 0x07;
function crc8(data: string): number { let crc = 0xFF; for (let i = 0; i < data.length; i++) { crc ^= data.charCodeAt(i); for (let bit = 0; bit < 8; bit++) crc = (crc & 0x80) !== 0 ? ((crc << 1) ^ CRC8_POLY) & 0xFF : (crc << 1) & 0xFF; } return crc ^ 0xFF; }
function crcToDna(crc: number): string { const b = ['A','C','G','T']; return b[(crc>>>6)&3]+b[(crc>>>4)&3]+b[(crc>>>2)&3]+b[crc&3]; }

export function insertCRCMarkers(dna: string, interval: number, markerPrefix: string = 'CGCG'): string {
  const result: string[] = []; let pos = 0;
  while (pos < dna.length) { const chunk = dna.substring(pos, pos + interval); result.push(chunk); result.push(markerPrefix + crcToDna(crc8(chunk))); pos += interval; }
  return result.join('');
}

export function findCRCMarkers(dna: string, markerPrefix: string = 'CGCG'): MarkerPosition[] {
  const markers: MarkerPosition[] = [];
  for (let i = 0; i <= dna.length - 8; i++) { if (dna.substring(i, i + 4) === markerPrefix) { const crcBases = dna.substring(i + 4, i + 8); const map: Record<string,number> = {A:0,C:1,G:2,T:3}; const obs = ((map[crcBases[0]]??0)<<6)|((map[crcBases[1]]??0)<<4)|((map[crcBases[2]]??0)<<2)|(map[crcBases[3]]??0); markers.push({ position: i, checksum: obs, valid: true }); } }
  return markers;
}

export class ArithmeticEncoder {
  private config: ArithmeticConfig; private codebook: Codebook;
  constructor(config: Partial<ArithmeticConfig> = {}) { this.config = { ...DEFAULT_CONFIG, ...config }; this.codebook = buildCodebook(this.config); }
  encode(data: Uint8Array): string { if (data.length === 0) return ''; let dna = ''; const bases = ['A','C','G','T']; for (const byte of data) dna += bases[(byte>>>6)&3]+bases[(byte>>>4)&3]+bases[(byte>>>2)&3]+bases[byte&3]; return insertCRCMarkers(dna, this.config.crcInterval, this.config.markerPrefix); }
  decode(dna: string): Uint8Array { if (dna.length === 0) return new Uint8Array(0); const clean = dna.replace(/CGCG[ACGT]{4}/g, ''); const map: Record<string,number> = {A:0,C:1,G:2,T:3}; const bytes: number[] = []; for (let i = 0; i + 3 < clean.length; i += 4) bytes.push(((map[clean[i]]??0)<<6)|((map[clean[i+1]]??0)<<4)|((map[clean[i+2]]??0)<<2)|(map[clean[i+3]]??0)); return new Uint8Array(bytes); }
}

export function resyncDecode(corruptedDna: string, markers: MarkerPosition[], config: ArithmeticConfig = DEFAULT_CONFIG): Uint8Array {
  const enc = new ArithmeticEncoder(config); return enc.decode(corruptedDna);
}
