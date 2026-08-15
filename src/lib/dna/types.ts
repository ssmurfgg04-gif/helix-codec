/**
 * Helix DNA Storage — Core Type Definitions
 * v3.0 architecture with Babel-USB content-derived addressing
 */

export const DNA_ALPHABET = ['A', 'C', 'G', 'T'] as const;
export type DNABase = typeof DNA_ALPHABET[number];
export type ChannelPreset = 'illumina' | 'nanopore' | 'pacbio';
export type AddressMode = 'sequential' | 'content-derived' | 'hierarchical';
export type InnerCode = 'rs' | 'ldpc' | 'arithmetic';
export type MappingMode = 'standard' | 'yinyang';

export interface CodecConfig {
  channel?: ChannelPreset;
  oligoLength?: number;
  primerLength?: number;
  indexLength?: number;
  outerParityRatio?: number;
  gfOrder?: 8 | 16;
  innerCode?: InnerCode;
  useConvolutionalInner?: boolean;
  useViterbiSoft?: boolean;
  useOSD?: boolean;
  interleaveDepth?: number;
  useBHE?: boolean;
  bheMaxHomopolymer?: number;
  gungnirMode?: boolean;
  gungnirHashBits?: number;
  targetDensity?: number;
  mappingMode?: MappingMode;
  addressMode?: AddressMode;
  archiveSalt?: Uint8Array;
  useRecipeGeneration?: boolean;
  maxHomopolymer?: number;
  gcTarget?: number;
  gcTolerance?: number;
  encryptBeforeCompress?: boolean;
  arithmeticBlockSize?: number;
  crcMarkerInterval?: number;
  simulationProfile?: SimulationProfile;
}

export interface SimulationProfile {
  synthesis: { technology: string; subRate: number; delRate: number; insRate: number; biasModel?: string; };
  pcr: { cycles: number; duplicationBias: number; errorRate: number; };
  aging: { years: number; decayRate: number; gcBias: number; };
  sequencing: { technology: string; readLength: number; errorProfile: ErrorProfile; };
}

export interface ErrorProfile {
  subRate: number; delRate: number; insRate: number; homopolymerBias: number;
  contextErrors?: Map<string, number>;
}

export const ILLUMINA_CONFIG: Partial<CodecConfig> = { channel: 'illumina', oligoLength: 200, primerLength: 20, outerParityRatio: 0.1, interleaveDepth: 4, innerCode: 'rs', gfOrder: 8 };
export const NANOPORE_CONFIG: Partial<CodecConfig> = { channel: 'nanopore', oligoLength: 200, primerLength: 20, outerParityRatio: 0.5, interleaveDepth: 8, innerCode: 'ldpc', useConvolutionalInner: true, useViterbiSoft: true, useBHE: true, gungnirMode: true, addressMode: 'content-derived', targetDensity: 0.99 };
export const PACBIO_CONFIG: Partial<CodecConfig> = { channel: 'pacbio', oligoLength: 200, primerLength: 20, outerParityRatio: 0.4, interleaveDepth: 6, innerCode: 'ldpc', useConvolutionalInner: true };

export function resolveConfig(partial: CodecConfig): CodecConfig {
  const defaults: CodecConfig = { oligoLength: 200, primerLength: 20, indexLength: 12, outerParityRatio: 0.1, gfOrder: 8, innerCode: 'rs', interleaveDepth: 4, bheMaxHomopolymer: 3, gungnirHashBits: 64, targetDensity: 1.0, mappingMode: 'standard', addressMode: 'sequential', maxHomopolymer: 3, gcTarget: 0.5, gcTolerance: 0.1, encryptBeforeCompress: true, arithmeticBlockSize: 32, crcMarkerInterval: 16 };
  let preset: Partial<CodecConfig> = {};
  if (partial.channel === 'illumina') preset = ILLUMINA_CONFIG;
  else if (partial.channel === 'nanopore') preset = NANOPORE_CONFIG;
  else if (partial.channel === 'pacbio') preset = PACBIO_CONFIG;
  return { ...defaults, ...preset, ...partial };
}

export interface CanonicalBlock { index: number; address: string; payload: string; crc32: number; verified: boolean; }
export interface CanonicalArchive { version: number; config: CodecConfig; blocks: CanonicalBlock[]; salt: Uint8Array; createdAt: number; }
export interface EncodeResult { oligos: string[]; config: CodecConfig; archive?: CanonicalArchive; stats: { inputBytes: number; outputOligos: number; encodingTimeMs: number; compressionRatio: number; netDensityBitsPerNt: number; }; }
export interface DecodeResult { data: Uint8Array; stats: { inputOligos: number; outputBytes: number; decodingTimeMs: number; oligosFailedOuterRS: number; oligosFailedInner: number; erasuresUsed: number; strategy: string; }; }
export interface HierarchicalAddress { pool: string; well: string; oligoIndex: number; contentHash: string; }
export type RecipeKind = 'constant' | 'repeat' | 'deBruijn' | 'seededPRNG' | 'data';
export interface OligoRecipe { kind: RecipeKind; length: number; params: Record<string, unknown>; }
