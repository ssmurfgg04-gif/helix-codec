/**
 * Native Rust module loader for the new (non-Viterbi) Rust modules.
 *
 * Provides lazy-loaded access to pack/bhe/compress/ecc/simulate native
 * functions, with graceful fallback to the existing TS implementations.
 *
 * Loaded from the same .so/.dylib as the Viterbi addon — all functions are
 * exported by the single helix-dna-napi crate.
 */
import * as path from 'node:path';
import * as fs from 'node:fs';

export interface NativeAddon {
  // Viterbi (already wired)
  napiVersion(): string;
  convK7Encode(data: Buffer): Buffer;
  convK9Encode(data: Buffer): Buffer;
  viterbiK7Decode(received: Buffer, config?: any): Buffer;
  viterbiK9Decode(received: Buffer, config?: any): Buffer;
  viterbiK9DecodeStandard(received: Buffer): Buffer;
  viterbiK9DecodeWithLlr(received: Buffer, llr: Float32Array, config?: any): Buffer;
  // Pack (new)
  packDnaToBits(dna: string): Uint8Array;
  unpackBitsToDna(bits: Uint8Array, numBases: number): string;
  complementPacked(bits: Uint8Array): Uint8Array;
  reverseComplementPacked(bits: Uint8Array, numBases: number): Uint8Array;
  bitParallelHamming(a: Uint8Array, b: Uint8Array): number;
  rollingHash(bits: Uint8Array, windowSize: number): Uint32Array;
  gcContent(dna: string): number;
  maxHomopolymerRun(dna: string): number;
  // BHE (new)
  bheEncode(data: Uint8Array, config?: { maxRun?: number; enforceGC?: boolean; gcMin?: number; gcMax?: number }): string;
  bheDecode(dna: string, expectedLen: number): Uint8Array;
  // Compress (new)
  compressZstd(data: Uint8Array, level?: number): Uint8Array;
  decompressZstd(data: Uint8Array): Uint8Array;
  isAlreadyCompressed(data: Uint8Array): boolean;
  // ECC (new)
  rsEncode(data: Uint8Array, nsym: number): Uint8Array;
  rsDecode(codeword: Uint8Array, nsym: number): Uint8Array;
  rsDecodeErasures(codeword: Uint8Array, nsym: number, erasurePositions: Uint32Array): Uint8Array;
  rsParity(data: Uint8Array, nsym: number): Uint8Array;
  rsVersion(): string;
  // Simulate (new)
  simulateOligoReads(oligoSeq: string, config?: any): Uint8Array;
  simulateBasic(seq: string, sub: number, ins: number, del: number, seed: number): string;
  readStats(reads: Uint8Array): Float64Array;
  simulateVersion(): string;
}

let _addon: NativeAddon | null = null;
let _loadAttempted = false;

function projectRootCandidates(): string[] {
  const roots = new Set<string>();
  try { roots.add(path.resolve(__dirname, '../../../../')); } catch { /* */ }
  try {
    const url = (import.meta as any)?.url;
    if (typeof url === 'string') {
      roots.add(path.resolve(new URL(url).pathname, '../../../../'));
    }
  } catch { /* */ }
  try {
    let dir = process.cwd();
    for (let i = 0; i < 8; i++) {
      if (fs.existsSync(path.join(dir, 'rust/helix-dna-napi'))) {
        roots.add(dir);
        break;
      }
      const parent = path.dirname(dir);
      if (parent === dir) break;
      dir = parent;
    }
  } catch { /* */ }
  return Array.from(roots);
}

export function tryLoadAddon(): NativeAddon | null {
  if (_loadAttempted) return _addon;
  _loadAttempted = true;
  try {
    const candidates: string[] = [];
    for (const root of projectRootCandidates()) {
      candidates.push(path.resolve(root, 'rust/helix-dna-napi/target/release/libhelix_dna_napi.so'));
      candidates.push(path.resolve(root, 'rust/helix-dna-napi/target/release/libhelix_dna_napi.dylib'));
    }
    for (const addonPath of candidates) {
      if (!fs.existsSync(addonPath)) continue;
      try {
        const mod = { exports: {} };
        (process as any).dlopen(mod, addonPath);
        const addon = mod.exports as NativeAddon;
        if (addon && typeof addon.packDnaToBits === 'function') {
          _addon = addon;
          return _addon;
        }
      } catch (e: any) {
        try { console.warn(`[helix-napi] dlopen failed for ${addonPath}: ${e.message}`); } catch { /* */ }
      }
    }
    return null;
  } catch { return null; }
}

export function getNativeAddon(): NativeAddon | null {
  return tryLoadAddon();
}

export function isNativeAvailable(): boolean {
  return tryLoadAddon() !== null;
}
