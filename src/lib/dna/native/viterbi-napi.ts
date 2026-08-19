/**
 * Native Viterbi Decoder — Hybrid native + JS approach
 *
 * v4.1: Fixed del_pen=1.5 (was 1.0, caused spurious D paths), drift_pen=0.5
 *
 * Strategy:
 *   - Clean channel (no indels): Use native standard Viterbi (~0.5ms)
 *   - Noisy channel (with indels): Use native indel-tolerant Viterbi (~100-200ms for K=9)
 *   - Fallback: JS IndelTolerantConvolutionalInnerCode
 */

import * as path from 'node:path';
import * as fs from 'node:fs';

export interface ViterbiNapiConfig {
  maxDrift?: number;
  insertionPenalty?: number;
  deletionPenalty?: number;
  useLlr?: boolean;
  expectedLength?: number;
  /** Number of information bits. If not provided, estimated from received length. */
  numInfoBits?: number;
}

interface NativeAddon {
  viterbiK9DecodeStandard(received: Buffer): Buffer;
  viterbiK9Decode(received: Buffer, config?: ViterbiNapiConfig): Buffer;
  viterbiK7Decode(received: Buffer, config?: ViterbiNapiConfig): Buffer;
  viterbiK9DecodeWithLlr(received: Buffer, llr: Float32Array, config?: ViterbiNapiConfig): Buffer;
  convK9Encode(data: Buffer): Buffer;
  convK7Encode(data: Buffer): Buffer;
  napiVersion(): string;
}

let _addon: NativeAddon | null = null;
let _loadAttempted = false;

// Resolve project root robustly. In tsx/ESM, `__dirname` is shimmed to the
// file's directory, so going up 4 levels from src/lib/dna/native/ lands at
// project root. We also try a few fallbacks to be safe.
function projectRootCandidates(): string[] {
  const roots = new Set<string>();
  try { roots.add(path.resolve(__dirname, '../../../../')); } catch { /* */ }
  try {
    // import.meta.url — works under native ESM (node --experimental-vm-modules etc.)
    const url = (import.meta as any)?.url;
    if (typeof url === 'string') {
      roots.add(path.resolve(new URL(url).pathname, '../../../../'));
    }
  } catch { /* */ }
  // Walk up from cwd looking for rust/helix-dna-napi
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

function tryLoadAddon(): NativeAddon | null {
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
        if (addon && typeof addon.viterbiK9DecodeStandard === 'function') {
          _addon = addon;
          return _addon;
        }
      } catch (e: any) {
        // Log the dlopen failure so silent misconfiguration is visible
        try { console.warn(`[viterbi-napi] dlopen failed for ${addonPath}: ${e.message}`); } catch { /* */ }
      }
    }
    return null;
  } catch { return null; }
}

export async function enableNativeViterbi(): Promise<boolean> {
  const addon = tryLoadAddon();
  if (addon) {
    try { console.log(`[viterbi-napi] ${addon.napiVersion()}`); } catch {}
  }
  return addon !== null;
}

export function isNativeViterbiActive(): boolean { return _addon !== null; }

/**
 * Standard K=9 Viterbi decode (no indel tolerance) — native, ~0.5ms.
 * Use when the channel is clean or after MSA consensus has removed indels.
 */
export function nativeViterbiK9DecodeStandard(received: Uint8Array | Buffer): Buffer {
  if (!_addon) throw new Error('Native Viterbi addon not loaded.');
  return _addon.viterbiK9DecodeStandard(Buffer.isBuffer(received) ? received : Buffer.from(received));
}

/**
 * Indel-tolerant K=9 Viterbi decode — native, ~100-200ms.
 * Production-hardened with I-chain propagation, zero-tail penalty, correct traceback.
 */
export function nativeViterbiK9Decode(received: Uint8Array | Buffer, config?: ViterbiNapiConfig): Buffer {
  if (!_addon) throw new Error('Native Viterbi addon not loaded.');
  return _addon.viterbiK9Decode(Buffer.isBuffer(received) ? received : Buffer.from(received), config);
}

/**
 * Indel-tolerant K=7 Viterbi decode — native.
 */
export function nativeViterbiK7Decode(received: Uint8Array | Buffer, config?: ViterbiNapiConfig): Buffer {
  if (!_addon) throw new Error('Native Viterbi addon not loaded.');
  return _addon.viterbiK7Decode(Buffer.isBuffer(received) ? received : Buffer.from(received), config);
}

/**
 * K=9 convolutional encode — native, <1ms.
 */
export function nativeConvK9Encode(data: Uint8Array | Buffer): Buffer {
  if (!_addon) throw new Error('Native Viterbi addon not loaded.');
  return _addon.convK9Encode(Buffer.isBuffer(data) ? data : Buffer.from(data));
}

/**
 * K=7 convolutional encode — native, <1ms.
 */
export function nativeConvK7Encode(data: Uint8Array | Buffer): Buffer {
  if (!_addon) throw new Error('Native Viterbi addon not loaded.');
  return _addon.convK7Encode(Buffer.isBuffer(data) ? data : Buffer.from(data));
}

/**
 * Hybrid Viterbi decode: use native standard decoder for near-clean channels,
 * native indel decoder for noisy channels, JS indel decoder as fallback.
 *
 * The decision is based on length comparison: if the received stream length
 * is close to the expected encoded length, use the fast native standard path.
 * Otherwise, use the native indel-tolerant decoder.
 *
 * @param received Received byte stream
 * @param expectedEncodedLen Expected encoded length (from conv encode)
 * @param numInfoBits Number of information bits (for indel decoder)
 * @param maxDrift Maximum drift to track (default 15)
 * @param jsIndelDecoder JS IndelTolerantConvolutionalInnerCode instance (fallback)
 * @returns Decoded byte stream
 */
export function hybridViterbiDecode(
  received: Uint8Array,
  expectedEncodedLen: number,
  numInfoBits: number,
  maxDrift: number = 15,
  jsIndelDecoder?: { decode(received: Uint8Array): { decoded: Uint8Array } },
): Uint8Array {
  const lenRatio = received.length / expectedEncodedLen;

  // If received length is within 2% of expected, use fast native standard decode
  if (_addon && lenRatio >= 0.98 && lenRatio <= 1.02) {
    try {
      const result = _addon.viterbiK9DecodeStandard(Buffer.from(received));
      return new Uint8Array(result);
    } catch {
      // Native standard decode failed, try indel
    }
  }

  // Use native indel-tolerant Viterbi for noisy channels
  if (_addon) {
    try {
      const result = _addon.viterbiK9Decode(Buffer.from(received), {
        maxDrift,
        insertionPenalty: 1.5,
        deletionPenalty: 1.5, // v4.1: MUST equal ins_pen (was 1.0 — caused spurious D paths)
        numInfoBits,
      });
      return new Uint8Array(result);
    } catch {
      // Native indel decode failed, try JS fallback
    }
  }

  // Use JS indel decoder as last resort
  if (jsIndelDecoder) {
    try {
      const { decoded } = jsIndelDecoder.decode(received);
      return decoded;
    } catch {
      // All decoders failed — return received as-is (best effort)
      return received;
    }
  }

  throw new Error('No Viterbi decoder available.');
}
