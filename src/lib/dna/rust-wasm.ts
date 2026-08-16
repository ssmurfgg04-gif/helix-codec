/**
 * Rust WASM acceleration layer — thin TypeScript wrapper.
 *
 * This module loads the compiled Rust WASM (helix-dna-wasm) and provides
 * a clean TS API that mirrors the original pure-TS modules. When WASM
 * is unavailable, it falls back to the original TS implementations.
 *
 * The Rust WASM contains the 9% hot-path code:
 *   - pack:    2-bit pack/unpack with SIMD
 *   - ecc:     RS GF(256) + LDPC belief propagation
 *   - compress: Arithmetic coding with context models
 *   - bhe:     u128 bit-parallel FSM encoding
 *   - simulate: Parallel per-oligo stochastic simulation
 */

import type * as WasmTypes from "./wasm-pkg-rust/helix_dna_wasm";

let wasmPkg: typeof WasmTypes | null = null;
let initPromise: Promise<boolean> | null = null;

/**
 * Initialize the Rust WASM module. Call once at startup.
 * Returns true if WASM loaded successfully, false if falling back to JS.
 */
export async function initRustWasm(): Promise<boolean> {
  if (wasmPkg) return true;
  if (initPromise) return initPromise;

  initPromise = (async () => {
    try {
      const pkg = await import("./wasm-pkg-rust/helix_dna_wasm.js");
      await pkg.default(); // WASM init
      wasmPkg = pkg as unknown as typeof WasmTypes;
      return true;
    } catch (e) {
      console.warn("helix-dna-wasm load failed, using JS fallback:", e);
      return false;
    }
  })();

  return initPromise;
}

/** Check if Rust WASM is available. */
export function isRustWasmReady(): boolean {
  return wasmPkg !== null;
}

// ---------------------------------------------------------------------------
// Pack / Unpack
// ---------------------------------------------------------------------------

/**
 * Pack DNA string to 2-bit bytes using Rust WASM (SIMD-accelerated).
 * Falls back to TS if WASM unavailable.
 */
export function rustPackDnaToBits(dna: string): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.pack_dna_to_bits(dna);
}

/**
 * Unpack 2-bit bytes to DNA ASCII using Rust WASM (SIMD-accelerated).
 * Falls back to TS if WASM unavailable.
 */
export function rustUnpackBitsToAscii(bits: Uint8Array, numBases: number): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.unpack_bits_to_ascii(bits, numBases);
}

/**
 * Unpack 2-bit bytes to DNA string using Rust WASM.
 */
export function rustUnpackBitsToDna(bits: Uint8Array, numBases: number): string | null {
  if (!wasmPkg) return null;
  return wasmPkg.unpack_bits_to_dna(bits, numBases);
}

// ---------------------------------------------------------------------------
// Hamming distance
// ---------------------------------------------------------------------------

export function rustBitParallelHamming(a: Uint8Array, b: Uint8Array): number | null {
  if (!wasmPkg) return null;
  return wasmPkg.bit_parallel_hamming(a, b);
}

// ---------------------------------------------------------------------------
// Complement
// ---------------------------------------------------------------------------

export function rustComplement(bits: Uint8Array): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.complement(bits);
}

// ---------------------------------------------------------------------------
// Rolling hash
// ---------------------------------------------------------------------------

export function rustRollingHash(bits: Uint8Array, windowSize: number): Uint32Array | null {
  if (!wasmPkg) return null;
  const result = wasmPkg.rolling_hash(bits, windowSize);
  return new Uint32Array(result.buffer, result.byteOffset, result.length);
}

// ---------------------------------------------------------------------------
// BHE encoding (u128 FSM)
// ---------------------------------------------------------------------------

export function rustBheEncodeK1(bytes: Uint8Array): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.bhe_encode_k1(bytes);
}

export function rustBheDecodeK1(dna: Uint8Array, numBytes: number): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.bhe_decode_k1(dna, numBytes);
}

export function rustBheEncodeFsm(bytes: Uint8Array, maxRun: number): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.bhe_encode_fsm(bytes, maxRun);
}

export function rustBheDecodeFsm(dna: Uint8Array, numBytes: number, maxRun: number): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.bhe_decode_fsm(dna, numBytes, maxRun);
}

export function rustBheValidate(dna: Uint8Array, maxRun: number): boolean | null {
  if (!wasmPkg) return null;
  return wasmPkg.bhe_validate(dna, maxRun);
}

// ---------------------------------------------------------------------------
// Arithmetic compression
// ---------------------------------------------------------------------------

export function rustArithCompress(data: Uint8Array): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.arith_compress(data);
}

export function rustArithDecompress(data: Uint8Array): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.arith_decompress(data);
}

export function rustDnaCompressOrder1(packed: Uint8Array, numBases: number): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.dna_compress_order1(packed, numBases);
}

export function rustDnaDecompressOrder1(data: Uint8Array): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.dna_decompress_order1(data);
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

export function rustSimulateSingle(
  oligo: Uint8Array,
  synthSub: number, synthIns: number, synthDel: number,
  seqSub: number, seqIns: number, seqDel: number,
  seed: number,
): Uint8Array | null {
  if (!wasmPkg) return null;
  return wasmPkg.simulate_single(oligo, synthSub, synthIns, synthDel, seqSub, seqIns, seqDel, seed);
}

// ---------------------------------------------------------------------------
// Version
// ---------------------------------------------------------------------------

export function rustWasmVersion(): string | null {
  if (!wasmPkg) return null;
  return wasmPkg.version();
}
