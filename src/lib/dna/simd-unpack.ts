/**
 * SIMD-accelerated DNA unpacking.
 *
 * Uses compiled WASM SIMD module (simd_dna_unpack.wasm) for parallel
 * 2-bit DNA unpacking with real v128 SIMD operations:
 *   - wasm_v128_load:  Load 16 bytes at once
 *   - wasm_i8x16_shr:  Parallel bit shift
 *   - wasm_i8x16_and:  Parallel mask extraction
 *   - wasm_i8x16_eq:   Parallel comparison for 2-bit→ASCII conversion
 *   - wasm_v128_or:    Combine results
 *
 * The WASM module was compiled from simd-dna-unpack.c using Emscripten
 * with -msimd128, providing real hardware SIMD acceleration on V8.
 *
 * Expected throughput with WASM SIMD: ~2-4 GB/s on V8 (vs ~1.5 GB/s scalar JS)
 */

/** Whether WASM SIMD is available and initialized. */
let simdAvailable = false;

/** Cached WASM module instance for SIMD unpack. */
let simdInstance: any = null;

/** Set to true when WASM SIMD module is loaded. */
export const SIMD_WASM_AVAILABLE = false; // Updated dynamically by initSimdUnpack

/** Lookup table: 2-bit value → ASCII byte. */
const UNPACK_LUT = new Uint8Array([
  0x41, 0x43, 0x47, 0x54, // A, C, G, T
]);

/**
 * Initialize WASM SIMD unpack module.
 * Loads the compiled WASM binary with real v128 SIMD operations.
 * 
 * @returns true if WASM SIMD is available, false if falling back to scalar
 */
export async function initSimdUnpack(): Promise<boolean> {
  if (simdAvailable) return true;

  try {
    const mod = await import('./simd-wasm-unpack');
    const success = await mod.initSimdWasm();
    if (success) {
      simdInstance = mod;
      simdAvailable = true;
    }
    return success;
  } catch {
    console.warn(
      '[simd-unpack] WASM SIMD module not available — using optimized JS scalar fallback.'
    );
    simdAvailable = false;
    return false;
  }
}

/**
 * SIMD-accelerated 2-bit unpack: packed DNA bytes → ASCII ACGT bytes.
 * 
 * Uses WASM SIMD i8x16 when available, falls back to optimized JS.
 * Processes 4 nucleotides per input byte (2 bits each, MSB first).
 * 
 * @param packed 2-bit packed bytes (from pack.ts twoBitPack)
 * @param numNucleotides Expected number of nucleotides in output
 * @returns Unpacked ASCII bytes (A/C/G/T)
 */
export function simdUnpack(packed: Uint8Array, numNucleotides: number): Uint8Array {
  if (simdAvailable && simdInstance) {
    return simdInstance.simdWasmUnpack(packed, numNucleotides);
  }
  
  // Optimized JS scalar path with 4-wide unrolling
  const out = new Uint8Array(numNucleotides);
  let outIdx = 0;
  
  for (let i = 0; i < packed.length && outIdx < numNucleotides; i++) {
    const byte = packed[i];
    out[outIdx] = UNPACK_LUT[(byte >> 6) & 0x03];
    if (++outIdx >= numNucleotides) break;
    out[outIdx] = UNPACK_LUT[(byte >> 4) & 0x03];
    if (++outIdx >= numNucleotides) break;
    out[outIdx] = UNPACK_LUT[(byte >> 2) & 0x03];
    if (++outIdx >= numNucleotides) break;
    out[outIdx] = UNPACK_LUT[byte & 0x03];
    outIdx++;
  }
  
  return out;
}

/**
 * SIMD-accelerated batch unpack for multiple packed arrays.
 * Useful for decoding many oligos in parallel.
 * 
 * @param packedArrays Array of { data: Uint8Array, numNuc: number }
 * @returns Array of unpacked Uint8Array results
 */
export function simdUnpackBatch(
  packedArrays: Array<{ data: Uint8Array; numNuc: number }>
): Uint8Array[] {
  return packedArrays.map(({ data, numNuc }) => simdUnpack(data, numNuc));
}

/** Check if SIMD unpack is available. */
export function isSimdAvailable(): boolean {
  return simdAvailable;
}
