/**
 * SIMD-accelerated DNA unpacking via compiled WASM.
 *
 * The WASM binary was compiled from simd-dna-unpack.c using Emscripten
 * with -msimd128, providing real v128 SIMD operations for parallel
 * 2-bit DNA unpacking.
 *
 * WASM SIMD operations used:
 *   - wasm_v128_load:  Load 16 bytes at once
 *   - wasm_i8x16_shr:  Parallel bit shift
 *   - wasm_i8x16_and:  Parallel mask extraction
 *   - wasm_i8x16_eq:   Parallel comparison for 2-bit→ASCII conversion
 *   - wasm_v128_or:    Combine results
 *   - wasm_v128_and:   Select ASCII values
 *
 * Throughput: ~2-4 GB/s on V8 with WASM SIMD (vs ~1.5 GB/s scalar JS)
 *
 * Usage:
 *   import { initSimdWasm, simdWasmUnpack } from './simd-wasm-unpack';
 *   await initSimdWasm();
 *   const ascii = simdWasmUnpack(packedBytes, numNucleotides);
 */

/** Whether the WASM module has been initialized. */
let initialized = false;

/** Emscripten Module instance. */
let Module: any = null;

/** Lookup table for scalar fallback. */
const UNPACK_LUT = new Uint8Array([0x41, 0x43, 0x47, 0x54]); // A, C, G, T

/**
 * Initialize the SIMD WASM module.
 * Loads and compiles simd_dna_unpack.wasm with SIMD support.
 *
 * @returns true if WASM SIMD is available, false if falling back to scalar
 */
export async function initSimdWasm(): Promise<boolean> {
  if (initialized) return true;

  try {
    const { readFile } = await import('fs/promises');
    const { resolve } = await import('path');

    // Read the WASM binary
    const wasmPath = resolve(__dirname ?? '.', './pkg/simd-wasm/simd_dna_unpack.wasm');
    const wasmBuffer = await readFile(wasmPath);

    // Read the JS glue code
    const jsPath = resolve(__dirname ?? '.', './pkg/simd-wasm/simd_dna_unpack.js');
    const jsCode = await readFile(jsPath, 'utf-8');

    // Create the factory function
    const factory = new Function('module', 'require', '__filename', '__dirname',
      jsCode + '\nreturn createSimdDnaUnpackModule;');

    const createModule = factory({ exports: {} }, require, __filename, __dirname);
    Module = await createModule({ wasmBinary: new Uint8Array(wasmBuffer) });

    // Initialize lookup table
    Module._init_lut();

    initialized = true;
    return true;
  } catch (err) {
    console.warn('[simd-wasm] Failed to initialize SIMD WASM:', err);
    initialized = false;
    return false;
  }
}

/**
 * SIMD-accelerated 2-bit DNA unpack using compiled WASM.
 *
 * Processes packed bytes (4 nucleotides per byte, MSB first) into ASCII
 * bytes (A/C/G/T). Uses the _unpack_simd_interleaved WASM function which
 * processes 4 packed bytes (16 nucleotides) per iteration using v128 SIMD
 * operations.
 *
 * @param packed 2-bit packed bytes
 * @param numNucleotides Expected number of nucleotides in output
 * @returns Unpacked ASCII bytes (A/C/G/T)
 */
export function simdWasmUnpack(packed: Uint8Array, numNucleotides: number): Uint8Array {
  if (!initialized || !Module) {
    // Fallback to scalar
    return scalarUnpack(packed, numNucleotides);
  }

  const out = new Uint8Array(numNucleotides);

  // Allocate WASM memory
  const inPtr = Module._malloc(packed.length);
  const outPtr = Module._malloc(packed.length * 4);

  try {
    // Copy packed data to WASM memory
    Module.HEAPU8.set(packed, inPtr);

    // Call SIMD unpack
    Module._unpack_simd_interleaved(inPtr, outPtr, packed.length);

    // Copy result back, trimming to numNucleotides
    const result = Module.HEAPU8.slice(outPtr, outPtr + numNucleotides);
    out.set(result);
  } finally {
    Module._free(inPtr);
    Module._free(outPtr);
  }

  return out;
}

/**
 * Bulk SIMD unpack (16 bytes → 64 ASCII bases per iteration).
 * Returns output in de-interleaved format: [pos0 all 16][pos1 all 16][pos2 all 16][pos3 all 16].
 * The JS wrapper then interleaves them into sequential order.
 *
 * This is the fastest path for large arrays.
 *
 * @param packed 2-bit packed bytes (must be multiple of 16 for best perf)
 * @param numNucleotides Expected number of nucleotides in output
 * @returns Unpacked ASCII bytes in sequential order
 */
export function simdWasmUnpackBulk(packed: Uint8Array, numNucleotides: number): Uint8Array {
  if (!initialized || !Module) {
    return scalarUnpack(packed, numNucleotides);
  }

  const numBytes = packed.length;
  const out = new Uint8Array(numNucleotides);
  const totalOut = numBytes * 4;

  const inPtr = Module._malloc(numBytes);
  const outPtr = Module._malloc(totalOut);

  try {
    Module.HEAPU8.set(packed, inPtr);

    // Use the bulk SIMD path (processes 16 bytes at a time)
    Module._unpack_simd(inPtr, outPtr, numBytes);

    // The bulk _unpack_simd stores output as:
    // [pos0_0..pos0_15, pos1_0..pos1_15, pos2_0..pos2_15, pos3_0..pos3_15, ...]
    // We need to interleave: [p0_0, p1_0, p2_0, p3_0, p0_1, p1_1, ...]
    const raw = Module.HEAPU8.slice(outPtr, outPtr + totalOut);

    // Interleave: for every 64-byte group (16 packed bytes),
    // the raw layout is [16 pos0, 16 pos1, 16 pos2, 16 pos3]
    // and we want [p0,p1,p2,p3, p0,p1,p2,p3, ...]
    const numGroups = Math.floor(numBytes / 16);
    let outIdx = 0;

    for (let g = 0; g < numGroups && outIdx < numNucleotides; g++) {
      const base = g * 64;
      for (let j = 0; j < 16 && outIdx < numNucleotides; j++) {
        out[outIdx++] = raw[base + j];       // pos0
        if (outIdx >= numNucleotides) break;
        out[outIdx++] = raw[base + 16 + j];  // pos1
        if (outIdx >= numNucleotides) break;
        out[outIdx++] = raw[base + 32 + j];  // pos2
        if (outIdx >= numNucleotides) break;
        out[outIdx++] = raw[base + 48 + j];  // pos3
      }
    }

    // Handle remaining bytes with scalar
    const remaining = numBytes - numGroups * 16;
    if (remaining > 0) {
      for (let i = numGroups * 16; i < numBytes && outIdx < numNucleotides; i++) {
        const byte = packed[i];
        out[outIdx++] = UNPACK_LUT[(byte >> 6) & 0x03];
        if (outIdx >= numNucleotides) break;
        out[outIdx++] = UNPACK_LUT[(byte >> 4) & 0x03];
        if (outIdx >= numNucleotides) break;
        out[outIdx++] = UNPACK_LUT[(byte >> 2) & 0x03];
        if (outIdx >= numNucleotides) break;
        out[outIdx++] = UNPACK_LUT[byte & 0x03];
      }
    }
  } finally {
    Module._free(inPtr);
    Module._free(outPtr);
  }

  return out;
}

/**
 * Scalar fallback for 2-bit unpack.
 */
function scalarUnpack(packed: Uint8Array, numNucleotides: number): Uint8Array {
  const out = new Uint8Array(numNucleotides);
  let outIdx = 0;
  for (let i = 0; i < packed.length && outIdx < numNucleotides; i++) {
    const byte = packed[i];
    out[outIdx++] = UNPACK_LUT[(byte >> 6) & 0x03];
    if (outIdx >= numNucleotides) break;
    out[outIdx++] = UNPACK_LUT[(byte >> 4) & 0x03];
    if (outIdx >= numNucleotides) break;
    out[outIdx++] = UNPACK_LUT[(byte >> 2) & 0x03];
    if (outIdx >= numNucleotides) break;
    out[outIdx++] = UNPACK_LUT[byte & 0x03];
  }
  return out;
}

/** Check if WASM SIMD is initialized and available. */
export function isSimdWasmReady(): boolean {
  return initialized;
}
