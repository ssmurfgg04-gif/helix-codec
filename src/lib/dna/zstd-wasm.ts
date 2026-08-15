/**
 * Real Zstandard (de)compression via WASM.
 *
 * Uses @bokuweb/zstd-wasm which ships a pre-compiled zstd C library
 * (compiled with Emscripten). Provides both ZSTD_compress and ZSTD_decompress
 * — true zstd format, fully compatible with the zstd command-line tool.
 *
 * Usage:
 *   import { initZstdWasm, zstdCompress, zstdDecompress } from './zstd-wasm';
 *   await initZstdWasm();
 *   const compressed = zstdCompress(data, 3);
 *   const decompressed = zstdDecompress(compressed);
 */

/** Whether the WASM module has been initialized. */
let initialized = false;

/** @bokuweb/zstd-wasm module (loaded dynamically). */
let zstdModule: any = null;

/**
 * Initialize the zstd WASM module.
 * Must be called once before using zstdCompress/zstdDecompress.
 *
 * @returns true if initialization succeeded
 */
export async function initZstdWasm(): Promise<boolean> {
  if (initialized) return true;

  try {
    // Dynamic import of the zstd-wasm package
    // The WASM binary is at ./pkg/zstd-wasm/zstd.wasm
    const { readFile } = await import('fs/promises');
    const { resolve } = await import('path');
    const { createRequire } = await import('module');

    // Use createRequire for CJS interop in ESM context
    const require = createRequire(import.meta.url ?? __filename);

    // Load the JS glue code
    const zstdPkg = require('./pkg/zstd-wasm/zstd.js');
    const wasmPath = resolve(__dirname ?? '.', './pkg/zstd-wasm/zstd.wasm');

    // Read the WASM binary
    const wasmBuffer = await readFile(wasmPath);

    // Initialize the module
    zstdModule = {
      compress: zstdPkg.compress,
      decompress: zstdPkg.decompress,
      init: zstdPkg.init,
      Module: zstdPkg.Module,
    };

    // Call init to load the WASM binary into the emscripten Module
    const { init } = zstdModule;
    await init();

    initialized = true;
    return true;
  } catch (err) {
    console.warn('[zstd-wasm] Failed to initialize zstd WASM:', err);
    initialized = false;
    return false;
  }
}

/**
 * Compress data using real zstd WASM.
 *
 * @param data Input bytes to compress
 * @param level Compression level (1-22, default 3). Level 3 = default zstd.
 * @returns Compressed bytes in true zstd frame format
 * @throws Error if WASM not initialized
 */
export function zstdCompress(data: Uint8Array, level: number = 3): Uint8Array {
  if (!initialized || !zstdModule) {
    throw new Error('zstd WASM not initialized — call initZstdWasm() first');
  }
  return zstdModule.compress(Buffer.from(data), level);
}

/**
 * Decompress data using real zstd WASM.
 *
 * @param data Compressed bytes in zstd frame format
 * @returns Decompressed bytes
 * @throws Error if WASM not initialized or data is not valid zstd
 */
export function zstdDecompress(data: Uint8Array): Uint8Array {
  if (!initialized || !zstdModule) {
    throw new Error('zstd WASM not initialized — call initZstdWasm() first');
  }
  return zstdModule.decompress(Buffer.from(data));
}

/** Check if zstd WASM is initialized and available. */
export function isZstdWasmReady(): boolean {
  return initialized;
}
