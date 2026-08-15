/**
 * Real Zstandard (de)compression via WASM.
 *
 * Uses @bokuweb/zstd-wasm which ships a pre-compiled zstd C library
 * (compiled with Emscripten). The `index.node.js` entry point handles
 * reading the WASM binary and initializing the Emscripten Module.
 *
 * Usage:
 *   import { initZstdWasm, zstdCompress, zstdDecompress } from './zstd-wasm';
 *   await initZstdWasm();
 *   const compressed = zstdCompress(data, 3);
 *   const decompressed = zstdDecompress(compressed);
 *
 * The compressed output starts with 0x28 0xB5 0x2F 0xFD — the zstd magic
 * number — making it fully compatible with the zstd command-line tool.
 */

/** Whether the WASM module has been initialized. */
let initialized = false;

/** Cached reference to the loaded zstd module compress/decompress. */
let zstdApi: {
  compress: (buf: Uint8Array, level?: number) => Uint8Array;
  decompress: (buf: Uint8Array, opts?: { defaultHeapSize?: number }) => Uint8Array;
} | null = null;

/**
 * Load the zstd-wasm index.node.js module.
 * Tries pkg/ first, then wasm-pkg/ as fallback.
 */
function loadZstdModule(): any | null {
  // Try pkg/zstd-wasm/index.node.js first
  try {
    // eslint-disable-next-line @typescript-eslint/no-require-imports
    const mod = require('./pkg/zstd-wasm/index.node.js');
    if (mod && typeof mod.init === 'function' && typeof mod.compress === 'function') {
      return mod;
    }
  } catch {
    // Not available from pkg/
  }

  // Try wasm-pkg/zstd-wasm/index.node.js
  try {
    // eslint-disable-next-line @typescript-eslint/no-require-imports
    const mod = require('./wasm-pkg/zstd-wasm/index.node.js');
    if (mod && typeof mod.init === 'function' && typeof mod.compress === 'function') {
      return mod;
    }
  } catch {
    // Not available from wasm-pkg/
  }

  return null;
}

/**
 * Initialize the zstd WASM module.
 * Must be called once before using zstdCompress/zstdDecompress.
 *
 * Tries both `pkg/zstd-wasm/` and `wasm-pkg/zstd-wasm/` directories.
 *
 * @returns true if initialization succeeded
 */
export async function initZstdWasm(): Promise<boolean> {
  if (initialized) return true;

  try {
    const mod = loadZstdModule();

    if (!mod) {
      throw new Error('Could not find zstd-wasm/index.node.js in pkg/ or wasm-pkg/');
    }

    // Call init() — this reads the WASM binary (zstd.wasm) and initializes
    // the Emscripten Module. After this, compress/decompress work.
    await mod.init();

    // Store references to compress and decompress
    zstdApi = {
      compress: mod.compress,
      decompress: mod.decompress,
    };

    initialized = true;
    return true;
  } catch (err) {
    console.warn('[zstd-wasm] Failed to initialize zstd WASM:', err);
    initialized = false;
    zstdApi = null;
    return false;
  }
}

/**
 * Compress data using real zstd WASM.
 *
 * @param data Input bytes to compress
 * @param level Compression level (1-22, default 3). Level 3 = default zstd.
 * @returns Compressed bytes in true zstd frame format (starts with 0x28 0xB5 0x2F 0xFD)
 * @throws Error if WASM not initialized
 */
export function zstdCompress(data: Uint8Array, level: number = 3): Uint8Array {
  if (!initialized || !zstdApi) {
    throw new Error('zstd WASM not initialized — call initZstdWasm() first');
  }
  return zstdApi.compress(data, level);
}

/**
 * Decompress data using real zstd WASM.
 *
 * @param data Compressed bytes in zstd frame format
 * @returns Decompressed bytes
 * @throws Error if WASM not initialized or data is not valid zstd
 */
export function zstdDecompress(data: Uint8Array): Uint8Array {
  if (!initialized || !zstdApi) {
    throw new Error('zstd WASM not initialized — call initZstdWasm() first');
  }
  return zstdApi.decompress(data);
}

/** Check if zstd WASM is initialized and available. */
export function isZstdWasmReady(): boolean {
  return initialized;
}
