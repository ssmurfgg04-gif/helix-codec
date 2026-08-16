/**
 * WASM DNA Compressor Bridge — registerDnaCompressorWasm()
 *
 * Loads the Emscripten-compiled C++ DNA compressors (GeCo2, DNA-Diff, DNA-QLS)
 * and registers them in the codec's compression pipeline.
 *
 * Usage:
 *   import { registerDnaCompressorWasm, getWasmCompressors } from './dna-compress-wasm';
 *   await registerDnaCompressorWasm();  // loads WASM module
 *   const compressors = getWasmCompressors();
 *   // Each compressor has .name, .compress(data), .decompress(data)
 */

import { readFileSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

// Type for a registered compressor
export interface DnaCompressor {
  /** Algorithm name */
  name: string;
  /** Algorithm index in WASM module */
  algoIndex: number;
  /** Compress DNA ASCII bytes → compressed bytes */
  compress(data: Uint8Array): Uint8Array;
  /** Decompress compressed bytes → DNA ASCII bytes */
  decompress(data: Uint8Array): Uint8Array;
}

// Module state
let wasmModule: any = null;
let wasmCompressors: DnaCompressor[] = [];
let wasmLoaded = false;

/**
 * Load the DNA compressors WASM module.
 *
 * Tries multiple paths to find the WASM file:
 *   1. Adjacent to this source file
 *   2. In wasm-src/ directory
 *   3. In the project root
 */
async function loadWasmModule(): Promise<any> {
  // Try to import the Emscripten-generated JS glue
  const __dirname = dirname(fileURLToPath(import.meta.url));

  const candidatePaths = [
    join(__dirname, '..', '..', '..', 'wasm-src', 'dna_compressors.js'),
    join(__dirname, 'dna_compressors.js'),
    join(process.cwd(), 'wasm-src', 'dna_compressors.js'),
  ];

  let mod: any = null;
  for (const p of candidatePaths) {
    try {
      const createModule = await import(p);
      mod = await createModule.default();
      break;
    } catch {
      continue;
    }
  }

  if (!mod) {
    // Fallback: load the WASM binary directly and instantiate
    const wasmPaths = [
      join(__dirname, '..', '..', '..', 'wasm-src', 'dna_compressors.wasm'),
      join(__dirname, 'dna_compressors.wasm'),
      join(process.cwd(), 'wasm-src', 'dna_compressors.wasm'),
    ];

    for (const p of wasmPaths) {
      try {
        const wasmBinary = readFileSync(p);
        // We need the JS glue for Emscripten modules, so this won't work
        // for a full Emscripten module. Let's use the Node.js approach.
        const { default: createModule } = await import(
          join(__dirname, '..', '..', '..', 'wasm-src', 'dna_compressors.js')
        );
        mod = await createModule({
          locateFile: (name: string) => {
            if (name.endsWith('.wasm')) return p;
            return name;
          },
        });
        break;
      } catch {
        continue;
      }
    }
  }

  return mod;
}

/**
 * Register WASM-compiled DNA compressors into the codec pipeline.
 *
 * This is the API that the user called to enable C++ DNA compressors
 * compiled to WASM via Emscripten, replacing the TypeScript arithmetic
 * coding implementations with native C++ performance.
 *
 * @returns Array of registered compressors
 */
export async function registerDnaCompressorWasm(): Promise<DnaCompressor[]> {
  if (wasmLoaded) return wasmCompressors;

  try {
    wasmModule = await loadWasmModule();

    if (!wasmModule) {
      console.warn('[dna-compress-wasm] Could not load WASM module. ' +
        'Run scripts/build-compressors-wasm.sh first.');
      return [];
    }

    const count = wasmModule._dna_compressor_count();
    wasmCompressors = [];

    for (let i = 0; i < count; i++) {
      const namePtr = wasmModule._dna_compressor_name(i);
      const name = wasmModule.UTF8ToString(namePtr);

      const compressor: DnaCompressor = {
        name,
        algoIndex: i,
        compress(data: Uint8Array): Uint8Array {
          const inputLen = data.length;
          const outputCap = inputLen + inputLen / 4 + 256;  // conservative upper bound

          const inputPtr = wasmModule._malloc(inputLen);
          const outputPtr = wasmModule._malloc(outputCap);

          try {
            // Copy input to WASM memory
            wasmModule.HEAPU8.set(data, inputPtr);

            const compressedLen = wasmModule._dna_compress(
              i, inputPtr, inputLen, outputPtr, outputCap
            );

            if (compressedLen < 0) {
              throw new Error(`WASM compression failed for algorithm ${name}`);
            }

            // Copy output from WASM memory
            return wasmModule.HEAPU8.slice(outputPtr, outputPtr + compressedLen);
          } finally {
            wasmModule._free(inputPtr);
            wasmModule._free(outputPtr);
          }
        },

        decompress(data: Uint8Array): Uint8Array {
          const inputLen = data.length;
          // Decompressed size is at most the original DNA length.
          // Read from the header (first 4 bytes are magic, next 4 are length).
          let outputCap = inputLen * 4 + 1024;  // conservative default
          if (inputLen >= 8) {
            const origLen =
              data[4] | (data[5] << 8) | (data[6] << 16) | (data[7] << 24);
            if (origLen > 0 && origLen < 1_000_000_000) {
              outputCap = origLen + 64;
            }
          }

          const inputPtr = wasmModule._malloc(inputLen);
          const outputPtr = wasmModule._malloc(outputCap);

          try {
            wasmModule.HEAPU8.set(data, inputPtr);

            const decompressedLen = wasmModule._dna_decompress(
              i, inputPtr, inputLen, outputPtr, outputCap
            );

            if (decompressedLen < 0) {
              throw new Error(`WASM decompression failed for algorithm ${name}`);
            }

            return wasmModule.HEAPU8.slice(outputPtr, outputPtr + decompressedLen);
          } finally {
            wasmModule._free(inputPtr);
            wasmModule._free(outputPtr);
          }
        },
      };

      wasmCompressors.push(compressor);
    }

    wasmLoaded = true;
    console.log(`[dna-compress-wasm] Registered ${count} WASM DNA compressors: ` +
      wasmCompressors.map(c => c.name).join(', '));

  } catch (err) {
    console.warn('[dna-compress-wasm] Failed to load WASM module:', (err as Error).message);
    console.warn('Run scripts/build-compressors-wasm.sh to compile the C++ compressors.');
  }

  return wasmCompressors;
}

/**
 * Get the currently registered WASM DNA compressors.
 * Returns empty array if registerDnaCompressorWasm() hasn't been called.
 */
export function getWasmCompressors(): DnaCompressor[] {
  return wasmCompressors;
}

/**
 * Check if WASM DNA compressors are loaded.
 */
export function isWasmCompressorsLoaded(): boolean {
  return wasmLoaded;
}

/**
 * Get the raw WASM module (for advanced usage).
 */
export function getWasmModule(): any {
  return wasmModule;
}
