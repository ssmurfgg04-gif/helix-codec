#!/bin/bash
# Build DNA compressors C++ -> WASM using Emscripten
#
# Compiles api.cpp (which includes geco.cpp, dna_diff.cpp, dna_qls.cpp)
# into a single WASM module with JS glue code.
#
# Output:
#   wasm-src/dna_compressors.js     - ES6 module glue code
#   wasm-src/dna_compressors.wasm   - WASM binary
#   wasm-src/dna_compressors.d.ts   - TypeScript declarations

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Source Emscripten
source /home/z/emsdk/emsdk_env.sh 2>/dev/null || true

# Verify emcc is available
if ! command -v emcc &>/dev/null; then
    echo "ERROR: emcc not found. Install Emscripten first:"
    echo "  git clone https://github.com/emscripten-core/emsdk.git /home/z/emsdk"
    echo "  cd /home/z/emsdk && ./emsdk install latest && ./emsdk activate latest"
    exit 1
fi

echo "=== Building DNA Compressors WASM ==="
echo "Emscripten version: $(emcc --version 2>/dev/null | head -1)"

CPP_SRC="$PROJECT_DIR/wasm-src/cpp-compressors/api.cpp"
OUTPUT_DIR="$PROJECT_DIR/wasm-src"
OUTPUT_BASE="$OUTPUT_DIR/dna_compressors"

# Compile with Emscripten
emcc "$CPP_SRC" \
    -o "$OUTPUT_BASE.js" \
    -O3 \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s INITIAL_MEMORY=16MB \
    -s MAXIMUM_MEMORY=256MB \
    -s EXPORTED_FUNCTIONS='["_dna_compress","_dna_decompress","_dna_compressor_name","_dna_compressor_count","_malloc","_free"]' \
    -s EXPORTED_RUNTIME_METHODS='["UTF8ToString","stringToUTF8","HEAPU8"]' \
    -s FILESYSTEM=0 \
    -s MINIFY_HTML=0 \
    -Wall \
    -Wno-unused-function \
    -std=c++17

echo "=== Build complete ==="
echo "Output files:"
ls -la "$OUTPUT_BASE.js" "$OUTPUT_BASE.wasm" 2>/dev/null

# Generate TypeScript declarations
cat > "$OUTPUT_BASE.d.ts" << 'DTS'
// DNA Compressors WASM — TypeScript declarations
export interface DnaCompressorsModule {
  _dna_compress(algo: number, input: number, input_len: number, output: number, output_cap: number): number;
  _dna_decompress(algo: number, input: number, input_len: number, output: number, output_cap: number): number;
  _dna_compressor_name(algo: number): number;
  _dna_compressor_count(): number;
  _malloc(size: number): number;
  _free(ptr: number): void;
  HEAPU8: Uint8Array;
  UTF8ToString(ptr: number): string;
}

export default function createModule(): Promise<DnaCompressorsModule>;
DTS

echo "TypeScript declarations written to $OUTPUT_BASE.d.ts"
echo "=== All done ==="
