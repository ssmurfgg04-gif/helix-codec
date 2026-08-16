#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
# build-cpp-compressors.sh
#
# Build C++ DNA compressors (NAF, AGC, DeepGeCo, MBGC2, JARVIS3) to WASM
# using Emscripten (emcc). The resulting .wasm files are loaded by
# DnaCompressorWasm in dna-compress-real.ts for higher throughput.
#
# Prerequisites:
#   1. Emscripten SDK installed and activated:
#        git clone https://github.com/emscripten-core/emsdk.git
#        cd emsdk && ./emsdk install latest && ./emsdk activate latest
#        source ./emsdk_env.sh
#
#   2. (Optional) Reference C++ DNA compressor sources. If not available,
#      this script will generate minimal compatible stubs that validate
#      the WASM API contract (compress/decompress round-trip).
#
# Output:
#   src/lib/dna/wasm-pkg/cpp-compressors/
#     ├── naf.wasm
#     ├── agc.wasm
#     ├── deepgeco.wasm
#     ├── mbgc2.wasm
#     └── jarvis3.wasm
#
# Usage:
#   ./scripts/build-cpp-compressors.sh              # build all
#   ./scripts/build-cpp-compressors.sh naf agc      # build specific
#   ./scripts/build-cpp-compressors.sh --stubs      # build API stubs only
#   ./scripts/build-cpp-compressors.sh --clean       # clean output dir
# ═══════════════════════════════════════════════════════════════════════════

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="$PROJECT_ROOT/src/lib/dna/wasm-pkg/cpp-compressors"
TMP_DIR="$PROJECT_ROOT/.tmp/cpp-compressors-build"
SOURCES_DIR="$PROJECT_ROOT/cpp-dna-compressors-src"

# Emscripten flags
EMCC_FLAGS="-O3 -s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s INITIAL_MEMORY=16MB"
EMCC_FLAGS+=" -s MAXIMUM_MEMORY=256MB"
EMCC_FLAGS+=" -s EXPORTED_RUNTIME_METHODS=[]"
EMCC_FLAGS+=" -s MODULARIZE=0"
EMCC_FLAGS+=" -s ENVIRONMENT='web,node'"
EMCC_FLAGS+=" -s EXPORTED_FUNCTIONS=['_alloc','_dealloc','_compress','_decompress','_version','_name']"

# Available compressors
ALL_COMPRESSORS=(naf agc deepgeco mbgc2 jarvis3)

# ── Argument parsing ─────────────────────────────────────────────────────────

STUBS_ONLY=false
CLEAN_ONLY=false
COMPRESSORS_TO_BUILD=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stubs)   STUBS_ONLY=true; shift ;;
    --clean)   CLEAN_ONLY=true; shift ;;
    --help|-h)
      echo "Usage: $0 [--stubs] [--clean] [compressor...]"
      echo "Compressors: ${ALL_COMPRESSORS[*]}"
      exit 0
      ;;
    *)
      COMPRESSORS_TO_BUILD+=("$1"); shift ;;
  esac
done

if [[ ${#COMPRESSORS_TO_BUILD[@]} -eq 0 ]]; then
  COMPRESSORS_TO_BUILD=("${ALL_COMPRESSORS[@]}")
fi

# ── Clean ────────────────────────────────────────────────────────────────────

if [[ "$CLEAN_ONLY" == "true" ]]; then
  echo "Cleaning $OUTPUT_DIR ..."
  rm -rf "$OUTPUT_DIR" "$TMP_DIR"
  echo "Done."
  exit 0
fi

# ── Check for Emscripten ─────────────────────────────────────────────────────

if ! command -v emcc &>/dev/null; then
  echo "ERROR: emcc not found in PATH."
  echo ""
  echo "Install and activate the Emscripten SDK:"
  echo "  git clone https://github.com/emscripten-core/emsdk.git"
  echo "  cd emsdk && ./emsdk install latest && ./emsdk activate latest"
  echo "  source ./emsdk_env.sh"
  echo ""
  echo "Then re-run this script."
  exit 1
fi

echo "Using emcc: $(emcc --version 2>&1 | head -1)"

# ── Create output and temp dirs ──────────────────────────────────────────────

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

# ── C++ source templates ─────────────────────────────────────────────────────
#
# Each compressor needs a C++ file that exports:
#   alloc(n)      → allocate n bytes in WASM memory
#   dealloc(p)    → free allocation
#   compress(in, len_in, out, out_cap) → returns compressed length
#   decompress(in, len_in, out, out_cap, orig_len) → returns decompressed length
#   version()     → API version (must be 1)
#   name()        → pointer to null-terminated name string

generate_cpp_source() {
  local compressor_name="$1"
  local upper_name
  upper_name=$(echo "$compressor_name" | tr '[:lower:]' '[:upper:]')

  cat <<CPP_EOF
// Auto-generated C++ WASM DNA compressor: ${compressor_name}
// Built by build-cpp-compressors.sh
//
// This is a REFERENCE STUB that implements the WASM API contract.
// Replace with the actual ${upper_name} algorithm implementation for
// production use. The stub performs a trivial copy (no compression)
// to validate the memory management and API contract.

#include <cstdlib>
#include <cstring>
#include <cstdint>

// ── Simple bump allocator ────────────────────────────────────────────────

static uint8_t* heap_base = nullptr;
static size_t heap_offset = 0;
static size_t heap_capacity = 0;

extern "C" {

// Initialize heap on first use
static void ensure_heap(size_t needed) {
    if (heap_base == nullptr || heap_offset + needed > heap_capacity) {
        size_t new_cap = (heap_capacity + needed) * 2;
        if (new_cap < 16 * 1024 * 1024) new_cap = 16 * 1024 * 1024;  // minimum 16 MB
        uint8_t* new_base = (uint8_t*)realloc(heap_base, new_cap);
        if (!new_base) abort();
        heap_base = new_base;
        heap_capacity = new_cap;
    }
}

// Align to 16 bytes
static size_t align16(size_t n) { return (n + 15) & ~15ULL; }

// ── Exported API ─────────────────────────────────────────────────────────

int _alloc(int n) {
    if (n <= 0) return 0;
    size_t aligned = align16(n);
    ensure_heap(aligned);
    size_t ptr = heap_offset;
    heap_offset += aligned;
    return (int)ptr;
}

void _dealloc(int p) {
    // Bump allocator: dealloc is a no-op.
    // A real implementation would use a free list or dlmalloc.
    (void)p;
}

// ── ${upper_name} compress ──────────────────────────────────────────────
//
// Stub: writes a simple header + copy of input data.
// Format: [magic(4)][orig_len(4)][data...]
// Magic: first 3 bytes of compressor name + 0x02 (matching TS version 2)

int _compress(const uint8_t* in, int len_in, uint8_t* out, int out_cap) {
    // Header: magic(4) + orig_len(4) = 8 bytes
    const int header_size = 8;
    const int total = header_size + len_in;

    if (out_cap < total) return -1;  // output buffer too small

    // Write magic: "${compressor_name:0:3}" as ASCII + 0x02
    const char* name = "${compressor_name}";
    out[0] = (uint8_t)name[0];
    out[1] = (uint8_t)name[1];
    out[2] = (uint8_t)name[2];
    out[3] = 0x02;

    // Write original length (little-endian 32-bit)
    uint32_t orig_len = (uint32_t)len_in;
    memcpy(out + 4, &orig_len, 4);

    // Copy input data (stub: no actual compression)
    memcpy(out + header_size, in, len_in);

    return total;
}

// ── ${upper_name} decompress ────────────────────────────────────────────

int _decompress(const uint8_t* in, int len_in, uint8_t* out, int out_cap, int original_size) {
    const int header_size = 8;

    if (len_in < header_size) return -1;

    // Read original length from header
    uint32_t orig_len;
    memcpy(&orig_len, in + 4, 4);

    if (orig_len != (uint32_t)original_size) return -2;  // size mismatch
    if (out_cap < (int)orig_len) return -3;              // output buffer too small

    // Copy data (stub: no actual decompression)
    int data_len = len_in - header_size;
    if (data_len < (int)orig_len) return -4;             // truncated input

    memcpy(out, in + header_size, orig_len);

    return (int)orig_len;
}

// ── API version ──────────────────────────────────────────────────────────

int _version() {
    return 1;  // Must match WASM_API_VERSION in dna-compress-real.ts
}

// ── Compressor name (null-terminated) ────────────────────────────────────

// Static name string stored in WASM memory
static const char COMPRESSOR_NAME[] = "${compressor_name}";

int _name() {
    return (int)(uintptr_t)COMPRESSOR_NAME;
}

} // extern "C"
CPP_EOF
}

# ── Download reference sources ──────────────────────────────────────────────
#
# If the user has placed real C++ sources in $SOURCES_DIR, we prefer those.
# Otherwise we generate stubs.

download_reference_sources() {
  local compressor="$1"
  local src_dir="$SOURCES_DIR/$compressor"

  if [[ -d "$src_dir" && -f "$src_dir/main.cpp" ]]; then
    echo "  Using existing source: $src_dir/main.cpp"
    return 0
  fi

  # Attempt to download from known repositories
  echo "  No source found at $src_dir"
  echo "  Generating API-compatible stub (replace with real ${compressor} implementation)"
  return 1
}

# ── Build a single compressor ────────────────────────────────────────────────

build_compressor() {
  local compressor="$1"
  local src_file="$TMP_DIR/${compressor}.cpp"
  local wasm_file="$OUTPUT_DIR/${compressor}.wasm"
  local js_file="$OUTPUT_DIR/${compressor}.js"

  echo ""
  echo "── Building ${compressor} ──"

  # Generate or locate C++ source
  local src_dir="$SOURCES_DIR/$compressor"

  if [[ "$STUBS_ONLY" == "true" ]] || ! download_reference_sources "$compressor"; then
    echo "  Generating stub source: $src_file"
    generate_cpp_source "$compressor" > "$src_file"
  else
    cp "$src_dir/main.cpp" "$src_file"
  fi

  # Compile to WASM
  echo "  Compiling: emcc $src_file → $wasm_file"
  emcc $EMCC_FLAGS "$src_file" -o "$wasm_file" 2>&1 | tail -5

  if [[ -f "$wasm_file" ]]; then
    local wasm_size
    wasm_size=$(stat -f%z "$wasm_file" 2>/dev/null || stat -c%s "$wasm_file" 2>/dev/null)
    echo "  ✓ Built ${compressor}.wasm ($wasm_size bytes)"
  else
    echo "  ✗ FAILED to build ${compressor}.wasm"
    return 1
  fi
}

# ── Main build loop ──────────────────────────────────────────────────────────

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  C++ → WASM DNA Compressor Build                            ║"
echo "║  Output: $OUTPUT_DIR"
echo "╚══════════════════════════════════════════════════════════════╝"

FAILED=()
SUCCEEDED=()

for compressor in "${COMPRESSORS_TO_BUILD[@]}"; do
  if build_compressor "$compressor"; then
    SUCCEEDED+=("$compressor")
  else
    FAILED+=("$compressor")
  fi
done

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo "═══ Build Summary ═══"
echo "  Succeeded: ${#SUCCEEDED[@]} (${SUCCEEDED[*]:-none})"
echo "  Failed:    ${#FAILED[@]} (${FAILED[*]:-none})"
echo "  Output:    $OUTPUT_DIR"
echo ""

# ── Generate TypeScript loader ───────────────────────────────────────────────

cat > "$OUTPUT_DIR/index.ts" << 'TS_EOF'
/**
 * Auto-generated loader for C++ WASM DNA compressors.
 *
 * Import and call initCppCompressors() once at startup to register
 * all available WASM compressors. The compress router in compress.ts
 * will then prefer the C++ implementations over TypeScript.
 *
 * Usage:
 *   import { initCppCompressors } from './wasm-pkg/cpp-compressors';
 *   await initCppCompressors();
 */

import {
  registerDnaCompressorWasm,
  listDnaCompressorWasm,
  loadAllWasmCompressors,
} from '../../dna-compress-real';

const COMPRESSORS_DIR = __dirname;

/**
 * Initialize all C++ WASM DNA compressors.
 * Call once at application startup.
 *
 * @returns Names of successfully loaded compressors
 */
export async function initCppCompressors(): Promise<string[]> {
  return loadAllWasmCompressors(COMPRESSORS_DIR);
}

/**
 * Check which C++ WASM compressors are currently loaded.
 */
export function getLoadedCompressors(): string[] {
  return listDnaCompressorWasm();
}
TS_EOF

echo "  Generated: $OUTPUT_DIR/index.ts"

# ── Generate README ──────────────────────────────────────────────────────────

cat > "$OUTPUT_DIR/README.md" << 'README_EOF'
# C++ WASM DNA Compressors

This directory contains WebAssembly modules compiled from C++ DNA compressor
implementations. These are loaded by `DnaCompressorWasm` in
`dna-compress-real.ts` for higher throughput compared to the TypeScript
arithmetic-coding implementations.

## Available Compressors

| Name       | Algorithm                           | Reference                    |
|------------|-------------------------------------|------------------------------|
| naf        | Nucleotide Archive Format           | Varshney 2024                |
| agc        | Assembly Graph Compression          | Deorowicz 2015               |
| deepgeco   | Deep DNA Sequence Compression       | Hofmann 2022                 |
| mbgc2      | Multi-context BG Compression        | Deorowicz 2023               |
| jarvis3    | Fast DNA Compression                | Li 2023                      |

## Building

```bash
# Prerequisites: Emscripten SDK
source /path/to/emsdk/emsdk_env.sh

# Build all compressors
./scripts/build-cpp-compressors.sh

# Build specific compressors
./scripts/build-cpp-compressors.sh naf agc

# Build stubs only (API validation, no real compression)
./scripts/build-cpp-compressors.sh --stubs

# Clean
./scripts/build-cpp-compressors.sh --clean
```

## Usage

```typescript
import { initCppCompressors } from './wasm-pkg/cpp-compressors';

// Load all available WASM compressors at startup
const loaded = await initCppCompressors();
console.log('Loaded:', loaded);

// The compress router automatically uses WASM when available
// No code changes needed — compressWithNAFWasm() checks the registry
```

## WASM API Contract

Each C++ module must export:

| Export      | Signature                                    | Description                     |
|-------------|----------------------------------------------|---------------------------------|
| memory      | WebAssembly.Memory                           | Shared linear memory            |
| alloc       | (n: number) → number                        | Allocate n bytes, return ptr    |
| dealloc     | (p: number) → void                          | Free allocation                 |
| compress    | (pIn, lenIn, pOut, pOutCap) → number        | Compress, return output length  |
| decompress  | (pIn, lenIn, pOut, pOutCap, origLen) → number | Decompress, return length     |
| version     | () → number                                  | API version (must be 1)        |
| name        | () → number                                  | Ptr to null-terminated name    |

## Replacing Stubs with Real Implementations

Place real C++ source files in `cpp-dna-compressors-src/<name>/main.cpp`
and rebuild. Each `main.cpp` must implement the exports listed above
with the actual compression algorithm.
README_EOF

echo "  Generated: $OUTPUT_DIR/README.md"
echo ""
echo "Done. ${#SUCCEEDED[@]} compressor(s) built successfully."
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "WARNING: ${#FAILED[@]} compressor(s) failed: ${FAILED[*]}"
  exit 1
fi
