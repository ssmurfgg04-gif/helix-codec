// Benchmark SIMD batch RS encode vs single encode.
import * as fs from "fs";
import * as path from "path";

async function main() {
  const wasmPath = path.join(__dirname, "..", "rust-dna", "pkg", "helix_dna_wasm_bg.wasm");
  const wasmBuffer = fs.readFileSync(wasmPath);
  let wasmExports: any;

  const importObject = {
    "./helix_dna_wasm_bg.js": {
      __wbg_now_8b265300afd5f2b9: function() { return Date.now(); },
      __wbg___wbindgen_throw_bb96b2010945f0bc: function(arg: number, arg1: number) {
        const mem = new Uint8Array(wasmExports.memory.buffer);
        throw new Error(new TextDecoder().decode(mem.slice(arg, arg + arg1)));
      },
      __wbg___wbindgen_copy_to_typed_array_c7f28e53671b41e8: function() {},
      __wbindgen_init_externref_table: function() {},
    },
  };

  const { instance } = await WebAssembly.instantiate(wasmBuffer, importObject);
  wasmExports = instance.exports;
  wasmExports.init_gf();

  console.log("=== RS Pipeline SIMD Benchmark ===\n");
  console.log(`WASM binary: ${(wasmBuffer.length / 1024).toFixed(0)} KB\n`);

  // Single RS encode
  const singleMbps = wasmExports.bench_rs_encode(10000, 40, 32);
  console.log(`Single RS(40,32):       ${singleMbps.toFixed(1)} MB/s`);

  // Batch RS encode (scalar)
  const batchMbps = wasmExports.bench_batch_rs_encode(64, 40, 32);
  console.log(`Batch RS x64 (scalar):  ${batchMbps.toFixed(1)} MB/s`);

  // SIMD batch RS encode
  const simdBatchMbps = wasmExports.bench_batch_rs_encode_simd(64, 40, 32);
  console.log(`Batch RS x64 (SIMD128): ${simdBatchMbps.toFixed(1)} MB/s`);

  // SIMD GF multiply (raw throughput ceiling)
  const simdMulMbps = wasmExports.bench_simd_mul(1024, 10000);
  console.log(`SIMD GF multiply:       ${simdMulMbps.toFixed(1)} MB/s`);

  // DNA mapping
  const dnaMbps = wasmExports.bench_dna_mapping(10240, 1000);
  console.log(`DNA mapping 10KB:       ${dnaMbps.toFixed(1)} MB/s`);

  // Larger batch
  const largeBatch = wasmExports.bench_batch_rs_encode_simd(256, 40, 32);
  console.log(`Batch RS x256 (SIMD128): ${largeBatch.toFixed(1)} MB/s`);

  console.log("\n=== Speedup Analysis ===");
  console.log(`SIMD batch vs single:   ${(simdBatchMbps / singleMbps).toFixed(1)}x`);
  console.log(`SIMD batch vs scalar:   ${(simdBatchMbps / batchMbps).toFixed(1)}x`);
  console.log(`SIMD batch vs GF ceil:  ${((simdMulMbps / simdBatchMbps)).toFixed(2)}x of raw GF throughput`);
  console.log(`\nGap to 3.4 GB/s ceiling: ${((simdMulMbps - simdBatchMbps) / simdMulMbps * 100).toFixed(0)}% (interleave overhead + XOR)`);

  console.log("\n=== Final Comparison ===");
  console.log("| Metric                | Pure-JS   | Rust/WASM SIMD128 | bioarc  |");
  console.log("|-----------------------|-----------|-------------------|---------|");
  console.log(`| RS(40,32) single      | 0.5 MB/s  | ${singleMbps.toFixed(0).padStart(17)} MB/s | 120     |`);
  console.log(`| RS(40,32) batch x64   | N/A       | ${simdBatchMbps.toFixed(0).padStart(17)} MB/s | N/A     |`);
  console.log(`| RS(40,32) batch x256  | N/A       | ${largeBatch.toFixed(0).padStart(17)} MB/s | N/A     |`);
  console.log(`| GF multiply (raw)     | N/A       | ${simdMulMbps.toFixed(0).padStart(17)} MB/s | ~1,000  |`);
  console.log(`| DNA mapping           | 18 MB/s   | ${dnaMbps.toFixed(0).padStart(17)} MB/s | 200     |`);
}

main().catch(console.error);
