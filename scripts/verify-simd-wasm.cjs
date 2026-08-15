#!/usr/bin/env node
/**
 * SIMD WASM Verification Test
 *
 * Loads the real SIMD WASM module, packs ACGTACGT as 2-bit,
 * unpacks with SIMD, verifies output, and benchmarks vs scalar.
 */

const path = require('path');

async function main() {
  console.log('=== SIMD WASM Verification Test ===\n');

  // ── Step 1: Load the WASM module ──
  console.log('[1] Loading SIMD WASM module...');

  const wasmJsPath = path.join(__dirname, '..', 'src', 'lib', 'dna', 'pkg', 'simd-wasm', 'simd_dna_unpack.js');
  const createModule = require(wasmJsPath);

  let Module;
  try {
    Module = await createModule({
      locateFile: (filename, scriptDir) => {
        return scriptDir + filename;
      },
    });
    Module._init_lut();
    console.log('    ✓ WASM module loaded successfully');
    console.log('    Available exports: _init_lut, _unpack_simd, _unpack_scalar, _unpack_simd_interleaved, _malloc, _free');
  } catch (err) {
    console.error('    ✗ Failed to load WASM:', err.message);
    process.exit(1);
  }

  // ── Step 2: Pack ACGTACGT as 2-bit ──
  // A=00, C=01, G=10, T=11
  // ACGT = 00_01_10_11 = 0x1B
  // ACGTACGT = 0x1B, 0x1B
  console.log('\n[2] Packing ACGTACGT as 2-bit...');
  const packed = new Uint8Array([0x1B, 0x1B]);
  console.log(`    packed = [0x${packed[0].toString(16).toUpperCase().padStart(2,'0')}, 0x${packed[1].toString(16).toUpperCase().padStart(2,'0')}]`);
  console.log(`    expected: [0x1B, 0x1B]`);

  // ── Step 3: Unpack with SIMD (interleaved) ──
  console.log('\n[3] Unpacking with SIMD (_unpack_simd_interleaved)...');
  const numNucleotides = 8;
  const inPtr = Module._malloc(packed.length);
  const outPtr = Module._malloc(packed.length * 4); // 8 bytes output

  try {
    Module.HEAPU8.set(packed, inPtr);
    Module._unpack_simd_interleaved(inPtr, outPtr, packed.length);

    const result = Module.HEAPU8.slice(outPtr, outPtr + numNucleotides);
    const resultStr = Array.from(result).map(b => String.fromCharCode(b)).join('');

    console.log(`    raw bytes: [${Array.from(result).map(b => '0x' + b.toString(16).toUpperCase()).join(', ')}]`);
    console.log(`    decoded:   "${resultStr}"`);
    console.log(`    expected:  "ACGTACGT"`);

    if (resultStr === 'ACGTACGT') {
      console.log('    ✓ SIMD interleaved output CORRECT');
    } else {
      console.error('    ✗ SIMD interleaved output MISMATCH');
      process.exit(1);
    }
  } finally {
    Module._free(inPtr);
    Module._free(outPtr);
  }

  // ── Step 4: Unpack with SIMD (bulk, 16 bytes to trigger SIMD path) ──
  console.log('\n[4] Unpacking with SIMD (_unpack_simd, 16-byte block)...');
  // Use 16 bytes of 0x1B (= ACGT repeated) to trigger the SIMD path
  const bulkPacked = new Uint8Array(16).fill(0x1B);
  const bulkNuc = 64; // 16 bytes * 4 nuc/byte
  const bulkExpected = 'ACGT'.repeat(16);
  const inPtr2 = Module._malloc(bulkPacked.length);
  const outPtr2 = Module._malloc(bulkPacked.length * 4);

  try {
    Module.HEAPU8.set(bulkPacked, inPtr2);
    Module._unpack_simd(inPtr2, outPtr2, bulkPacked.length);

    // _unpack_simd with 16 bytes uses SIMD path:
    // stores as [pos0 all 16][pos1 all 16][pos2 all 16][pos3 all 16]
    const raw = Module.HEAPU8.slice(outPtr2, outPtr2 + bulkPacked.length * 4);

    // Interleave to get sequential order
    const numBytes = bulkPacked.length;
    const interleaved = new Uint8Array(bulkNuc);
    let outIdx = 0;
    const numGroups = Math.floor(numBytes / 16);
    for (let g = 0; g < numGroups && outIdx < bulkNuc; g++) {
      const base = g * 64;
      for (let j = 0; j < 16 && outIdx < bulkNuc; j++) {
        interleaved[outIdx++] = raw[base + j];       // pos0
        if (outIdx >= bulkNuc) break;
        interleaved[outIdx++] = raw[base + 16 + j];  // pos1
        if (outIdx >= bulkNuc) break;
        interleaved[outIdx++] = raw[base + 32 + j];  // pos2
        if (outIdx >= bulkNuc) break;
        interleaved[outIdx++] = raw[base + 48 + j];  // pos3
      }
    }

    const bulkStr = Array.from(interleaved).map(b => String.fromCharCode(b)).join('');
    console.log(`    decoded:  "${bulkStr.substring(0, 32)}..."`);
    console.log(`    expected: "${bulkExpected.substring(0, 32)}..."`);

    if (bulkStr === bulkExpected) {
      console.log('    ✓ SIMD bulk output CORRECT');
    } else {
      console.error('    ✗ SIMD bulk output MISMATCH');
      process.exit(1);
    }
  } finally {
    Module._free(inPtr2);
    Module._free(outPtr2);
  }

  // ── Step 5: Unpack with scalar ──
  console.log('\n[5] Unpacking with scalar (_unpack_scalar)...');
  const inPtr3 = Module._malloc(packed.length);
  const outPtr3 = Module._malloc(packed.length * 4);

  try {
    Module.HEAPU8.set(packed, inPtr3);
    Module._unpack_scalar(inPtr3, outPtr3, packed.length);

    const scalarResult = Module.HEAPU8.slice(outPtr3, outPtr3 + numNucleotides);
    const scalarStr = Array.from(scalarResult).map(b => String.fromCharCode(b)).join('');
    console.log(`    decoded: "${scalarStr}"`);

    if (scalarStr === 'ACGTACGT') {
      console.log('    ✓ Scalar output CORRECT');
    } else {
      console.error('    ✗ Scalar output MISMATCH');
      process.exit(1);
    }
  } finally {
    Module._free(inPtr3);
    Module._free(outPtr3);
  }

  // ── Step 6: Benchmark SIMD vs Scalar ──
  console.log('\n[6] Benchmark: SIMD vs Scalar...');

  // Generate a large test array (1MB of packed data = 4M nucleotides)
  const benchSize = 1 << 20; // 1 MB of packed bytes
  const benchPacked = new Uint8Array(benchSize);
  for (let i = 0; i < benchSize; i++) {
    benchPacked[i] = (Math.random() * 256) | 0;
  }
  const benchNumNuc = benchSize * 4;

  // Warmup
  const warmIn = Module._malloc(benchSize);
  const warmOut = Module._malloc(benchSize * 4);
  Module.HEAPU8.set(benchPacked, warmIn);
  Module._unpack_simd_interleaved(warmIn, warmOut, benchSize);
  Module._unpack_scalar(warmIn, warmOut, benchSize);
  Module._free(warmIn);
  Module._free(warmOut);

  // Benchmark SIMD interleaved
  const simdIters = 20;
  const bIn1 = Module._malloc(benchSize);
  const bOut1 = Module._malloc(benchSize * 4);
  Module.HEAPU8.set(benchPacked, bIn1);

  const tSimdStart = performance.now();
  for (let i = 0; i < simdIters; i++) {
    Module._unpack_simd_interleaved(bIn1, bOut1, benchSize);
  }
  const tSimdEnd = performance.now();
  Module._free(bIn1);
  Module._free(bOut1);

  // Benchmark scalar
  const bIn2 = Module._malloc(benchSize);
  const bOut2 = Module._malloc(benchSize * 4);
  Module.HEAPU8.set(benchPacked, bIn2);

  const tScalarStart = performance.now();
  for (let i = 0; i < simdIters; i++) {
    Module._unpack_scalar(bIn2, bOut2, benchSize);
  }
  const tScalarEnd = performance.now();
  Module._free(bIn2);
  Module._free(bOut2);

  const simdMs = tSimdEnd - tSimdStart;
  const scalarMs = tScalarEnd - tScalarStart;
  const simdThroughput = (benchSize * simdIters / simdMs / 1000).toFixed(1);  // MB/ms = GB/s
  const scalarThroughput = (benchSize * simdIters / scalarMs / 1000).toFixed(1);
  const speedup = (scalarMs / simdMs).toFixed(2);

  console.log(`    SIMD interleaved: ${simdMs.toFixed(1)} ms (${simdThroughput} GB/s)`);
  console.log(`    Scalar:           ${scalarMs.toFixed(1)} ms (${scalarThroughput} GB/s)`);
  console.log(`    Speedup:          ${speedup}x`);

  // ── Step 7: Larger correctness test ──
  console.log('\n[7] Extended correctness test (256 patterns)...');
  let passCount = 0;
  const UNPACK_LUT = [0x41, 0x43, 0x47, 0x54]; // A, C, G, T
  for (let b0 = 0; b0 < 256; b0++) {
    const testPacked = new Uint8Array([b0]);
    const expected = String.fromCharCode(
      UNPACK_LUT[(b0 >> 6) & 3],
      UNPACK_LUT[(b0 >> 4) & 3],
      UNPACK_LUT[(b0 >> 2) & 3],
      UNPACK_LUT[b0 & 3]
    );

    const tIn = Module._malloc(1);
    const tOut = Module._malloc(4);
    Module.HEAPU8.set(testPacked, tIn);
    Module._unpack_simd_interleaved(tIn, tOut, 1);
    const actual = String.fromCharCode(
      Module.HEAPU8[tOut], Module.HEAPU8[tOut+1],
      Module.HEAPU8[tOut+2], Module.HEAPU8[tOut+3]
    );
    Module._free(tIn);
    Module._free(tOut);

    if (actual === expected) {
      passCount++;
    } else {
      console.error(`    ✗ Mismatch at 0x${b0.toString(16).padStart(2,'0')}: expected "${expected}", got "${actual}"`);
    }
  }
  console.log(`    ${passCount}/256 patterns passed`);
  if (passCount !== 256) {
    process.exit(1);
  }

  // ── Final ──
  console.log('\n' + '='.repeat(40));
  console.log('SIMD WASM: REAL ✓');
  console.log('='.repeat(40));
}

main().catch(err => {
  console.error('Fatal error:', err);
  process.exit(1);
});
