#!/usr/bin/env node
/**
 * Verification test for zstd WASM module.
 *
 * Loads the real zstd WASM, compresses test data, verifies the output
 * starts with the zstd magic number (0x28 0xB5 0x2F 0xFD), then
 * decompresses and verifies roundtrip integrity.
 *
 * Prints "ZSTD WASM: REAL ✓" on success.
 */

'use strict';

const path = require('path');

async function main() {
  console.log('=== zstd WASM Verification Test ===\n');

  // Step 1: Load the index.node.js entry point
  const pkgDir = path.resolve(__dirname, '../src/lib/dna/pkg/zstd-wasm');
  const wasmPkgDir = path.resolve(__dirname, '../src/lib/dna/wasm-pkg/zstd-wasm');

  let mod;
  let loadDir;

  try {
    mod = require(path.join(pkgDir, 'index.node.js'));
    loadDir = 'pkg/zstd-wasm';
    console.log(`[1] Loaded module from: ${loadDir}`);
  } catch (e) {
    try {
      mod = require(path.join(wasmPkgDir, 'index.node.js'));
      loadDir = 'wasm-pkg/zstd-wasm';
      console.log(`[1] Loaded module from: ${loadDir}`);
    } catch (e2) {
      console.error('[1] FAIL: Could not load zstd-wasm from either pkg/ or wasm-pkg/');
      console.error('   pkg error:', e.message);
      console.error('   wasm-pkg error:', e2.message);
      process.exit(1);
    }
  }

  // Step 2: Verify module has expected exports
  if (typeof mod.init !== 'function') {
    console.error('[2] FAIL: Module does not export init()');
    process.exit(1);
  }
  if (typeof mod.compress !== 'function') {
    console.error('[2] FAIL: Module does not export compress()');
    process.exit(1);
  }
  if (typeof mod.decompress !== 'function') {
    console.error('[2] FAIL: Module does not export decompress()');
    process.exit(1);
  }
  console.log('[2] Module exports: init, compress, decompress ✓');

  // Step 3: Initialize the WASM module
  console.log('[3] Calling init()...');
  try {
    await mod.init();
    console.log('[3] init() completed ✓');
  } catch (e) {
    console.error('[3] FAIL: init() threw:', e.message);
    process.exit(1);
  }

  // Step 4: Compress test data
  const testData = new Uint8Array(
    Array.from('Hello, zstd WASM! This is a verification test for real zstd compression.'.split(''),
      c => c.charCodeAt(0))
  );
  console.log(`[4] Input data: ${testData.length} bytes`);

  let compressed;
  try {
    compressed = mod.compress(testData, 3);
    console.log(`[4] Compressed: ${compressed.length} bytes (ratio: ${(testData.length / compressed.length).toFixed(2)}x)`);
  } catch (e) {
    console.error('[4] FAIL: compress() threw:', e.message);
    process.exit(1);
  }

  // Step 5: Verify zstd magic number (0x28 0xB5 0x2F 0xFD)
  const ZSTD_MAGIC = [0x28, 0xB5, 0x2F, 0xFD];
  const header = Array.from(compressed.slice(0, 4));
  const magicMatch = header[0] === ZSTD_MAGIC[0] &&
                     header[1] === ZSTD_MAGIC[1] &&
                     header[2] === ZSTD_MAGIC[2] &&
                     header[3] === ZSTD_MAGIC[3];

  if (magicMatch) {
    console.log(`[5] Magic bytes: 0x${header.map(b => b.toString(16).padStart(2, '0')).join(' ')} — zstd magic ✓`);
  } else {
    console.error(`[5] FAIL: Expected magic 0x28 b5 2f fd, got 0x${header.map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    process.exit(1);
  }

  // Step 6: Decompress and verify roundtrip
  let decompressed;
  try {
    decompressed = mod.decompress(compressed);
    console.log(`[6] Decompressed: ${decompressed.length} bytes`);
  } catch (e) {
    console.error('[6] FAIL: decompress() threw:', e.message);
    process.exit(1);
  }

  // Step 7: Verify data integrity
  let match = decompressed.length === testData.length;
  if (match) {
    for (let i = 0; i < testData.length; i++) {
      if (decompressed[i] !== testData[i]) {
        match = false;
        console.error(`[7] FAIL: Mismatch at byte ${i}: expected ${testData[i]}, got ${decompressed[i]}`);
        process.exit(1);
      }
    }
  } else {
    console.error(`[7] FAIL: Length mismatch: expected ${testData.length}, got ${decompressed.length}`);
    process.exit(1);
  }
  console.log('[7] Roundtrip integrity: original === decompressed ✓');

  // Step 8: Test with larger data (1KB)
  const largeData = new Uint8Array(1024);
  for (let i = 0; i < 1024; i++) largeData[i] = i % 256;
  const largeCompressed = mod.compress(largeData, 3);
  const largeDecompressed = mod.decompress(largeCompressed);
  const largeHeader = Array.from(largeCompressed.slice(0, 4));
  const largeMagicOk = largeHeader[0] === 0x28 && largeHeader[1] === 0xB5 && largeHeader[2] === 0x2F && largeHeader[3] === 0xFD;
  const largeRoundtripOk = largeDecompressed.length === largeData.length &&
    largeDecompressed.every((b, i) => b === largeData[i]);

  console.log(`[8] 1KB test: magic=${largeMagicOk ? '✓' : '✗'}, roundtrip=${largeRoundtripOk ? '✓' : '✗'}, ratio=${(largeData.length / largeCompressed.length).toFixed(2)}x`);

  if (!largeMagicOk || !largeRoundtripOk) {
    console.error('[8] FAIL: 1KB test failed');
    process.exit(1);
  }

  // Step 9: Test different compression levels
  for (const level of [1, 3, 9, 19]) {
    const c = mod.compress(testData, level);
    const d = mod.decompress(c);
    const h = [c[0], c[1], c[2], c[3]];
    const magicOk = h[0] === 0x28 && h[1] === 0xB5 && h[2] === 0x2F && h[3] === 0xFD;
    const roundtripOk = d.length === testData.length && d.every((b, i) => b === testData[i]);
    console.log(`[9] Level ${level.toString().padStart(2)}: ${c.length} bytes, magic=${magicOk ? '✓' : '✗'}, roundtrip=${roundtripOk ? '✓' : '✗'}`);
    if (!magicOk || !roundtripOk) {
      console.error(`[9] FAIL: Level ${level} test failed`);
      process.exit(1);
    }
  }

  // Final
  console.log('\n=== All checks passed ===');
  console.log('ZSTD WASM: REAL ✓');
}

main().catch(e => {
  console.error('Unhandled error:', e);
  process.exit(1);
});
