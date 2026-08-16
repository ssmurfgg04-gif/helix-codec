/**
 * Debug WASM DNA Compressors — check round-trip with small example
 */

import * as path from 'node:path';

async function debugRoundTrip() {
  const wasmPath = path.join(process.cwd(), 'wasm-src', 'dna_compressors.js');
  const createModule = await import(wasmPath);
  const mod = await createModule.default();

  // Test with very simple DNA: "ACGTACGT"
  const testDna = "-ACGTACGT-"; // we'll use just ACGTACGT
  const input = Buffer.from("ACGTACGT", 'ascii');
  
  console.log(`Input: "${input.toString('ascii')}" (${input.length} bytes)`);
  console.log(`Input bytes: [${Array.from(input).join(', ')}]`);

  for (let algo = 0; algo < 3; algo++) {
    const namePtr = mod._dna_compressor_name(algo);
    const name = mod.UTF8ToString(namePtr);
    console.log(`\n--- ${name} ---`);

    // Compress
    const inputLen = input.length;
    const outputCap = inputLen * 2 + 1024;
    const inputPtr = mod._malloc(inputLen);
    const outputPtr = mod._malloc(outputCap);
    mod.HEAPU8.set(input, inputPtr);

    const compressedLen = mod._dna_compress(algo, inputPtr, inputLen, outputPtr, outputCap);
    console.log(`Compressed: ${compressedLen} bytes`);

    if (compressedLen < 0) {
      console.log('Compression failed!');
      mod._free(inputPtr);
      mod._free(outputPtr);
      continue;
    }

    const compressed = Buffer.from(mod.HEAPU8.slice(outputPtr, outputPtr + compressedLen));
    console.log(`Compressed bytes: [${Array.from(compressed.slice(0, Math.min(32, compressedLen))).map(b => b.toString(16).padStart(2,'0')).join(', ')}...]`);
    mod._free(inputPtr);
    mod._free(outputPtr);

    // Decompress
    const decompressCap = inputLen * 2 + 1024;
    const compPtr = mod._malloc(compressedLen);
    const decompPtr = mod._malloc(decompressCap);
    mod.HEAPU8.set(compressed, compPtr);

    const decompressedLen = mod._dna_decompress(algo, compPtr, compressedLen, decompPtr, decompressCap);
    console.log(`Decompressed: ${decompressedLen} bytes`);

    if (decompressedLen < 0) {
      console.log('Decompression failed!');
      mod._free(compPtr);
      mod._free(decompPtr);
      continue;
    }

    const decompressed = Buffer.from(mod.HEAPU8.slice(decompPtr, decompPtr + decompressedLen));
    console.log(`Decompressed bytes: [${Array.from(decompressed2).join(', ')}]`);
    console.log(`Decompressed string: "${decompressed2.toString('ascii')}"`);
    console.log(`Round-trip: ${input.toString('ascii') === decompressed2.toString('ascii') ? 'OK' : 'FAIL'}`);
    
    // Show where they differ
    if (input.toString('ascii') !== decompressed2.toString('ascii')) {
      for (let i = 0; i < Math.max(input.length, decompressedLen); i++) {
        const a = i < input.length ? input[i] : -1;
        const b = i < decompressedLen ? decompressed2[i] : -1;
        if (a !== b) {
          console.log(`  DiffGpos ${i}: orig=${a} (${String.fromCharCode(a)}) decomp=${b} (${b >= 0 ? String.fromCharCode(b) : 'N/A'})`);
        }
      }
    }

    mod._free(compPtr);
    mod._free(decompPtr);
  }
}

debugRoundTrip().catch(err => {
  console.error('Debug failed:', err);
  process.exit(1);
});
