/**
 * Test WASM DNA Compressors — Round-trip validation against synthetic DNA data.
 *
 * Loads the Emscripten-compiled C++ compressors, runs compress→decompress
 * round-trips on various DNA sequences, and measures compression ratios.
 *
 * Usage:
 *   npx tsx scripts/test-dna-compressors.ts
 */

import * as fs from 'node:fs';
import * as path from 'node:path';

// ---------------------------------------------------------------------------
// DNA test data generators
// ---------------------------------------------------------------------------

function generateRandomDna(length: number, rng: () => number): string {
  const bases = 'ACGT';
  let dna = '';
  for (let i = 0; i < length; i++) {
    dna += bases[Math.floor(rng() * 4)];
  }
  return dna;
}

function generateRepeatDna(unitLength: number, repeats: number, rng: () => number): string {
  const unit = generateRandomDna(unitLength, rng);
  return unit.repeat(repeats);
}

function generateHumanLikeDna(length: number, rng: () => number): string {
  // Simulate human genome: 41% GC, with some low-complexity regions
  // GC bias: P(G) = P(C) = 0.205, P(A) = P(T) = 0.145
  const cumFreq = [0.145, 0.35, 0.555, 1.0];  // A, C, G, T
  const bases = 'ACGT';
  let dna = '';
  for (let i = 0; i < length; i++) {
    const r = rng();
    let idx = 0;
    while (idx < 3 && r > cumFreq[idx]) idx++;
    dna += bases[idx];
  }
  return dna;
}

// ---------------------------------------------------------------------------
// Simple PRNG
// ---------------------------------------------------------------------------

class Rng {
  private state: number;
  constructor(seed: number) { this.state = (seed >>> 0) || 1; }
  next(): number {
    this.state ^= this.state << 13;
    this.state ^= this.state >>> 17;
    this.state ^= this.state << 5;
    this.state = this.state >>> 0;
    return this.state / 0x100000000;
  }
}

// ---------------------------------------------------------------------------
// Test using the WASM module directly (without the TS bridge)
// ---------------------------------------------------------------------------

async function testWithWasmModule() {
  const wasmPath = path.join(process.cwd(), 'wasm-src', 'dna_compressors.js');
  const wasmBinaryPath = path.join(process.cwd(), 'wasm-src', 'dna_compressors.wasm');

  if (!fs.existsSync(wasmPath) || !fs.existsSync(wasmBinaryPath)) {
    console.error('WASM module not found. Run scripts/build-compressors-wasm.sh first.');
    process.exit(1);
  }

  console.log('=== WASM DNA Compressor Test ===\n');

  // Load the WASM module
  const createModule = await import(wasmPath);
  const mod = await createModule.default();

  const count = mod._dna_compressor_count();
  console.log(`Available compressors: ${count}`);
  for (let i = 0; i < count; i++) {
    const namePtr = mod._dna_compressor_name(i);
    const name = mod.UTF8ToString(namePtr);
    console.log(`  [${i}] ${name}`);
  }
  console.log('');

  // Generate test data
  const rng = new Rng(42);
  const testSequences: { name: string; dna: string }[] = [
    { name: 'Random 1KB', dna: generateRandomDna(1000, () => rng.next()) },
    { name: 'Random 10KB', dna: generateRandomDna(10000, () => rng.next()) },
    { name: 'Random 100KB', dna: generateRandomDna(100000, () => rng.next()) },
    { name: 'Repeat 100×10', dna: generateRepeatDna(100, 10, () => rng.next()) },
    { name: 'Human-like 10KB', dna: generateHumanLikeDna(10000, () => rng.next()) },
  ];

  // Also load E. coli genome if available
  const ecoliPath = path.join(process.cwd(), 'test-data', 'ecoli.fasta');
  if (fs.existsSync(ecoliPath)) {
    const fasta = fs.readFileSync(ecoliPath, 'utf-8');
    const lines = fasta.split('\n').filter(l => !l.startsWith('>'));
    const dna = lines.join('');
    if (dna.length > 0) {
      testSequences.push({ name: `E.coli ${Math.round(dna.length / 1000)}KB`, dna: dna.substring(0, 100000) });
    }
  }

  // Test each compressor on each sequence
  console.log('┌──────────────┬────────────────┬──────────┬──────────┬───────────┬────────────┐');
  console.log('│ Compressor   │ Sequence       │ Orig (B) │ Comp (B) │ Ratio     │ Round-trip │');
  console.log('├──────────────┼────────────────┼──────────┼──────────┼───────────┼────────────┤');

  for (let algo = 0; algo < count; algo++) {
    const namePtr = mod._dna_compressor_name(algo);
    const name = mod.UTF8ToString(namePtr);

    for (const { name: seqName, dna } of testSequences) {
      const inputBytes = Buffer.from(dna, 'ascii');
      const inputLen = inputBytes.length;
      const outputCap = inputLen + inputLen / 2 + 1024;

      // Compress
      const inputPtr = mod._malloc(inputLen);
      const outputPtr = mod._malloc(outputCap);
      mod.HEAPU8.set(inputBytes, inputPtr);

      const compressedLen = mod._dna_compress(algo, inputPtr, inputLen, outputPtr, outputCap);
      const compressed = mod.HEAPU8.slice(outputPtr, outputPtr + compressedLen);
      mod._free(inputPtr);
      mod._free(outputPtr);

      if (compressedLen < 0) {
        console.log(`│ ${name.padEnd(12)} │ ${seqName.padEnd(14)} │ ERROR    │          │           │            │`);
        continue;
      }

      // Decompress
      const decompressCap = inputLen + 64;
      const compPtr = mod._malloc(compressedLen);
      const decompPtr = mod._malloc(decompressCap);
      mod.HEAPU8.set(compressed, compPtr);

      const decompressedLen = mod._dna_decompress(algo, compPtr, compressedLen, decompPtr, decompressCap);
      const decompressed = mod.HEAPU8.slice(decompPtr, decompPtr + decompressedLen);
      mod._free(compPtr);
      mod._free(decompPtr);

      // Check round-trip
      const originalStr = inputBytes.toString('ascii');
      const decompressedStr = Buffer.from(decompressed).toString('ascii');
      const roundTripOk = originalStr === decompressedStr;

      const ratio = compressedLen / inputLen;

      console.log(
        `│ ${name.padEnd(12)} │ ${seqName.padEnd(14)} │ ` +
        `${String(inputLen).padStart(8)} │ ` +
        `${String(compressedLen).padStart(8)} │ ` +
        `${ratio.toFixed(3).padStart(9)} │ ` +
        `${(roundTripOk ? 'OK' : 'FAIL').padStart(10)} │`
      );
    }
  }

  console.log('└──────────────┴────────────────┴──────────┴──────────┴───────────┴────────────┘');
  console.log('\n=== WASM DNA Compressor Test Complete ===');
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

testWithWasmModule().catch(err => {
  console.error('Test failed:', err);
  process.exit(1);
});
