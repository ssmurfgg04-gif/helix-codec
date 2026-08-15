const fs = require('fs');
const path = require('path');

async function test() {
  // Read the JS glue code
  const jsCode = fs.readFileSync(path.join(__dirname, 'simd_dna_unpack_mod.js'), 'utf-8');
  
  // The file defines createSimdDnaUnpackModule as a global variable
  // We need to eval it in the right context
  const wasmBinary = fs.readFileSync(path.join(__dirname, 'simd_dna_unpack_mod.wasm'));
  
  // Use Function constructor to create the factory
  const factory = new Function('module', 'require', '__filename', '__dirname', 
    jsCode + '\nreturn createSimdDnaUnpackModule;');
  const createModule = factory({ exports: {} }, require, __filename, __dirname);
  
  console.log('createModule type:', typeof createModule);
  
  const Module = await createModule({ wasmBinary });
  Module._init_lut();
  
  const packed = new Uint8Array([0x1B, 0x1B, 0x1B, 0x1B]);
  const numBytes = packed.length;
  const outSize = numBytes * 4;
  const inPtr = Module._malloc(numBytes);
  const outPtr = Module._malloc(outSize);
  Module.HEAPU8.set(packed, inPtr);
  
  Module._unpack_simd_interleaved(inPtr, outPtr, numBytes);
  const result = Array.from(Module.HEAPU8.slice(outPtr, outPtr + outSize));
  const decoded = String.fromCharCode(...result);
  console.log('SIMD interleaved:', decoded);
  console.log('Match:', decoded === 'ACGTACGTACGTACGT');
  
  Module._unpack_scalar(inPtr, outPtr, numBytes);
  const result2 = Array.from(Module.HEAPU8.slice(outPtr, outPtr + outSize));
  const decoded2 = String.fromCharCode(...result2);
  console.log('Scalar:', decoded2, 'Match:', decoded2 === 'ACGTACGTACGTACGT');
  
  Module._free(inPtr);
  Module._free(outPtr);
  
  // Benchmark
  const bigPacked = new Uint8Array(1000000);
  bigPacked.fill(0x1B);
  const bigInPtr = Module._malloc(bigPacked.length);
  const bigOutPtr = Module._malloc(bigPacked.length * 4);
  Module.HEAPU8.set(bigPacked, bigInPtr);
  
  const t1 = Date.now();
  for (let i = 0; i < 10; i++) Module._unpack_simd_interleaved(bigInPtr, bigOutPtr, bigPacked.length);
  const t2 = Date.now();
  for (let i = 0; i < 10; i++) Module._unpack_scalar(bigInPtr, bigOutPtr, bigPacked.length);
  const t3 = Date.now();
  
  console.log('SIMD interleaved: 10x 1M bytes in', t2-t1, 'ms');
  console.log('Scalar: 10x 1M bytes in', t3-t2, 'ms');
  if (t2-t1 > 0) console.log('Speedup:', ((t3-t2)/(t2-t1)).toFixed(2) + 'x');
  
  Module._free(bigInPtr);
  Module._free(bigOutPtr);
}

test().catch(console.error);
