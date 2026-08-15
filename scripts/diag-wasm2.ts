// Test WASM helper functions to see if WASM is fundamentally working
const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");

// Test 1: bytes_to_dna / dna_to_bytes roundtrip
const testBytes = new Uint8Array([0x1b, 0x4b, 0x24, 0x6d, 0x00, 0xff]);
console.log("input bytes:", Array.from(testBytes).map(b=>b.toString(16).padStart(2,'0')).join(' '));

// Try test_arithmetic_decode (a simple function to verify WASM works)
try {
  // simple DNA: ACGT repeated 16 times = 16 bytes of 0x1B
  const dna = "ACGT".repeat(16);
  const dnaBytes = new Uint8Array(dna.length);
  for (let i=0; i<dna.length; i++) dnaBytes[i] = dna.charCodeAt(i);
  const out = wasm.test_arithmetic_decode(dnaBytes, 3, 16);
  console.log("test_arithmetic_decode result:", Array.from(out).map(b=>b.toString(16).padStart(2,'0')).join(' '));
} catch (e) {
  console.log("test_arithmetic_decode err:", e.message);
}

// Try test_rs216_decode
try {
  // simple RS test: encode then decode
  const nsym = 4;
  const n = 8, k = 4;
  // Make a codeword where the first 4 are data
  // For simplicity, just feed zeros and see what happens
  const cw = new Uint8Array(n*2);
  for (let i=0; i<n*2; i++) cw[i] = 0;
  // erased = positions 0, 1 (each u32 LE)
  const erased = new Uint8Array(8);
  erased[0]=0; erased[4]=1;
  const out = wasm.test_rs216_decode(cw, erased, n, k);
  console.log("test_rs216_decode (all-zero) result:", Array.from(out).map(b=>b.toString(16).padStart(2,'0')).join(' '));
} catch (e) {
  console.log("test_rs216_decode err:", e.message);
}

// Look at what full_decode returns when given simple input
try {
  // Build a trivial decode: 1 oligo, no reads at all
  const allReads = new Uint8Array(0);
  const readOffsets = new Uint8Array(0);
  const readLengths = new Uint8Array(0);
  const fwdPrimer = new Uint8Array(0);
  const revPrimer = new Uint8Array(0);
  const out = wasm.full_decode(allReads, readOffsets, readLengths, fwdPrimer, revPrimer, 0, 0, 0, 0, 0, 0, 0, 0, false);
  console.log("full_decode (empty) returned len:", out.length);
} catch (e) {
  console.log("full_decode empty err:", e.message);
}
