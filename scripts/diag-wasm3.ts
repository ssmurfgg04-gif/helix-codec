// Test WASM LDPC encode/decode roundtrip
const wasm = require("../src/lib/dna/wasm-pkg/helix_dna_wasm.js");

// Build LDPC code with same dimensions as full_decode uses
// innerK=58, innerN=62 (default config)
const innerK = 58;
const innerN = 62;
const code = new wasm.LdpcCode(innerN, innerK);

// Encode some test data
const data = new Uint8Array(innerK);
for (let i=0; i<innerK; i++) data[i] = (i * 31) & 0xff;
console.log("input data:", Array.from(data.slice(0,8)).map(b=>b.toString(16).padStart(2,'0')).join(' '), "...");

const encoded = code.encode(data);
console.log("encoded length:", encoded.length, "first 8:", Array.from(encoded.slice(0,8)).map(b=>b.toString(16).padStart(2,'0')).join(' '));

// Verify encode: data bytes should match input
const dataRoundtrip = encoded.slice(0, innerK);
let match = true;
for (let i=0; i<innerK; i++) {
  if (dataRoundtrip[i] !== data[i]) { match = false; console.log("MISMATCH at", i, data[i], dataRoundtrip[i]); break; }
}
console.log("encode systematic data preserved:", match);

// Decode (no errors) — should return same data
const decoded = code.decode(encoded);
console.log("decoded length:", decoded.length, "first 8:", Array.from(decoded.slice(0,8)).map(b=>b.toString(16).padStart(2,'0')).join(' '));

let decMatch = true;
for (let i=0; i<innerK; i++) {
  if (decoded[i] !== data[i]) { decMatch = false; console.log("DECODE MISMATCH at", i, "expected", data[i], "got", decoded[i]); break; }
}
console.log("decode roundtrip match:", decMatch);

// Now test with 1 bit flipped
const corrupted = new Uint8Array(innerN);
for (let i=0; i<innerN; i++) corrupted[i] = encoded[i];
corrupted[20] ^= 0x40; // flip one bit
const dec2 = code.decode(corrupted);
let dec2Match = true;
for (let i=0; i<innerK; i++) {
  if (dec2[i] !== data[i]) { dec2Match = false; console.log("DECODE2 MISMATCH at", i, "expected", data[i], "got", dec2[i]); break; }
}
console.log("decode with 1 bit flip match:", dec2Match);
