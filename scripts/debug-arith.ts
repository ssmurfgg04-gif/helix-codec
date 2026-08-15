import {
  ArithmeticEncoder,
  ArithmeticDecoder,
  AdaptiveFrequencyModel,
  adaptiveArithEncode,
  adaptiveArithDecode,
} from '../src/lib/dna/arithmetic-coder';

// Test 1: Simple 4-symbol roundtrip
console.log("=== 4-symbolE-symbol test ===");
const syms = new Uint8Array([0, 1, 2, 3, 0, 0, 0, 1, 2, 3]);
const comp = adaptiveArithEncode(syms, 4, 0);
const dec = adaptiveArithDecode(comp, 4, 10, 0);
console.log("Match:", JSON.stringify(Array.from(syms)) === JSON.stringify(Array.from(dec)));

// Test 2: 256-symbol roundtrip
console.log("\n=== 256-symbol test ===");
{
  const bytes = new Uint8Array(100);
  for (let i = 0; i < 100; i++) bytes[i] = Math.floor(Math.random() * 256);
  const comp2 = adaptiveArithEncode(bytes, 256, 0);
  const dec2 = adaptiveArithDecode(comp2, 256, 100, 0);
  let match = true;
  for (let i = 0; i < 100; i++) { if (bytes[i] !== dec2[i]) { match = false; break; } }
  console.log("100-byte 256-symbol roundtrip:", match);
  if (!match) {
    console.log("First 20 input:", Array.from(bytes.slice(0, 20)));
    console.log("First 20 decoded:", Array.from(dec2.slice(0, 20)));
  }
}

// Test 3: NAF
console.log("\n=== NAF test ===");
import {
  compressWithNAF,
  decompressWithNAF,
} from '../src/lib/dna/dna-compress-real';

{
  const dna = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  const dnaBytes = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) dnaBytes[i] = dna.charCodeAt(i);

  try {
    const compressed = compressWithNAF(dnaBytes, 6);
    console.log("Compressed:", compressed.length, "bytes, magic:", compressed[0].toString(16), compressed[1].toString(16));
    const decompressed = decompressWithNAF(compressed);
    let match = true;
    for (let i = 0; i < dnaBytes.length; i++) {
      if (dnaBytes[i] !== decompressed[i]) { match = false; break; }
    }
    console.log("Roundtrip:", match);
  } catch (e: any) {
    console.error("Error:", e.message);
  }
}
