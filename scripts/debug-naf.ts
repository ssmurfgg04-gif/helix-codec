import {
  compressWithNAF,
  decompressWithNAF,
  compressWithAGC,
  decompressWithAGC,
} from '../src/lib/dna/dna-compress-real';

// Test NAF
console.log("=== NAF test ===");
{
  const dna = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  const dnaBytes = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) dnaBytes[i] = dna.charCodeAt(i);

  console.log("DNA length:", dnaBytes.length);

  try {
    const compressed = compressWithNAF(dnaBytes, 6);
    console.log("Compressed:", compressed.length, "bytes");
    console.log("Header bytes:", Array.from(compressed.slice(0, 20)).map(b => b.toString(16).padStart(2, '0')).join(' '));
    
    // Parse header manually
    const v = new DataView(compressed.buffer, compressed.byteOffset);
    const seqLen = v.getUint32(4, true);
    const rleLen = v.getUint32(8, true);
    const modelLen = v.getUint32(12, true);
    console.log("seqLen:", seqLen, "rleLen:", rleLen, "modelLen:", modelLen);
    console.log("Expected data after header:", 16 + modelLen, "to", compressed.length);
    
    const decompressed = decompressWithNAF(compressed);
    let match = true;
    for (let i = 0; i < dnaBytes.length; i++) {
      if (dnaBytes[i] !== decompressed[i]) { match = false; console.log(`Mismatch at ${i}: expected ${dnaBytes[i]} got ${decompressed[i]}`); break; }
    }
    console.log("Roundtrip:", match);
  } catch (e: any) {
    console.error("Error:", e.message);
    console.error(e.stack?.split('\n').slice(0, 5).join('\n'));
  }
}

// Test AGC
console.log("\n=== AGC test ===");
{
  const dna = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  const dnaBytes = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) dnaBytes[i] = dna.charCodeAt(i);

  try {
    const compressed = compressWithAGC(dnaBytes, 6);
    console.log("Compressed:", compressed.length, "bytes");
    const decompressed = decompressWithAGC(compressed);
    let match = true;
    for (let i = 0; i < dnaBytes.length; i++) {
      if (dnaBytes[i] !== decompressed[i]) { match = false; console.log(`AGC Mismatch at ${i}`); break; }
    }
    console.log("AGC Roundtrip:", match);
  } catch (e: any) {
    console.error("AGC Error:", e.message);
  }
}
