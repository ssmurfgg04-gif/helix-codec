/**
 * Quick unit test for constrained mapping roundtrip with GC codebooks.
 */
import { bytesToConstrainedDna, constrainedDnaToBytesWithErasure, bytesToSplitConstrainedDna, splitConstrainedDnaToBytesWithErasure } from "../src/lib/dna/constrained-mapping";

// Test 1: Simple constrained mapping roundtrip
console.log("=== Test 1: bytesToConstrainedDna roundtrip ===");
for (let len = 1; len <= 20; len++) {
  const data = new Uint8Array(len);
  for (let i = 0; i < len; i++) data[i] = (i * 37 + 42) & 0xff;

  const encoded = bytesToConstrainedDna(data, 3, 0.4, 0.6);
  const decoded = constrainedDnaToBytesWithErasure(encoded.dna, 3, len, encoded.codebookSequence);

  let match = decoded.data.length === data.length;
  if (match) {
    for (let i = 0; i < data.length; i++) {
      if (decoded.data[i] !== data[i]) { match = false; break; }
    }
  }
  if (!match) {
    console.log(`FAIL: len=${len}, encoded=${encoded.dna.slice(0, 40)}...`);
    console.log(`  expected: [${Array.from(data).join(",")}]`);
    console.log(`  got:      [${Array.from(decoded.data).join(",")}]`);
    console.log(`  erasures: ${decoded.erasures.filter(e => e).length}`);
    process.exit(1);
  }
}
console.log("PASS: All lengths 1-20");

// Test 2: Split constrained mapping roundtrip
console.log("\n=== Test 2: bytesToSplitConstrainedDna roundtrip ===");
for (let len = 5; len <= 30; len++) {
  const data = new Uint8Array(len);
  for (let i = 0; i < len; i++) data[i] = (i * 37 + 42) & 0xff;

  const encoded = bytesToSplitConstrainedDna(data, 3, 4, 0.4, 0.6);
  const decoded = splitConstrainedDnaToBytesWithErasure(encoded.dna, 3, 4, len, encoded.codebookSequence);

  let match = decoded.data.length === data.length;
  if (match) {
    for (let i = 0; i < data.length; i++) {
      if (decoded.data[i] !== data[i]) { match = false; break; }
    }
  }
  if (!match) {
    console.log(`FAIL: len=${len}`);
    console.log(`  expected: [${Array.from(data).join(",")}]`);
    console.log(`  got:      [${Array.from(decoded.data).join(",")}]`);
    process.exit(1);
  }
}
console.log("PASS: All lengths 5-30");

// Test 3: Split constrained WITHOUT codebookSequence (decoder reconstructs)
console.log("\n=== Test 3: Split constrained roundtrip WITHOUT codebookSequence ===");
let fail3 = 0;
for (let len = 5; len <= 30; len++) {
  const data = new Uint8Array(len);
  for (let i = 0; i < len; i++) data[i] = (i * 37 + 42) & 0xff;

  const encoded = bytesToSplitConstrainedDna(data, 3, 4, 0.4, 0.6);
  // Decode WITHOUT codebookSequence — decoder must reconstruct from position
  const decoded = splitConstrainedDnaToBytesWithErasure(encoded.dna, 3, 4, len);

  let match = decoded.data.length === data.length;
  if (match) {
    for (let i = 0; i < data.length; i++) {
      if (decoded.data[i] !== data[i]) { match = false; break; }
    }
  }
  if (!match) {
    fail3++;
    if (fail3 <= 3) {
      console.log(`FAIL: len=${len}`);
      console.log(`  expected: [${Array.from(data).join(",")}]`);
      console.log(`  got:      [${Array.from(decoded.data).join(",")}]`);
    }
  }
}
if (fail3 > 0) {
  console.log(`FAIL: ${fail3} lengths failed without codebookSequence`);
} else {
  console.log("PASS: All lengths 5-30 (without codebookSequence)");
}

// Test 4: Check GC distribution of encoded DNA
console.log("\n=== Test 4: GC distribution of encoded DNA ===");
const bigData = new Uint8Array(1000);
for (let i = 0; i < 1000; i++) bigData[i] = (i * 7 + 13) & 0xff;
const bigEncoded = bytesToSplitConstrainedDna(bigData, 3, 4, 0.4, 0.6);
let gc = 0;
for (const c of bigEncoded.dna) {
  if (c === 'G' || c === 'C') gc++;
}
const gcFrac = gc / bigEncoded.dna.length;
console.log(`GC fraction: ${gcFrac.toFixed(4)} (target: 0.40-0.60)`);

// Check max homopolymer
let maxRun = 1, curRun = 1;
for (let i = 1; i < bigEncoded.dna.length; i++) {
  if (bigEncoded.dna[i] === bigEncoded.dna[i-1]) {
    curRun++;
    if (curRun > maxRun) maxRun = curRun;
  } else {
    curRun = 1;
  }
}
console.log(`Max homopolymer: ${maxRun} (target: <= 3)`);

console.log("\nDone.");
