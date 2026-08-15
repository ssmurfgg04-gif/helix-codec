// Test fixed-rate constrained mapping (2.0 bits/nt, homopolymer-free, with erasure).
import { bytesToConstrainedDna, constrainedDnaToBytesWithErasure, satisfiesHomopolymer } from "../src/lib/dna/constrained-mapping";
import { gcContent, maxHomopolymerRun } from "../src/lib/dna/mapping";

function randomBytes(n: number, seed: number): Uint8Array {
  let s = seed >>> 0;
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    s ^= s << 13; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
    out[i] = s & 0xff;
  }
  return out;
}

async function main() {
  console.log("=== Fixed-Rate Constrained Mapping Test ===\n");

  // Test 1: Homopolymer constraint
  let hpOk = 0, gcOk = 0, lengthOk = 0;
  let totalHp = 0, totalGc = 0, totalErasures = 0;
  for (let trial = 0; trial < 100; trial++) {
    const data = randomBytes(40, trial);
    const dna = bytesToConstrainedDna(data, 3);
    const hp = maxHomopolymerRun(dna);
    const gc = gcContent(dna);
    if (hp <= 3) hpOk++;
    if (gc >= 0.4 && gc <= 0.6) gcOk++;
    if (dna.length === data.length * 4) lengthOk++;
    totalHp += hp;
    totalGc += gc;

    // Check erasures
    const { erasures } = constrainedDnaToBytesWithErasure(dna, 3, data.length);
    totalErasures += erasures.filter(e => e).length;
  }
  console.log(`Homopolymer constraint (≤3): ${hpOk}/100 pass`);
  console.log(`  Avg max homopolymer: ${(totalHp / 100).toFixed(2)}`);
  console.log(`GC constraint (40-60%): ${gcOk}/100 pass`);
  console.log(`  Avg GC: ${(totalGc / 100 * 100).toFixed(1)}%`);
  console.log(`Fixed length (4:1): ${lengthOk}/100 pass`);
  console.log(`Avg erasures per 40-byte block: ${(totalErasures / 100).toFixed(1)} bits (${(totalErasures / 100 / 320 * 100).toFixed(1)}% erasure rate)\n`);

  // Test 2: Density
  const data = randomBytes(1000, 42);
  const dna = bytesToConstrainedDna(data, 3);
  const density = (data.length * 8) / dna.length;
  console.log(`Density: ${density.toFixed(3)} bits/nt (theoretical 2.0)`);
  console.log(`DNA length: ${dna.length} nt for ${data.length} bytes (4 nt/byte ✓)`);
  console.log(`Max homopolymer: ${maxHomopolymerRun(dna)}`);
  console.log(`GC: ${(gcContent(dna) * 100).toFixed(1)}%\n`);

  // Test 3: Round-trip with erasure info
  let fullPass = 0, partialPass = 0;
  for (let trial = 0; trial < 100; trial++) {
    const data = randomBytes(40, trial);
    const dna = bytesToConstrainedDna(data, 3);
    const { data: recovered, erasures } = constrainedDnaToBytesWithErasure(dna, 3, data.length);

    // Check how many bits match (ignoring erased bits)
    let totalBits = 0, correctBits = 0, erasedBits = 0;
    for (let i = 0; i < data.length; i++) {
      for (let b = 7; b >= 0; b--) {
        const origBit = (data[i] >> b) & 1;
        const recBit = (recovered[i] >> b) & 1;
        const bitIdx = i * 8 + (7 - b);
        totalBits++;
        if (erasures[bitIdx]) {
          erasedBits++;
        } else if (origBit === recBit) {
          correctBits++;
        }
      }
    }
    const nonErasedCorrect = correctBits / (totalBits - erasedBits);
    if (nonErasedCorrect === 1.0) fullPass++;
    if (nonErasedCorrect >= 0.99) partialPass++;
  }
  console.log(`Round-trip (non-erased bits match): ${fullPass}/100 fully correct, ${partialPass}/100 ≥99% correct`);
  console.log(`  (Erased bits need BP decoder to recover)\n`);

  // Test 4: Worst case — all same base (maximum homopolymer pressure)
  console.log("--- Worst case: all 0xFF bytes ---");
  const worstData = new Uint8Array(40).fill(0xFF);
  const worstDna = bytesToConstrainedDna(worstData, 3);
  console.log(`  DNA: ${worstDna.slice(0, 40)}...`);
  console.log(`  Max homopolymer: ${maxHomopolymerRun(worstDna)}`);
  console.log(`  GC: ${(gcContent(worstDna) * 100).toFixed(1)}%`);
  console.log(`  Length: ${worstDna.length} (expected ${40 * 4})`);
  const { erasures: worstErasures } = constrainedDnaToBytesWithErasure(worstDna, 3, 40);
  console.log(`  Erasures: ${worstErasures.filter(e => e).length} / ${40 * 8} bits (${(worstErasures.filter(e => e).length / (40 * 8) * 100).toFixed(1)}%)`);
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
