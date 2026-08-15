// Test modified-SRT constrained coding.
import { bytesToSrtDna, srtDnaToBytes, satisfiesHomopolymer } from "../src/lib/dna/srt-constrained";
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
  console.log("=== Modified-SRT Constrained Coding Test ===\n");

  // Test 1: Homopolymer constraint
  let hpPass = 0, gcPass = 0;
  let totalHp = 0, totalGc = 0;
  for (let trial = 0; trial < 100; trial++) {
    const data = randomBytes(40, trial);
    const dna = bytesToSrtDna(data, 3, 160); // 40 bytes → 160 nt target
    const hp = maxHomopolymerRun(dna);
    const gc = gcContent(dna);
    if (hp <= 3) hpPass++;
    if (gc >= 0.4 && gc <= 0.6) gcPass++;
    totalHp += hp;
    totalGc += gc;
  }
  console.log(`Homopolymer constraint (≤3): ${hpPass}/100 pass`);
  console.log(`  Avg max homopolymer: ${(totalHp / 100).toFixed(2)}`);
  console.log(`GC constraint (40-60%): ${gcPass}/100 pass`);
  console.log(`  Avg GC: ${(totalGc / 100 * 100).toFixed(1)}%\n`);

  // Test 2: Density
  const data = randomBytes(1000, 42);
  const dna = bytesToSrtDna(data, 3, 4000); // Target 4000 nt (same as direct)
  const density = (data.length * 8) / dna.length;
  console.log(`Density: ${density.toFixed(3)} bits/nt (theoretical 2.0, ℓ=3 capacity 1.982)`);
  console.log(`DNA length: ${dna.length} nt for ${data.length} bytes`);
  console.log(`Max homopolymer: ${maxHomopolymerRun(dna)}`);
  console.log(`GC: ${(gcContent(dna) * 100).toFixed(1)}%\n`);

  // Test 3: Round-trip
  let roundTripPass = 0;
  let partialPass = 0;
  for (let trial = 0; trial < 100; trial++) {
    const data = randomBytes(40, trial);
    const dna = bytesToSrtDna(data, 3, 160);
    const recovered = srtDnaToBytes(dna, 3, 40);

    let matching = 0;
    for (let i = 0; i < Math.min(data.length, recovered.length); i++) {
      if (recovered[i] === data[i]) matching++;
    }
    if (matching === data.length) roundTripPass++;
    if (matching >= data.length * 0.95) partialPass++;
  }
  console.log(`Round-trip (exact): ${roundTripPass}/100 pass`);
  console.log(`Round-trip (≥95%): ${partialPass}/100 pass\n`);

  // Test 4: Worst case (all same byte → maximum homopolymer pressure)
  console.log("--- Worst case: all 0xFF bytes ---");
  const worstData = new Uint8Array(40).fill(0xFF);
  const worstDna = bytesToSrtDna(worstData, 3, 160);
  console.log(`  Max homopolymer: ${maxHomopolymerRun(worstDna)}`);
  console.log(`  GC: ${(gcContent(worstDna) * 100).toFixed(1)}%`);
  console.log(`  Length: ${worstDna.length}`);

  // Test 5: All zeros (maximum homopolymer)
  console.log("\n--- All 0x00 bytes ---");
  const zeroData = new Uint8Array(40).fill(0x00);
  const zeroDna = bytesToSrtDna(zeroData, 3, 160);
  console.log(`  Max homopolymer: ${maxHomopolymerRun(zeroDna)}`);
  console.log(`  GC: ${(gcContent(zeroDna) * 100).toFixed(1)}%`);
  console.log(`  Length: ${zeroDna.length}`);
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
