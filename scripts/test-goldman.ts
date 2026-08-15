// Test the existing Goldman mapping implementation.
import { bytesToGoldmanDna, goldmanDnaToBytes, hasHomopolymer } from "../src/lib/dna/goldman";
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
  console.log("=== Goldman Mapping Test ===\n");

  // Test round-trip
  let pass = 0, fail = 0;
  for (let trial = 0; trial < 100; trial++) {
    const data = randomBytes(40, trial);
    const dna = bytesToGoldmanDna(data);
    const recovered = goldmanDnaToBytes(dna);
    const ok = recovered.every((b, i) => b === data[i]);
    if (ok) pass++; else fail++;
    if (trial < 3) {
      console.log(`Trial ${trial}: ${data.length} bytes → ${dna.length} nt, gc=${gcContent(dna).toFixed(2)}, maxHP=${maxHomopolymerRun(dna)}, hasHP=${hasHomopolymer(dna)}`);
      console.log(`  DNA: ${dna.slice(0, 40)}...`);
    }
  }
  console.log(`\nRound-trip: ${pass}/100 pass\n`);

  // Check density
  const data = randomBytes(1000, 42);
  const dna = bytesToGoldmanDna(data);
  const density = (data.length * 8) / dna.length;
  console.log(`Density: ${density.toFixed(3)} bits/nt (theoretical log2(3)/2 = ${(Math.log2(3)/2).toFixed(3)})`);
  console.log(`Homopolymer-free: ${!hasHomopolymer(dna) ? "✅ YES" : "❌ NO"}`);
  console.log(`Max homopolymer run: ${maxHomopolymerRun(dna)}`);
  console.log(`GC content: ${gcContent(dna).toFixed(3)}`);
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
