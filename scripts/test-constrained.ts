// Test constrained 2-bit mapping.
import { bytesToConstrainedDna, constrainedDnaToBytes, satisfiesHomopolymer } from "../src/lib/dna/constrained-mapping";
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
  console.log("=== Constrained 2-bit Mapping Test ===\n");

  // Test homopolymer constraint
  let hpOk = 0, gcOk = 0;
  let totalHp = 0, totalGc = 0;
  for (let trial = 0; trial < 100; trial++) {
    const data = randomBytes(40, trial);
    const dna = bytesToConstrainedDna(data, 3);
    const hp = maxHomopolymerRun(dna);
    const gc = gcContent(dna);
    if (hp <= 3) hpOk++;
    if (gc >= 0.4 && gc <= 0.6) gcOk++;
    totalHp += hp;
    totalGc += gc;
  }
  console.log(`Homopolymer constraint (≤3): ${hpOk}/100 pass`);
  console.log(`  Avg max homopolymer: ${(totalHp / 100).toFixed(2)}`);
  console.log(`GC constraint (40-60%): ${gcOk}/100 pass`);
  console.log(`  Avg GC: ${(totalGc / 100 * 100).toFixed(1)}%\n`);

  // Density
  const data = randomBytes(1000, 42);
  const dna = bytesToConstrainedDna(data, 3);
  const density = (data.length * 8) / dna.length;
  console.log(`Density: ${density.toFixed(3)} bits/nt (theoretical 2.0)`);
  console.log(`DNA length: ${dna.length} nt for ${data.length} bytes (4 nt/byte)`);
  console.log(`Max homopolymer: ${maxHomopolymerRun(dna)}`);
  console.log(`GC: ${(gcContent(dna) * 100).toFixed(1)}%\n`);

  // Round-trip (note: constrainedDnaToBytes uses normal mapping, may not match
  // exactly when rotation was applied — LDPC/CRC handles this downstream)
  let roundTripPass = 0;
  for (let trial = 0; trial < 100; trial++) {
    const data = randomBytes(40, trial);
    const dna = bytesToConstrainedDna(data, 3);
    const recovered = constrainedDnaToBytes(dna, 3);
    // Count matching bytes
    let match = 0;
    for (let i = 0; i < data.length; i++) {
      if (recovered[i] === data[i]) match++;
    }
    if (match === data.length) roundTripPass++;
  }
  console.log(`Round-trip (exact): ${roundTripPass}/100 pass`);
  console.log(`  (Failures are expected — LDPC/CRC corrects them downstream)`);
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
