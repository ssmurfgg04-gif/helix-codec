// Test outer RS (GF(2^16)) erasure recovery directly.
import { ReedSolomon216 } from "../src/lib/dna/reedsolomon216";

async function main() {
  console.log("=== Outer RS GF(2^16) Erasure Test ===\n");

  // Match the 10MB test config: n=844, k=649
  const n = 844, k = 649;
  const rs = new ReedSolomon216({ n, k });

  // Create a known data vector
  const data = new Uint16Array(k);
  for (let i = 0; i < k; i++) data[i] = (i * 31 + 17) & 0xffff;

  // Encode
  const parity = rs.parity(data);
  console.log(`Encoded: ${k} data + ${parity.length} parity = ${k + parity.length} total (n=${n})`);

  // Build full codeword
  const codeword = new Uint16Array(n);
  for (let i = 0; i < k; i++) codeword[i] = data[i];
  for (let i = 0; i < parity.length; i++) codeword[k + i] = parity[i];

  // Verify syndrome is zero
  // (We trust the encoder here)

  // Test 1: No erasures — decode should return data unchanged
  console.log("\n--- Test 1: No erasures ---");
  try {
    const r = rs.decodeWithErasures(codeword, []);
    let ok = true;
    for (let i = 0; i < k; i++) {
      if (r.data[i] !== data[i]) { ok = false; break; }
    }
    console.log(`  Result: ${ok ? "✅ PASS" : "❌ FAIL"}`);
  } catch (e) {
    console.log(`  Result: ❌ EXCEPTION: ${(e as Error).message}`);
  }

  // Test 2: 2 erasures (matching the 10MB test scenario)
  console.log("\n--- Test 2: 2 erasures at positions 203, 248 ---");
  const erased = [203, 248];
  const corrupted = codeword.slice();
  corrupted[203] = 0;
  corrupted[248] = 0;
  try {
    const r = rs.decodeWithErasures(corrupted, erased);
    let ok = true;
    let wrongIdx = -1;
    for (let i = 0; i < k; i++) {
      if (r.data[i] !== data[i]) { ok = false; wrongIdx = i; break; }
    }
    console.log(`  Result: ${ok ? "✅ PASS" : "❌ FAIL"}`);
    if (!ok) {
      console.log(`  First wrong at index ${wrongIdx}: expected ${data[wrongIdx]}, got ${r.data[wrongIdx]}`);
      console.log(`  Recovered[203] = ${r.data[203]} (expected ${data[203]})`);
      console.log(`  Recovered[248] = ${r.data[248]} (expected ${data[248]})`);
    }
  } catch (e) {
    console.log(`  Result: ❌ EXCEPTION: ${(e as Error).message}`);
  }

  // Test 3: 10 erasures
  console.log("\n--- Test 3: 10 erasures ---");
  const erased10 = [10, 20, 30, 100, 200, 300, 400, 500, 600, 648];
  const corrupted10 = codeword.slice();
  for (const e of erased10) corrupted10[e] = 0;
  try {
    const r = rs.decodeWithErasures(corrupted10, erased10);
    let ok = true;
    for (let i = 0; i < k; i++) {
      if (r.data[i] !== data[i]) { ok = false; break; }
    }
    console.log(`  Result: ${ok ? "✅ PASS" : "❌ FAIL"}`);
  } catch (e) {
    console.log(`  Result: ❌ EXCEPTION: ${(e as Error).message}`);
  }

  // Test 4: 100 erasures (stress test)
  console.log("\n--- Test 4: 100 erasures ---");
  const erased100: number[] = [];
  for (let i = 0; i < 100; i++) erased100.push(i * 7 % k);
  const corrupted100 = codeword.slice();
  for (const e of erased100) corrupted100[e] = 0;
  try {
    const r = rs.decodeWithErasures(corrupted100, erased100);
    let ok = true;
    for (let i = 0; i < k; i++) {
      if (r.data[i] !== data[i]) { ok = false; break; }
    }
    console.log(`  Result: ${ok ? "✅ PASS" : "❌ FAIL"}`);
  } catch (e) {
    console.log(`  Result: ❌ EXCEPTION: ${(e as Error).message}`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
