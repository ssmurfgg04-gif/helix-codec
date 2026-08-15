// Quick LDPC codec round-trip test.
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { crc16Bytes, verifyCrc16 } from "../src/lib/dna/crc16";

function randomBytes(n: number, seed: number): Uint8Array {
  let s = seed >>> 0;
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    s ^= s << 13; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
    out[i] = s & 0xff;
  }
  return out;
}

function injectErrors(bits: Uint8Array, numErrors: number, seed: number): Uint8Array {
  // Inject bit errors at random positions
  let s = seed >>> 0;
  const positions = new Set<number>();
  while (positions.size < numErrors) {
    s ^= s << 13; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
    positions.add(s % bits.length);
  }
  const out = bits.slice();
  for (const p of positions) out[p] ^= 1;
  return out;
}

function bytesToBits(data: Uint8Array): Uint8Array {
  const bits = new Uint8Array(data.length * 8);
  for (let i = 0; i < data.length; i++) {
    for (let bit = 0; bit < 8; bit++) {
      bits[i * 8 + bit] = (data[i] >> (7 - bit)) & 1;
    }
  }
  return bits;
}

function bitsToBytes(bits: Uint8Array): Uint8Array {
  const out = new Uint8Array(bits.length / 8);
  for (let i = 0; i < out.length; i++) {
    let b = 0;
    for (let bit = 0; bit < 8; bit++) b |= bits[i * 8 + bit] << (7 - bit);
    out[i] = b;
  }
  return out;
}

async function main() {
  console.log("=== LDPC Inner Code Round-Trip Test ===\n");

  // Test with the default Helix parameters:
  // 200 nt oligo → 40 bytes inner → 4B addr + 30B payload + 4B parity + 2B CRC
  // LDPC codeword (excluding CRC): 38 bytes = 304 bits
  // LDPC info: 34 bytes = 272 bits (addr + payload)
  // LDPC parity: 4 bytes = 32 bits
  const k = 34; // 4B addr + 30B payload
  const n = 38; // 34 + 4 parity
  const ldpc = new LDPCInnerCode({ n, k });

  console.log(`LDPC(${n * 8} bits, ${k * 8} bits, ${(n - k) * 8} bits parity)`);
  console.log(`  Rate: ${(k / n).toFixed(4)} (theoretical ${(k * 8 / (n * 8)).toFixed(4)})`);
  console.log(`  Theoretical density (k/n): ${(k * 8 / (n - 4) / 2).toFixed(3)} bits/nt (excluding CRC bytes)`);
  console.log();

  // Test 1: No errors
  console.log("--- Test 1: No errors ---");
  let pass = 0, fail = 0;
  for (let trial = 0; trial < 100; trial++) {
    const info = randomBytes(k, trial);
    const codeword = ldpc.encode(info);
    const decoded = ldpc.decode(codeword);
    const ok = decoded.data.every((b, i) => b === info[i]);
    if (ok) pass++; else fail++;
  }
  console.log(`  100 trials: ${pass} pass, ${fail} fail\n`);

  // Test 2: 1-bit errors
  console.log("--- Test 2: 1-bit errors ---");
  pass = 0; fail = 0;
  for (let trial = 0; trial < 100; trial++) {
    const info = randomBytes(k, trial);
    const codeword = ldpc.encode(info);
    const bits = bytesToBits(codeword);
    const noisyBits = injectErrors(bits, 1, trial + 1000);
    const noisyBytes = bitsToBytes(noisyBits);
    try {
      const decoded = ldpc.decode(noisyBytes);
      const ok = decoded.data.every((b, i) => b === info[i]);
      if (ok) pass++; else fail++;
    } catch {
      fail++;
    }
  }
  console.log(`  100 trials: ${pass} pass, ${fail} fail\n`);

  // Test 3: 2-bit errors
  console.log("--- Test 3: 2-bit errors ---");
  pass = 0; fail = 0;
  for (let trial = 0; trial < 100; trial++) {
    const info = randomBytes(k, trial);
    const codeword = ldpc.encode(info);
    const bits = bytesToBits(codeword);
    const noisyBits = injectErrors(bits, 2, trial + 2000);
    const noisyBytes = bitsToBytes(noisyBits);
    try {
      const decoded = ldpc.decode(noisyBytes);
      const ok = decoded.data.every((b, i) => b === info[i]);
      if (ok) pass++; else fail++;
    } catch {
      fail++;
    }
  }
  console.log(`  100 trials: ${pass} pass, ${fail} fail\n`);

  // Test 4: 4-bit errors
  console.log("--- Test 4: 4-bit errors ---");
  pass = 0; fail = 0;
  for (let trial = 0; trial < 100; trial++) {
    const info = randomBytes(k, trial);
    const codeword = ldpc.encode(info);
    const bits = bytesToBits(codeword);
    const noisyBits = injectErrors(bits, 4, trial + 4000);
    const noisyBytes = bitsToBytes(noisyBits);
    try {
      const decoded = ldpc.decode(noisyBytes);
      const ok = decoded.data.every((b, i) => b === info[i]);
      if (ok) pass++; else fail++;
    } catch {
      fail++;
    }
  }
  console.log(`  100 trials: ${pass} pass, ${fail} fail\n`);

  // Test 5: 8-bit errors (should mostly fail — only 16 correctable)
  console.log("--- Test 5: 8-bit errors ---");
  pass = 0; fail = 0;
  for (let trial = 0; trial < 100; trial++) {
    const info = randomBytes(k, trial);
    const codeword = ldpc.encode(info);
    const bits = bytesToBits(codeword);
    const noisyBits = injectErrors(bits, 8, trial + 8000);
    const noisyBytes = bitsToBytes(noisyBits);
    try {
      const decoded = ldpc.decode(noisyBytes);
      const ok = decoded.data.every((b, i) => b === info[i]);
      if (ok) pass++; else fail++;
    } catch {
      fail++;
    }
  }
  console.log(`  100 trials: ${pass} pass, ${fail} fail\n`);

  // Test 6: Timing — 1000 encodes + 1000 decodes
  console.log("--- Test 6: Performance ---");
  const N = 1000;
  const testInfo = randomBytes(k, 42);
  const t0 = Date.now();
  for (let i = 0; i < N; i++) {
    ldpc.encode(testInfo);
  }
  const encMs = Date.now() - t0;
  const codeword = ldpc.encode(testInfo);
  const t1 = Date.now();
  for (let i = 0; i < N; i++) {
    try { ldpc.decode(codeword); } catch {}
  }
  const decMs = Date.now() - t1;
  console.log(`  Encode: ${N} ops in ${encMs}ms (${(N / (encMs / 1000)).toFixed(0)} ops/sec)`);
  console.log(`  Decode: ${N} ops in ${decMs}ms (${(N / (decMs / 1000)).toFixed(0)} ops/sec)\n`);

  console.log("=== LDPC Test Complete ===");
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
