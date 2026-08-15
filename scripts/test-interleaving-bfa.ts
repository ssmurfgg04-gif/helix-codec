// Test interleaving + BFA decoder
import { interleaveCodewords, deinterleaveCodewords, maxErrorsPerCodeword, minDepthForBurst } from "../src/lib/dna/interleaving";
import { BFADecoder } from "../src/lib/dna/bfa-decoder";
import { makeLDPCInner } from "../src/lib/dna/ldpc-codec";

function testInterleaving() {
  console.log("=== Interleaving Test ===\n");

  // Test 1: Basic interleave/deinterleave round-trip
  const cw1 = new Uint8Array([0x01, 0x02, 0x03, 0x04, 0x05]);
  const cw2 = new Uint8Array([0x10, 0x20, 0x30, 0x40, 0x50]);
  const cw3 = new Uint8Array([0x11, 0x22, 0x33, 0x44, 0x55]);

  const interleaved = interleaveCodewords([cw1, cw2, cw3]);
  console.log(`Interleaved ${cw1.length * 3} bytes from 3 codewords`);
  console.log(`  Interleaved: ${Array.from(interleaved).map(b => b.toString(16).padStart(2, "0")).join(" ")}`);

  const deinterleaved = deinterleaveCodewords(interleaved, 3);
  let match = true;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 5; j++) {
      if (deinterleaved[i][j] !== [cw1, cw2, cw3][i][j]) match = false;
    }
  }
  console.log(`  Round-trip: ${match ? "PASS" : "FAIL"}\n`);

  // Test 2: Burst error analysis
  console.log("Burst error analysis:");
  console.log(`  Burst of 6 bytes, depth=3: max ${maxErrorsPerCodeword(6, 3)} errors/cw (LDPC corrects 4)`);
  console.log(`  Burst of 12 bytes, depth=3: max ${maxErrorsPerCodeword(12, 3)} errors/cw`);
  console.log(`  Burst of 8 bytes, depth=4: max ${maxErrorsPerCodeword(8, 4)} errors/cw`);
  console.log(`  Min depth for burst=8, LDPC cap=4: ${minDepthForBurst(8, 4)}`);
  console.log(`  Min depth for burst=12, LDPC cap=4: ${minDepthForBurst(12, 4)}\n`);

  // Test 3: Burst error spread
  console.log("Burst error spread simulation:");
  const depth = 4;
  const cwLen = 10;
  const codewords: Uint8Array[] = [];
  for (let i = 0; i < depth; i++) {
    const cw = new Uint8Array(cwLen);
    for (let j = 0; j < cwLen; j++) cw[j] = i * 16 + j;
    codewords.push(cw);
  }
  const interleaved2 = interleaveCodewords(codewords);

  // Inject a burst of 8 bytes at position 10
  for (let i = 10; i < 18; i++) {
    interleaved2[i] ^= 0xFF;
  }

  const deint2 = deinterleaveCodewords(interleaved2, depth);
  let errorsPerCw = [0, 0, 0, 0];
  for (let i = 0; i < depth; i++) {
    for (let j = 0; j < cwLen; j++) {
      if (deint2[i][j] !== codewords[i][j]) errorsPerCw[i]++;
    }
  }
  console.log(`  Burst of 8 bytes in interleaved stream → errors per codeword: [${errorsPerCw.join(", ")}]`);
  console.log(`  Max errors in any codeword: ${Math.max(...errorsPerCw)} (LDPC can correct 4)\n`);
}

function testBFA() {
  console.log("=== BFA Decoder Test ===\n");

  // Create LDPC code with standard config (4 parity bytes = 32 parity bits)
  const ldpc = makeLDPCInner(4, 50, 4); // 4 parity bytes, 50 payload, 4 address
  const bfa = new BFADecoder(ldpc, { maxFlipBits: 2, candidateCount: 15, maxIter: 30 });

  // Test data
  const data = new Uint8Array(ldpc.k);
  for (let i = 0; i < ldpc.k; i++) data[i] = (i * 37 + 11) & 0xff;
  const codeword = ldpc.encode(data);

  // Test 1: No errors
  try {
    const result = bfa.decode(codeword);
    let match = true;
    for (let i = 0; i < ldpc.k; i++) {
      if (data[i] !== result.data[i]) { match = false; break; }
    }
    console.log(`  No errors: ${match ? "PASS" : "FAIL"}`);
  } catch (e: any) {
    console.log(`  No errors: FAIL (${e.message})`);
  }

  // Test 2: 1-bit error
  let passCount = 0;
  let totalCount = 10;
  for (let trial = 0; trial < totalCount; trial++) {
    const recv = codeword.slice();
    const bytePos = trial % ldpc.n;
    const bitPos = trial % 8;
    recv[bytePos] ^= (1 << bitPos);
    try {
      const result = bfa.decode(recv);
      let match = true;
      for (let i = 0; i < ldpc.k; i++) {
        if (data[i] !== result.data[i]) { match = false; break; }
      }
      if (match) passCount++;
    } catch {}
  }
  console.log(`  1-bit error: ${passCount}/${totalCount} PASS`);

  // Test 3: 2-bit errors
  passCount = 0;
  totalCount = 10;
  for (let trial = 0; trial < totalCount; trial++) {
    const recv = codeword.slice();
    recv[trial % ldpc.n] ^= 0x01;
    recv[(trial + 5) % ldpc.n] ^= 0x80;
    try {
      const result = bfa.decode(recv);
      let match = true;
      for (let i = 0; i < ldpc.k; i++) {
        if (data[i] !== result.data[i]) { match = false; break; }
      }
      if (match) passCount++;
    } catch {}
  }
  console.log(`  2-bit errors: ${passCount}/${totalCount} PASS (BFA enhanced)`);

  console.log();
}

testInterleaving();
testBFA();
