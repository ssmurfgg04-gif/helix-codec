/**
 * v52 Convolutional Inner Code Unit Test
 *
 * Verifies that ConvolutionalInnerCode encode → decode is the identity
 * for an arbitrary byte input (zero-noise BSC channel).
 */

import { ConvolutionalInnerCode } from "../src/lib/dna/convolutional";

const REPORT = (msg: string) => console.log(`[v52-conv-unit] ${msg}`);

function testRoundTrip(N: number) {
  const conv = new ConvolutionalInnerCode(N);
  // Random input
  const input = new Uint8Array(N);
  for (let i = 0; i < N; i++) input[i] = (i * 37 + 13) & 0xff;

  const encoded = conv.encode(input);
  const decoded = conv.decode(encoded);

  let errors = 0;
  for (let i = 0; i < N; i++) {
    if (input[i] !== decoded[i]) errors++;
  }
  REPORT(`N=${N}: encoded=${encoded.length}B, decoded=${decoded.length}B, errors=${errors}/${N}`);
  if (errors === 0) {
    REPORT(`  ✅ identity round-trip`);
  } else {
    REPORT(`  ❌ round-trip FAILED`);
    // Print first few mismatches
    let shown = 0;
    for (let i = 0; i < N && shown < 5; i++) {
      if (input[i] !== decoded[i]) {
        REPORT(`    byte[${i}]: expected=${input[i]}, got=${decoded[i]}`);
        shown++;
      }
    }
  }
  return errors === 0;
}

REPORT("=== ConvolutionalInnerCode Unit Test ===");
let allPass = true;
for (const N of [1, 4, 8, 16, 24, 28, 29, 30, 32, 50]) {
  if (!testRoundTrip(N)) allPass = false;
}
REPORT(allPass ? "\n=== ALL PASSED ===" : "\n=== FAILED ===");
process.exit(allPass ? 0 : 1);
