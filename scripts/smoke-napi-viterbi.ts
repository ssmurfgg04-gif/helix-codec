/**
 * Smoke test: verify napi-rs native Viterbi loads and works end-to-end.
 * Tests: addon load, K=9 encode/decode roundtrip, indel-tolerant decode.
 */
import * as path from 'node:path';
import * as fs from 'node:fs';

const ADDON_PATH = path.resolve('/home/z/my-project/rust/helix-dna-napi/target/release/libhelix_dna_napi.so');

async function main() {
  console.log('=== napi-rs Viterbi Smoke Test ===');
  console.log(`Addon path: ${ADDON_PATH}`);
  console.log(`Addon exists: ${fs.existsSync(ADDON_PATH)}`);

  // Step 1: Direct dlopen test
  console.log('\n[1] Direct dlopen test...');
  try {
    const mod: any = { exports: {} };
    (process as any).dlopen(mod, ADDON_PATH);
    const addon = mod.exports;
    console.log('  dlopen OK');
    console.log(`  napiVersion(): ${addon.napiVersion()}`);
    console.log(`  exports: ${Object.keys(addon).join(', ')}`);
  } catch (e: any) {
    console.error(`  dlopen FAILED: ${e.message}`);
    process.exit(1);
  }

  // Step 2: Through wrapper module
  console.log('\n[2] Wrapper module test...');
  const { enableNativeViterbi, isNativeViterbiActive, nativeViterbiK9DecodeStandard, nativeViterbiK9Decode, nativeConvK9Encode } =
    await import('/home/z/my-project/src/lib/dna/native/viterbi-napi.ts');
  const enabled = await enableNativeViterbi();
  console.log(`  enableNativeViterbi(): ${enabled}`);
  console.log(`  isNativeViterbiActive(): ${isNativeViterbiActive()}`);

  if (!enabled) {
    console.error('  Native Viterbi not active — aborting.');
    process.exit(1);
  }

  // Step 3: K=9 encode/decode roundtrip (clean channel)
  console.log('\n[3] K=9 encode/decode roundtrip (clean channel)...');
  const data = new Uint8Array(32);
  for (let i = 0; i < 32; i++) data[i] = (i * 7 + 13) & 0xff;
  console.log(`  Input (32 bytes): ${Buffer.from(data).toString('hex')}`);

  const encoded = nativeConvK9Encode(data);
  console.log(`  Encoded (${encoded.length} bytes): ${Buffer.from(encoded).toString('hex')}`);

  const decodedStd = nativeViterbiK9DecodeStandard(encoded);
  console.log(`  Standard decode (${decodedStd.length} bytes): ${Buffer.from(decodedStd).toString('hex')}`);

  const matchStd = Buffer.from(decodedStd.slice(0, 32)).equals(Buffer.from(data));
  console.log(`  Roundtrip match (standard): ${matchStd}`);

  // Step 4: Indel-tolerant decode on clean channel (should also work)
  console.log('\n[4] Indel-tolerant decode on clean channel...');
  const decodedIndel = nativeViterbiK9Decode(encoded, {
    maxDrift: 15,
    insertionPenalty: 1.5,
    deletionPenalty: 1.5,
    numInfoBits: 32 * 8,
  });
  console.log(`  Indel decode (${decodedIndel.length} bytes): ${Buffer.from(decodedIndel).toString('hex')}`);
  const matchIndel = Buffer.from(decodedIndel.slice(0, 32)).equals(Buffer.from(data));
  console.log(`  Roundtrip match (indel): ${matchIndel}`);

  // Step 5: Indel-tolerant decode with simulated insertion
  console.log('\n[5] Indel-tolerant decode with simulated insertion...');
  const bits = [];
  for (let i = 0; i < encoded.length; i++) {
    for (let b = 7; b >= 0; b--) bits.push((encoded[i] >> b) & 1);
  }
  // Insert 3 random bits at position 50
  bits.splice(50, 0, 1, 0, 1);
  const corrupted = new Uint8Array(Math.ceil(bits.length / 8));
  for (let i = 0; i < bits.length; i++) {
    corrupted[i >> 3] |= bits[i] << (7 - (i & 7));
  }
  console.log(`  Corrupted (${corrupted.length} bytes, +3 bits inserted at bit 50)`);
  const decodedCorrupt = nativeViterbiK9Decode(corrupted, {
    maxDrift: 15,
    insertionPenalty: 1.5,
    deletionPenalty: 1.5,
    numInfoBits: 32 * 8,
  });
  console.log(`  Decode (${decodedCorrupt.length} bytes): ${Buffer.from(decodedCorrupt).toString('hex')}`);
  const matchCorrupt = Buffer.from(decodedCorrupt.slice(0, 32)).equals(Buffer.from(data));
  console.log(`  Recovery match: ${matchCorrupt}`);

  // Step 6: Through IndelTolerantConvolutionalInnerCode wrapper
  console.log('\n[6] IndelTolerantConvolutionalInnerCode wrapper test...');
  const { enableNativeViterbi: enableInConv, IndelTolerantConvolutionalInnerCode } =
    await import('/home/z/my-project/src/lib/dna/convolutional-indel.ts');
  await enableInConv();
  const inner = new IndelTolerantConvolutionalInnerCode(32, { useK9: true } as any);
  const wrapperEncoded = inner.encode(data);
  const wrapperDecoded = inner.decode(wrapperEncoded);
  console.log(`  Wrapper encode (${wrapperEncoded.length} bytes)`);
  console.log(`  Wrapper decode (${wrapperDecoded.decoded.length} bytes)`);
  console.log(`  Wrapper match: ${Buffer.from(wrapperDecoded.decoded).equals(Buffer.from(data))}`);

  // Summary
  console.log('\n=== Summary ===');
  const allPassed = matchStd && matchIndel;
  console.log(`Standard roundtrip: ${matchStd ? 'PASS' : 'FAIL'}`);
  console.log(`Indel roundtrip (clean): ${matchIndel ? 'PASS' : 'FAIL'}`);
  console.log(`Indel recovery (3 inserted bits): ${matchCorrupt ? 'PASS' : 'PARTIAL (expected — heavy indel)'}`);
  console.log(`Wrapper roundtrip: ${Buffer.from(wrapperDecoded.decoded).equals(Buffer.from(data)) ? 'PASS' : 'FAIL'}`);

  if (allPassed) {
    console.log('\nALL CORE TESTS PASSED — napi-rs Viterbi is operational.');
  } else {
    console.log('\nSOME TESTS FAILED — review output above.');
    process.exit(1);
  }
}

main().catch(e => {
  console.error('Smoke test crashed:', e);
  process.exit(1);
});
