/**
 * Quick K=9 penalty diagnostic: test ONE config to see if it works at all.
 */
import { enableNativeViterbi, nativeConvK9Encode, nativeViterbiK9Decode } from '../src/lib/dna/native/viterbi-napi';

function bytesToBits(b: Uint8Array): number[] {
  const out: number[] = [];
  for (let i = 0; i < b.length; i++) for (let k = 7; k >= 0; k--) out.push((b[i] >> k) & 1);
  return out;
}
function bitsToBytes(bits: number[]): Uint8Array {
  const n = Math.floor(bits.length / 8);
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    let v = 0;
    for (let k = 7; k >= 0; k--) v |= (bits[i * 8 + (7 - k)] & 1) << k;
    out[i] = v;
  }
  return out;
}

async function main() {
  console.log('=== Quick K=9 diagnostic ===');
  await enableNativeViterbi();

  const p = new Uint8Array(32);
  for (let i = 0; i < 32; i++) p[i] = (i * 7 + 13) & 0xff;
  console.log(`Payload: ${Buffer.from(p).toString('hex')}`);

  const enc = nativeConvK9Encode(p);
  console.log(`Encoded: ${enc.length} bytes`);

  // Test 1: clean
  console.log('\n[1] Clean decode...');
  const t0 = Date.now();
  const dec1 = nativeViterbiK9Decode(enc, { maxDrift: 15, insertionPenalty: 1.5, deletionPenalty: 1.5, numInfoBits: 256 });
  console.log(`  Decoded in ${Date.now() - t0}ms, len=${dec1.length}, match=${Buffer.from(dec1.slice(0, 32)).equals(Buffer.from(p))}`);

  // Test 2: 1 ins + 1 del
  console.log('\n[2] 1 insertion + 1 deletion...');
  const bits = bytesToBits(enc);
  bits.splice(50, 0, 1);
  bits.splice(100, 1);
  const c2 = bitsToBytes(bits);
  const t1 = Date.now();
  try {
    const dec2 = nativeViterbiK9Decode(c2, { maxDrift: 15, insertionPenalty: 1.5, deletionPenalty: 1.5, numInfoBits: 256 });
    console.log(`  Decoded in ${Date.now() - t1}ms, len=${dec2.length}, match=${Buffer.from(dec2.slice(0, 32)).equals(Buffer.from(p))}`);
  } catch (e: any) {
    console.log(`  FAILED in ${Date.now() - t1}ms: ${e.message}`);
  }

  // Test 3: 3 ins + 2 del
  console.log('\n[3] 3 insertions + 2 deletions...');
  const bits3 = bytesToBits(enc);
  bits3.splice(40, 0, 1, 0, 1);
  bits3.splice(120, 1);
  bits3.splice(180, 1);
  const c3 = bitsToBytes(bits3);
  const t2 = Date.now();
  try {
    const dec3 = nativeViterbiK9Decode(c3, { maxDrift: 15, insertionPenalty: 1.5, deletionPenalty: 1.5, numInfoBits: 256 });
    console.log(`  Decoded in ${Date.now() - t2}ms, len=${dec3.length}, match=${Buffer.from(dec3.slice(0, 32)).equals(Buffer.from(p))}`);
  } catch (e: any) {
    console.log(`  FAILED in ${Date.now() - t2}ms: ${e.message}`);
  }

  // Test 4: large drift=20, 3i+2d
  console.log('\n[4] Same scenario with drift=20...');
  const t3 = Date.now();
  try {
    const dec4 = nativeViterbiK9Decode(c3, { maxDrift: 20, insertionPenalty: 1.5, deletionPenalty: 1.5, numInfoBits: 256 });
    console.log(`  Decoded in ${Date.now() - t3}ms, len=${dec4.length}, match=${Buffer.from(dec4.slice(0, 32)).equals(Buffer.from(p))}`);
  } catch (e: any) {
    console.log(`  FAILED in ${Date.now() - t3}ms: ${e.message}`);
  }
}

main().catch(e => { console.error(e); process.exit(1); });
