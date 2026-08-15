/**
 * v63 Single-read K=9 Viterbi test — verify decoder works, measure timing
 */
import { IndelTolerantConvolutionalInnerCode } from "../src/lib/dna/convolutional-indel";
import { NASA_K9_CONFIG } from "../src/lib/dna/convolutional-k9";

async function main() {
  const inputBytes = 63; // matches 300nt nanopore config
  const innerCode = new IndelTolerantConvolutionalInnerCode(inputBytes, { conv: NASA_K9_CONFIG });

  // Encode random data
  const data = new Uint8Array(inputBytes);
  for (let i = 0; i < inputBytes; i++) data[i] = Math.floor(Math.random() * 256);
  const encoded = innerCode.encode(data);
  console.log(`K=9 indel Viterbi: inputBytes=${inputBytes}, encoded=${encoded.length}B`);

  // Decode clean (no errors)
  const t0 = performance.now();
  const decoded = innerCode.decode(encoded);
  const decMs = performance.now() - t0;
  console.log(`clean decode: ${decMs.toFixed(0)}ms`);
  const match = Buffer.compare(Buffer.from(data), Buffer.from(decoded)) === 0;
  console.log(`  data match: ${match ? "OK ✅" : "FAIL ❌"}`);

  // Test with 1 substitution
  const noisy = encoded.slice();
  noisy[10] ^= 0x40; // flip 1 bit
  const t1 = performance.now();
  const decoded1 = innerCode.decode(noisy);
  const dec1Ms = performance.now() - t1;
  const match1 = Buffer.compare(Buffer.from(data), Buffer.from(decoded1)) === 0;
  console.log(`1-sub decode: ${dec1Ms.toFixed(0)}ms, match: ${match1 ? "OK ✅" : "FAIL ❌"}`);

  // Project full nanopore decode time
  const perDecodeMs = decMs;
  const totalReads = 250; // 25 oligos × 10× coverage
  const projectedSec = (perDecodeMs * totalReads) / 1000;
  console.log(`\nProjected full nanopore decode (250 reads): ${projectedSec.toFixed(0)}s`);
}

main().catch(e => { console.error(e); process.exit(1); });
