/**
 * Debug: dump the actual LDPC codeword produced by the codec.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { dnaToBytes, unwhitenAddress } from "../src/lib/dna/mapping";
import { crc16Bytes } from "../src/lib/dna/crc16";
import { computeLayout } from "../src/lib/dna/types";

async function main() {
  const payload = new Uint8Array(64);
  for (let i = 0; i < 64; i++) payload[i] = i + 1;

  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "t.bin", contentType: "application/octet-stream" });
  const cfg = ULTIMATE_V55_DENSITY_CONFIG;
  const layout = computeLayout(cfg);
  console.log("layout:", JSON.stringify({
    addressBytes: layout.addressBytes,
    payloadBytes: layout.payloadBytes,
    innerParityBytes: layout.innerParityBytes,
    crcBytes: layout.crcBytes,
    totalInnerBytes: layout.totalInnerBytes,
  }));
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  console.log(`innerK=${innerK}, innerN=${innerN}`);

  // Get oligo 0
  const oligo0 = enc.encoded.oligos[0];
  const fwd = enc.encoded.forwardPrimer;
  const rev = enc.encoded.reversePrimer;
  const inner = oligo0.sequence.slice(fwd.length, oligo0.sequence.length - rev.length);
  console.log(`inner DNA len: ${inner.length}`);
  const innerBlock = dnaToBytes(inner);
  console.log(`innerBlock len: ${innerBlock.length}`);
  console.log(`first 16 bytes: ${Array.from(innerBlock.slice(0, 16)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  console.log(`bytes 149-165 (parity+CRC+pad): ${Array.from(innerBlock.slice(innerK)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);

  // Now manually encode with LDPC and compare
  const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
  const rsData = innerBlock.slice(0, innerK);  // address + payload
  console.log(`\nrsData (info): ${Array.from(rsData.slice(0, 16)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  const codeword = ldpc.encode(rsData);
  console.log(`ldpc.encode parity: ${Array.from(codeword.slice(innerK)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  console.log(`stored parity:      ${Array.from(innerBlock.slice(innerK, innerN)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  
  // Now decode and check
  const decoded = ldpc.decode(innerBlock.slice(0, innerN));
  console.log(`\nldpc.decode: corrected=${decoded.corrected}`);
  console.log(`decoded info matches stored info: ${
    decoded.data.every((b, i) => b === innerBlock[i]) ? "YES" : "NO"}`);

  // Check syndrome
  // The syndrome = H * recv. For a valid codeword, syndrome = 0.
  // If stored parity is zero but should be non-zero, syndrome != 0.
  // But the decoder's fast path (mBits > 32) doesn't compute syndrome!
  console.log(`\n=== BUG CONFIRMATION ===`);
  console.log(`LDPC mBits=${ldpc.mBits} (>32: ${ldpc.mBits > 32})`);
  console.log(`If mBits>32, decoder fast path SKIPS syndrome computation → always returns corrected=0`);
  console.log(`This means: errors in the parity region go UNDETECTED, and the wrong codeword is silently accepted.`);
}

main().catch(e => { console.error(e); process.exit(1); });
