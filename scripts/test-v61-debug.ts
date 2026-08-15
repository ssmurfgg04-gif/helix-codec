/**
 * v61 debug: Find where the encoder/decoder disagree (zero-noise hash FAIL).
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { crc16Bytes } from "../src/lib/dna/crc16";
import { dnaToBytes, bytesToDna, whitenAddress, unwhitenAddress } from "../src/lib/dna/mapping";
import * as fs from "fs";

async function main() {
  // Tiny payload (fits in 1 oligo after RS)
  const payload = new Uint8Array(64);
  for (let i = 0; i < 64; i++) payload[i] = i + 1;

  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "t.bin", contentType: "application/octet-stream" });
  console.log(`oligos=${enc.encoded.oligos.length}, payloadBytes per oligo=${enc.encoded.metadata.outerRS.k > 0 ? "?" : "?"}`);
  console.log(`metadata:`, JSON.stringify({
    oligoCount: enc.encoded.metadata.oligoCount,
    outerRS: enc.encoded.metadata.outerRS,
    rawSize: enc.encoded.metadata.rawSize,
    fileSize: enc.encoded.metadata.fileSize,
    fileHash: enc.encoded.metadata.fileHash,
    mappingMode: enc.encoded.metadata.mappingMode,
    innerCode: enc.encoded.metadata.innerCode,
    useConvolutionalInner: enc.encoded.metadata.useConvolutionalInner,
  }, null, 2));

  // Test 1: Direct LDPC encode/decode roundtrip with the actual config
  console.log("\n=== Test 1: Direct LDPC roundtrip (no codec wrapper) ===");
  const layout = { addressBytes: 4, payloadBytes: 145, innerParityBytes: 8, crcBytes: 2 };
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
  console.log(`LDPC: n=${ldpc.n}, k=${ldpc.k}, nsym=${ldpc.nsym}, nBits=${ldpc.nBits}, kBits=${ldpc.kBits}, mBits=${ldpc.mBits}`);

  // Test data
  const testData = new Uint8Array(innerK);
  for (let i = 0; i < innerK; i++) testData[i] = (i * 7 + 3) & 0xff;
  const codeword = ldpc.encode(testData);
  console.log(`encoded ${testData.length}B → ${codeword.length}B`);

  // Decode
  const decoded = ldpc.decode(codeword);
  console.log(`decoded: corrected=${decoded.corrected}, success=${
    decoded.data.length === testData.length &&
    decoded.data.every((b, i) => b === testData[i]) ? "YES" : "NO"}`);

  // Test 2: Add 1-bit error and check if decoder finds it (mBits=64 case)
  console.log("\n=== Test 2: 1-bit error detection (mBits=64) ===");
  const corrupted = codeword.slice();
  corrupted[10] ^= 0x80; // flip MSB of byte 10 (which is bit 80 = info bit 80)
  const dec2 = ldpc.decode(corrupted);
  console.log(`1-bit error: corrected=${dec2.corrected}, success=${
    dec2.data.length === testData.length &&
    dec2.data.every((b, i) => b === testData[i]) ? "YES" : "NO"}`);

  // Test 3: Full codec roundtrip with 0 noise
  console.log("\n=== Test 3: Full codec roundtrip (0 noise) ===");
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 5, seed: 42, subRate: 0, insRate: 0, delRate: 0 });
  console.log(`simulated ${sim.reads.length} reads`);

  const dec = await decodeReads(sim.reads, enc.encoded.metadata, ULTIMATE_V55_DENSITY_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`hash matches: ${dec.hashMatches}`);
  console.log(`stats:`, JSON.stringify(dec.stats));

  // Dump per-oligo info for first 3 oligos
  for (let i = 0; i < Math.min(3, dec.perOligo.length); i++) {
    const p = dec.perOligo[i];
    console.log(`oligo ${i}: idx=${p.index}, crcPassed=${p.crcPassed}, corrected=${p.innerRS.corrected}, payloadLen=${p.payloadBytes.length}, first8=${Array.from(p.payloadBytes.slice(0, 8)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  }

  // Test 4: Manually decode the FIRST oligo and verify
  console.log("\n=== Test 4: Manual decode of oligo 0 ===");
  const oligo0 = enc.encoded.oligos[0];
  const fwd = enc.encoded.forwardPrimer;
  const rev = enc.encoded.reversePrimer;
  // Strip primers
  const inner = oligo0.sequence.slice(fwd.length, oligo0.sequence.length - rev.length);
  console.log(`inner DNA length: ${inner.length} (expected ${660})`);
  const innerBlock = dnaToBytes(inner);
  console.log(`innerBlock length: ${innerBlock.length}`);
  const rsCodeword = innerBlock.slice(0, innerN);
  const crcBytes = innerBlock.slice(innerN);
  // Compute CRC of rsCodeword
  const computedCrc = crc16Bytes(rsCodeword);
  console.log(`CRC: stored=${crcBytes[0].toString(16).padStart(2,'0')}${crcBytes[1].toString(16).padStart(2,'0')}, computed=${computedCrc[0].toString(16).padStart(2,'0')}${computedCrc[1].toString(16).padStart(2,'0')}, match=${crcBytes[0] === computedCrc[0] && crcBytes[1] === computedCrc[1]}`);
  // LDPC decode
  const dec4 = ldpc.decode(rsCodeword);
  console.log(`LDPC decode: corrected=${dec4.corrected}`);
  const whitenedAddr = dec4.data.slice(0, 4);
  const addr = unwhitenAddress(whitenedAddr);
  const decodedIdx = (addr[0] << 16) | (addr[1] << 8) | addr[2];
  console.log(`decoded address: idx=${decodedIdx}, seed=${addr[3]} (expected idx=0, seed=0)`);
  console.log(`first 8 payload bytes: ${Array.from(dec4.data.slice(4, 12)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);

  // Test 5: Re-encode the decoded data and verify it matches the original codeword
  console.log("\n=== Test 5: Re-encode decoded data, compare to original ===");
  const reEncoded = ldpc.encode(dec4.data);
  let mismatches = 0;
  for (let i = 0; i < innerN; i++) {
    if (reEncoded[i] !== rsCodeword[i]) {
      if (mismatches < 5) console.log(`  mismatch at byte ${i}: stored=${rsCodeword[i].toString(16).padStart(2,'0')}, re-encoded=${reEncoded[i].toString(16).padStart(2,'0')}`);
      mismatches++;
    }
  }
  console.log(`total mismatches: ${mismatches} / ${innerN}`);
}

main().catch(e => { console.error(e); process.exit(1); });
