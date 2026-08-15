/**
 * Debug: investigate oligo 65's mismatch.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { dnaToBytes, unwhitenAddress, whitenAddress, xorWithSeed } from "../src/lib/dna/mapping";
import { crc16Bytes } from "../src/lib/dna/crc16";
import { computeLayout } from "../src/lib/dna/types";
import { deflate } from "pako";

async function main() {
  const payload = new Uint8Array(256 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  
  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "t.bin", contentType: "application/octet-stream" });
  const cfg = ULTIMATE_V55_DENSITY_CONFIG;
  const layout = computeLayout(cfg);
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  
  // Get oligo 65
  const oligo65 = enc.encoded.oligos[65];
  const fwd = enc.encoded.forwardPrimer;
  const rev = enc.encoded.reversePrimer;
  const inner = oligo65.sequence.slice(fwd.length, oligo65.sequence.length - rev.length);
  const innerBlock = dnaToBytes(inner);
  const rsCodeword = innerBlock.slice(0, innerN);
  const crcBytes = innerBlock.slice(innerN, innerN + 2);
  
  console.log(`oligo 65: seed=${oligo65.seed}, gc=${oligo65.gc.toFixed(3)}, maxHp=${oligo65.maxHomopolymer}`);
  
  // Decode
  const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
  const decoded = ldpc.decode(rsCodeword);
  console.log(`LDPC decode: corrected=${decoded.corrected}`);
  
  const whitenedAddr = decoded.data.slice(0, 4);
  const addr = unwhitenAddress(whitenedAddr);
  const decodedIdx = (addr[0] << 16) | (addr[1] << 8) | addr[2];
  const seed = addr[3];
  console.log(`decoded: idx=${decodedIdx}, seed=${seed}`);
  
  // Extract payload (with XOR if seed != 0)
  let payloadBytes = decoded.data.slice(4, 4 + layout.payloadBytes);
  if (seed !== 0) payloadBytes = xorWithSeed(payloadBytes, seed);
  
  // Check CRC
  const computedCrc = crc16Bytes(rsCodeword);
  console.log(`CRC: stored=${crcBytes[0].toString(16).padStart(2,'0')}${crcBytes[1].toString(16).padStart(2,'0')}, computed=${computedCrc[0].toString(16).padStart(2,'0')}${computedCrc[1].toString(16).padStart(2,'0')}`);
  
  // Re-encode and check
  const reEncoded = ldpc.encode(decoded.data);
  let parityMatch = true;
  for (let i = innerK; i < innerN; i++) {
    if (reEncoded[i] !== rsCodeword[i]) {
      parityMatch = false;
      break;
    }
  }
  console.log(`parity match after re-encode: ${parityMatch}`);
  
  // Compare to expected payload
  const compressed = deflate(payload, { level: 9 });
  const chunkSize = layout.payloadBytes;
  const padded = new Uint8Array(enc.encoded.metadata.outerRS.k * chunkSize);
  padded.set(compressed, 0);
  const expectedPayload = padded.slice(65 * chunkSize, 66 * chunkSize);
  
  let mismatch = 0;
  let firstMismatch = -1;
  for (let i = 0; i < chunkSize; i++) {
    if (expectedPayload[i] !== payloadBytes[i]) {
      mismatch++;
      if (firstMismatch === -1) firstMismatch = i;
    }
  }
  console.log(`payload mismatch: ${mismatch}/${chunkSize} bytes, first at byte ${firstMismatch}`);
  if (firstMismatch >= 0) {
    console.log(`  expected[${firstMismatch}] = ${expectedPayload[firstMismatch].toString(16).padStart(2,'0')}`);
    console.log(`  got[${firstMismatch}]      = ${payloadBytes[firstMismatch].toString(16).padStart(2,'0')}`);
  }
  
  // Check if seed=0 was the original (try re-encoding with seed=0)
  console.log(`\n=== Trying seed=0 manually ===`);
  const addr0 = new Uint8Array(4);
  addr0[0] = (65 >> 16) & 0xff;
  addr0[1] = (65 >> 8) & 0xff;
  addr0[2] = 65 & 0xff;
  addr0[3] = 0;
  const whitenedAddr0 = whitenAddress(addr0);
  const rsData0 = new Uint8Array(innerK);
  rsData0.set(whitenedAddr0, 0);
  rsData0.set(expectedPayload, 4);  // seed=0, no XOR
  const codeword0 = ldpc.encode(rsData0);
  const crc0 = crc16Bytes(codeword0);
  console.log(`seed=0 codeword parity: ${Array.from(codeword0.slice(innerK)).map(b => b.toString(16).padStart(2,'0')).join(' ')}`);
  console.log(`stored parity:          ${Array.from(rsCodeword.slice(innerK)).map(b => b.toString(16).padStart(2,'0')).join(' ')}`);
  console.log(`seed=0 CRC: ${crc0[0].toString(16).padStart(2,'0')}${crc0[1].toString(16).padStart(2,'0')}`);
  console.log(`stored CRC: ${crcBytes[0].toString(16).padStart(2,'0')}${crcBytes[1].toString(16).padStart(2,'0')}`);
}

main().catch(e => { console.error(e); process.exit(1); });
