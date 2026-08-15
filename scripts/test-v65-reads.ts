/**
 * Debug: find which read was accepted for oligo 65.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { dnaToBytes, unwhitenAddress, xorWithSeed } from "../src/lib/dna/mapping";
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
  
  // Get expected payload for oligo 65
  const compressed = deflate(payload, { level: 9 });
  const chunkSize = layout.payloadBytes;
  const padded = new Uint8Array(enc.encoded.metadata.outerRS.k * chunkSize);
  padded.set(compressed, 0);
  const expectedPayload = padded.slice(65 * chunkSize, 66 * chunkSize);
  
  // Simulate
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  
  // Find reads clustered to oligo 65
  const fwd = enc.encoded.forwardPrimer;
  const rev = enc.encoded.reversePrimer;
  const addressNt = layout.addressBytes * 4;
  const expectedDnaLen = layout.totalInnerBytes * 4;
  
  const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
  
  let readIdx = 0;
  for (const read of sim.reads) {
    // Trim primers (simplified)
    if (read.sequence.length < fwd.length + rev.length) continue;
    const fwdMatch = read.sequence.slice(0, fwd.length);
    const revMatch = read.sequence.slice(read.sequence.length - rev.length);
    if (fwdMatch !== fwd && revMatch !== rev) continue;
    const inner = read.sequence.slice(fwd.length, read.sequence.length - rev.length);
    if (inner.length < addressNt) continue;
    
    // Extract address
    const addressDna = inner.slice(0, addressNt);
    let addressBytes: Uint8Array;
    try {
      addressBytes = dnaToBytes(addressDna);
    } catch { continue; }
    const unwhitened = unwhitenAddress(addressBytes);
    const idx = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
    
    if (idx !== 65) continue;
    
    // This read is clustered to oligo 65
    let dna = inner;
    if (dna.length > expectedDnaLen) dna = dna.slice(0, expectedDnaLen);
    else if (dna.length < expectedDnaLen) dna = dna + "A".repeat(expectedDnaLen - dna.length);
    
    const innerBlock = dnaToBytes(dna);
    const rsCodeword = innerBlock.slice(0, innerN);
    const crcBytes = innerBlock.slice(innerN, innerN + 2);
    
    // LDPC decode
    let decoded;
    try {
      decoded = ldpc.decode(rsCodeword);
    } catch {
      continue;
    }
    
    // CRC check (fast path)
    if (decoded.corrected === 0) {
      const expectedCrc = crc16Bytes(rsCodeword);
      if (expectedCrc[0] !== crcBytes[0] || expectedCrc[1] !== crcBytes[1]) continue;
    } else {
      const reEncoded = ldpc.encode(decoded.data);
      const expectedCrc = crc16Bytes(reEncoded);
      if (expectedCrc[0] !== crcBytes[0] || expectedCrc[1] !== crcBytes[1]) continue;
    }
    
    // Address verification
    const whitenedAddr = decoded.data.slice(0, layout.addressBytes);
    const addr = unwhitenAddress(whitenedAddr);
    const decodedIdx = (addr[0] << 16) | (addr[1] << 8) | addr[2];
    if (decodedIdx !== 65) continue;
    
    // This read would be ACCEPTED
    const seed = addr[3];
    let payloadBytes = decoded.data.slice(4, 4 + layout.payloadBytes);
    if (seed !== 0) payloadBytes = xorWithSeed(payloadBytes, seed);
    
    let mismatch = 0;
    let firstMismatch = -1;
    for (let i = 0; i < chunkSize; i++) {
      if (expectedPayload[i] !== payloadBytes[i]) {
        mismatch++;
        if (firstMismatch === -1) firstMismatch = i;
      }
    }
    
    console.log(`read ${readIdx}: corrected=${decoded.corrected}, seed=${seed}, mismatch=${mismatch}/${chunkSize}` +
      (firstMismatch >= 0 ? `, first at byte ${firstMismatch} (expected ${expectedPayload[firstMismatch].toString(16)}, got ${payloadBytes[firstMismatch].toString(16)})` : ""));
    
    readIdx++;
    if (readIdx >= 15) break;
  }
}

main().catch(e => { console.error(e); process.exit(1); });
