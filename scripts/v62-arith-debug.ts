/**
 * v62: Debug script for arithmetic-v2 decode.
 * Encodes a single oligo and traces the decode path.
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { ULTIMATE_V61_ARITHMETIC_CONFIG } from "../src/lib/dna/presets";
import { simulate } from "../src/lib/dna/simulate";
import { computeLayoutAuto } from "../src/lib/dna/types";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { arithmeticDnaToBytesCrc } from "../src/lib/dna/markov-arithmetic";
import { dnaToBytes, unwhitenAddress } from "../src/lib/dna/mapping";

async function main() {
  console.log("=== v62 Arithmetic-v2 Debug ===\n");

  // Small payload
  const data = new Uint8Array(256);
  for (let i = 0; i < data.length; i++) data[i] = (i * 37 + 13) & 0xff;

  const cfg = { ...ULTIMATE_V61_ARITHMETIC_CONFIG };
  const layout = computeLayoutAuto(cfg);
  console.log("Layout:", {
    addressBytes: layout.addressBytes,
    payloadBytes: layout.payloadBytes,
    innerParityBytes: layout.innerParityBytes,
    crcBytes: layout.crcBytes,
    totalInnerBytes: layout.totalInnerBytes,
  });

  const encodeResult = await encodeFile(data, cfg, {
    fileName: "test.bin",
    contentType: "application/octet-stream",
  });
  console.log(`Encoded: ${encodeResult.encoded.oligos.length} oligos`);
  console.log(`Metadata innerRS: ${JSON.stringify(encodeResult.encoded.metadata.innerRS)}`);

  // Check first oligo
  const oligo = encodeResult.encoded.oligos[0];
  const fwd = encodeResult.encoded.forwardPrimer;
  const rev = encodeResult.encoded.reversePrimer;
  const inner = oligo.sequence.slice(fwd.length, oligo.sequence.length - rev.length);
  console.log(`\nOligo 0 inner DNA length: ${inner.length}`);
  console.log(`  Address DNA (first 16): ${inner.slice(0, 16)}`);
  console.log(`  Arithmetic DNA (next 80): ${inner.slice(16, 96)}`);

  // Simulate reads (noiseless)
  const simResult = simulate(encodeResult.encoded.oligos, {
    coverage: 10,
    substitutionRate: 0,
    insertionRate: 0,
    deletionRate: 0,
    dropoutRate: 0,
    seed: 42,
  });
  console.log(`\nSimulated ${simResult.reads.length} reads`);

  // Try to decode oligo 0 manually
  const read = simResult.reads[0];
  const readInner = read.sequence; // already trimmed? Let me check
  console.log(`Read 0 sequence length: ${read.sequence.length}`);
  // Trim primers
  const fwdMatch = read.sequence.startsWith(fwd);
  const revMatch = read.sequence.endsWith(rev);
  console.log(`  FWD primer match: ${fwdMatch}, REV primer match: ${revMatch}`);

  const innerDna = read.sequence.slice(fwd.length, read.sequence.length - rev.length);
  console.log(`  Inner DNA length: ${innerDna.length}`);

  // Decode address
  const addressDna = innerDna.slice(0, 16);
  const addressBytes = dnaToBytes(addressDna);
  const unwhitened = unwhitenAddress(addressBytes);
  const idx = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
  console.log(`  Decoded address index: ${idx} (expected: ${oligo.index})`);

  // Decode arithmetic stream
  const arithmeticDna = innerDna.slice(16);
  console.log(`  Arithmetic DNA length: ${arithmeticDna.length}`);

  const blockSize = 80;
  const innerNArith = layout.payloadBytes + layout.innerParityBytes;
  console.log(`  Expected LDPC codeword bytes: ${innerNArith}`);

  const result = arithmeticDnaToBytesCrc(arithmeticDna, 3, innerNArith, blockSize);
  console.log(`  Decoded ${result.data.length} bytes from arithmetic stream`);
  const eraseCount = result.erasures.filter((e) => e).length;
  console.log(`  Erasures: ${eraseCount} bytes`);

  // LDPC decode
  const innerK = layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  console.log(`  LDPC: k=${innerK}, n=${innerN}`);

  const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
  const rsCodeword = result.data; // LDPC codeword = payload + parity
  console.log(`  rsCodeword length: ${rsCodeword.length}`);

  try {
    if (eraseCount > 0) {
      const erasePos: number[] = [];
      for (let i = 0; i < result.erasures.length && i < innerN; i++) {
        if (result.erasures[i]) {
          for (let bit = 0; bit < 8; bit++) {
            erasePos.push(i * 8 + bit);
          }
        }
      }
      console.log(`  Erase positions (bits): ${erasePos.length}`);
      const r = ldpc.decodeWithErasures(rsCodeword, erasePos);
      console.log(`  LDPC erasure decode SUCCESS! corrected=${r.corrected}, erased=${r.erased}`);
      console.log(`  Decoded data length: ${r.data.length}`);
    } else {
      const r = ldpc.decode(rsCodeword);
      console.log(`  LDPC hard decode SUCCESS! corrected=${r.corrected}`);
    }
  } catch (e) {
    console.log(`  LDPC decode FAILED: ${e}`);
  }

  // Now try the full decode
  console.log("\n--- Full Decode ---");
  const decodeResult = await decodeReads(
    simResult.reads,
    encodeResult.encoded.metadata,
    cfg,
    fwd,
    rev,
  );
  console.log(`Decoded: ${decodeResult.data.length} bytes (expected ${data.length})`);
  console.log(`Hash match: ${decodeResult.hashMatches ? "✅ YES" : "❌ NO"}`);

  // Compare bytes
  if (decodeResult.data.length === data.length) {
    let mismatches = 0;
    let firstMismatch = -1;
    for (let i = 0; i < data.length; i++) {
      if (decodeResult.data[i] !== data[i]) {
        mismatches++;
        if (firstMismatch === -1) firstMismatch = i;
      }
    }
    console.log(`Byte mismatches: ${mismatches}/${data.length}, first at offset ${firstMismatch}`);
  }

  // Check hash
  const { createHash } = await import("crypto");
  const decodedHash = createHash("sha256").update(decodeResult.data).digest("hex");
  const originalHash = createHash("sha256").update(data).digest("hex");
  console.log(`Decoded hash:  ${decodedHash}`);
  console.log(`Original hash: ${originalHash}`);
  console.log(`Metadata hash: ${encodeResult.encoded.metadata.fileHash}`);
  console.log(`Hash match (manual): ${decodedHash === originalHash ? "✅ YES" : "❌ NO"}`);
  console.log(`Hash match (metadata): ${decodedHash === encodeResult.encoded.metadata.fileHash ? "✅ YES" : "❌ NO"}`);
}

main().catch(console.error);
