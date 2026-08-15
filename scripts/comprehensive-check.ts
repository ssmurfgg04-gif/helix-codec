// Comprehensive check: compare ALL decoded oligos (including outer RS recovered) with originals.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG, computeLayout } from "../src/lib/dna/types";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { ReedSolomon } from "../src/lib/dna/reedsolomon";
import { dnaToBytes, unwhitenAddress, xorWithSeed } from "../src/lib/dna/mapping";

async function main() {
  console.log("=== Comprehensive oligo check ===\n");

  // Use text payload (same as erlich-end-to-end.ts)
  const textPattern = "The Helix Codec encodes digital data into synthetic DNA oligos for archival storage. ";
  const repeatCount = Math.ceil((10 * 1024 * 1024) / textPattern.length);
  const bigPayload = Buffer.from(textPattern.repeat(repeatCount).slice(0, 10 * 1024 * 1024));

  const bigConfig = { ...DEFAULT_CONFIG, oligoLength: 300, primerLength: 20, outerParityRatio: 0.3 };
  const encBig = await encodeFile(bigPayload, bigConfig, { fileName: "10mb.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${encBig.encoded.oligos.length} oligos, ${encBig.stats.compressedSize} compressed`);
  console.log(`  Outer RS: n=${encBig.encoded.metadata.outerRS.n}, k=${encBig.encoded.metadata.outerRS.k}`);
  console.log(`  Inner: n=${encBig.encoded.metadata.innerRS.n}, k=${encBig.encoded.metadata.innerRS.k}, code=${encBig.encoded.metadata.innerCode}`);

  // Extract original per-oligo payloads
  const layout = computeLayout(bigConfig);
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  const useLDPC = (bigConfig.innerCode ?? "rs") === "ldpc";
  const innerLdpc = useLDPC ? new LDPCInnerCode({ n: innerN, k: innerK }) : null;
  const innerRs = new ReedSolomon({ n: innerN, k: innerK });

  const originalPayloads = new Map<number, Uint8Array>();
  for (const oligo of encBig.encoded.oligos) {
    const inner = oligo.sequence.slice(bigConfig.primerLength, oligo.sequence.length - bigConfig.primerLength);
    const bytes: number[] = [];
    for (let i = 0; i < inner.length; i += 4) {
      const bits = ["A", "C", "G", "T"].indexOf(inner[i]) << 6 |
                   ["A", "C", "G", "T"].indexOf(inner[i+1]) << 4 |
                   ["A", "C", "G", "T"].indexOf(inner[i+2]) << 2 |
                   ["A", "C", "G", "T"].indexOf(inner[i+3]);
      bytes.push(bits);
    }
    const innerBlock = new Uint8Array(bytes);
    const codeword = innerBlock.slice(0, innerN);
    let info: Uint8Array;
    try {
      if (useLDPC && innerLdpc) {
        info = innerLdpc.decode(codeword).data;
      } else {
        info = innerRs.decode(codeword).data;
      }
    } catch (e) {
      console.log(`  Oligo ${oligo.index}: FAILED to decode original — ${(e as Error).message.slice(0, 80)}`);
      continue;
    }
    const whitenedAddr = info.slice(0, layout.addressBytes);
    const addr = unwhitenAddress(whitenedAddr);
    const seed = addr[3];
    let payload = info.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
    if (seed !== 0) payload = xorWithSeed(payload, seed);
    originalPayloads.set(oligo.index, payload);
  }

  console.log(`Original payloads extracted: ${originalPayloads.size}`);

  // Simulate and decode
  const sim = simulate(encBig.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const dec = await decodeReads(sim.reads, encBig.encoded.metadata, bigConfig, encBig.encoded.forwardPrimer, encBig.encoded.reversePrimer, true);

  console.log(`\nDecode: hashMatches=${dec.hashMatches}, oligosRecovered=${dec.stats.oligosRecovered}/${encBig.encoded.oligos.length}`);

  // Check constraint violations in encoded oligos
  const violatedOligos = encBig.encoded.oligos.filter(o => 
    o.gc < bigConfig.constraints.gcMin || o.gc > bigConfig.constraints.gcMax || o.maxHomopolymer > bigConfig.constraints.maxHomopolymer
  );
  console.log(`\nConstraint-violating oligos: ${violatedOligos.length}/${encBig.encoded.oligos.length}`);
  for (const v of violatedOligos.slice(0, 10)) {
    console.log(`  Oligo ${v.index}: gc=${v.gc.toFixed(3)}, maxHomopolymer=${v.maxHomopolymer}, seed=${v.seed}`);
  }

  // Check which oligos were erased
  const erasedSet = new Set<number>();
  for (const p of dec.perOligo) {
    if (!p.crcPassed) erasedSet.add(p.index);
  }
  console.log(`\nErased oligos (in decode output): ${erasedSet.size}`);

  // Now check the FINAL recovered data against originals
  // The decode output has `data` which is the decompressed file. Let's check the per-oligo payloads
  // by re-encoding the recovered data.
  // Actually, let's check the perOligo array for wrong data.
  let wrongDataCount = 0;
  let wrongDataIndices: number[] = [];
  for (const p of dec.perOligo) {
    if (!p.crcPassed) continue; // skip erased
    const original = originalPayloads.get(p.index);
    if (!original) continue;
    let mismatch = false;
    for (let i = 0; i < original.length; i++) {
      if (original[i] !== p.payloadBytes[i]) { mismatch = true; break; }
    }
    if (mismatch) {
      wrongDataCount++;
      wrongDataIndices.push(p.index);
    }
  }
  console.log(`\nOligos with WRONG data (false positives): ${wrongDataCount}`);
  if (wrongDataIndices.length > 0) {
    console.log(`  Indices: ${wrongDataIndices.slice(0, 20).join(", ")}`);
  }

  // Check if the erased oligos overlap with constraint-violating oligos
  const erasedArr = Array.from(erasedSet);
  const violatedSet = new Set(violatedOligos.map(o => o.index));
  const overlap = erasedArr.filter(i => violatedSet.has(i));
  console.log(`\nErased oligos that also violate constraints: ${overlap.length}/${erasedArr.length}`);
  if (erasedArr.length > 0) {
    console.log(`  Erased indices: ${erasedArr.join(", ")}`);
  }

  // Now check the FINAL recovered data by re-encoding and comparing per-oligo
  // The decode result has `data` which is the decompressed file. Let's re-encode it
  // and compare per-oligo payloads.
  if (dec.data) {
    console.log(`\nRecovered data size: ${dec.data.length} bytes (expected ${bigPayload.length})`);
    // Re-encode the recovered data to get per-oligo payloads
    const reEnc = await encodeFile(dec.data, bigConfig, { fileName: "recovered.bin", contentType: "application/octet-stream" });
    console.log(`Re-encoded: ${reEnc.encoded.oligos.length} oligos`);

    // Compare per-oligo payloads
    let wrongPayloads = 0;
    for (const oligo of reEnc.encoded.oligos) {
      const original = originalPayloads.get(oligo.index);
      if (!original) continue;

      // Extract re-encoded payload
      const inner = oligo.sequence.slice(bigConfig.primerLength, oligo.sequence.length - bigConfig.primerLength);
      const bytes: number[] = [];
      for (let i = 0; i < inner.length; i += 4) {
        const bits = ["A", "C", "G", "T"].indexOf(inner[i]) << 6 |
                     ["A", "C", "G", "T"].indexOf(inner[i+1]) << 4 |
                     ["A", "C", "G", "T"].indexOf(inner[i+2]) << 2 |
                     ["A", "C", "G", "T"].indexOf(inner[i+3]);
        bytes.push(bits);
      }
      const innerBlock = new Uint8Array(bytes);
      const codeword = innerBlock.slice(0, innerN);
      let info: Uint8Array;
      try {
        if (useLDPC && innerLdpc) {
          info = innerLdpc.decode(codeword).data;
        } else {
          info = innerRs.decode(codeword).data;
        }
      } catch { continue; }
      const whitenedAddr = info.slice(0, layout.addressBytes);
      const addr = unwhitenAddress(whitenedAddr);
      const seed = addr[3];
      let payload = info.slice(layout.addressBytes, layout.addressBytes + layout.payloadBytes);
      if (seed !== 0) payload = xorWithSeed(payload, seed);

      let mismatch = false;
      for (let i = 0; i < original.length; i++) {
        if (original[i] !== payload[i]) { mismatch = true; break; }
      }
      if (mismatch) {
        wrongPayloads++;
        if (wrongPayloads <= 5) {
          console.log(`  Oligo ${oligo.index}: RECOVERED PAYLOAD MISMATCH`);
          console.log(`    Original: ${Array.from(original.slice(0, 16)).map(b => b.toString(16).padStart(2, "0")).join(" ")}`);
          console.log(`    Recovered: ${Array.from(payload.slice(0, 16)).map(b => b.toString(16).padStart(2, "0")).join(" ")}`);
        }
      }
    }
    console.log(`\nOligos with wrong recovered payload: ${wrongPayloads}`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
