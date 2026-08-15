// Find which oligo has wrong data (false positive).
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import { computeLayout } from "../src/lib/dna/types";
import { unwhitenAddress, xorWithSeed } from "../src/lib/dna/mapping";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { ReedSolomon } from "../src/lib/dna/reedsolomon";

async function main() {
  console.log("=== Find false positive oligo ===\n");

  const bigPayload = new Uint8Array(10 * 1024 * 1024);
  for (let i = 0; i < bigPayload.length; i++) bigPayload[i] = (i * 31 + 17) & 0xff;

  const bigConfig = { ...DEFAULT_CONFIG, oligoLength: 300, primerLength: 20, outerParityRatio: 0.3 };
  const encBig = await encodeFile(bigPayload, bigConfig, { fileName: "10mb.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${encBig.encoded.oligos.length} oligos`);

  // Get the original per-oligo payloads (from the encoder)
  const layout = computeLayout(bigConfig);
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  const useLDPC = (bigConfig.innerCode ?? "rs") === "ldpc";
  const innerLdpc = useLDPC ? new LDPCInnerCode({ n: innerN, k: innerK }) : null;
  const innerRs = new ReedSolomon({ n: innerN, k: innerK });

  const originalPayloads = new Map<number, Uint8Array>();
  for (const oligo of encBig.encoded.oligos) {
    // Re-derive the payload from the encoded oligo
    const inner = oligo.sequence.slice(bigConfig.primerLength, oligo.sequence.length - bigConfig.primerLength);
    // Convert DNA to bytes
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
    // Decode to get the info
    let info: Uint8Array;
    if (useLDPC && innerLdpc) {
      // LDPC codeword is already valid (no errors), so decode just returns it
      info = innerLdpc.decode(codeword).data;
    } else {
      info = innerRs.decode(codeword).data;
    }
    // Extract payload
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

  // Check each recovered oligo's payload against the original
  let mismatches = 0;
  for (const p of dec.perOligo) {
    if (!p.crcPassed) continue; // skip failed oligos
    const original = originalPayloads.get(p.index);
    if (!original) continue;
    let mismatch = false;
    for (let i = 0; i < original.length; i++) {
      if (original[i] !== p.payloadBytes[i]) {
        mismatch = true;
        break;
      }
    }
    if (mismatch) {
      mismatches++;
      if (mismatches <= 5) {
        console.log(`\n  Oligo ${p.index}: PAYLOAD MISMATCH (false positive!)`);
        console.log(`    readCount: ${p.readCount}, corrected: ${p.innerRS.corrected}`);
        console.log(`    Original: ${Array.from(original.slice(0, 16)).map(b => b.toString(16).padStart(2, "0")).join(" ")}`);
        console.log(`    Decoded:  ${Array.from(p.payloadBytes.slice(0, 16)).map(b => b.toString(16).padStart(2, "0")).join(" ")}`);
      }
    }
  }
  console.log(`\nTotal false positives: ${mismatches}`);
  console.log(`Hash matches: ${dec.hashMatches}`);
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
