/**
 * Debug: find which oligo's payload is wrong.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { computeLayout } from "../src/lib/dna/types";
import { deflate } from "pako";

async function main() {
  const payload = new Uint8Array(256 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  
  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "t.bin", contentType: "application/octet-stream" });
  const cfg = ULTIMATE_V55_DENSITY_CONFIG;
  const layout = computeLayout(cfg);
  console.log(`layout: payloadBytes=${layout.payloadBytes}, innerParityBytes=${layout.innerParityBytes}, totalInnerBytes=${layout.totalInnerBytes}`);
  console.log(`metadata: outerRS(n=${enc.encoded.metadata.outerRS.n}, k=${enc.encoded.metadata.outerRS.k}), oligoCount=${enc.encoded.metadata.oligoCount}`);
  
  // Re-compute what the compressed payload should be
  const compressed = deflate(payload, { level: 9 });
  console.log(`compressed length: ${compressed.length}`);
  console.log(`expected data oligos: ${enc.encoded.metadata.outerRS.k}, payloadBytes each: ${layout.payloadBytes}`);
  console.log(`expected total payload bytes: ${enc.encoded.metadata.outerRS.k * layout.payloadBytes}`);
  
  // Simulate with noise
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`\nhash matches: ${dec.hashMatches}`);
  
  // Reconstruct the expected per-oligo payloads
  // The encoder pads the compressed data to k * payloadBytes, then chunks it.
  const k = enc.encoded.metadata.outerRS.k;
  const chunkSize = layout.payloadBytes;
  const padded = new Uint8Array(k * chunkSize);
  padded.set(compressed, 0);
  
  // Compare each oligo's recovered payload to the expected
  let mismatchCount = 0;
  for (let i = 0; i < k; i++) {
    const expected = padded.slice(i * chunkSize, (i + 1) * chunkSize);
    const recovered = dec.perOligo.find(p => p.index === i)?.payloadBytes;
    if (!recovered) {
      console.log(`oligo ${i}: NO RECOVERED PAYLOAD`);
      mismatchCount++;
      continue;
    }
    let mismatch = false;
    let firstMismatch = -1;
    for (let j = 0; j < chunkSize; j++) {
      if (expected[j] !== recovered[j]) {
        mismatch = true;
        if (firstMismatch === -1) firstMismatch = j;
      }
    }
    if (mismatch) {
      mismatchCount++;
      if (mismatchCount <= 5) {
        console.log(`oligo ${i}: MISMATCH at byte ${firstMismatch}, expected=${expected[firstMismatch].toString(16)}, got=${recovered[firstMismatch].toString(16)}`);
      }
    }
  }
  console.log(`\ntotal mismatched oligos: ${mismatchCount}/${k}`);
}

main().catch(e => { console.error(e); process.exit(1); });
