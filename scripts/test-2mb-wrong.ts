import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { computeLayout } from "../src/lib/dna/types";
import { deflate } from "pako";

async function main() {
  const payload = new Uint8Array(2.1 * 1024 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  
  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "t.bin", contentType: "application/octet-stream" });
  const cfg = ULTIMATE_V55_DENSITY_CONFIG;
  const layout = computeLayout(cfg);
  const compressed = deflate(payload, { level: 9 });
  const chunkSize = layout.payloadBytes;
  const padded = new Uint8Array(enc.encoded.metadata.outerRS.k * chunkSize);
  padded.set(compressed, 0);
  
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
  
  let wrongCount = 0;
  for (let i = 0; i < enc.encoded.metadata.outerRS.k; i++) {
    const expected = padded.slice(i * chunkSize, (i + 1) * chunkSize);
    const recovered = dec.perOligo.find(p => p.index === i)?.payloadBytes;
    if (!recovered) continue;
    for (let j = 0; j < chunkSize; j++) {
      if (expected[j] !== recovered[j]) {
        wrongCount++;
        if (wrongCount <= 10) {
          console.log(`oligo ${i}: mismatch at byte ${j}, expected=${expected[j].toString(16)}, got=${recovered[j].toString(16)}, corrected=${dec.perOligo.find(p=>p.index===i)?.innerRS.corrected}`);
        }
        break;
      }
    }
  }
  console.log(`\ntotal wrong oligos: ${wrongCount}/${enc.encoded.metadata.outerRS.k}`);
}
main().catch(e => { console.error(e); process.exit(1); });
