import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";

async function main() {
  const payload = new Uint8Array(2.1 * 1024 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  
  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "t.bin", contentType: "application/octet-stream" });
  console.log(`oligos=${enc.encoded.oligos.length}, outerRS(n=${enc.encoded.metadata.outerRS.n}, k=${enc.encoded.metadata.outerRS.k})`);
  
  // ZERO noise
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 5, seed: 42, substitutionRate: 0, insertionRate: 0, deletionRate: 0 });
  console.log(`simulated ${sim.reads.length} reads (ZERO noise)`);
  
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, ULTIMATE_V55_DENSITY_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
  console.log(`stats:`, JSON.stringify(dec.stats));
}
main().catch(e => { console.error(e); process.exit(1); });
