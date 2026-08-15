import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";

async function main() {
  // Small payload, perfect conditions
  const payload = new Uint8Array(64 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  
  const t0 = Date.now();
  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`encode: ${Date.now() - t0}ms, oligos=${enc.encoded.oligos.length}`);
  console.log(`metadata: oligoCount=${enc.encoded.metadata.oligoCount}, outerRS(n=${enc.encoded.metadata.outerRS.n},k=${enc.encoded.metadata.outerRS.k})`);
  console.log(`mappingMode=${enc.encoded.metadata.mappingMode}, innerCode=${enc.encoded.metadata.innerCode}`);
  
  // NO NOISE — perfect reads
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42, subRate: 0, insRate: 0, delRate: 0 });
  console.log(`simulated ${sim.reads.length} reads (no noise)`);
  
  const t1 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, ULTIMATE_V55_DENSITY_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`decode: ${Date.now() - t1}ms`);
  console.log(`hash matches: ${dec.hashMatches}`);
  console.log(`expected: ${enc.encoded.metadata.fileHash}`);
  console.log(`got:      ${dec.hash}`);
  console.log(`stats:`, JSON.stringify(dec.stats, null, 2));
  
  // Check first few oligos
  let failures = 0;
  for (const p of dec.perOligo) {
    if (!p.crcPassed) {
      failures++;
      if (failures <= 3) {
        console.log(`oligo ${p.index}: crcPassed=${p.crcPassed}, readCount=${p.readCount}, corrected=${p.innerRS.corrected}, success=${p.innerRS.success}`);
      }
    }
  }
  console.log(`oligos with crcPassed=false: ${failures}/${dec.perOligo.length}`);
}

main().catch(e => { console.error(e); process.exit(1); });
