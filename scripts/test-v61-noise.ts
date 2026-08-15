/**
 * Debug: test with realistic Illumina noise.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";

async function main() {
  const payload = new Uint8Array(2.1 * 1024 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  
  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "t.bin", contentType: "application/octet-stream" });
  console.log(`oligos=${enc.encoded.oligos.length}`);
  
  // Use PRESET_ILLUMINA (sub=0.001, ins=0.0005, del=0.001, cov=20)
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  console.log(`simulated ${sim.reads.length} reads (PRESET_ILLUMINA cov=10)`);
  
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, ULTIMATE_V55_DENSITY_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`hash matches: ${dec.hashMatches}`);
  console.log(`stats:`, JSON.stringify(dec.stats));
  
  // Check per-oligo failures
  let crcFailCount = 0;
  let zeroReadCount = 0;
  let totalCorrected = 0;
  for (const p of dec.perOligo) {
    if (!p.crcPassed) crcFailCount++;
    if (p.readCount === 0) zeroReadCount++;
    totalCorrected += p.innerRS.corrected;
  }
  console.log(`oligos CRC failed: ${crcFailCount}/${dec.perOligo.length}`);
  console.log(`oligos with 0 reads: ${zeroReadCount}/${dec.perOligo.length}`);
  console.log(`total bits corrected: ${totalCorrected}`);
  
  // Find first failed oligo (crcPassed=false but readCount>0)
  for (const p of dec.perOligo) {
    if (!p.crcPassed && p.readCount > 0) {
      console.log(`\nfirst failed oligo: idx=${p.index}, readCount=${p.readCount}, corrected=${p.innerRS.corrected}`);
      break;
    }
  }
  
  // Find first oligo with 0 reads
  for (const p of dec.perOligo) {
    if (p.readCount === 0) {
      console.log(`\nfirst oligo with 0 reads: idx=${p.index}`);
      break;
    }
  }
}

main().catch(e => { console.error(e); process.exit(1); });
