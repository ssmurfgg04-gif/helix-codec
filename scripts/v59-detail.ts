/**
 * v59 Detailed diagnostic — check LDPC decode success rate per cluster.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { ULTIMATE_NANOPORE_V52_CONFIG } from "../src/lib/dna/presets";

const TAG = "[v59-detail]";

async function main() {
  const payload = new Uint8Array(16 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;

  const cfg = ULTIMATE_NANOPORE_V52_CONFIG;
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`${TAG} Encoded ${enc.encoded.oligos.length} oligos`);

  const sim = simulate(enc.encoded.oligos, { ...PRESET_NANOPORE, coverage: 10, seed: 42 });
  console.log(`${TAG} Simulated ${sim.reads.length} reads`);

  // Decode with debug
  process.env.HELIX_DEBUG = "1";
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);

  console.log(`\n${TAG} Decode stats:`);
  console.log(`${TAG}   hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
  console.log(`${TAG}   totalReads: ${dec.stats.totalReads}`);
  console.log(`${TAG}   readsUsed: ${dec.stats.readsUsed}`);
  console.log(`${TAG}   clustersFormed: ${dec.stats.clustersFormed}`);
  console.log(`${TAG}   oligosRecovered: ${dec.stats.oligosRecovered}`);
  console.log(`${TAG}   oligosErased: ${dec.stats.oligosErased}`);
  console.log(`${TAG}   oligosFailedInnerRS: ${dec.stats.oligosFailedInnerRS}`);
  console.log(`${TAG}   oligosFailedOuterRS: ${dec.stats.oligosFailedOuterRS}`);
  console.log(`${TAG}   consensusSuccessRate: ${(dec.stats.consensusSuccessRate * 100).toFixed(1)}%`);
  console.log(`${TAG}   decodeTimeMs: ${dec.stats.decodeTimeMs}`);

  // Per-oligo breakdown
  console.log(`\n${TAG} Per-oligo recovery (first 10):`);
  for (let i = 0; i < Math.min(10, dec.perOligo.length); i++) {
    const p = dec.perOligo[i];
    console.log(`${TAG}   oligo[${i}]: reads=${p.readCount}, crcPassed=${p.crcPassed}, innerRS=${p.innerRS.success ? "OK" : "FAIL"}(${p.innerRS.corrected})`);
  }

  // Count oligos with 0 reads
  const emptyOligos = dec.perOligo.filter(p => p.readCount === 0).length;
  const failedCrc = dec.perOligo.filter(p => p.readCount > 0 && !p.crcPassed).length;
  const passedCrc = dec.perOligo.filter(p => p.crcPassed).length;
  console.log(`\n${TAG} Summary:`);
  console.log(`${TAG}   Oligos with 0 reads: ${emptyOligos}`);
  console.log(`${TAG}   Oligos with reads but CRC failed: ${failedCrc}`);
  console.log(`${TAG}   Oligos with CRC passed: ${passedCrc}`);
}
main().catch(e => { console.error(e); process.exit(1); });
