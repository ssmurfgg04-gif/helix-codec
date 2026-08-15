import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG, computeDensity } from "../src/lib/dna/presets";

async function bench(name: string, cfg: any, size: number, cov: number) {
  const payload = new Uint8Array(size);
  for (let i = 0; i < size; i++) payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "b.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: cov, seed: 42 });
  const t1 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;
  const encMBs = (size / 1024 / 1024) / (encMs / 1000);
  const decMBs = (size / 1024 / 1024) / (decMs / 1000);
  console.log(`${name}: ${enc.encoded.oligos.length} oligos, ${cov}x cov`);
  console.log(`  encode: ${encMs}ms = ${encMBs.toFixed(2)} MB/s`);
  console.log(`  decode: ${decMs}ms = ${decMBs.toFixed(2)} MB/s`);
  console.log(`  hash: ${dec.hashMatches ? "OK" : "FAIL"}, recovered: ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}, erased: ${dec.stats.oligosErased}`);
  console.log(`  density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
}

async function main() {
  console.log("=== v61 Honest Benchmark ===\n");
  await bench("v55-density 256KB 10x", ULTIMATE_V55_DENSITY_CONFIG, 256 * 1024, 10);
  console.log();
  await bench("v55-density 2.1MB 10x", ULTIMATE_V55_DENSITY_CONFIG, 2.1 * 1024 * 1024, 10);
}
main().catch(e => { console.error(e); process.exit(1); });
