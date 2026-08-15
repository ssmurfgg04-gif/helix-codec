/**
 * v60: Quick regression test — encode + decode with Illumina channel.
 * Verifies the HMM + indel Viterbi changes didn't break the basic path.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { DEFAULT_CONFIG } from "../src/lib/dna/codec";

async function testConfig(name: string, cfg: any, payloadSize: number, coverage: number) {
  const payload = new Uint8Array(payloadSize);
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;
  
  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage, seed: 42 });
  
  const t1 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;
  
  console.log(`${name}:`);
  console.log(`  ${cfg.oligoLength}nt, ${enc.encoded.oligos.length} oligos, ${coverage}x cov`);
  console.log(`  encode: ${(encMs/1000).toFixed(2)}s, decode: ${(decMs/1000).toFixed(2)}s`);
  console.log(`  hash: ${dec.hashMatches ? "OK" : "FAIL"}, recovered: ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
  console.log(`  density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
  console.log();
}

async function main() {
  console.log("=== v60 Regression Test ===\n");
  
  // Test 1: Default config (RS inner, direct mapping, illumina)
  await testConfig("DEFAULT_CONFIG", DEFAULT_CONFIG, 256 * 1024, 10);
  
  // Test 2: v55-density (LDPC inner, direct mapping, illumina)
  await testConfig("ULTIMATE_V55_DENSITY", ULTIMATE_V55_DENSITY_CONFIG, 256 * 1024, 10);
  
  // Test 3: v55-density at 2x coverage
  await testConfig("V55 @ 2x cov", ULTIMATE_V55_DENSITY_CONFIG, 256 * 1024, 2);
  
  // Test 4: v55-density at 3x coverage
  await testConfig("V55 @ 3x cov", ULTIMATE_V55_DENSITY_CONFIG, 256 * 1024, 3);
}

main().catch(e => { console.error(e); process.exit(1); });
