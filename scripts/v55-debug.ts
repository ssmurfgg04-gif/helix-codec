/**
 * v59 v55-density debug — test with different outer RS parity.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";

const TAG = "[v55-debug]";

async function testConfig(label: string, cfgOverrides: any) {
  const cfg = { ...ULTIMATE_V55_DENSITY_CONFIG, ...cfgOverrides };
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 262144);

  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);

  console.log(`${TAG} ${label}:`);
  console.log(`${TAG}   oligos: ${enc.encoded.oligos.length}, density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
  console.log(`${TAG}   outerRS: n=${enc.encoded.metadata.outerRS.n}, k=${enc.encoded.metadata.outerRS.k}`);
  console.log(`${TAG}   hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
  console.log();
}

async function main() {
  console.log(`${TAG} === v55-density hash FAIL investigation ===\n`);

  // Original config (3% RS)
  await testConfig("v55-density (3% RS)", {});

  // With 5% RS
  await testConfig("v55-density (5% RS)", { outerParityRatio: 0.05 });

  // With 10% RS
  await testConfig("v55-density (10% RS)", { outerParityRatio: 0.10 });

  // With 8B inner parity (instead of 4B)
  await testConfig("v55-density (4B→8B parity)", { innerParityBytes: 8 });

  // With both stronger
  await testConfig("v55-density (8B parity + 10% RS)", { innerParityBytes: 8, outerParityRatio: 0.10 });
}
main().catch(e => { console.error(e); process.exit(1); });
