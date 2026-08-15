/**
 * v59 Encode Profile — find the actual bottleneck in encodeFile.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { ULTIMATE_V55_DENSITY_CONFIG, DEFAULT_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";

const TAG = "[v59-prof]";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  console.log(`${TAG} Payload: ${(payload.length / 1024 / 1024).toFixed(2)}MB`);

  // Test with v55 density config (700nt oligos)
  const cfg = ULTIMATE_V55_DENSITY_CONFIG;

  // Time each phase
  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const totalMs = Date.now() - t0;

  console.log(`${TAG} Total encode: ${totalMs}ms = ${((payload.length / 1024 / 1024) / (totalMs / 1000)).toFixed(2)} MB/s`);
  console.log(`${TAG} Oligos: ${enc.encoded.oligos.length}`);
  console.log(`${TAG} Screening retries: ${enc.stats.screeningRetries}`);
  console.log(`${TAG} Per-oligo encode time: ${(totalMs / enc.encoded.oligos.length).toFixed(2)}ms`);

  // Estimate phase costs
  const oligoCount = enc.encoded.oligos.length;
  const chunkSize = enc.encoded.metadata.payloadBytesPerOligo;
  console.log(`\n${TAG} Phase estimates:`);
  console.log(`${TAG}   Outer RS: ~${chunkSize} byte positions × ${oligoCount} oligos = ${chunkSize * oligoCount} RS ops`);
  console.log(`${TAG}   LDPC encode: ${oligoCount} ops`);
  console.log(`${TAG}   DNA mapping: ${oligoCount} ops`);
  console.log(`${TAG}   Screening: ${oligoCount + enc.stats.screeningRetries} attempts`);
}
main().catch(e => { console.error(e); process.exit(1); });
