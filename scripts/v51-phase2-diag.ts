/**
 * Diagnostic: low-coverage recovery at 4x with detailed stats.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, MutationConfig } from "../src/lib/dna/simulate";
import { ULTIMATE_LOW_COVERAGE_CONFIG } from "../src/lib/dna/presets";
import * as crypto from "crypto";

async function main() {
  const payload = crypto.randomBytes(64 * 1024);
  const config = {
    ...ULTIMATE_LOW_COVERAGE_CONFIG,
    mappingMode: "direct" as const,
    oligoLength: 300,
    innerParityBytes: 8,
    outerParityRatio: 0.25,
  };

  const enc = await encodeFile(
    Buffer.from(payload),
    config,
    { fileName: "test.bin", contentType: "application/octet-stream" },
  );

  console.log(`Encoded: ${enc.encoded.oligos.length} oligos`);
  console.log(`Outer RS: n=${enc.encoded.metadata.outerRS.n}, k=${enc.encoded.metadata.outerRS.k}`);
  console.log(`  → Can correct up to ${enc.encoded.metadata.outerRS.n - enc.encoded.metadata.outerRS.k} erasures`);

  for (const cov of [3, 4, 5]) {
    const baseNoise: MutationConfig = {
      substitutionRate: 0.001,
      insertionRate: 0.0005,
      deletionRate: 0.001,
      dropoutRate: 0.0,
      coverage: cov,
      seed: 42,
    };
    const sim = simulate(enc.encoded.oligos, baseNoise);
    const dec = await decodeReads(
      sim.reads,
      enc.encoded.metadata,
      config,
      enc.encoded.forwardPrimer,
      enc.encoded.reversePrimer,
      true,
    );
    console.log(`\n${cov}x coverage:`);
    console.log(`  Hash matches: ${dec.hashMatches}`);
    console.log(`  Stats:`, JSON.stringify(dec.stats, null, 2));
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
