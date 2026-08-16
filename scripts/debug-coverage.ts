/**
 * Test coverage effect on round-trip.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { readFile } from "fs/promises";
import { join } from "path";

const TEST_DATA_DIR = join(import.meta.dirname ?? ".", "..", "test-data");

async function main() {
  const buf = await readFile(join(TEST_DATA_DIR, "random_100kb.bin"));
  const data = new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
  console.log(`Input: ${data.length} bytes\n`);

  const cfg = {
    oligoLength: 300, primerLength: 12,
    innerCode: "ldpc" as const, ldpcDecoder: "auto" as const,
    mappingMode: "constrained" as const,
    innerParityBytes: 8, outerParityRatio: 0.1,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true, maxRetries: 1, interleaveDepth: 0,
    channel: "illumina" as const, lowCoverageTrigger: 5,
  };

  const enc = await encodeFile(data, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${enc.stats.oligoCount} oligos, density ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt\n`);

  for (const coverage of [1, 5, 10, 20]) {
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage, simulator: "basic" });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    console.log(`coverage=${coverage}: hash=${dec.hashMatches} rec=${dec.stats.oligosRecovered}/${enc.stats.oligoCount} erased=${dec.stats.oligosErased} failInner=${dec.stats.oligosFailedInnerRS} failOuter=${dec.stats.oligosFailedOuterRS} decMs=${dec.stats.decodeTimeMs}`);
  }

  // Also try with lowCoverageTrigger=999 (WASM path)
  console.log("\nWith lowCoverageTrigger=999:");
  const cfg2 = { ...cfg, lowCoverageTrigger: 999 };
  const enc2 = await encodeFile(data, cfg2, { fileName: "test.bin", contentType: "application/octet-stream" });
  for (const coverage of [1, 5, 10, 20]) {
    const sim = simulate(enc2.encoded.oligos, { ...PRESET_CLEAN, coverage, simulator: "basic" });
    const dec = await decodeReads(sim.reads, enc2.encoded.metadata, cfg2, enc2.encoded.forwardPrimer, enc2.encoded.reversePrimer);
    console.log(`coverage=${coverage}: hash=${dec.hashMatches} rec=${dec.stats.oligosRecovered}/${enc2.stats.oligoCount} erased=${dec.stats.oligosErased} failInner=${dec.stats.oligosFailedInnerRS} failOuter=${dec.stats.oligosFailedOuterRS}`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
