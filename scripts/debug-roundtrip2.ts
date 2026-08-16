/**
 * Test with different lowCoverageTrigger values and mapping modes.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import type { CodecConfig } from "../src/lib/dna/types";
import { readFile } from "fs/promises";
import { join } from "path";

const TEST_DATA_DIR = join(import.meta.dirname ?? ".", "..", "test-data");

async function testRoundtrip(data: Uint8Array, cfg: CodecConfig, label: string) {
  try {
    const enc = await encodeFile(data, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const status = dec.hashMatches ? "✓" : "✗";
    const strats = dec.perOligo.slice(0, 5).map(p => p.strategy ?? "?").join(",");
    console.log(`${status} ${label} | enc ${enc.stats.encodeTimeMs}ms | dec ${dec.stats.decodeTimeMs}ms | ${enc.stats.oligoCount} oligos | rec ${dec.stats.oligosRecovered} | erased ${dec.stats.oligosErased} | failInner ${dec.stats.oligosFailedInnerRS} | failOuter ${dec.stats.oligosFailedOuterRS} | ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt | strats: ${strats}`);
  } catch (e: any) {
    console.log(`✗ ${label} | ERROR: ${e.message?.slice(0, 100)}`);
  }
}

async function main() {
  const buf = await readFile(join(TEST_DATA_DIR, "random_1kb.bin"));
  const data = new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
  console.log(`Input: ${data.length} bytes\n`);

  const base = {
    innerCode: "ldpc" as const,
    ldpcDecoder: "auto" as const,
    innerParityBytes: 4,
    outerParityRatio: 0.1,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 1,
    interleaveDepth: 0,
    channel: "illumina" as const,
  };

  // Test different configs
  await testRoundtrip(data, { ...base, oligoLength: 300, primerLength: 12, mappingMode: "constrained", lowCoverageTrigger: 5 }, "300nt/constrained/lc5");
  await testRoundtrip(data, { ...base, oligoLength: 300, primerLength: 12, mappingMode: "constrained", lowCoverageTrigger: 999 }, "300nt/constrained/lc999");
  await testRoundtrip(data, { ...base, oligoLength: 300, primerLength: 12, mappingMode: "constrained", lowCoverageTrigger: 0 }, "300nt/constrained/lc0");
  await testRoundtrip(data, { ...base, oligoLength: 300, primerLength: 12, mappingMode: "goldman", goldmanMode: "fast", lowCoverageTrigger: 5 }, "300nt/goldman-fast/lc5");
  await testRoundtrip(data, { ...base, oligoLength: 200, primerLength: 20, mappingMode: "goldman", goldmanMode: "dense", lowCoverageTrigger: 5 }, "200nt/goldman-dense/lc5");
  await testRoundtrip(data, { ...base, oligoLength: 300, primerLength: 12, mappingMode: "srt", lowCoverageTrigger: 5 }, "300nt/srt/lc5");

  // v55 density config variations
  await testRoundtrip(data, { ...base, oligoLength: 700, primerLength: 12, mappingMode: "constrained", innerParityBytes: 8, outerParityRatio: 0.03, lowCoverageTrigger: 5 }, "700nt/constrained/lc5");
  await testRoundtrip(data, { ...base, oligoLength: 700, primerLength: 12, mappingMode: "constrained", innerParityBytes: 8, outerParityRatio: 0.03, lowCoverageTrigger: 999 }, "700nt/constrained/lc999");

  // Nanopore config with corrected oligo length (128 divisible by 4)
  await testRoundtrip(data, {
    oligoLength: 152, primerLength: 12,
    innerCode: "ldpc", ldpcDecoder: "auto", mappingMode: "constrained",
    innerParityBytes: 8, outerParityRatio: 0.4,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true, maxRetries: 1, interleaveDepth: 0,
    channel: "nanopore", lowCoverageTrigger: 5,
    useConvolutionalInner: true,
  }, "152nt/nanopore-conv/lc5");
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
