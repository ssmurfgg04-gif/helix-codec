/**
 * v60: Test nanopore 9% IDS WITHOUT conv inner code.
 * The MSA consensus should handle indels; LDPC handles residual substitutions.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { ULTIMATE_NANOPORE_V52_CONFIG, ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";

async function testNanoporeNoConv() {
  console.log("=== Test: Nanopore 9% IDS WITHOUT conv inner code ===\n");

  // Use v55-density config (no conv) but set channel to nanopore
  const cfg = {
    ...ULTIMATE_V55_DENSITY_CONFIG,
    channel: "nanopore" as const,
    useConvolutionalInner: false,
    outerParityRatio: 0.25,
  };

  const payload = new Uint8Array(16 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;

  console.log(`Config: ${cfg.oligoLength}nt, channel=${cfg.channel}, conv=${cfg.useConvolutionalInner}`);
  console.log(`Inner parity: ${cfg.innerParityBytes}B, outer RS: ${(cfg.outerParityRatio * 100).toFixed(0)}%`);

  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos, ${((Date.now() - t0) / 1000).toFixed(2)}s`);

  const coverage = 10;
  const sim = simulate(enc.encoded.oligos, { ...PRESET_NANOPORE, coverage, seed: 42 });
  console.log(`Simulated ${sim.reads.length} reads at ${coverage}x coverage, 9% IDS\n`);

  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;

  console.log(`Result:`);
  console.log(`  decode: ${(decMs / 1000).toFixed(2)}s`);
  console.log(`  hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
  console.log(`  oligos recovered: ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
  console.log();
}

async function testNanoporeWithConv() {
  console.log("=== Test: Nanopore 9% IDS WITH conv inner code (v59 baseline) ===\n");

  const cfg = ULTIMATE_NANOPORE_V52_CONFIG;

  const payload = new Uint8Array(16 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;

  console.log(`Config: ${cfg.oligoLength}nt, channel=${cfg.channel}, conv=${cfg.useConvolutionalInner}`);

  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos, ${((Date.now() - t0) / 1000).toFixed(2)}s`);

  const coverage = 10;
  const sim = simulate(enc.encoded.oligos, { ...PRESET_NANOPORE, coverage, seed: 42 });
  console.log(`Simulated ${sim.reads.length} reads at ${coverage}x coverage, 9% IDS\n`);

  const t1 = Date.now();
  const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - t1;

  console.log(`Result:`);
  console.log(`  decode: ${(decMs / 1000).toFixed(2)}s`);
  console.log(`  hash: ${dec.hashMatches ? "OK" : "FAIL"}`);
  console.log(`  oligos recovered: ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);
  console.log();
}

async function main() {
  await testNanoporeWithConv();
  await testNanoporeNoConv();
}

main().catch(e => { console.error(e); process.exit(1); });
