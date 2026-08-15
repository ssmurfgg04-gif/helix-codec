/**
 * v58: Test progressive MSA for Nanopore 9% IDS tolerance.
 * Uses LDPC inner (so MSA path triggers) and decodes via decodeReads.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import * as crypto from "crypto";

async function test(label: string, payload: Uint8Array, cfg: any, noise: any) {
  const enc = await encodeFile(Buffer.from(payload), cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const sim = simulate(enc.encoded.oligos, noise);
  const t0 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
  const ms = Date.now() - t0;
  const pct = ((dec.stats.oligosRecovered / enc.encoded.oligos.length) * 100).toFixed(1);
  console.log(`${label.padEnd(35)} | ${dec.hashMatches?"PASS":"FAIL"} | oligos=${dec.stats.oligosRecovered}/${enc.encoded.oligos.length} (${pct}%) | failInner=${dec.stats.oligosFailedInnerRS} failOuter=${dec.stats.oligosFailedOuterRS} | ${ms}ms`);
}

async function main() {
  console.log("=== v58 Progressive MSA — Nanopore 9% IDS Tolerance ===\n");
  const payload = crypto.randomBytes(8 * 1024); // 8KB
  console.log(`Payload: 8KB random, sub=2% ins=3% del=4% (9% total IDS)\n`);

  // LDPC inner (so MSA path triggers), 8 parity bytes, 30% outer RS
  const cfg = {
    oligoLength: 300,
    primerLength: 20,
    innerCode: "ldpc",
    ldpcDecoder: "auto",
    mappingMode: "direct",
    innerParityBytes: 8,
    outerParityRatio: 0.30,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 1,
    interleaveDepth: 0,
    channel: "nanopore",
    lowCoverageTrigger: 999, // force JS path (so MSA triggers)
    useConvolutionalInner: false,
  };

  console.log("Cov | Test                                | Hash | Oligos                     | Fails          | Time");
  console.log("----|-------------------------------------|------|----------------------------|----------------|-----");

  // Test 1: 9% IDS at 10x coverage
  await test(`9% IDS @ 10x`, payload, cfg, { ...PRESET_NANOPORE, coverage: 10, seed: 42 });

  // Test 2: 9% IDS at 15x
  await test(`9% IDS @ 15x`, payload, cfg, { ...PRESET_NANOPORE, coverage: 15, seed: 42 });

  // Test 3: 9% IDS at 25x
  await test(`9% IDS @ 25x`, payload, cfg, { ...PRESET_NANOPORE, coverage: 25, seed: 42 });

  // Test 4: 5% IDS at 10x (easier)
  await test(`5% IDS @ 10x`, payload, cfg, { ...PRESET_NANOPORE, substitutionRate: 0.01, insertionRate: 0.015, deletionRate: 0.025, coverage: 10, seed: 42 });

  // Test 5: 3% IDS at 10x (easy)
  await test(`3% IDS @ 10x`, payload, cfg, { ...PRESET_NANOPORE, substitutionRate: 0.005, insertionRate: 0.01, deletionRate: 0.015, coverage: 10, seed: 42 });

  // Test 6: pure substitutions at 9% (no indels — should work great with LDPC)
  await test(`9% sub-only @ 10x`, payload, cfg, { ...PRESET_NANOPORE, substitutionRate: 0.09, insertionRate: 0, deletionRate: 0, coverage: 10, seed: 42 });
}
main().catch(e => { console.error(e); process.exit(1); });
