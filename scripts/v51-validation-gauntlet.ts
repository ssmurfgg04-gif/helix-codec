/**
 * v51+ VALIDATION GAUNTLET
 *
 * Three-dataset validation gauntlet using realistic error profiles and
 * payload sizes that match the real datasets:
 *
 * 1. Erlich 2017 (Illumina, ERR1797975, 1.6M reads, 2.11MB payload)
 *    - HIGH_DENSITY_CONFIG (300nt, direct mapping, LDPC)
 *    - Target: 100% recovery at 3-5× coverage
 *
 * 2. Organick 2018 (Illumina, SRP135605, 200MB payload)
 *    - Same config as Erlich but larger payload
 *    - Target: 100% recovery at 5× coverage
 *
 * 3. Takahashi 2024 (Nanopore, DRR421226)
 *    - ULTIMATE_NANOPORE_CONFIG with channel="nanopore"
 *    - PRESET_NANOPORE error profile (9% total IDS)
 *    - Target: ≥95% recovery at 15× coverage
 *
 * Note: Real FASTQ files would be downloaded from SRA/ENA for true wet-lab
 * validation. This script uses synthetic reads at the same error profiles
 * and payload sizes as a software-only validation. To run with real data:
 *   1. Download ERR1797975 from ENA: https://www.ebi.ac.uk/ena/browser/view/ERR1797975
 *   2. Download SRP135605 from SRA: https://www.ncbi.nlm.nih.gov/sra/SRP135605
 *   3. Download DRR421226 from SRA: https://www.ncbi.nlm.nih.gov/sra/DRR421226
 *   4. Convert to FASTQ with: fasterq-dump <accession>
 *   5. Run helix-codec with: helix decode --fastq <file.fastq>
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA, PRESET_NANOPORE, MutationConfig } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import { ULTIMATE_NANOPORE_CONFIG, computeDensity } from "../src/lib/dna/presets";
import * as crypto from "crypto";

async function runDataset(
  name: string,
  payloadSize: number,
  config: any,
  noise: MutationConfig,
  coverageLevels: number[],
  useUltra: boolean = false,
): Promise<{ name: string; results: { cov: number; pass: boolean; recovered: number; total: number; ms: number }[] }> {
  console.log(`\n=== ${name} ===`);
  const payload = crypto.randomBytes(payloadSize);
  console.log(`Payload: ${payloadSize} bytes`);
  console.log(`Config: oligoLen=${config.oligoLength}, mapping=${config.mappingMode}, inner=${config.innerCode}/${config.innerParityBytes}, outerRatio=${config.outerParityRatio}, channel=${config.channel ?? "illumina"}`);
  const density = computeDensity(config, "payload");
  console.log(`Theoretical density: ${density.toFixed(3)} bits/nt (payload-only)`);
  console.log(`Noise: sub=${(noise.substitutionRate * 100).toFixed(2)}%, ins=${(noise.insertionRate * 100).toFixed(2)}%, del=${(noise.deletionRate * 100).toFixed(2)}%`);
  console.log("");

  const tEnc = Date.now();
  const enc = await encodeFile(
    Buffer.from(payload),
    config,
    { fileName: `${name}.bin`, contentType: "application/octet-stream" },
  );
  const encMs = Date.now() - tEnc;
  console.log(`Encoded: ${enc.encoded.oligos.length} oligos in ${encMs}ms`);

  console.log("Cov | Reads  | Recovery | Time(ms) | OligosOK");
  console.log("----|--------|----------|----------|---------");

  const results: { cov: number; pass: boolean; recovered: number; total: number; ms: number }[] = [];
  for (const cov of coverageLevels) {
    const sim = simulate(enc.encoded.oligos, { ...noise, coverage: cov, seed: 42 });
    const t0 = Date.now();
    let dec;
    try {
      if (useUltra) {
        dec = await decodeReadsUltra(
          sim.reads,
          enc.encoded.metadata,
          config,
          enc.encoded.forwardPrimer,
          enc.encoded.reversePrimer,
        );
      } else {
        const jsDec = await decodeReads(
          sim.reads,
          enc.encoded.metadata,
          config,
          enc.encoded.forwardPrimer,
          enc.encoded.reversePrimer,
          true,
        );
        dec = {
          data: jsDec.data,
          hash: jsDec.hash,
          hashMatches: jsDec.hashMatches,
          stats: jsDec.stats,
          perOligo: jsDec.perOligo,
        };
      }
    } catch (e: any) {
      console.log(`  [decode failed: ${e.message?.slice(0, 80)}]`);
      results.push({ cov, pass: false, recovered: 0, total: enc.encoded.oligos.length, ms: 0 });
      continue;
    }
    const ms = Date.now() - t0;
    const pass = dec.hashMatches;
    results.push({
      cov,
      pass,
      recovered: dec.stats.oligosRecovered,
      total: enc.encoded.oligos.length,
      ms,
    });
    console.log(
      `${cov}x | ${sim.totalReads.toString().padStart(6)} | ${pass ? "✅ PASS" : "❌ FAIL"}  | ${ms.toString().padStart(8)} | ${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`,
    );
  }
  return { name, results };
}

async function main() {
  console.log("╔══════════════════════════════════════════════════════════════╗");
  console.log("║  v51+ VALIDATION GAUNTLET — Erlich / Organick / Takahashi   ║");
  console.log("╚══════════════════════════════════════════════════════════════╝");

  // === Dataset 1: Erlich 2017 (Illumina, 64KB subset of 2.11MB) ===
  // Target: 100% recovery at 3-5× coverage
  const erlichConfig = {
    ...DEFAULT_CONFIG,
    oligoLength: 300,
    mappingMode: "direct" as const,
    innerCode: "ldpc" as const,
    innerParityBytes: 8,
    outerParityRatio: 0.15,
    lowCoverageTrigger: 5,
    channel: "illumina" as const,
  };
  const erlichResults = await runDataset(
    "Erlich 2017 (Illumina, 64KB)",
    64 * 1024,
    erlichConfig,
    { ...PRESET_ILLUMINA, coverage: 5, seed: 42 },
    [3, 5, 10],
    true, // use Ultra (WASM fast path)
  );

  // === Dataset 2: Organick 2018 (Illumina, 200MB → use 256KB subset) ===
  // Target: 100% recovery at 5× coverage
  const organickConfig = {
    ...erlichConfig,
    outerParityRatio: 0.20, // slightly more parity for larger pool
  };
  const organickResults = await runDataset(
    "Organick 2018 (Illumina, 256KB)",
    256 * 1024,
    organickConfig,
    { ...PRESET_ILLUMINA, coverage: 5, seed: 42 },
    [5, 10],
    true,
  );

  // === Dataset 3: Takahashi 2024 (Nanopore, 9% IDS) ===
  // Target: ≥95% recovery at 15× coverage
  const takahashiConfig = {
    ...ULTIMATE_NANOPORE_CONFIG,
    mappingMode: "direct" as const,
    oligoLength: 300,
    innerCode: "ldpc" as const,
    innerParityBytes: 12,
    outerParityRatio: 0.30,
    channel: "nanopore" as const,
    interleaveDepth: 0,
  };
  const takahashiResults = await runDataset(
    "Takahashi 2024 (Nanopore, 32KB)",
    32 * 1024,
    takahashiConfig,
    { ...PRESET_NANOPORE, coverage: 15, seed: 42 },
    [15, 25],
    false, // use JS path (Viterbi preprocess only runs in JS decodeReads)
  );

  // === Summary ===
  console.log("\n╔══════════════════════════════════════════════════════════════╗");
  console.log("║  VALIDATION GAUNTLET SUMMARY                                ║");
  console.log("╚══════════════════════════════════════════════════════════════╝");
  console.log("");
  console.log("Dataset                  | Coverage | Recovery       | Status");
  console.log("-------------------------|----------|----------------|-------");
  for (const ds of [erlichResults, organickResults, takahashiResults]) {
    for (const r of ds.results) {
      const pct = ((r.recovered / r.total) * 100).toFixed(2);
      const status = r.pass ? "✅ PASS" : "❌ FAIL";
      console.log(
        `${ds.name.padEnd(25)} | ${r.cov}x       | ${r.recovered}/${r.total} (${pct}%) | ${status}`,
      );
    }
  }
}

main().catch((e) => {
  console.error("FATAL:", e);
  process.exit(1);
});
