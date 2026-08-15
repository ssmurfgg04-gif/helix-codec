// Run the full adversarial fuzzing + scale testing + empirical density.
import { runAdversarialFuzz, measureEmpiricalDensity, runScaleTest } from "../src/lib/dna/adversarial-fuzz";
import * as fs from "fs";

async function main() {
  console.log("=".repeat(70));
  console.log("HELIX CODEC v13.0 — RIGOROUS TESTING & HARDENING");
  console.log("=".repeat(70));

  // 1. Adversarial Fuzzing
  console.log("\n\n");
  const fuzzPayload = new TextEncoder().encode(
    "The quick brown fox jumps over the lazy dog. ".repeat(30),
  );
  const fuzzResult = await runAdversarialFuzz(fuzzPayload);

  console.log("\n--- Cliff Edge Analysis ---\n");
  for (const c of fuzzResult.cliffEdge) {
    console.log(`  ${c.metric}: breaks at ${c.breakingPoint}, 50% recovery: ${c.recoveryAt50pct}`);
  }

  fs.writeFileSync("benchmarks/adversarial_fuzz_results.md", fuzzResult.markdownTable);
  console.log("\n✅ Fuzz results saved to benchmarks/adversarial_fuzz_results.md");

  // 2. Empirical Density (true, measured)
  console.log("\n\n=== True Empirical Density Measurement ===\n");
  const densityPayload = new Uint8Array(1024);
  for (let i = 0; i < 1024; i++) densityPayload[i] = (i * 31 + 17) & 0xff;

  const density = await measureEmpiricalDensity(densityPayload);
  console.log(`  Payload: ${density.payloadBytes} bytes (${density.payloadBits} bits)`);
  console.log(`  Total nucleotides generated: ${density.totalNucleotides.toLocaleString()}`);
  console.log(`  Empirical density: ${density.empiricalDensity.toFixed(4)} bits/nt`);
  console.log(`  Theoretical density: ${density.theoreticalDensity.toFixed(4)} bits/nt`);
  console.log(`  Gap: ${((density.theoreticalDensity - density.empiricalDensity) / density.theoreticalDensity * 100).toFixed(1)}% (padding + parity oligos)`);
  console.log("\n  Overhead breakdown:");
  for (const o of density.overheadBreakdown) {
    console.log(`    ${o.component}: ${o.nucleotides.toLocaleString()} nt (${o.percentage.toFixed(1)}%)`);
  }

  // 3. Scale Test
  console.log("\n\n");
  const scaleResult = await runScaleTest([256, 1024, 4096, 16384, 65536]);
  fs.writeFileSync("benchmarks/scale_test_results.md", scaleResult.markdownTable);
  console.log("\n✅ Scale test results saved to benchmarks/scale_test_results.md");

  // 4. Summary
  console.log("\n\n" + "=".repeat(70));
  console.log("SUMMARY");
  console.log("=".repeat(70));

  const passCount = fuzzResult.results.filter(r => r.hashMatch).length;
  const failCount = fuzzResult.results.filter(r => !r.hashMatch).length;
  const crashCount = fuzzResult.results.filter(r => r.error).length;

  console.log(`\nAdversarial Fuzzing: ${fuzzResult.results.length} tests`);
  console.log(`  ✅ Pass: ${passCount}`);
  console.log(`  ❌ Fail: ${failCount} (expected at high error rates)`);
  console.log(`  💥 Crash: ${crashCount}`);

  console.log(`\nEmpirical Density: ${density.empiricalDensity.toFixed(4)} bits/nt`);
  console.log(`  (theoretical: ${density.theoreticalDensity.toFixed(4)}, gap: ${((density.theoreticalDensity - density.empiricalDensity) / density.theoreticalDensity * 100).toFixed(1)}%)`);

  console.log(`\nScale Test: ${scaleResult.results.length} sizes tested`);
  const maxScale = scaleResult.results[scaleResult.results.length - 1];
  console.log(`  Largest: ${maxScale.size.toLocaleString()} bytes → ${maxScale.oligos} oligos`);
  console.log(`  Encode: ${maxScale.encodeMBps.toFixed(2)} MB/s, Decode: ${maxScale.decodeMBps.toFixed(2)} MB/s`);
  console.log(`  Errors: ${scaleResult.results.filter(r => r.error).length}`);

  // Save full report
  const fullReport = [
    "# Helix Codec v13.0 — Rigorous Testing Report",
    "",
    `**Date:** ${new Date().toISOString()}`,
    "",
    "## 1. Adversarial Fuzzing",
    "",
    fuzzResult.markdownTable,
    "",
    "## 2. Empirical Density",
    "",
    `| Metric | Value |`,
    `|--------|-------|`,
    `| Payload | ${density.payloadBytes} bytes (${density.payloadBits} bits) |`,
    `| Total nucleotides | ${density.totalNucleotides.toLocaleString()} |`,
    `| **Empirical density** | **${density.empiricalDensity.toFixed(4)} bits/nt** |`,
    `| Theoretical density | ${density.theoreticalDensity.toFixed(4)} bits/nt |`,
    "",
    "### Overhead Breakdown",
    "",
    "| Component | Nucleotides | % of total |",
    "|-----------|-------------|-----------|",
    ...density.overheadBreakdown.map(o => `| ${o.component} | ${o.nucleotides.toLocaleString()} | ${o.percentage.toFixed(1)}% |`),
    "",
    "## 3. Scale Test",
    "",
    scaleResult.markdownTable,
    "",
    "## 4. Conclusions",
    "",
    `1. **Adversarial fuzzing**: ${passCount}/${fuzzResult.results.length} tests passed. ${crashCount} crashes.`,
    `2. **Empirical density**: ${density.empiricalDensity.toFixed(4)} bits/nt (measured, not theoretical).`,
    `3. **Scale**: tested up to ${maxScale.size.toLocaleString()} bytes at ${maxScale.encodeMBps.toFixed(2)} MB/s encode.`,
    `4. **Cliff edges**: documented for substitution, deletion, and dropout rates.`,
  ].join("\n");

  fs.writeFileSync("benchmarks/v13_testing_report.md", fullReport);
  console.log("\n✅ Full report saved to benchmarks/v13_testing_report.md");
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
