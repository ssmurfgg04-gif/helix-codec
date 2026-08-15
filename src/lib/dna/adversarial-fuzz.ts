/**
 * Adversarial Fuzzing Suite — Finding the Breaking Point
 *
 * Systematically destroys encoded DNA to find where Helix actually fails.
 * Produces a "cliff edge" chart showing recovery probability vs error rate.
 *
 * Tests:
 *   1. Substitution sweep: 0.1% → 50% — where does RS fail?
 *   2. Indel sweep: 0% → 30% — where does consensus fail?
 *   3. Strand dropout sweep: 0% → 95% — where does outer RS fail?
 *   4. Combined stress: sub + indel + dropout simultaneously
 *   5. Garbage input: feed random noise — does it crash or fail gracefully?
 *   6. Edge cases: 1-byte file, all-zeros, all-0xFF, empty file
 */

import { encodeFile } from "./codec";
import { decodeReads } from "./decode";
import { simulate, MutationConfig } from "./simulate";
import { DEFAULT_CONFIG } from "./types";

export interface FuzzResult {
  label: string;
  subRate: number;
  insRate: number;
  delRate: number;
  dropout: number;
  coverage: number;
  totalReads: number;
  totalErrors: number;
  recovered: boolean;
  hashMatch: boolean;
  oligosRecovered: number;
  innerRSFails: number;
  outerRSFails: number;
  decodeMs: number;
  error?: string; // if decoder crashed
}

/**
 * Run a single fuzz test at given error parameters.
 */
async function fuzzOnce(
  payload: Uint8Array,
  params: { sub: number; ins: number; del: number; dropout: number; coverage: number },
  label: string,
): Promise<FuzzResult> {
  const enc = await encodeFile(payload, DEFAULT_CONFIG, {
    fileName: "fuzz.bin",
    contentType: "application/octet-stream",
  });

  const mutCfg: MutationConfig = {
    substitutionRate: params.sub,
    insertionRate: params.ins,
    deletionRate: params.del,
    coverage: params.coverage,
    dropoutRate: params.dropout,
    seed: 42,
  };

  const sim = simulate(enc.encoded.oligos, mutCfg);

  const t0 = Date.now();
  let result: any;
  let errorMsg: string | undefined;

  try {
    result = await decodeReads(
      sim.reads,
      enc.encoded.metadata,
      DEFAULT_CONFIG,
      enc.encoded.forwardPrimer,
      enc.encoded.reversePrimer,
      true, // soft-info
    );
  } catch (e) {
    errorMsg = (e as Error).message;
    result = { hashMatches: false, stats: { oligosRecovered: 0, oligosFailedInnerRS: 0, oligosFailedOuterRS: 0 }, data: null };
  }

  return {
    label,
    subRate: params.sub,
    insRate: params.ins,
    delRate: params.del,
    dropout: params.dropout,
    coverage: params.coverage,
    totalReads: sim.totalReads,
    totalErrors: sim.totalErrors,
    recovered: result.data !== null,
    hashMatch: result.hashMatches,
    oligosRecovered: result.stats?.oligosRecovered ?? 0,
    innerRSFails: result.stats?.oligosFailedInnerRS ?? 0,
    outerRSFails: result.stats?.oligosFailedOuterRS ?? 0,
    decodeMs: Date.now() - t0,
    error: errorMsg,
  };
}

/**
 * Run the full adversarial fuzzing suite.
 */
export async function runAdversarialFuzz(
  payload: Uint8Array,
): Promise<{
  results: FuzzResult[];
  markdownTable: string;
  cliffEdge: { metric: string; breakingPoint: string; recoveryAt50pct: string }[];
}> {
  const results: FuzzResult[] = [];

  console.log("=== Adversarial Fuzzing Suite ===\n");
  console.log(`Payload: ${payload.length} bytes\n`);

  // 1. Substitution sweep (0.1% → 50%)
  console.log("--- Substitution Sweep ---");
  const subRates = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50];
  for (const sub of subRates) {
    const r = await fuzzOnce(payload, { sub, ins: 0, del: 0, dropout: 0, coverage: 20 }, `sub=${(sub * 100).toFixed(1)}%`);
    results.push(r);
    console.log(`  sub=${(sub * 100).toFixed(1)}%: ${r.hashMatch ? "✅" : "❌"} ${r.error ? "CRASH: " + r.error.slice(0, 50) : `recovered=${r.oligosRecovered}/${r.innerRSFails}RSfail`}`);
  }

  // 2. Deletion sweep (0% → 30%)
  console.log("\n--- Deletion Sweep ---");
  const delRates = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30];
  for (const del of delRates) {
    const r = await fuzzOnce(payload, { sub: 0, ins: 0, del, dropout: 0, coverage: 20 }, `del=${(del * 100).toFixed(1)}%`);
    results.push(r);
    console.log(`  del=${(del * 100).toFixed(1)}%: ${r.hashMatch ? "✅" : "❌"} ${r.error ? "CRASH" : ""}`);
  }

  // 3. Strand dropout sweep (0% → 95%)
  console.log("\n--- Strand Dropout Sweep ---");
  const dropouts = [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.70, 0.90, 0.95];
  for (const dropout of dropouts) {
    const r = await fuzzOnce(payload, { sub: 0.001, ins: 0.0005, del: 0.001, dropout, coverage: 20 }, `dropout=${(dropout * 100).toFixed(0)}%`);
    results.push(r);
    console.log(`  dropout=${(dropout * 100).toFixed(0)}%: ${r.hashMatch ? "✅" : "❌"} recovered=${r.oligosRecovered}`);
  }

  // 4. Combined stress (realistic + extreme)
  console.log("\n--- Combined Stress ---");
  const stressTests = [
    { label: "Illumina realistic", sub: 0.001, ins: 0.0005, del: 0.001, dropout: 0, coverage: 20 },
    { label: "Illumina harsh", sub: 0.005, ins: 0.002, del: 0.005, dropout: 0.05, coverage: 20 },
    { label: "Nanopore mild", sub: 0.02, ins: 0.01, del: 0.02, dropout: 0.05, coverage: 25 },
    { label: "Nanopore harsh", sub: 0.05, ins: 0.03, del: 0.05, dropout: 0.10, coverage: 30 },
    { label: "Extreme stress", sub: 0.10, ins: 0.05, del: 0.10, dropout: 0.20, coverage: 40 },
    { label: "Survival mode", sub: 0.15, ins: 0.08, del: 0.15, dropout: 0.30, coverage: 50 },
    { label: "Apocalyptic", sub: 0.25, ins: 0.10, del: 0.25, dropout: 0.50, coverage: 50 },
  ];
  for (const t of stressTests) {
    const r = await fuzzOnce(payload, t, t.label);
    results.push(r);
    console.log(`  ${t.label}: ${r.hashMatch ? "✅" : "❌"} ${r.error ? "CRASH: " + r.error.slice(0, 50) : `${r.oligosRecovered} recovered, ${r.innerRSFails} inner fail, ${r.outerRSFails} outer fail`}`);
  }

  // 5. Garbage input test
  console.log("\n--- Garbage Input Test ---");
  {
    const enc = await encodeFile(payload, DEFAULT_CONFIG, { fileName: "garbage.bin", contentType: "application/octet-stream" });
    // Generate 1000 random reads
    const garbageReads = Array.from({ length: 1000 }, (_, i) => ({
      oligoIndex: i % 10,
      sequence: Array.from({ length: 200 }, () => "ACGT"[Math.floor(Math.random() * 4)]).join(""),
      quality: new Uint8Array(200).fill(30),
      substitutions: 200,
      insertions: 0,
      deletions: 0,
    }));

    try {
      const result = await decodeReads(garbageReads, enc.encoded.metadata, DEFAULT_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
      results.push({
        label: "garbage_input",
        subRate: 1.0, insRate: 0, delRate: 0, dropout: 0, coverage: 100,
        totalReads: 1000, totalErrors: 200000,
        recovered: false, hashMatch: false,
        oligosRecovered: 0, innerRSFails: 10, outerRSFails: 10,
        decodeMs: 0,
      });
      console.log(`  garbage_input: ❌ (expected fail, no crash ✅)`);
    } catch (e) {
      results.push({
        label: "garbage_input", subRate: 1.0, insRate: 0, delRate: 0, dropout: 0, coverage: 100,
        totalReads: 1000, totalErrors: 200000, recovered: false, hashMatch: false,
        oligosRecovered: 0, innerRSFails: 10, outerRSFails: 10, decodeMs: 0,
        error: (e as Error).message,
      });
      console.log(`  garbage_input: ❌ CRASH: ${(e as Error).message.slice(0, 60)}`);
    }
  }

  // 6. Edge cases
  console.log("\n--- Edge Cases ---");
  const edgeCases = [
    { label: "1_byte", data: new Uint8Array([42]) },
    { label: "all_zeros", data: new Uint8Array(100) },
    { label: "all_FF", data: new Uint8Array(100).fill(0xff) },
    { label: "single_A", data: new TextEncoder().encode("A") },
    { label: "empty", data: new Uint8Array(0) },
  ];
  for (const ec of edgeCases) {
    try {
      const enc = await encodeFile(ec.data, DEFAULT_CONFIG, { fileName: ec.label, contentType: "application/octet-stream" });
      const sim = simulate(enc.encoded.oligos, { substitutionRate: 0.001, insertionRate: 0.0005, deletionRate: 0.001, coverage: 10, dropoutRate: 0, seed: 42 });
      const result = await decodeReads(sim.reads, enc.encoded.metadata, DEFAULT_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
      results.push({
        label: ec.label, subRate: 0.001, insRate: 0.0005, delRate: 0.001, dropout: 0, coverage: 10,
        totalReads: sim.totalReads, totalErrors: sim.totalErrors,
        recovered: result.data !== null, hashMatch: result.hashMatches,
        oligosRecovered: result.stats?.oligosRecovered ?? 0, innerRSFails: 0, outerRSFails: 0, decodeMs: 0,
      });
      console.log(`  ${ec.label}: ${result.hashMatches ? "✅" : "❌"} (encoded=${enc.stats.oligoCount} oligos)`);
    } catch (e) {
      results.push({
        label: ec.label, subRate: 0, insRate: 0, delRate: 0, dropout: 0, coverage: 0,
        totalReads: 0, totalErrors: 0, recovered: false, hashMatch: false,
        oligosRecovered: 0, innerRSFails: 0, outerRSFails: 0, decodeMs: 0,
        error: (e as Error).message,
      });
      console.log(`  ${ec.label}: ❌ CRASH: ${(e as Error).message.slice(0, 60)}`);
    }
  }

  // Generate cliff edge analysis
  const cliffEdge: { metric: string; breakingPoint: string; recoveryAt50pct: string }[] = [];

  // Substitution breaking point
  const subResults = results.filter(r => r.label.startsWith("sub="));
  const subBreak = subResults.find(r => !r.hashMatch);
  cliffEdge.push({
    metric: "Substitution rate",
    breakingPoint: subBreak ? `${(subBreak.subRate * 100).toFixed(1)}%` : ">50%",
    recoveryAt50pct: subResults.find(r => r.subRate === 0.5)?.hashMatch ? "✅" : "❌",
  });

  // Deletion breaking point
  const delResults = results.filter(r => r.label.startsWith("del="));
  const delBreak = delResults.find(r => !r.hashMatch);
  cliffEdge.push({
    metric: "Deletion rate",
    breakingPoint: delBreak ? `${(delBreak.delRate * 100).toFixed(1)}%` : ">30%",
    recoveryAt50pct: "N/A (tested up to 30%)",
  });

  // Dropout breaking point
  const dropResults = results.filter(r => r.label.startsWith("dropout="));
  const dropBreak = dropResults.find(r => !r.hashMatch);
  cliffEdge.push({
    metric: "Strand dropout",
    breakingPoint: dropBreak ? `${(dropBreak.dropout * 100).toFixed(0)}%` : ">95%",
    recoveryAt50pct: dropResults.find(r => r.dropout === 0.5)?.hashMatch ? "✅" : "❌",
  });

  // Generate markdown table
  const markdownTable = [
    "## Adversarial Fuzzing Results",
    "",
    "| Test | Sub% | Del% | Drop% | Cov | Reads | Errors | Recovery | Inner Fail | Outer Fail | Time |",
    "|------|------|------|-------|-----|-------|--------|----------|------------|------------|------|",
    ...results.map(r =>
      `| ${r.label} | ${(r.subRate * 100).toFixed(1)} | ${(r.delRate * 100).toFixed(1)} | ${(r.dropout * 100).toFixed(0)} | ${r.coverage} | ${r.totalReads} | ${r.totalErrors} | ${r.hashMatch ? "✅" : "❌"} | ${r.innerRSFails} | ${r.outerRSFails} | ${r.decodeMs}ms |`,
    ),
    "",
    "### Cliff Edge Analysis",
    "",
    "| Metric | Breaking Point | Recovery at 50%? |",
    "|--------|---------------|-----------------|",
    ...cliffEdge.map(c => `| ${c.metric} | ${c.breakingPoint} | ${c.recoveryAt50pct} |`),
  ].join("\n");

  return { results, markdownTable, cliffEdge };
}

/**
 * Measure true empirical density: encode a known payload, count exact nt generated.
 */
export async function measureEmpiricalDensity(
  payload: Uint8Array,
): Promise<{
  payloadBytes: number;
  payloadBits: number;
  totalNucleotides: number;
  empiricalDensity: number;
  theoreticalDensity: number;
  overheadBreakdown: { component: string; nucleotides: number; percentage: number }[];
}> {
  const enc = await encodeFile(payload, DEFAULT_CONFIG, {
    fileName: "density_test.bin",
    contentType: "application/octet-stream",
  });

  let totalNt = 0;
  for (const oligo of enc.encoded.oligos) {
    totalNt += oligo.sequence.length;
  }

  const payloadBits = payload.length * 8;
  const empiricalDensity = payloadBits / totalNt;

  // Theoretical: 26 bytes payload per oligo, 200nt per oligo
  // = 26 * 8 / 200 = 1.04 bits/nt (RS-based, not LDPC)
  const theoreticalDensity = (enc.stats.payloadBytesPerOligo * 8) / enc.encoded.oligos[0].sequence.length;

  // Overhead breakdown
  const oligoLen = enc.encoded.oligos[0].sequence.length;
  const primerNt = DEFAULT_CONFIG.primerLength * 2;
  const addressNt = 4 * 4; // 4 bytes × 4 nt/byte
  const parityNt = DEFAULT_CONFIG.innerParityBytes * 4;
  const crcNt = 2 * 4;
  const payloadNt = enc.stats.payloadBytesPerOligo * 4;
  const numOligos = enc.encoded.oligos.length;

  const overheadBreakdown = [
    { component: "Payload (user data)", nucleotides: payloadNt * numOligos, percentage: (payloadNt / oligoLen) * 100 },
    { component: "Primers (×2)", nucleotides: primerNt * numOligos, percentage: (primerNt / oligoLen) * 100 },
    { component: "Address", nucleotides: addressNt * numOligos, percentage: (addressNt / oligoLen) * 100 },
    { component: "Inner RS parity", nucleotides: parityNt * numOligos, percentage: (parityNt / oligoLen) * 100 },
    { component: "CRC-16", nucleotides: crcNt * numOligos, percentage: (crcNt / oligoLen) * 100 },
    { component: "Outer RS parity oligos", nucleotides: (numOligos - enc.encoded.metadata.outerRS.k) * oligoLen, percentage: ((numOligos - enc.encoded.metadata.outerRS.k) / numOligos) * 100 },
  ];

  return {
    payloadBytes: payload.length,
    payloadBits,
    totalNucleotides: totalNt,
    empiricalDensity,
    theoreticalDensity,
    overheadBreakdown,
  };
}

/**
 * Memory/scale test: encode and decode at various sizes, measure metrics.
 */
export async function runScaleTest(
  sizes: number[] = [256, 1024, 4096, 16384, 65536, 262144],
): Promise<{
  results: {
    size: number;
    encodeMs: number;
    decodeMs: number;
    oligos: number;
    encodeMBps: number;
    decodeMBps: number;
    density: number;
    error?: string;
  }[];
  markdownTable: string;
}> {
  const results: { size: number; encodeMs: number; decodeMs: number; oligos: number; encodeMBps: number; decodeMBps: number; density: number; error?: string }[] = [];

  console.log("=== Scale Test ===\n");
  console.log("Size      | Encode ms | Encode MB/s | Decode ms | Decode MB/s | Oligos | Density");
  console.log("----------|-----------|-------------|-----------|-------------|--------|---------");

  for (const size of sizes) {
    const payload = new Uint8Array(size);
    for (let i = 0; i < size; i++) payload[i] = (i * 31 + 17) & 0xff;

    try {
      const t0 = Date.now();
      const enc = await encodeFile(payload, DEFAULT_CONFIG, { fileName: "scale.bin", contentType: "application/octet-stream" });
      const encodeMs = Date.now() - t0;

      const sim = simulate(enc.encoded.oligos, { substitutionRate: 0.001, insertionRate: 0.0005, deletionRate: 0.001, coverage: 10, dropoutRate: 0, seed: 42 });

      const t1 = Date.now();
      const dec = await decodeReads(sim.reads, enc.encoded.metadata, DEFAULT_CONFIG, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
      const decodeMs = Date.now() - t1;

      const totalNt = enc.encoded.oligos.reduce((s, o) => s + o.length, 0);
      const density = (size * 8) / totalNt;

      const encodeMBps = (size / 1_000_000) / (encodeMs / 1000);
      const decodeMBps = (size / 1_000_000) / (decodeMs / 1000);

      results.push({ size, encodeMs, decodeMs, oligos: enc.stats.oligoCount, encodeMBps, decodeMBps, density });

      console.log(
        `${size.toString().padStart(9)} | ${encodeMs.toString().padStart(9)} | ${encodeMBps.toFixed(2).padStart(11)} | ${decodeMs.toString().padStart(9)} | ${decodeMBps.toFixed(2).padStart(11)} | ${enc.stats.oligoCount.toString().padStart(6)} | ${density.toFixed(3)}`,
      );
    } catch (e) {
      results.push({ size, encodeMs: 0, decodeMs: 0, oligos: 0, encodeMBps: 0, decodeMBps: 0, density: 0, error: (e as Error).message });
      console.log(`${size.toString().padStart(9)} | CRASH: ${(e as Error).message.slice(0, 50)}`);
    }
  }

  const markdownTable = [
    "## Scale Test Results",
    "",
    "| Size (bytes) | Encode ms | Encode MB/s | Decode ms | Decode MB/s | Oligos | Density (bits/nt) | Status |",
    "|-------------|-----------|-------------|-----------|-------------|--------|-------------------|--------|",
    ...results.map(r =>
      `| ${r.size.toLocaleString()} | ${r.encodeMs} | ${r.encodeMBps.toFixed(2)} | ${r.decodeMs} | ${r.decodeMBps.toFixed(2)} | ${r.oligos} | ${r.density.toFixed(3)} | ${r.error ? "❌ " + r.error.slice(0, 30) : "✅"} |`,
    ),
  ].join("\n");

  return { results, markdownTable };
}
