// Test: Deletion tolerance with MSA-integrated decode.
// This should show improvement from 0.1% → 5%+ after the progressive MSA integration.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";

async function main() {
  console.log("=== Deletion Tolerance Test (Post-MSA Integration) ===\n");

  const payload = new TextEncoder().encode(
    "The quick brown fox jumps over the lazy dog. DNA storage is the future. ".repeat(10),
  );
  console.log(`Payload: ${payload.length} bytes\n`);

  const enc = await encodeFile(payload, DEFAULT_CONFIG, {
    fileName: "del_test.bin",
    contentType: "application/octet-stream",
  });
  console.log(`Encoded: ${enc.stats.oligoCount} oligos\n`);

  // Test deletion rates from 0.1% to 10%
  const delRates = [0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20];
  const coverage = 20;

  console.log("Del%   | Coverage | Reads  | Errors  | Recovery | Inner Fail | Time");
  console.log("-------|----------|--------|---------|----------|------------|-----");

  for (const delRate of delRates) {
    const sim = simulate(enc.encoded.oligos, {
      substitutionRate: 0.001, // small sub rate
      insertionRate: 0.0005,
      deletionRate: delRate,
      coverage,
      dropoutRate: 0,
      seed: 42,
    });

    const t0 = Date.now();
    let result: any;
    let error: string | undefined;
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
      error = (e as Error).message;
      result = { hashMatches: false, stats: { oligosRecovered: 0, oligosFailedInnerRS: 0 }, data: null };
    }
    const ms = Date.now() - t0;

    console.log(
      `${(delRate * 100).toFixed(1)}%  | ${coverage}x      | ${sim.totalReads.toString().padStart(6)} | ${sim.totalErrors.toString().padStart(7)} | ${result.hashMatches ? "✅ PASS" : "❌ FAIL"}  | ${result.stats?.oligosFailedInnerRS ?? 0}          | ${ms}ms ${error ? "CRASH" : ""}`,
    );
  }

  // Also test with higher coverage for high deletion rates
  console.log("\n--- Higher Coverage for High Deletion Rates ---\n");

  // Use higher outer parity for these tests
  const highParityConfig = { ...DEFAULT_CONFIG, outerParityRatio: 0.5 };
  const encHP = await encodeFile(payload, highParityConfig, {
    fileName: "del_test_hp.bin",
    contentType: "application/octet-stream",
  });
  console.log(`  High-parity encoded: ${encHP.stats.oligoCount} oligos (50% outer parity)\n`);

  for (const delRate of [0.03, 0.05, 0.08, 0.10]) {
    for (const cov of [30, 50]) {
      const sim = simulate(encHP.encoded.oligos, {
        substitutionRate: 0.001,
        insertionRate: 0.0005,
        deletionRate: delRate,
        coverage: cov,
        dropoutRate: 0,
        seed: 42,
      });

      const t0 = Date.now();
      let result: any;
      try {
        result = await decodeReads(
          sim.reads,
          encHP.encoded.metadata,
          highParityConfig,
          encHP.encoded.forwardPrimer,
          encHP.encoded.reversePrimer,
          true,
        );
      } catch (e) {
        result = { hashMatches: false, stats: { oligosRecovered: 0, oligosFailedInnerRS: 0 }, data: null };
      }
      const ms = Date.now() - t0;

      console.log(
        `del=${(delRate * 100).toFixed(0)}% cov=${cov}x parity=50%: ${result.hashMatches ? "✅ PASS" : "❌ FAIL"} (${result.stats?.oligosRecovered ?? 0} recovered, ${result.stats?.oligosFailedInnerRS ?? 0} inner fail, ${ms}ms)`,
      );
    }
  }

  // True empirical density with incompressible data
  console.log("\n=== True Empirical Density (Incompressible Data) ===\n");
  {
    const randomData = new Uint8Array(1024);
    for (let i = 0; i < 1024; i++) randomData[i] = Math.floor(Math.random() * 256);

    const encRand = await encodeFile(randomData, DEFAULT_CONFIG, {
      fileName: "random.bin",
      contentType: "application/octet-stream",
    });

    let totalNt = 0;
    for (const oligo of encRand.encoded.oligos) totalNt += oligo.sequence.length;

    const density = (randomData.length * 8) / totalNt;
    console.log(`  Payload: ${randomData.length} bytes (incompressible random)`);
    console.log(`  Total nucleotides: ${totalNt.toLocaleString()}`);
    console.log(`  Empirical density: ${density.toFixed(4)} bits/nt`);
    console.log(`  Oligos: ${encRand.encoded.oligos.length}`);
    console.log(`  Compressed: ${encRand.stats.compressedSize} bytes (ratio: ${(encRand.stats.compressedSize / randomData.length * 100).toFixed(0)}%)`);
  }

  // Scale test: 1MB
  console.log("\n=== Scale Test: 1MB ===\n");
  {
    const bigPayload = new Uint8Array(1048576); // 1MB
    for (let i = 0; i < bigPayload.length; i++) bigPayload[i] = (i * 31 + 17) & 0xff;

    const t0 = Date.now();
    const encBig = await encodeFile(bigPayload, DEFAULT_CONFIG, {
      fileName: "1mb.bin",
      contentType: "application/octet-stream",
    });
    const encodeMs = Date.now() - t0;

    let totalNt = 0;
    for (const oligo of encBig.encoded.oligos) totalNt += oligo.sequence.length;

    console.log(`  Payload: ${bigPayload.length.toLocaleString()} bytes (1 MB)`);
    console.log(`  Encoded: ${encBig.encoded.oligos.length.toLocaleString()} oligos`);
    console.log(`  Total nucleotides: ${totalNt.toLocaleString()}`);
    console.log(`  Encode time: ${encodeMs}ms (${(bigPayload.length / 1_000_000 / (encodeMs / 1000)).toFixed(2)} MB/s)`);
    console.log(`  Density: ${((bigPayload.length * 8) / totalNt).toFixed(4)} bits/nt`);

    // Simulate + decode at 10x coverage
    const sim = simulate(encBig.encoded.oligos, {
      substitutionRate: 0.001,
      insertionRate: 0.0005,
      deletionRate: 0.001,
      coverage: 10,
      dropoutRate: 0,
      seed: 42,
    });

    const t1 = Date.now();
    const dec = await decodeReads(
      sim.reads,
      encBig.encoded.metadata,
      DEFAULT_CONFIG,
      encBig.encoded.forwardPrimer,
      encBig.encoded.reversePrimer,
      true,
    );
    const decodeMs = Date.now() - t1;

    console.log(`  Decode time: ${decodeMs}ms (${(bigPayload.length / 1_000_000 / (decodeMs / 1000)).toFixed(2)} MB/s)`);
    console.log(`  Recovery: ${dec.hashMatches ? "✅ PASS" : "❌ FAIL"}`);
    console.log(`  Oligos recovered: ${dec.stats.oligosRecovered}/${encBig.encoded.oligos.length}`);
    console.log(`  Memory: ${(process.memoryUsage().heapUsed / 1024 / 1024).toFixed(0)} MB heap`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
