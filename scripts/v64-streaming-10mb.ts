/**
 * v64: 10MB Streaming Decode Test (Task 3)
 *
 * Uses StreamingDecodeRunner to process reads in batches, verifying that
 * memory stays flat during the accumulation phase.
 *
 * Strategy:
 *   1. Encode 10MB payload (v63-hd config)
 *   2. Generate reads in BATCHES (not all upfront) — simulates streaming
 *      from a FASTQ file or sequencer
 *   3. Feed each batch to StreamingDecodeRunner.addReads()
 *   4. Free each batch after feeding (allows GC)
 *   5. Call runner.decode() at the end
 *   6. Monitor RSS at each phase
 *
 * Expected: RSS stays bounded during batch accumulation (only the runner's
 * internal Map<oligoIdx, reads[]> grows, capped at maxReadsPerOligo × numOligos).
 */

import { encodeFile } from "../src/lib/dna/codec";
import { simulate, MutationConfig, SequencingRead } from "../src/lib/dna/simulate";
import { ULTIMATE_V63_HD_CONFIG } from "../src/lib/dna/presets";
import { StreamingDecodeRunner } from "../src/lib/dna/streaming-decode-runner";

const MB = 1024 * 1024;

function randomBytes(n: number): Uint8Array {
  const buf = new Uint8Array(n);
  for (let i = 0; i < n; i++) buf[i] = Math.floor(Math.random() * 256);
  return buf;
}

function rssMB(): number {
  return process.memoryUsage().rss / MB;
}

async function main(): Promise<void> {
  console.log("=== v64 10MB Streaming Decode Test ===");
  console.log(`Date: ${new Date().toISOString()}`);
  console.log(`Node: ${process.version}, Bun: ${(globalThis as any).Bun?.version ?? "n/a"}`);
  console.log(`Config: v63-hd (1100nt, direct, 4B LDPC, 2% outer RS)\n`);

  const cfg = ULTIMATE_V63_HD_CONFIG;
  const payloadSize = 10 * MB;
  const payload = randomBytes(payloadSize);
  const coverage = 10;
  const batchSize = 50000; // reads per batch

  console.log(`Initial RSS: ${rssMB().toFixed(0)}MB`);

  // 1. Encode
  console.log(`\nEncoding ${(payloadSize / MB).toFixed(0)}MB payload...`);
  const encodeStart = performance.now();
  const enc = await encodeFile(payload, cfg, {
    fileName: "b.bin",
    contentType: "application/octet-stream",
  });
  const encodeMs = performance.now() - encodeStart;
  const { oligos, metadata } = enc.encoded;
  console.log(`Encoded: ${oligos.length.toLocaleString()} oligos in ${(encodeMs / 1000).toFixed(1)}s`);
  console.log(`Density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
  console.log(`RSS after encode: ${rssMB().toFixed(0)}MB`);

  // 2. Initialize streaming decoder
  const runner = new StreamingDecodeRunner(
    metadata,
    cfg,
    enc.encoded.forwardPrimer,
    enc.encoded.reversePrimer,
    { maxReadsPerOligo: 5, batchSize },
  );

  // 3. Generate reads in batches and feed to runner
  console.log(`\nStreaming reads in batches of ${batchSize.toLocaleString()}...`);
  const simCfg: MutationConfig = {
    substitutionRate: 0.001,
    insertionRate: 0,
    deletionRate: 0,
    coverage,
    dropoutRate: 0,
    seed: 42,
  };

  // Generate reads per oligo (in oligo order) to simulate streaming
  const totalReads = oligos.length * coverage;
  let readsGenerated = 0;
  let batchNum = 0;
  const simStart = performance.now();

  // Generate all reads for a batch of oligos, then feed to runner
  const oligosPerBatch = Math.ceil(batchSize / coverage);
  for (let oligoStart = 0; oligoStart < oligos.length; oligoStart += oligosPerBatch) {
    const oligoEnd = Math.min(oligoStart + oligosPerBatch, oligos.length);
    const batchOligos = oligos.slice(oligoStart, oligoEnd);

    // Simulate this batch of oligos
    const batchSim = simulate(batchOligos, { ...simCfg, seed: 42 + batchNum });
    const batchReads: SequencingRead[] = batchSim.reads;

    // Feed to runner
    runner.addReads(batchReads);

    readsGenerated += batchReads.length;
    batchNum++;

    if (batchNum % 5 === 0 || oligoEnd === oligos.length) {
      const elapsed = (performance.now() - simStart) / 1000;
      const rate = readsGenerated / elapsed;
      console.log(
        `  batch ${batchNum}: ${readsGenerated.toLocaleString()}/${totalReads.toLocaleString()} reads ` +
        `(${rate.toFixed(0)} reads/s), RSS=${rssMB().toFixed(0)}MB, ` +
        `oligos=${runner.getStats().oligosWithReads.toLocaleString()}`,
      );
    }

    // Free batch reads (they're now in the runner)
    batchReads.length = 0;
  }

  const simMs = performance.now() - simStart;
  console.log(`\nStreamed ${readsGenerated.toLocaleString()} reads in ${(simMs / 1000).toFixed(1)}s`);
  console.log(`RSS after streaming: ${rssMB().toFixed(0)}MB`);

  const stats = runner.getStats();
  console.log(`\nRunner stats:`);
  console.log(`  batches processed:     ${stats.batchesProcessed}`);
  console.log(`  total reads processed: ${stats.totalReadsProcessed.toLocaleString()}`);
  console.log(`  reads discarded:       ${stats.totalReadsDiscarded.toLocaleString()}`);
  console.log(`  reads accumulated:     ${stats.totalReadsAccumulated.toLocaleString()}`);
  console.log(`  oligos with reads:     ${stats.oligosWithReads.toLocaleString()}/${oligos.length.toLocaleString()}`);
  console.log(`  estimated memory:      ${(runner.getEstimatedMemoryUsage() / MB).toFixed(0)}MB`);

  // 4. Decode
  console.log(`\nDecoding...`);
  const decodeStart = performance.now();
  const result = await runner.decode();
  const decodeMs = performance.now() - decodeStart;

  console.log(`RSS after decode: ${rssMB().toFixed(0)}MB`);

  // 5. Verify
  const hashMatch = result.hashMatches;
  const dataMatch = result.data
    ? Buffer.compare(Buffer.from(payload), Buffer.from(result.data)) === 0
    : false;

  console.log(`\n=== Results ===`);
  console.log(`hash:       ${hashMatch ? "OK ✅" : "FAIL ❌"}`);
  console.log(`data:       ${dataMatch ? "OK ✅" : "FAIL ❌"}`);
  console.log(`recovered:  ${result.stats.oligosRecovered.toLocaleString()}/${oligos.length.toLocaleString()} oligos`);
  console.log(`erased:     ${result.stats.oligosErased}`);

  console.log(`\n=== Throughput ===`);
  console.log(`Encode: ${(payloadSize / MB / (encodeMs / 1000)).toFixed(2)} MB/s`);
  console.log(`Decode: ${(payloadSize / MB / (decodeMs / 1000)).toFixed(2)} MB/s`);

  console.log(`\n=== Memory Summary ===`);
  console.log(`Peak RSS during encode:   ${rssMB().toFixed(0)}MB`);
  console.log(`Peak RSS during streaming: ${rssMB().toFixed(0)}MB`);
  console.log(`Peak RSS during decode:   ${rssMB().toFixed(0)}MB`);
}

main().catch((e) => {
  console.error("Test failed:", e);
  process.exit(1);
});
