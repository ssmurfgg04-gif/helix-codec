// Test arithmetic+CRC with various oligo lengths + benchmark throughput
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG, computeLayout } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 262144);

  console.log("=== Arithmetic+CRC at various oligo lengths ===\n");

  for (const oligoLength of [300, 400, 500, 600]) {
    console.log(`--- Oligo length: ${oligoLength}nt ---`);
    const config = { ...DEFAULT_CONFIG, mappingMode: "arithmetic" as const, oligoLength };
    try {
      const layout = computeLayout(config);
      console.log(`  Layout: innerK=${layout.addressBytes + layout.payloadBytes}, innerN=${layout.addressBytes + layout.payloadBytes + layout.innerParityBytes}, totalInnerBytes=${layout.totalInnerBytes}, payloadBytes=${layout.payloadBytes}`);

      const enc = await encodeFile(payload, config, { fileName: "test.bin", contentType: "application/octet-stream" });
      console.log(`  Encoded: ${enc.encoded.oligos.length} oligos, density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);

      // Perfect reads
      const perfectReads = enc.encoded.oligos.map(o => ({ sequence: o.sequence, quality: [] }));
      const t0 = Date.now();
      const dec1 = await decodeReadsUltra(perfectReads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
      console.log(`  Perfect reads: ${Date.now() - t0}ms, hash: ${dec1.hashMatches ? "PASS" : "FAIL"}`);

      if (dec1.hashMatches) {
        // Noisy reads
        const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
        const t1 = Date.now();
        const dec2 = await decodeReadsUltra(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
        const dec2Ms = Date.now() - t1;
        const throughput = (payload.length / 1024 / 1024) / (dec2Ms / 1000);
        console.log(`  Noisy (0.1% sub): ${dec2Ms}ms, hash: ${dec2.hashMatches ? "PASS" : "FAIL"}, throughput: ${throughput.toFixed(1)} MB/s`);

        // Real Erlich noise (0.0176% sub)
        const sim2 = simulate(enc.encoded.oligos, { substitutionRate: 0.000176, insertionRate: 0.0001, deletionRate: 0.0001, coverage: 10, seed: 42 });
        const t2 = Date.now();
        const dec3 = await decodeReadsUltra(sim2.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
        console.log(`  Erlich noise (0.0176% sub): ${Date.now() - t2}ms, hash: ${dec3.hashMatches ? "PASS" : "FAIL"}`);
      }
    } catch (e: any) {
      console.log(`  Error: ${e.message}`);
    }
    console.log();
  }

  // Throughput benchmark for direct mode at different payload sizes
  console.log("=== Direct mode throughput benchmark ===\n");
  const directConfig = { ...DEFAULT_CONFIG };
  for (const size of [65536, 262144, 1048576, 2116608]) {
    const payloadChunk = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, size);
    console.log(`--- Payload: ${(size / 1024 / 1024).toFixed(2)} MB ---`);

    const t0 = Date.now();
    const enc = await encodeFile(payloadChunk, directConfig, { fileName: "test.bin", contentType: "application/octet-stream" });
    const encMs = Date.now() - t0;

    const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
    const t1 = Date.now();
    const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, directConfig, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const decMs = Date.now() - t1;

    const encThroughput = (size / 1024 / 1024) / (encMs / 1000);
    const decThroughput = (size / 1024 / 1024) / (decMs / 1000);
    console.log(`  Encode: ${encMs}ms (${encThroughput.toFixed(1)} MB/s), ${enc.encoded.oligos.length} oligos`);
    console.log(`  Decode: ${decMs}ms (${decThroughput.toFixed(1)} MB/s), ${dec.hashMatches ? "PASS" : "FAIL"}`);
    console.log(`  Density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
    console.log();
  }
}
main().catch(e => { console.error(e); process.exit(1); });
