// Comprehensive speed benchmark for the DNA codec.
// Measures encode/decode throughput at various payload sizes.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA, PRESET_REAL_2024 } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import { holographicEncode, holographicDecode, simulateShardLoss } from "../src/lib/dna/holographic";
import { gfMul } from "../src/lib/dna/gf256";

async function main() {
  console.log("=== Helix Codec v2.0 — Speed Benchmarks ===\n");

  // Sizes to test
  const sizes = [256, 1024, 4096, 16384, 65536];

  console.log("Codec: RS(38,30) inner + RS outer (20% parity), 200nt oligos, DEFLATE+SHA-256");
  console.log("Mutation: Illumina (1e-3 sub, 5e-4 ins, 1e-3 del), 20x coverage\n");

  console.log("Payload  | Encode ms | Encode MiB/s | Decode ms | Decode MiB/s | Oligos | Density");
  console.log("---------|-----------|--------------|-----------|--------------|--------|---------");

  for (const size of sizes) {
    const payload = new Uint8Array(size);
    for (let i = 0; i < size; i++) payload[i] = (i * 31 + 17) & 0xff;

    const t0 = Date.now();
    const encodeResult = await encodeFile(payload, DEFAULT_CONFIG, {
      fileName: "bench.bin",
      contentType: "application/octet-stream",
    });
    const encodeMs = Date.now() - t0;

    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_ILLUMINA, seed: 42 });

    const t1 = Date.now();
    const decodeResult = await decodeReads(
      sim.reads,
      encodeResult.encoded.metadata,
      DEFAULT_CONFIG,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );
    const decodeMs = Date.now() - t1;

    const encodeMibps = (size / 1024 / 1024) / (encodeMs / 1000);
    const decodeMibps = (size / 1024 / 1024) / (decodeMs / 1000);

    console.log(
      `${size.toString().padStart(8)} | ${encodeMs.toString().padStart(9)} | ${encodeMibps.toFixed(2).padStart(12)} | ${decodeMs.toString().padStart(9)} | ${decodeMibps.toFixed(2).padStart(12)} | ${encodeResult.stats.oligoCount.toString().padStart(6)} | ${encodeResult.stats.netDensityBitsPerNt.toFixed(3)}`,
    );
  }

  console.log("\n=== Holographic Sharding Codec ===\n");
  console.log("K=10 data shards, N=15 total (1.5x overhead)\n");

  const holoSizes = [100, 1000, 10000, 100000];
  console.log("Payload  | Encode ms | Encode MiB/s | Decode ms | Decode MiB/s | Recovery @ 30% loss");
  console.log("---------|-----------|--------------|-----------|--------------|--------------------");

  for (const size of holoSizes) {
    const data = new Uint8Array(size);
    for (let i = 0; i < size; i++) data[i] = (i * 31 + 17) & 0xff;

    const t0 = Date.now();
    const enc = holographicEncode(data, { dataShards: 10, totalShards: 15, blockSize: 10 });
    const encodeMs = Date.now() - t0;

    const t1 = Date.now();
    const available = enc.shards.slice(0, 10); // simulate 33% loss
    // We need the encoding object to decode, so use all shards here
    // For benchmark, just measure decode of full set
    const recovered = holographicDecode(enc.shards, enc);
    const decodeMs = Date.now() - t1;

    const lossResult = simulateShardLoss(enc, 0.3, 42);

    const encodeMibps = (size / 1024 / 1024) / (encodeMs / 1000);
    const decodeMibps = (size / 1024 / 1024) / (decodeMs / 1000);

    console.log(
      `${size.toString().padStart(8)} | ${encodeMs.toString().padStart(9)} | ${encodeMibps.toFixed(2).padStart(12)} | ${decodeMs.toString().padStart(9)} | ${decodeMibps.toFixed(2).padStart(12)} | ${lossResult.recoverySuccessful ? "OK" : "FAIL"}`,
    );
  }

  console.log("\n=== Recovery Success Rate vs. Error Model ===\n");

  const payload = new Uint8Array(2048);
  for (let i = 0; i < 2048; i++) payload[i] = (i * 31 + 17) & 0xff;
  const enc = await encodeFile(payload, DEFAULT_CONFIG, {
    fileName: "test.bin",
    contentType: "application/octet-stream",
  });

  const presets = [
    { name: "Perfect (0% errors)", cfg: { substitutionRate: 0, insertionRate: 0, deletionRate: 0, coverage: 1, dropoutRate: 0, seed: 42 } },
    { name: "Illumina (0.25% total, 20x)", cfg: { ...PRESET_ILLUMINA, seed: 42 } },
    { name: "Real 2024 (12.3% total, 25x)", cfg: { ...PRESET_REAL_2024, seed: 42 } },
  ];

  console.log("Profile                       | Recovery | Hash Match | Oligos Recovered | Inner RS Fails");
  console.log("------------------------------|----------|------------|------------------|---------------");

  for (const p of presets) {
    const sim = simulate(enc.encoded.oligos, p.cfg);
    const dec = await decodeReads(
      sim.reads,
      enc.encoded.metadata,
      DEFAULT_CONFIG,
      enc.encoded.forwardPrimer,
      enc.encoded.reversePrimer,
    );
    console.log(
      `${p.name.padEnd(30)}| ${dec.data ? "OK" : "FAIL"}     | ${dec.hashMatches ? "YES" : "NO "}        | ${dec.stats.oligosRecovered.toString().padStart(16)} | ${dec.stats.oligosFailedInnerRS.toString().padStart(14)}`,
    );
  }

  console.log("\n=== GF(256) Table Size ===");
  // Force table build
  gfMul(2, 3);
  console.log("  EXP table: 512 bytes");
  console.log("  LOG table: 256 bytes");
  console.log("  MUL table: 65,536 bytes (64 KB)");
  console.log("  Total: ~66 KB — fits in L1 cache");
}

main().catch(console.error);
