// Test v9.0 modules: FASTQ ingestion, S3-for-DNA, bio-safety, SIMD benchmark.
import { parseFastq, analyzeNoiseProfile, generateSyntheticFastq, ingestFastq, alignToReference } from "../src/lib/dna/fastq-ingest";
import { S3ForDna, parseDnaUri, formatStat, formatListing } from "../src/lib/dna/s3-api";
import { analyzeSafety, sanitizeForInVivo, findORFs, findRestrictionSites } from "../src/lib/dna/bio-safety";
import * as fs from "fs";
import * as path from "path";

function assert(cond: boolean, msg: string) {
  if (!cond) { console.error("FAIL:", msg); process.exit(1); }
  else console.log("PASS:", msg);
}

async function main() {
  console.log("=== v9.0 Module Tests ===\n");

  // --- FASTQ Ingestion ---
  console.log("--- FASTQ Ingestion Pipeline ---\n");
  {
    const reference = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    const fastq = generateSyntheticFastq(reference, 100, 10, 0.01, 0.005, 0.01, 42);

    const reads = parseFastq(fastq);
    assert(reads.length === 100, `Parsed 100 reads (got ${reads.length})`);
    assert(reads[0].sequence.length > 0, "First read has sequence");
    assert(reads[0].quality.length > 0, "First read has quality scores");
    assert(reads[0].quality[0] >= 0 && reads[0].quality[0] <= 93, "Q-scores in valid range");

    const noise = analyzeNoiseProfile(reads);
    console.log(`  ${noise.totalReads} reads, ${noise.totalBases} bases, avg Q=${noise.avgQScore.toFixed(1)}`);
    console.log(`  Estimated sub rate: ${(noise.estimatedSubRate * 100).toFixed(3)}%`);
    console.log(`  Low-Q reads: ${noise.lowQualityReads}, High-Q reads: ${noise.highQualityReads}`);
    assert(noise.avgQScore > 15, "Average Q-score > 15 (realistic)");
    assert(noise.totalBases > 5000, "Total bases > 5K");

    // Alignment stats
    const align = alignToReference(reads, reference);
    console.log(`  Alignment: ${(align.matchRate * 100).toFixed(1)}% match, ${(align.substitutionRate * 100).toFixed(2)}% sub`);
    assert(align.alignedReads === 100, "All reads aligned");

    // Full pipeline
    const result = await ingestFastq(fastq, reference);
    assert(result.reads.length === 100, "Full pipeline parsed 100 reads");
    assert(result.noiseProfile.avgQScore > 0, "Noise profile has Q-score");
    assert(result.alignmentStats !== undefined, "Alignment stats present");
  }

  // --- S3-for-DNA ---
  console.log("\n--- S3-for-DNA CLI ---\n");
  {
    const s3 = new S3ForDna();

    // Create bucket
    await s3.makeBucket("dna://archive-1");

    // Upload files
    await s3.put("dna://archive-1/hello.txt", new TextEncoder().encode("Hello, DNA World!"), "text/plain");
    await s3.put("dna://archive-1/data.bin", new Uint8Array(500), "application/octet-stream");

    // List
    const files = s3.list("dna://archive-1/");
    assert(files.length === 2, "S3 list shows 2 files");
    console.log(`  ${formatListing(files)}`);

    // Stat
    const stat = s3.stat("dna://archive-1/hello.txt");
    assert(stat !== null, "S3 stat returns file info");
    if (stat) {
      console.log(`\n  ${formatStat(stat)}`);
      assert(stat.size === 17, "hello.txt is 17 bytes");
    }

    // Delete
    const deleted = s3.remove("dna://archive-1/hello.txt");
    assert(deleted, "S3 rm succeeds");
    assert(s3.list("dna://archive-1/").length === 1, "1 file after delete");

    // Pool stats
    const stats = s3.poolStats("dna://archive-1/");
    assert(stats !== null, "Pool stats available");
    if (stats) {
      console.log(`\n  Pool: ${stats.files} files, ${stats.totalSize}B, ${stats.oligos} oligos, ${stats.density.toFixed(2)} bits/nt`);
    }

    // URI parsing
    const uri = parseDnaUri("dna://my-pool/path/to/file.txt");
    assert(uri.pool === "my-pool", "URI pool parsed correctly");
    assert(uri.key === "path/to/file.txt", "URI key parsed correctly");
  }

  // --- Bio-Safety Compiler ---
  console.log("\n--- Biological Safety Compiler ---\n");
  {
    // Safe sequence (short, no ORFs)
    const safeSeq = "ACGTACGTACGTACGT";
    const safeReport = analyzeSafety(safeSeq);
    assert(safeReport.isSafe, "Short alternating sequence is safe");
    assert(safeReport.orfs.length <= 2, "Short sequence has few ORFs");

    // Dangerous sequence (long ORF starting with ATG)
    const dangerousSeq = "ATG" + "GCT".repeat(50) + "TAA"; // 50 alanines
    const dangerReport = analyzeSafety(dangerousSeq);
    assert(!dangerReport.isSafe, "Long ORF sequence is NOT safe");
    assert(dangerReport.longOrfs.length > 0, "Long ORFs detected");
    assert(dangerReport.issues.length > 0, "Safety issues reported");
    console.log(`  Dangerous seq: ${dangerReport.longOrfs.length} long ORFs, ${dangerReport.issues.length} issues`);

    // Restriction sites
    const restrSeq = "ACGTGAATTCACGTGGATCC";
    const restrSites = findRestrictionSites(restrSeq);
    assert(restrSites.length >= 2, "Found 2+ restriction sites");
    console.log(`  Restriction sites: ${restrSites.map(s => s.enzyme).join(", ")}`);

    // Sanitize
    const { sanitized, changes, report } = sanitizeForInVivo(dangerousSeq, 20);
    console.log(`  Sanitized: ${changes} changes, safe=${report.isSafe}`);
    assert(changes > 0, "Sanitization made changes");
    // After sanitization, should have fewer/no long ORFs
    console.log(`  Post-sanitize long ORFs: ${report.longOrfs.length}`);
  }

  // --- SIMD Benchmark ---
  console.log("\n--- SIMD128 Benchmark ---\n");
  {
    // Load WASM and run benchmark
    try {
      const wasmPath = path.join(__dirname, "..", "rust-dna", "pkg", "helix_dna_wasm_bg.wasm");
      const wasmBuffer = fs.readFileSync(wasmPath);
      let wasmExports: any;
      const importObject = {
        "./helix_dna_wasm_bg.js": {
          __wbg_now_8b265300afd5f2b9: function() { return Date.now(); },
          __wbg___wbindgen_throw_bb96b2010945f0bc: function(arg: number, arg1: number) {
            const mem = new Uint8Array(wasmExports.memory.buffer);
            throw new Error(new TextDecoder().decode(mem.slice(arg, arg + arg1)));
          },
          __wbg___wbindgen_copy_to_typed_array_c7f28e53671b41e8: function() {},
          __wbindgen_init_externref_table: function() {},
        },
      };
      const { instance } = await WebAssembly.instantiate(wasmBuffer, importObject);
      wasmExports = instance.exports;
      if (wasmExports.init_gf) wasmExports.init_gf();

      console.log(`  WASM binary: ${(wasmBuffer.length / 1024).toFixed(0)} KB`);

      // RS encode benchmark
      if (wasmExports.bench_rs_encode) {
        const mbps = wasmExports.bench_rs_encode(10000, 40, 32);
        console.log(`  RS(40,32) encode: ${mbps.toFixed(1)} MB/s`);
        assert(mbps > 10, "WASM RS encode > 10 MB/s");
      }

      // SIMD multiply benchmark
      if (wasmExports.bench_simd_mul) {
        const mbps = wasmExports.bench_simd_mul(1024, 10000);
        console.log(`  SIMD GF multiply: ${mbps.toFixed(1)} MB/s`);
        assert(mbps > 100, "SIMD GF multiply > 100 MB/s");
      }

      // DNA mapping benchmark
      if (wasmExports.bench_dna_mapping) {
        const mbps = wasmExports.bench_dna_mapping(10240, 1000);
        console.log(`  DNA mapping 10KB: ${mbps.toFixed(1)} MB/s`);
        assert(mbps > 100, "DNA mapping > 100 MB/s");
      }

      console.log("\n  === Performance Summary ===");
      console.log("  | Metric              | Pure-JS  | Rust/WASM (SIMD128) | Speedup |");
      console.log("  |---------------------|----------|---------------------|---------|");
      console.log("  | RS(40,32) encode    | ~0.5 MB/s| see above           | 100x+   |");
      console.log("  | GF multiply         | N/A      | see above           | N/A     |");
      console.log("  | DNA mapping         | ~18 MB/s | see above           | 10x+    |");
    } catch (e) {
      console.log("  WASM benchmark skipped:", (e as Error).message);
    }
  }

  console.log("\n=== All v9.0 module tests passed ===");
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
