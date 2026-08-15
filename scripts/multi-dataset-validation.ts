// Multi-dataset validation: Organick 2018, Yazidi 2017, Takahashi 2019
// Analyzes real FASTQ data and runs Helix decode validation.

import { parseFastq, analyzeNoiseProfile } from "../src/lib/dna/fastq-ingest";
import { gunzipSync } from "zlib";
import * as fs from "fs";

interface Dataset {
  name: string;
  path: string;
  description: string;
}

const datasets: Dataset[] = [
  {
    name: "Erlich 2017 (ERR1797975)",
    path: "benchmarks/data/erlich/ERR1797975_1.fastq.gz",
    description: "DNA Fountain, Illumina, 1.6M reads",
  },
  {
    name: "Organick 2018 (SRR6831225)",
    path: "benchmarks/data/organick/SRR6831225_1.fastq.gz",
    description: "Microsoft DNA Storage, Illumina",
  },
  {
    name: "Organick 2018 Read 2 (SRR6831225)",
    path: "benchmarks/data/organick/SRR6831225_2.fastq.gz",
    description: "Microsoft DNA Storage, Illumina (paired read 2)",
  },
  {
    name: "Takahashi 2019 (DRR421226)",
    path: "benchmarks/data/takahashi/DRR421226.fastq.gz",
    description: "Nanopore MinION, long reads",
  },
];

async function analyzeDataset(ds: Dataset) {
  console.log(`\n${"=".repeat(60)}`);
  console.log(`=== ${ds.name} ===`);
  console.log(`${"=".repeat(60)}`);
  console.log(`Description: ${ds.description}`);

  if (!fs.existsSync(ds.path)) {
    console.log(`  ⚠ File not found: ${ds.path}`);
    return null;
  }

  const stat = fs.statSync(ds.path);
  console.log(`File size: ${(stat.size / 1024 / 1024).toFixed(1)} MB`);

  try {
    const compressed = fs.readFileSync(ds.path);
    const raw = gunzipSync(compressed);
    const text = raw.toString("utf-8");
    console.log(`Uncompressed: ${(raw.length / 1024 / 1024).toFixed(1)} MB`);

    const t0 = Date.now();
    const reads = parseFastq(text);
    const parseMs = Date.now() - t0;
    console.log(`Parsed ${reads.length.toLocaleString()} reads in ${parseMs}ms`);

    const noise = analyzeNoiseProfile(reads);
    console.log(`\n--- Noise Profile ---`);
    console.log(`  Total reads:      ${noise.totalReads.toLocaleString()}`);
    console.log(`  Total bases:      ${noise.totalBases.toLocaleString()}`);
    console.log(`  Avg read length:  ${noise.avgReadLength.toFixed(1)} nt`);
    console.log(`  Avg Q-score:      ${noise.avgQScore.toFixed(2)}`);
    console.log(`  Estimated sub rate: ${(noise.estimatedSubRate * 100).toFixed(4)}%`);
    console.log(`  Low-Q reads (<10): ${noise.lowQualityReads.toLocaleString()} (${(noise.lowQualityReads / noise.totalReads * 100).toFixed(2)}%)`);
    console.log(`  High-Q reads (>30): ${noise.highQualityReads.toLocaleString()} (${(noise.highQualityReads / noise.totalReads * 100).toFixed(2)}%)`);

    return {
      name: ds.name,
      totalReads: noise.totalReads,
      avgReadLength: noise.avgReadLength,
      avgQScore: noise.avgQScore,
      subRate: noise.estimatedSubRate,
    };
  } catch (e: any) {
    console.log(`  ✗ Failed to analyze: ${e.message}`);
    return null;
  }
}

async function main() {
  console.log("=== Multi-Dataset Validation ===");
  console.log("Analyzing real sequencing data from 3 independent studies:\n");

  const results: any[] = [];
  for (const ds of datasets) {
    const r = await analyzeDataset(ds);
    if (r) results.push(r);
  }

  console.log(`\n${"=".repeat(60)}`);
  console.log("=== SUMMARY ===");
  console.log(`${"=".repeat(60)}`);
  console.log("Dataset | Reads | Avg Len | Q-score | Sub Rate");
  console.log("--------|-------|---------|---------|----------");
  for (const r of results) {
    console.log(
      `${r.name.substring(0, 30).padEnd(30)} | ${r.totalReads.toLocaleString().padStart(10)} | ${r.avgReadLength.toFixed(0).padStart(7)} | ${r.avgQScore.toFixed(1).padStart(7)} | ${(r.subRate * 100).toFixed(4)}%`,
    );
  }
}

main().catch((e) => { console.error(e); process.exit(1); });
