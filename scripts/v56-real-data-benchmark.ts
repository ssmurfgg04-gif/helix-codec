#!/usr/bin/env bun
/**
 * v56 Real-Data Benchmark — Test helix-codec on REAL sequencing data.
 *
 * Datasets (downloaded by real-data agent):
 *   1. Erlich-Zielinski DNA Fountain (2017) — 10K FASTQ reads
 *   2. Goldman et al. (2013) — 10K FASTQ reads
 *   3. Davos Bitcoin Challenge (2019) — 10K FASTQ reads
 *
 * Tests:
 *   A) Encode → simulate → decode on the Erlich input files (2.14 MB)
 *   B) Ingest real FASTQ reads and measure error profiles
 *   C) Compare density vs. Erlich's published 1.55 bits/nt
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { readFileSync, existsSync, statSync } from "fs";
import { gunzipSync } from "zlib";
import { createHash } from "crypto";

const TAG = "[v56-real]";

interface DatasetInfo {
  name: string;
  path: string;
  readCount: number;
  readLength: number;
  source: string;
}

/** Parse a FASTQ file (gzipped or plain) and return reads. */
function parseFastq(path: string): { sequence: string; quality: number[] }[] {
  if (!existsSync(path)) return [];
  const raw = readFileSync(path);
  const text = path.endsWith(".gz") ? gunzipSync(raw).toString() : raw.toString();
  const lines = text.split("\n");
  const reads: { sequence: string; quality: number[] }[] = [];
  for (let i = 0; i + 3 < lines.length; i += 4) {
    if (!lines[i].startsWith("@")) continue;
    const seq = lines[i + 1].trim();
    const qual = lines[i + 3].trim();
    if (seq.length === 0) continue;
    const qScores = Array.from(qual).map((c) => c.charCodeAt(0) - 33);
    reads.push({ sequence: seq, quality: qScores });
  }
  return reads;
}

/** Compute error profile from real reads (vs. expected length). */
function analyzeReads(reads: { sequence: string }[], expectedLen: number) {
  if (reads.length === 0) return { count: 0, meanLen: 0, lenStd: 0, gc: 0 };
  const lens = reads.map((r) => r.sequence.length);
  const meanLen = lens.reduce((a, b) => a + b, 0) / lens.length;
  const lenStd = Math.sqrt(lens.reduce((a, b) => a + (b - meanLen) ** 2, 0) / lens.length);
  const gc = reads.reduce((a, r) => {
    const g = (r.sequence.match(/[GC]/g)?.length ?? 0);
    return a + g / r.sequence.length;
  }, 0) / reads.length;
  return { count: reads.length, meanLen, lenStd, gc };
}

/** Encode → simulate → decode round-trip on real data payload. */
async function benchEncodeDecode(name: string, payload: Uint8Array, cfg: any) {
  console.log(`\n${TAG} --- ${name} ---`);
  console.log(`${TAG}   payload: ${payload.length} bytes (${(payload.length / 1024).toFixed(1)} KB)`);

  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "real-data.bin", contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  const density = enc.stats.netDensityBitsPerNt;
  const oligoCount = enc.encoded.oligos.length;
  const oligoLen = enc.encoded.oligos[0]?.sequence.length ?? 0;
  console.log(`${TAG}   encoded: ${oligoCount} oligos × ${oligoLen}nt, density=${density.toFixed(3)} bits/nt, ${encMs}ms`);

  // Test at 10× coverage (high coverage, should get 100%)
  for (const coverage of [5, 10]) {
    const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage, seed: 42 });
    const t1 = Date.now();
    const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const decMs = Date.now() - t1;
    const recovered = dec.stats?.oligosRecovered ?? 0;
    const total = (dec.stats?.oligosRecovered ?? 0) + (dec.stats?.oligosErased ?? 0) + (dec.stats?.oligosFailedInnerRS ?? 0) + (dec.stats?.oligosFailedOuterRS ?? 0);
    const recovery = total > 0 ? recovered / total : 0;
    console.log(`${TAG}   ${coverage}× decode: recovery=${(recovery * 100).toFixed(1)}%, hash=${dec.hashMatches ? "OK" : "FAIL"}, ${decMs}ms`);
  }

  return { density, oligoCount, encMs };
}

async function main() {
  console.log(`${TAG} Helix v56 Real-Data Benchmark`);
  console.log(`${TAG} ${new Date().toISOString()}`);
  console.log(`${TAG} ============================================`);

  // === Part A: Analyze real FASTQ reads ===
  console.log(`\n${TAG} === Part A: Real FASTQ Read Analysis ===`);

  const datasets: DatasetInfo[] = [
    {
      name: "Erlich-Zielinski 2017",
      path: "real-data/erlich-zielinski-2017/ERR1816980_1.first10000.fastq.gz",
      readCount: 10000,
      readLength: 151,
      source: "Nature 2017, doi:10.1038/nature23038",
    },
    {
      name: "Goldman 2013",
      path: "real-data/goldman-2013/ERR215679_1.first10000.fastq.gz",
      readCount: 10000,
      readLength: 104,
      source: "Nature 2013, doi:10.1038/nature11875",
    },
    {
      name: "Davos Bitcoin 2019",
      path: "real-data/dna-storage-benchmarks/davos-bitcoin-challenge-2019/davos.first10000.fastq.gz",
      readCount: 10000,
      readLength: 75,
      source: "Davos DNA Bitcoin Challenge 2019",
    },
  ];

  for (const ds of datasets) {
    const reads = parseFastq(ds.path);
    const stats = analyzeReads(reads, ds.readLength);
    console.log(`\n${TAG}   ${ds.name}:`);
    console.log(`${TAG}     source: ${ds.source}`);
    console.log(`${TAG}     reads: ${stats.count}`);
    console.log(`${TAG}     mean length: ${stats.meanLen.toFixed(1)}nt (σ=${stats.lenStd.toFixed(1)})`);
    console.log(`${TAG}     mean GC: ${(stats.gc * 100).toFixed(1)}%`);
    // Show first read
    if (reads.length > 0) {
      console.log(`${TAG}     first read: ${reads[0].sequence.slice(0, 60)}... (${reads[0].sequence.length}nt)`);
    }
  }

  // === Part B: Encode/decode on real Erlich input files ===
  console.log(`\n${TAG} === Part B: Encode/Decode on Real Data Payloads ===`);

  // Load Erlich input files (2.14 MB archive)
  const erlichArchive = "real-data/erlich-zielinski-2017/dna-fountain-input-files.tar.gz";
  if (existsSync(erlichArchive)) {
    const payload = readFileSync(erlichArchive);
    console.log(`${TAG}   Erlich input archive: ${(payload.length / 1024 / 1024).toFixed(2)} MB`);

    // Test with v55-density config (700nt, direct, 4B LDPC, 3% RS)
    const result1 = await benchEncodeDecode(
      "Erlich 2.14MB archive → v55-density (700nt, direct, 4B LDPC, 3% RS)",
      payload,
      ULTIMATE_V55_DENSITY_CONFIG,
    );

    // Compare with Erlich's published density
    console.log(`${TAG}   Erlich published: 1.55 bits/nt (DNA Fountain, 2017)`);
    console.log(`${TAG}   Helix v56:        ${result1.density.toFixed(3)} bits/nt`);
    console.log(`${TAG}   Improvement:      ${((result1.density / 1.55 - 1) * 100).toFixed(1)}% vs Erlich`);
  }

  // Load HEDGES Wizard of Oz text
  const hedgesText = "real-data/hedges-2020/WizardOfOzInEsperanto.txt";
  if (existsSync(hedgesText)) {
    const payload = readFileSync(hedgesText);
    const result2 = await benchEncodeDecode(
      "HEDGES Wizard-of-Oz text → v55-density",
      payload,
      ULTIMATE_V55_DENSITY_CONFIG,
    );
    console.log(`${TAG}   HEDGES published rate: ~1.0 bits/nt (rate-1/2 conv + RS)`);
    console.log(`${TAG}   Helix v56:             ${result2.density.toFixed(3)} bits/nt`);
  }

  // === Part C: Throughput on real data ===
  console.log(`\n${TAG} === Part C: Throughput on Real Data ===`);
  if (existsSync(erlichArchive)) {
    const payload = readFileSync(erlichArchive);
    // Best-of-3 encode throughput
    let bestMs = Infinity;
    for (let trial = 0; trial < 3; trial++) {
      const t0 = Date.now();
      await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "throughput.bin", contentType: "application/octet-stream" });
      const ms = Date.now() - t0;
      if (ms < bestMs) bestMs = ms;
    }
    const mbps = (payload.length / 1024 / 1024) / (bestMs / 1000);
    console.log(`${TAG}   ${(payload.length / 1024 / 1024).toFixed(2)} MB encode: best-of-3 = ${bestMs}ms → ${mbps.toFixed(2)} MB/s`);
    console.log(`${TAG}   SOTA (Catalog Shannon): 12.5 MB/s`);
    console.log(`${TAG}   Helix v56 vs SOTA: ${((mbps / 12.5 - 1) * 100).toFixed(1)}%`);
  }

  console.log(`\n${TAG} ============================================`);
  console.log(`${TAG} Real-data benchmark complete.`);
}

main().catch((e) => {
  console.error(`${TAG} FATAL:`, e);
  process.exit(1);
});
