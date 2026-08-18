/**
 * Combined Medium + Large Tier Real-Dataset Validation
 * 
 * Tests E. coli K-12, Yeast, and Human chr21 against v51-default.
 * Each dataset is processed one at a time with progress tracking.
 * Saves results incrementally for resume capability.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(content: string): string {
  return content.split("\n").filter(l => !l.startsWith(">") && l.trim()).join("");
}

interface Result {
  file: string;
  preset: string;
  pass: boolean;
  oligoCount: number;
  density: number;
  encMs: number;
  decMs: number;
  gcV: number;
  hpV: number;
  gcMin: number;
  gcMax: number;
  maxHp: number;
  roundtrip: boolean;
  hashOk: boolean;
  dataSize: number;
  screeningRetries: number;
  error?: string;
}

async function testOne(data: Uint8Array, name: string): Promise<Result> {
  const cfg = V51_DEFAULT_CONFIG;
  console.log(`\n  === ${name} ===`);
  console.log(`  Size: ${(data.length/1024/1024).toFixed(2)} MB (${data.length} bytes)`);
  console.log(`  Start: ${new Date().toISOString()}`);

  const t0 = Date.now();
  let enc: any, stats: any;
  try {
    const result = await encodeFile(data, cfg, { fileName: name, contentType: "application/octet-stream" });
    enc = result.encoded;
    stats = result.stats;
  } catch (e: any) {
    return { file: name, preset: "v51-default", pass: false, oligoCount: 0, density: 0, encMs: Date.now() - t0, decMs: 0, gcV: 0, hpV: 0, gcMin: 0, gcMax: 0, maxHp: 0, roundtrip: false, hashOk: false, dataSize: data.length, screeningRetries: 0, error: e.message?.slice(0, 300) };
  }
  const encMs = Date.now() - t0;

  let gcV = 0, hpV = 0, gcMin = 1, gcMax = 0, maxHp = 0;
  for (const o of enc.oligos) {
    if (o.gc < cfg.constraints.gcMin || o.gc > cfg.constraints.gcMax) gcV++;
    if (o.maxHomopolymer > cfg.constraints.maxHomopolymer) hpV++;
    if (o.gc < gcMin) gcMin = o.gc;
    if (o.gc > gcMax) gcMax = o.gc;
    if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
  }

  console.log(`  Encode: ${encMs}ms, ${stats.oligoCount} oligos, density=${stats.netDensityBitsPerNt.toFixed(3)} b/nt, retries=${stats.screeningRetries}`);
  console.log(`  Constraints: GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcViol=${gcV}, hpViol=${hpV}`);

  // 10× coverage for faster decode
  console.log(`  Simulating reads (10×)...`);
  const sim = simulate(enc.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: "basic" });
  console.log(`  ${sim.totalReads} reads, decoding...`);

  const d0 = Date.now();
  let dec: any;
  try {
    dec = await decodeReads(sim.reads, enc.metadata, cfg, enc.forwardPrimer, enc.reversePrimer);
  } catch (e: any) {
    return { file: name, preset: "v51-default", pass: false, oligoCount: stats.oligoCount, density: stats.netDensityBitsPerNt, encMs, decMs: Date.now() - d0, gcV, hpV, gcMin, gcMax, maxHp, roundtrip: false, hashOk: false, dataSize: data.length, screeningRetries: stats.screeningRetries, error: e.message?.slice(0, 300) };
  }
  const decMs = Date.now() - d0;

  let roundtrip = false;
  if (dec.data && dec.data.length === data.length) {
    roundtrip = true;
    for (let i = 0; i < data.length; i++) { if (dec.data[i] !== data[i]) { roundtrip = false; break; } }
  }
  const hashOk = dec.hashMatches;
  const pass = roundtrip && hashOk;
  console.log(`  Decode: ${decMs}ms, roundtrip=${roundtrip}, hash=${hashOk}`);
  console.log(`  [${pass ? "PASS" : "FAIL"}] ${name}`);

  return { file: name, preset: "v51-default", pass, oligoCount: stats.oligoCount, density: stats.netDensityBitsPerNt, encMs, decMs, gcV, hpV, gcMin, gcMax, maxHp, roundtrip, hashOk, dataSize: data.length, screeningRetries: stats.screeningRetries };
}

async function main() {
  console.log("=".repeat(70));
  console.log("  helix-codec Medium + Large Tier Validation");
  console.log("=".repeat(70));

  const results: Result[] = [];
  const outPath = path.join(__dirname, "..", "datasets", "all-tier-results.json");

  // Resume from existing results
  try {
    const existing = JSON.parse(fs.readFileSync(outPath, "utf-8"));
    if (Array.isArray(existing) && existing.length > 0) {
      results.push(...existing);
      console.log(`  Resuming: ${results.length} tests already completed`);
    }
  } catch {}

  const completed = new Set(results.map(r => r.file));
  const datasets = [
    { file: "ecoli-k12.fa", dir: "medium" },
    { file: "yeast.fa", dir: "medium" },
    { file: "chr21.fa", dir: "large" },
  ];

  for (const ds of datasets) {
    if (completed.has(ds.file)) {
      console.log(`  SKIP ${ds.file} (already completed)`);
      continue;
    }

    const p = path.join(__dirname, "..", "datasets", ds.dir, ds.file);
    if (!fs.existsSync(p)) {
      console.log(`  SKIP: ${ds.file} not found`);
      continue;
    }

    console.log(`\n  Loading ${ds.file}...`);
    const seq = parseFasta(fs.readFileSync(p, "utf-8"));
    const data = new Uint8Array(Buffer.from(seq, "utf-8"));

    const result = await testOne(data, ds.file);
    results.push(result);

    // Save incrementally
    fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
    console.log(`  Saved results`);
  }

  // Summary
  console.log("\n" + "=".repeat(70));
  console.log("  MEDIUM + LARGE TIER SUMMARY");
  console.log("=".repeat(70));
  for (const r of results) {
    const status = r.pass ? "PASS" : "FAIL";
    console.log(`  [${status}] ${r.file}: ${r.oligoCount} oligos, density=${r.density.toFixed(3)} b/nt, enc=${r.encMs}ms, dec=${r.decMs}ms, roundtrip=${r.roundtrip}, hash=${r.hashOk}`);
    if (r.gcV > 0 || r.hpV > 0) console.log(`    Constraints: gcViol=${r.gcV}, hpViol=${r.hpV}, maxHp=${r.maxHp}`);
  }
  const passCount = results.filter(r => r.pass).length;
  console.log(`\n  Total: ${passCount}/${results.length} passed`);
}

main().catch(e => { console.error(e); process.exit(1); });
