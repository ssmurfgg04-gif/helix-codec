/**
 * Robust medium-tier test: E. coli K-12 + Yeast S288C
 * 
 * Processes one dataset at a time, saves incrementally, 
 * works within long timeouts.
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
  error?: string;
}

async function testOne(data: Uint8Array, name: string): Promise<Result> {
  const cfg = V51_DEFAULT_CONFIG;
  console.log(`\n  Encoding ${name} (${(data.length/1024/1024).toFixed(2)} MB, ${data.length} bytes)...`);
  console.log(`  Start: ${new Date().toISOString()}`);
  
  const t0 = Date.now();
  let enc: any, stats: any;
  try {
    const result = await encodeFile(data, cfg, { fileName: name, contentType: "application/octet-stream" });
    enc = result.encoded;
    stats = result.stats;
  } catch (e: any) {
    console.log(`  ENCODE ERROR: ${e.message?.slice(0, 300)}`);
    return { file: name, preset: "v51-default", pass: false, oligoCount: 0, density: 0, encMs: Date.now() - t0, decMs: 0, gcV: 0, hpV: 0, gcMin: 0, gcMax: 0, maxHp: 0, roundtrip: false, hashOk: false, dataSize: data.length, error: e.message?.slice(0, 300) };
  }
  const encMs = Date.now() - t0;
  console.log(`  Encode done: ${encMs}ms`);

  let gcV = 0, hpV = 0, gcMin = 1, gcMax = 0, maxHp = 0;
  for (const o of enc.oligos) {
    if (o.gc < cfg.constraints.gcMin || o.gc > cfg.constraints.gcMax) gcV++;
    if (o.maxHomopolymer > cfg.constraints.maxHomopolymer) hpV++;
    if (o.gc < gcMin) gcMin = o.gc;
    if (o.gc > gcMax) gcMax = o.gc;
    if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
  }

  console.log(`  Encoded: ${stats.oligoCount} oligos, density=${stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
  console.log(`  Constraints: GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcViol=${gcV}, hpViol=${hpV}`);

  // Use 10× coverage for faster decode
  console.log(`  Simulating reads (10× coverage)...`);
  const sim = simulate(enc.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: "basic" });
  console.log(`  Simulated ${sim.totalReads} reads, decoding...`);
  
  const d0 = Date.now();
  let dec: any;
  try {
    dec = await decodeReads(sim.reads, enc.metadata, cfg, enc.forwardPrimer, enc.reversePrimer);
  } catch (e: any) {
    console.log(`  DECODE ERROR: ${e.message?.slice(0, 300)}`);
    return { file: name, preset: "v51-default", pass: false, oligoCount: stats.oligoCount, density: stats.netDensityBitsPerNt, encMs, decMs: Date.now() - d0, gcV, hpV, gcMin, gcMax, maxHp, roundtrip: false, hashOk: false, dataSize: data.length, error: e.message?.slice(0, 300) };
  }
  const decMs = Date.now() - d0;
  
  let roundtrip = false;
  if (dec.data && dec.data.length === data.length) {
    roundtrip = true;
    for (let i = 0; i < data.length; i++) {
      if (dec.data[i] !== data[i]) { roundtrip = false; break; }
    }
  }
  
  const hashOk = dec.hashMatches;
  const pass = roundtrip && hashOk && gcV === 0 && hpV === 0;
  console.log(`  Decode done: ${decMs}ms, roundtrip=${roundtrip}, hash=${hashOk}`);
  console.log(`  [${pass ? "PASS" : "FAIL"}] ${name}`);
  
  return { file: name, preset: "v51-default", pass, oligoCount: stats.oligoCount, density: stats.netDensityBitsPerNt, encMs, decMs, gcV, hpV, gcMin, gcMax, maxHp, roundtrip, hashOk, dataSize: data.length };
}

async function main() {
  console.log("=".repeat(70));
  console.log("  helix-codec Medium-Tier Real-Dataset Validation");
  console.log("=".repeat(70));

  const results: Result[] = [];
  const dir = path.join(__dirname, "..", "datasets", "medium");

  // Load any existing results for resume
  const outPath = path.join(__dirname, "..", "datasets", "medium-results.json");
  const existingResults: Result[] = (() => {
    try { return JSON.parse(fs.readFileSync(outPath, "utf-8")); } catch { return []; }
  })();
  const completed = new Set(existingResults.map(r => `${r.file}:${r.preset}`));
  if (completed.size > 0) {
    console.log(`  Resuming: ${completed.size} tests already completed`);
    results.push(...existingResults);
  }

  for (const f of ["ecoli-k12.fa", "yeast.fa"]) {
    const key = `${f}:v51-default`;
    if (completed.has(key)) {
      console.log(`  SKIP ${f} (already completed)`);
      continue;
    }
    
    const p = path.join(dir, f);
    if (!fs.existsSync(p)) {
      console.log(`  SKIP: ${f} not found at ${p}`);
      continue;
    }
    
    console.log(`\n  Loading ${f}...`);
    const seq = parseFasta(fs.readFileSync(p, "utf-8"));
    const data = new Uint8Array(Buffer.from(seq, "utf-8"));
    console.log(`  Loaded: ${data.length} bytes`);
    
    const result = await testOne(data, f);
    results.push(result);

    // Save incrementally
    fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
    console.log(`  Saved incremental results`);
  }

  // Summary
  console.log("\n" + "=".repeat(70));
  console.log("  MEDIUM-TIER SUMMARY");
  console.log("=".repeat(70));
  for (const r of results) {
    const status = r.pass ? "PASS" : "FAIL";
    console.log(`  [${status}] ${r.file}: ${r.oligoCount} oligos, density=${r.density.toFixed(3)} b/nt, enc=${r.encMs}ms, dec=${r.decMs}ms, roundtrip=${r.roundtrip}, hash=${r.hashOk}, gcViol=${r.gcV}, hpViol=${r.hpV}`);
  }
  const passCount = results.filter(r => r.pass).length;
  console.log(`\n  Total: ${passCount}/${results.length} passed`);
}

main().catch(e => { console.error(e); process.exit(1); });
