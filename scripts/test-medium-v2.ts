/**
 * Medium-tier real-dataset validation — E. coli K-12 + Yeast S288C
 *
 * v2: Extended timeout, incremental save, both v51 and v59 presets.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import {
  V51_DEFAULT_CONFIG,
  ULTIMATE_V59_FAST_ENCODE_CONFIG,
} from "../src/lib/dna/presets";
import type { CodecConfig } from "../src/lib/dna/types";
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
  error?: string;
}

async function testOne(data: Uint8Array, name: string, presetName: string, cfg: CodecConfig): Promise<Result> {
  console.log(`  Encoding ${name} with ${presetName} (${(data.length/1024/1024).toFixed(2)} MB, ${data.length} bytes)...`);
  const t0 = Date.now();
  let enc: any, stats: any;
  try {
    const result = await encodeFile(data, cfg, { fileName: name, contentType: "application/octet-stream" });
    enc = result.encoded;
    stats = result.stats;
  } catch (e: any) {
    console.log(`  ENCODE ERROR: ${e.message?.slice(0, 200)}`);
    return { file: name, preset: presetName, pass: false, oligoCount: 0, density: 0, encMs: Date.now() - t0, decMs: 0, gcV: 0, hpV: 0, gcMin: 0, gcMax: 0, maxHp: 0, roundtrip: false, hashOk: false, error: e.message?.slice(0, 200) };
  }
  const encMs = Date.now() - t0;

  let gcV = 0, hpV = 0, gcMin = 1, gcMax = 0, maxHp = 0;
  for (const o of enc.oligos) {
    if (o.gc < 0.4 || o.gc > 0.6) gcV++;
    if (o.maxHomopolymer > 3) hpV++;
    if (o.gc < gcMin) gcMin = o.gc;
    if (o.gc > gcMax) gcMax = o.gc;
    if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
  }

  console.log(`  Encoded: ${stats.oligoCount} oligos, density=${stats.netDensityBitsPerNt.toFixed(3)} b/nt, ${encMs}ms`);
  console.log(`  Constraints: GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcViol=${gcV}, hpViol=${hpV}`);

  // Simulate clean reads at lower coverage (10×) for faster decode
  const sim = simulate(enc.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: "basic" });
  console.log(`  Simulated ${sim.totalReads} reads (10× coverage), decoding...`);
  const d0 = Date.now();
  let dec: any;
  try {
    dec = await decodeReads(sim.reads, enc.metadata, cfg, enc.forwardPrimer, enc.reversePrimer);
  } catch (e: any) {
    console.log(`  DECODE ERROR: ${e.message?.slice(0, 200)}`);
    return { file: name, preset: presetName, pass: false, oligoCount: stats.oligoCount, density: stats.netDensityBitsPerNt, encMs, decMs: Date.now() - d0, gcV, hpV, gcMin, gcMax, maxHp, roundtrip: false, hashOk: false, error: e.message?.slice(0, 200) };
  }
  const decMs = Date.now() - d0;
  const ok = dec.data && dec.data.length === data.length && dec.data.every((b: number, i: number) => b === data[i]);

  console.log(`  ${ok && dec.hashMatches ? "PASS" : "FAIL"}: dec=${decMs}ms, roundtrip=${ok}, hash=${dec.hashMatches}`);
  return { file: name, preset: presetName, pass: ok && dec.hashMatches, oligoCount: stats.oligoCount, density: stats.netDensityBitsPerNt, encMs, decMs, gcV, hpV, gcMin, gcMax, maxHp, roundtrip: ok, hashOk: dec.hashMatches };
}

async function main() {
  console.log("=".repeat(70));
  console.log("  helix-codec Real-Dataset Validation — Medium Tier v2");
  console.log("=".repeat(70));

  const results: Result[] = [];
  const dir = path.join(__dirname, "..", "datasets", "medium");

  const presets: { name: string; config: CodecConfig }[] = [
    { name: "v51-default", config: V51_DEFAULT_CONFIG },
    { name: "v59-fast", config: ULTIMATE_V59_FAST_ENCODE_CONFIG },
  ];

  for (const f of ["ecoli-k12.fa", "yeast.fa"]) {
    const p = path.join(dir, f);
    if (!fs.existsSync(p)) {
      console.log(`  SKIP: ${f} not found`);
      continue;
    }
    const seq = parseFasta(fs.readFileSync(p, "utf-8"));
    const data = new Uint8Array(Buffer.from(seq, "utf-8"));

    for (const preset of presets) {
      const result = await testOne(data, f, preset.name, preset.config);
      results.push(result);

      // Incremental save after each test
      const outPath = path.join(__dirname, "..", "datasets", "medium-results.json");
      fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
    }
  }

  console.log("\n" + "=".repeat(70));
  console.log("  SUMMARY");
  console.log("=".repeat(70));
  for (const r of results) {
    const status = r.pass ? "PASS" : "FAIL";
    console.log(`  [${status}] ${r.file} × ${r.preset}: density=${r.density.toFixed(3)} b/nt, enc=${r.encMs}ms, dec=${r.decMs}ms, roundtrip=${r.roundtrip}, hash=${r.hashOk}${r.error ? `, err=${r.error}` : ""}`);
  }
  const passCount = results.filter(r => r.pass).length;
  console.log(`\n  Total: ${passCount} passed, ${results.length - passCount} failed`);
  console.log(`  Results saved to datasets/medium-results.json`);
}

main().catch(e => { console.error(e); process.exit(1); });
