/**
 * Focused real-dataset test: SARS-CoV-2 + UniProt with v51-default only.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { V51_DEFAULT_CONFIG, ULTIMATE_V59_FAST_ENCODE_CONFIG } from "../src/lib/dna/presets";
import type { CodecConfig } from "../src/lib/dna/types";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(content: string): string {
  return content.split("\n").filter(l => !l.startsWith(">") && l.trim()).join("");
}

async function testOne(data: Uint8Array, name: string, cfg: CodecConfig, label: string) {
  console.log(`\n  ${label}: encoding ${data.length} bytes...`);
  const enc = await encodeFile(data, cfg, { fileName: name, contentType: "application/octet-stream" });
  const stats = enc.stats;

  // Check constraints
  let gcV = 0, hpV = 0, gcMin = 1, gcMax = 0, maxHp = 0;
  for (const o of enc.encoded.oligos) {
    if (o.gc < cfg.constraints.gcMin || o.gc > cfg.constraints.gcMax) gcV++;
    if (o.maxHomopolymer > cfg.constraints.maxHomopolymer) hpV++;
    if (o.gc < gcMin) gcMin = o.gc;
    if (o.gc > gcMax) gcMax = o.gc;
    if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
  }

  // Simulate + decode
  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const roundtrip = dec.data && dec.data.length === data.length && dec.data.every((b, i) => b === data[i]);

  const status = roundtrip && dec.hashMatches ? "PASS" : "FAIL";
  console.log(`  [${status}] ${label}: ${stats.oligoCount} oligos, density=${stats.netDensityBitsPerNt.toFixed(3)} b/nt, ` +
    `enc=${stats.encodeTimeMs}ms, GC=[${gcMin.toFixed(2)},${gcMax.toFixed(2)}], maxHp=${maxHp}, ` +
    `gcViol=${gcV}, hpViol=${hpV}, roundtrip=${roundtrip}, hash=${dec.hashMatches}`);
}

async function main() {
  console.log("helix-codec Real-Dataset Validation — Small Tier");

  const smallDir = path.join(__dirname, "..", "datasets", "small");
  const presets = [
    { name: "v51-default", cfg: V51_DEFAULT_CONFIG },
    { name: "v59-fast", cfg: ULTIMATE_V59_FAST_ENCODE_CONFIG },
  ];

  const datasets: { name: string; bytes: Uint8Array }[] = [];
  for (const f of ["sars-cov-2.fa", "uniprot-p00533.fa"]) {
    const p = path.join(smallDir, f);
    if (fs.existsSync(p)) {
      const seq = parseFasta(fs.readFileSync(p, "utf-8"));
      datasets.push({ name: f, bytes: new Uint8Array(Buffer.from(seq, "utf-8")) });
    }
  }

  for (const ds of datasets) {
    console.log(`\n=== ${ds.name} (${(ds.bytes.length / 1024).toFixed(1)} KB) ===`);
    for (const preset of presets) {
      await testOne(ds.bytes, ds.name, preset.cfg, preset.name);
    }
  }

  // Medium tier — just v51-default
  const mediumDir = path.join(__dirname, "..", "datasets", "medium");
  for (const f of ["ecoli-k12.fa", "yeast.fa"]) {
    const p = path.join(mediumDir, f);
    if (fs.existsSync(p)) {
      const seq = parseFasta(fs.readFileSync(p, "utf-8"));
      const bytes = new Uint8Array(Buffer.from(seq, "utf-8"));
      console.log(`\n=== ${f} (${(bytes.length / 1024 / 1024).toFixed(1)} MB) ===`);
      await testOne(bytes, f, V51_DEFAULT_CONFIG, "v51-default");
    }
  }
}

main().catch(e => { console.error(e); process.exit(1); });
