/**
 * Encode-only test for medium/large datasets.
 * Tests encoding + constraint validation (no decode to save time).
 */
import { encodeFile } from "../src/lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(content: string): string {
  return content.split("\n").filter(l => !l.startsWith(">") && l.trim()).join("");
}

async function testEncode(data: Uint8Array, name: string) {
  const cfg = V51_DEFAULT_CONFIG;
  const sizeMB = (data.length / 1024 / 1024).toFixed(2);
  console.log(`  Encoding ${name} (${sizeMB} MB)...`);
  const t0 = Date.now();
  const enc = await encodeFile(data, cfg, { fileName: name, contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  const stats = enc.stats;

  let gcV = 0, hpV = 0, gcMin = 1, gcMax = 0, maxHp = 0;
  for (const o of enc.encoded.oligos) {
    if (o.gc < 0.4 || o.gc > 0.6) gcV++;
    if (o.maxHomopolymer > 3) hpV++;
    if (o.gc < gcMin) gcMin = o.gc;
    if (o.gc > gcMax) gcMax = o.gc;
    if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
  }

  const totalV = gcV + hpV;
  const violPct = (totalV / stats.oligoCount * 100).toFixed(1);
  console.log(`  ${totalV === 0 ? "PASS" : "WARN"}: ${stats.oligoCount} oligos, density=${stats.netDensityBitsPerNt.toFixed(3)} b/nt, enc=${encMs}ms`);
  console.log(`    GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcViol=${gcV}, hpViol=${hpV} (${violPct}% of oligos)`);

  return {
    file: name, sizeBytes: data.length, oligoCount: stats.oligoCount,
    density: stats.netDensityBitsPerNt, encMs,
    gcRange: { min: gcMin, max: gcMax }, maxHp,
    gcViolations: gcV, hpViolations: hpV,
  };
}

async function main() {
  console.log("helix-codec Real-Dataset Encode Validation\n");
  const results: any[] = [];

  // Small
  const smallDir = path.join(__dirname, "..", "datasets", "small");
  console.log("=== SMALL ===");
  for (const f of ["sars-cov-2.fa", "uniprot-p00533.fa"]) {
    const p = path.join(smallDir, f);
    if (fs.existsSync(p)) {
      const seq = parseFasta(fs.readFileSync(p, "utf-8"));
      const r = await testEncode(new Uint8Array(Buffer.from(seq, "utf-8")), f);
      results.push(r);
    }
  }

  // Medium
  const mediumDir = path.join(__dirname, "..", "datasets", "medium");
  console.log("\n=== MEDIUM ===");
  for (const f of ["ecoli-k12.fa", "yeast.fa"]) {
    const p = path.join(mediumDir, f);
    if (fs.existsSync(p)) {
      const seq = parseFasta(fs.readFileSync(p, "utf-8"));
      const r = await testEncode(new Uint8Array(Buffer.from(seq, "utf-8")), f);
      results.push(r);
    }
  }

  // Large
  const largeDir = path.join(__dirname, "..", "datasets", "large");
  console.log("\n=== LARGE ===");
  const chr21 = path.join(largeDir, "chr21.fa");
  if (fs.existsSync(chr21)) {
    const seq = parseFasta(fs.readFileSync(chr21, "utf-8"));
    const r = await testEncode(new Uint8Array(Buffer.from(seq, "utf-8")), "chr21.fa");
    results.push(r);
  }

  // Save
  const outPath = path.join(__dirname, "..", "datasets", "encode-results.json");
  fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
  console.log(`\nResults saved to ${outPath}`);
}

main().catch(e => { console.error(e); process.exit(1); });
