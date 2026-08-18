/**
 * Medium tier test: E. coli + Yeast with v51-default.
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

async function testOne(data: Uint8Array, name: string) {
  const cfg = V51_DEFAULT_CONFIG;
  console.log(`  Encoding ${name} (${(data.length/1024/1024).toFixed(2)} MB, ${data.length} bytes)...`);
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

  console.log(`  Encoded: ${stats.oligoCount} oligos, density=${stats.netDensityBitsPerNt.toFixed(3)} b/nt, ${encMs}ms`);
  console.log(`  Constraints: GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcViol=${gcV}, hpViol=${hpV}`);

  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });
  console.log(`  Simulated ${sim.totalReads} reads, decoding...`);
  const d0 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - d0;
  const ok = dec.data && dec.data.length === data.length && dec.data.every((b, i) => b === data[i]);

  console.log(`  ${ok && dec.hashMatches ? "PASS" : "FAIL"}: dec=${decMs}ms, roundtrip=${ok}, hash=${dec.hashMatches}`);
  return { pass: ok && dec.hashMatches, oligoCount: stats.oligoCount, density: stats.netDensityBitsPerNt, encMs, decMs, gcV, hpV, gcMin, gcMax, maxHp };
}

async function main() {
  console.log("helix-codec Real-Dataset Validation — Medium Tier\n");
  const results: anyGany[] = [];

  const dir = path.join(__dirname, "..", "datasets", "medium");
  for (const f of ["ecoli-k12.fa", "yeast.fa"]) {
    const p = path.join(dir, f);
    if (fs.existsSync(p)) {
      const seq = parseFasta(fs.readFileSync(p, "utf-8"));
      const result = await testOne(new Uint8Array(Buffer.from(seq, "utf-8")), f);
      results.push({ file: f, ...result });
    }
  }

  // Save results
  const outPath = path.join(__dirname, "..", "datasets", "medium-results.json");
  fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
  console.log(`\nResults saved to ${outPath}`);
}

main().catch(e => { console.error(e); process.exit(1); });
