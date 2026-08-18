/**
 * Quick test of v67 fixes: constraint retry + auto-sharding.
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

async function testEncode(filePath: string, label: string) {
  const content = fs.readFileSync(filePath, "utf-8");
  const seq = parseFasta(content);
  const data = new Uint8Array(Buffer.from(seq, "utf-8"));
  console.log(`\n${label}: ${data.length} bytes (${(data.length/1024).toFixed(1)} KB)`);
  
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const t0 = Date.now();
  const enc = await encodeFile(data, cfg, { fileName: label, contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  
  let gcV = 0, hpV = 0, gcMin = 1, gcMax = 0, maxHp = 0;
  for (const o of enc.encoded.oligos) {
    if (o.gc < 0.4 || o.gc > 0.6) gcV++;
    if (o.maxHomopolymer > 3) hpV++;
    if (o.gc < gcMin) gcMin = o.gc;
    if (o.gc > gcMax) gcMax = o.gc;
    if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
  }
  
  const ok = gcV === 0 && hpV === 0;
  console.log(`  [${ok ? 'PASS' : 'FAIL'}] ${enc.stats.oligoCount} oligos, density=${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt, enc=${encMs}ms`);
  console.log(`    GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcViol=${gcV}, hpViol=${hpV} (${((gcV+hpV)/enc.stats.oligoCount*100).toFixed(1)}%)`);
  console.log(`    retries=${enc.stats.screeningRetries}, shards=${enc.encoded.metadata.shardCount ?? 1}`);
  return ok;
}

async function main() {
  console.log("helix-codec v67 Real-Dataset Quick Test\n");
  const dsDir = path.join(__dirname, "..", "datasets");
  let allOk = true;
  
  // Small
  for (const [f, label] of [["small/sars-cov-2.fa", "SARS-CoV-2"], ["small/uniprot-p00533.fa", "UniProt EGFR"]]) {
    const ok = await testEncode(path.join(dsDir, f), label);
    if (!ok) allOk = false;
  }
  
  // Medium
  for (const [f, label] of [["medium/ecoli-k12.fa", "E.coli K-12"], ["medium/yeast.fa", "Yeast S288C"]]) {
    const ok = await testEncode(path.join(dsDir, f), label);
    if (!ok) allOk = false;
  }
  
  console.log(`\n${allOk ? 'ALL PASSED' : 'SOME FAILED'}`);
  if (!allOk) process.exit(1);
}

main().catch(e => { console.error(e); process.exit(1); });
