/**
 * Small-tier only test for v67 constraint fixes.
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

async function main() {
  console.log("helix-codec v67 Small-Tier Test\n");
  const dsDir = path.join(__dirname, "..", "datasets");
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  let allOk = true;

  // SARS-CoV-2
  {
    const content = fs.readFileSync(path.join(dsDir, "small/sars-cov-2.fa"), "utf-8");
    const seq = parseFasta(content);
    const data = new Uint8Array(Buffer.from(seq, "utf-8"));
    console.log(`SARS-CoV-2: ${data.length} bytes`);
    const enc = await encodeFile(data, cfg, { fileName: "sars-cov-2", contentType: "application/octet-stream" });
    let gcV=0, hpV=0, gcMin=1, gcMax=0, maxHp=0;
    for (const o of enc.encoded.oligos) {
      if (o.gc < 0.4 || o.gc > 0.6) gcV++;
      if (o.maxHomopolymer > 3) hpV++;
      if (o.gc < gcMin) gcMin = o.gc;
      if (o.gc > gcMax) gcMax = o.gc;
      if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
    }
    const ok = gcV === 0 && hpV === 0;
    console.log(`  [${ok?'PASS':'FAIL'}] ${enc.stats.oligoCount} oligos, density=${enc.stats.netDensityBitsPerNt.toFixed(3)}, enc=${enc.stats.encodeTimeMs}ms`);
    console.log(`    GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcV=${gcV}, hpV=${hpV}`);
    if (!ok) allOk = false;
  }

  // UniProt
  {
    const content = fs.readFileSync(path.join(dsDir, "small/uniprot-p00533.fa"), "utf-8");
    const seq = parseFasta(content);
    const data = new Uint8Array(Buffer.from(seq, "utf-8"));
    console.log(`UniProt EGFR: ${data.length} bytes`);
    const enc = await encodeFile(data, cfg, { fileName: "uniprot-egfr", contentType: "application/octet-stream" });
    let gcV=0, hpV=0, gcMin=1, gcMax=0, maxHp=0;
    for (const o of enc.encoded.oligos) {
      if (o.gc < 0.4 || o.gc > 0.6) gcV++;
      if (o.maxHomopolymer > 3) hpV++;
      if (o.gc < gcMin) gcMin = o.gc;
      if (o.gc > gcMax) gcMax = o.gc;
      if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
    }
    const ok = gcV === 0 && hpV === 0;
    console.log(`  [${ok?'PASS':'FAIL'}] ${enc.stats.oligoCount} oligos, density=${enc.stats.netDensityBitsPerNt.toFixed(3)}, enc=${enc.stats.encodeTimeMs}ms`);
    console.log(`    GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcV=${gcV}, hpV=${hpV}`);
    if (!ok) allOk = false;
  }

  console.log(`\n${allOk ? 'ALL PASSED' : 'SOME FAILED'}`);
  if (!allOk) process.exit(1);
}

main().catch(e => { console.error(e); process.exit(1); });
