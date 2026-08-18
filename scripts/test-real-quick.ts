/**
 * Quick test: SARS-CoV-2 + UniProt with v51-default.
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
  console.log(`  Encoding ${name} (${data.length} bytes)...`);
  const t0 = Date.now();
  const enc = await encodeFile(data, cfg, { fileName: name, contentType: "application/octet-stream" });
  const encMs = Date.now() - t0;
  const stats = enc.stats;

  let gcV = 0, hpV = 0;
  for (const o of enc.encoded.oligos) {
    if (o.gc < 0.4 || o.gc > 0.6) gcV++;
    if (o.maxHomopolymer > 3) hpV++;
  }

  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });
  const d0 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decMs = Date.now() - d0;
  const ok = dec.data && dec.data.length === data.length && dec.data.every((b, i) => b === data[i]);

  console.log(`  ${ok && dec.hashMatches ? "PASS" : "FAIL"}: ${stats.oligoCount} oligos, ` +
    `density=${stats.netDensityBitsPerNt.toFixed(3)} b/nt, enc=${encMs}ms, dec=${decMs}ms, ` +
    `gcViol=${gcV}, hpViol=${hpV}, roundtrip=${ok}, hash=${dec.hashMatches}`);
  return ok && dec.hashMatches;
}

async function main() {
  console.log("helix-codec Real-Dataset Validation — Small Tier\n");
  let allPass = true;

  const dir = path.join(__dirname, "..", "datasets", "small");
  for (const f of ["sars-cov-2.fa", "uniprot-p00533.fa"]) {
    const p = path.join(dir, f);
    if (fs.existsSync(p)) {
      const seq = parseFasta(fs.readFileSync(p, "utf-8"));
      const ok = await testOne(new Uint8Array(Buffer.from(seq, "utf-8")), f);
      if (!ok) allPass = false;
    }
  }
  if (!allPass) process.exit(1);
}

main().catch(e => { console.error(e); process.exit(1); });
