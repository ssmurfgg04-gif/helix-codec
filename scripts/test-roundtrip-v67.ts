/**
 * Full roundtrip test: encode → simulate → decode → verify
 * Tests against real biological datasets only.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import type { CodecConfig } from "../src/lib/dna/types";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(content: string): string {
  return content.split("\n").filter(l => !l.startsWith(">") && l.trim()).join("");
}

async function testRoundtrip(filePath: string, label: string, cfg: CodecConfig) {
  const content = fs.readFileSync(filePath, "utf-8");
  const seq = parseFasta(content);
  const data = new Uint8Array(Buffer.from(seq, "utf-8"));
  console.log(`\n  ${label}: ${data.length} bytes (${(data.length/1024).toFixed(1)} KB)`);

  // Encode
  const enc = await encodeFile(data, cfg, { fileName: label, contentType: "application/octet-stream" });

  // Constraint check
  let gcV=0, hpV=0, gcMin=1, gcMax=0, maxHp=0;
  for (const o of enc.encoded.oligos) {
    if (o.gc < 0.4 || o.gc > 0.6) gcV++;
    if (o.maxHomopolymer > 3) hpV++;
    if (o.gc < gcMin) gcMin = o.gc;
    if (o.gc > gcMax) gcMax = o.gc;
    if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
  }
  const constraintsOk = gcV === 0 && hpV === 0;
  console.log(`    Encode: ${enc.stats.oligoCount} oligos, density=${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt, ${enc.stats.encodeTimeMs}ms`);
  console.log(`    Constraints: ${constraintsOk ? 'OK' : 'FAIL'} GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}] maxHp=${maxHp}`);

  // Simulate clean reads (30x coverage)
  const simResult = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });

  // Decode
  const decodeStart = Date.now();
  const dec = await decodeReads(simResult.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decodeMs = Date.now() - decodeStart;

  // Verify roundtrip
  let roundtripOk = false;
  if (dec.data && dec.data.length === data.length) {
    roundtripOk = true;
    for (let i = 0; i < data.length; i++) {
      if (dec.data[i] !== data[i]) { roundtripOk = false; break; }
    }
  }
  const hashOk = dec.hashMatches;
  
  console.log(`    Decode: ${decodeMs}ms, roundtrip=${roundtripOk ? 'OK' : 'FAIL'}, hash=${hashOk ? 'OK' : 'FAIL'}`);

  return { constraintsOk, roundtripOk, hashOk };
}

async function main() {
  console.log("=".repeat(60));
  console.log("  helix-codec v67 Real-Dataset Roundtrip Test");
  console.log("=".repeat(60));

  const dsDir = path.join(__dirname, "..", "datasets");
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  let allOk = true;

  console.log("\n=== SMALL TIER ===");
  for (const [f, label] of [["small/sars-cov-2.fa", "SARS-CoV-2"], ["small/uniprot-p00533.fa", "UniProt EGFR"]]) {
    const r = await testRoundtrip(path.join(dsDir, f), label, cfg);
    if (!r.constraintsOk || !r.roundtripOk || !r.hashOk) allOk = false;
  }

  console.log(`\n${"=".repeat(60)}`);
  console.log(`  ${allOk ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED'}`);
  if (!allOk) process.exit(1);
}

main().catch(e => { console.error(e); process.exit(1); });
