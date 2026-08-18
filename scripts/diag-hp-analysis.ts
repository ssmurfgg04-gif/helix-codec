/**
 * Analyze WHERE homopolymer violations occur in encoded oligos
 */
import { encodeFile } from "../src/lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(c: string) { return c.split("\n").filter(l => !l.startsWith(">") && l.trim()).join(""); }

function findHomopolymers(dna: string, maxAllowed: number) {
  const violations: { pos: number; base: string; run: number }[] = [];
  let run = 1;
  for (let i = 1; i < dna.length; i++) {
    if (dna[i] === dna[i-1]) {
      run++;
      if (run === maxAllowed + 1) {
        violations.push({ pos: i - maxAllowed, base: dna[i], run });
      }
    } else {
      run = 1;
    }
  }
  return violations;
}

async function main() {
  const dsDir = path.join(__dirname, "..", "datasets");
  const content = fs.readFileSync(path.join(dsDir, "small/sars-cov-2.fa"), "utf-8");
  const seq = parseFasta(content);
  const data = new Uint8Array(Buffer.from(seq, "utf-8"));
  
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const enc = await encodeFile(data, cfg, { fileName: "sars-cov-2", contentType: "application/octet-stream" });
  
  const primerLen = cfg.primerLength;
  const addressNt = 16; // 4 bytes * 4 nt/byte
  
  let hpInAddr = 0, hpInPayload = 0, hpAtBoundary = 0, hpTotal = 0;
  
  for (let i = 0; i < enc.encoded.oligos.length; i++) {
    const oligo = enc.encoded.oligos[i];
    const innerDna = oligo.sequence.slice(primerLen, oligo.sequence.length - primerLen);
    const viols = findHomopolymers(innerDna, 3);
    
    if (viols.length > 0) {
      hpTotal += viols.length;
      for (const v of viols) {
        if (v.pos < addressNt - 1) hpInAddr++;
        else if (v.pos >= addressNt - 1 && v.pos <= addressNt + 1) hpAtBoundary++;
        else hpInPayload++;
      }
      if (i < 5 || viols.length > 2) {
        console.log(`  Oligo ${i}: ${viols.length} hp violations at positions:`, viols.map(v => `${v.base}×${v.run}@${v.pos}`).join(', '));
        console.log(`    addr region: "${innerDna.slice(0, addressNt)}" payload start: "${innerDna.slice(addressNt, addressNt + 20)}..."`);
      }
    }
  }
  
  console.log(`\nTotal HP violations: ${hpTotal}`);
  console.log(`  In address region: ${hpInAddr}`);
  console.log(`  At boundary: ${hpAtBoundary}`);
  console.log(`  In payload region: ${hpInPayload}`);
}

main().catch(e => { console.error(e); process.exit(1); });
