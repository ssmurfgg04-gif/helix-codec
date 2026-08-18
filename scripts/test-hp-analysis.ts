/**
 * Analyze where homopolymer violations occur.
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
  const fastaPath = path.join(__dirname, "..", "datasets", "small", "sars-cov-2.fa");
  const content = fs.readFileSync(fastaPath, "utf-8");
  const seq = parseFasta(content);
  const data = new Uint8Array(Buffer.from(seq, "utf-8"));

  const cfg = V51_DEFAULT_CONFIG;
  const result = await encodeFile(data, cfg, { fileName: "sars-cov-2", contentType: "application/octet-stream" });
  const primerLen = cfg.primerLength;

  let addrV = 0, constV = 0, boundV = 0;
  const addressBp = 16; // 4 address bytes = 16 nt

  for (const oligo of result.encoded.oligos) {
    if (oligo.maxHomopolymer <= 3) continue;
    const inner = oligo.sequence.slice(primerLen, oligo.sequence.length - primerLen);

    // Find all runs > 3
    let curRun = 1;
    for (let i = 1; i <= inner.length; i++) {
      if (i < inner.length && inner[i] === inner[i-1]) {
        curRun++;
      } else {
        if (curRun > 3) {
          const start = i - curRun;
          const end = i;
          if (end <= addressBp) addrV++;
          else if (start >= addressBp) constV++;
          else boundV++;
        }
        curRun = 1;
      }
    }
  }

  console.log(`Address region violations: ${addrV}`);
  console.log(`Constrained region violations: ${constV}`);
  console.log(`Boundary-spanning violations: ${boundV}`);

  // Show details for first 3
  let shown = 0;
  for (const oligo of result.encoded.oligos) {
    if (oligo.maxHomopolymer <= 3 || shown >= 3) continue;
    const inner = oligo.sequence.slice(primerLen, oligo.sequence.length - primerLen);
    console.log(`\nOligo ${oligo.index}: maxHp=${oligo.maxHomopolymer}, GC=${oligo.gc.toFixed(3)}`);
    console.log(`  First 20nt: ${inner.slice(0, 20)}`);

    let maxStart = 0, maxLen = 1, curLen = 1, curStart = 0;
    for (let i = 1; i <= inner.length; i++) {
      if (i < inner.length && inner[i] === inner[i-1]) { curLen++; }
      else {
        if (curLen > maxLen) { maxLen = curLen; maxStart = curStart; }
        curLen = 1; curStart = i;
      }
    }
    console.log(`  Longest run: "${inner.slice(maxStart, maxStart + maxLen)}" at pos ${maxStart} (len=${maxLen})`);
    shown++;
  }
}

main().catch(console.error);
