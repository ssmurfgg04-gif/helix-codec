/**
 * Analyze yinyang constraint violations
 */
import { encodeFile } from "../src/lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(c: string) { return c.split("\n").filter(l => !l.startsWith(">") && l.trim()).join(""); }

async function main() {
  const dsDir = path.join(__dirname, "..", "datasets");
  const content = fs.readFileSync(path.join(dsDir, "small/sars-cov-2.fa"), "utf-8");
  const seq = parseFasta(content);
  const data = new Uint8Array(Buffer.from(seq, "utf-8"));
  
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const enc = await encodeFile(data, cfg, { fileName: "Cov-2", contentType: "application/octet-stream" });
  
  const primerLen = cfg.primerLength;
  let count = 0;
  
  for (const oligo of enc.encoded.oligos) {
    if (oligo.maxHomopolymer > 3 || oligo.gc < 0.4 || oligo.gc > 0.6) {
      if (count < 3) {
        const seq = oligo.sequence;
        const inner = seq.slice(primerLen, seq.length - primerLen);
        // Find homopolymers
        let maxRun = 1, curRun = 1, maxPos = 0;
        for (let i = 1; i < inner.length; i++) {
          if (inner[i] === inner[i-1]) { curRun++; if (curRun > maxRun) { maxRun = curRun; maxPos = i - curRun + 1; } }
          else curRun = 1;
        }
        console.log(`Oligo ${oligo.index}: gc=${oligo.gc.toFixed(3)} maxHp=${oligo.maxHomopolymer} hp@${maxPos}="${inner.slice(Math.max(0,maxPos-2),maxPos+6)}" innerLen=${inner.length}`);
      }
      count++;
    }
  }
  console.log(`Total violations: ${count} / ${enc.encoded.oligos.length}`);
}

main().catch(e => { console.error(e); process.exit(1); });
