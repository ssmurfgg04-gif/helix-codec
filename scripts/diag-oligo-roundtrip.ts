/**
 * Single oligo roundtrip test: encode → DNA → decode → verify
 */
import { encodeFile } from "../src/lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { bytesToSplitConstrainedDna, splitConstrainedDnaToBytesWithErasure, constrainedDnaToBytes } from "../src/lib/dna/constrained-mapping";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(c: string) { return c.split("\n").filter(l => !l.startsWith(">") && l.trim()).join(""); }

async function main() {
  const dsDir = path.join(__dirname, "..", "datasets");
  const content = fs.readFileSync(path.join(dsDir, "small/uniprot-p00533.fa"), "utf-8");
  const seq = parseFasta(content);
  const data = new Uint8Array(Buffer.from(seq, "utf-8"));
  
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const enc = await encodeFile(data, cfg, { fileName: "test", contentType: "application/octet-stream" });
  
  const primerLen = cfg.primerLength;
  const totalInnerBytes = enc.encoded.metadata.innerRS.n + 2; // innerN + CRC
  
  console.log("Testing oligo-level roundtrip for", enc.encoded.oligos.length, "oligos");
  console.log("totalInnerBytes:", totalInnerBytes, "innerRS:", JSON.stringify(enc.encoded.metadata.innerRS));
  
  let passCount = 0;
  let failCount = 0;
  
  for (let i = 0; i < Math.min(enc.encoded.oligos.length, 5); i++) {
    const oligo = enc.encoded.oligos[i];
    const innerDna = oligo.sequence.slice(primerLen, oligo.sequence.length - primerLen);
    
    // Try split constrained decode with directBytes=0
    try {
      const result = splitConstrainedDnaToBytesWithErasure(innerDna, 3, 0, totalInnerBytes);
      const decoded = result.data;
      const erasures = result.erasures;
      
      // Also encode the decoded bytes back and compare
      const reEncoded = bytesToSplitConstrainedDna(decoded, 3, 0, 0.4, 0.6);
      
      const roundtripOk = reEncoded.dna === innerDna;
      console.log(`  Oligo ${i}: decoded ${decoded.length}B, erasures=${erasures.length}, re-encode match=${roundtripOk}`);
      if (!roundtripOk) {
        // Find first mismatch
        for (let j = 0; j < innerDna.length; j++) {
          if (reEncoded.dna[j] !== innerDna[j]) {
            console.log(`    First mismatch at nt ${j}: encoded='${innerDna[j]}' re-encoded='${reEncoded.dna[j]}'`);
            break;
          }
        }
      }
      if (roundtripOk) passCount++; else failCount++;
    } catch (e: any) {
      console.log(`  Oligo ${i}: DECODE ERROR: ${e.message}`);
      failCount++;
    }
  }
  
  console.log(`\nPass: ${passCount}, Fail: ${failCount}`);
}

main().catch(e => { console.error(e); process.exit(1); });
