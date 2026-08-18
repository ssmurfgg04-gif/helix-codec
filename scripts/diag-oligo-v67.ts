/**
 * Oligo-level encode → decode roundtrip diagnostic
 */
import { encodeFile } from "../src/lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { bytesToSplitConstrainedDna, splitConstrainedDnaToBytesWithErasure } from "../src/lib/dna/constrained-mapping";
import { addressToHomopolymerSafeDna, homopolymerSafeDnaToAddress, unwhitenAddress, whitenAddress } from "../src/lib/dna/mapping";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(c: string) { return c.split("\n").filter(l => !l.startsWith(">") && l.trim()).join(""); }

async function main() {
  // Test homopolymer-safe address encode/decode roundtrip first
  console.log("=== Address encode/decode roundtrip ===");
  for (const idx of [0, 1, 100, 1000, 50000, 255*256*256-1]) {
    const addr = new Uint8Array(4);
    addr[0] = (idx >> 16) & 0xff;
    addr[1] = (idx >> 8) & 0xff;
    addr[2] = idx & 0xff;
    addr[3] = 0;
    const whitened = whitenAddress(addr);
    const dna = addressToHomopolymerSafeDna(whitened);
    const decoded = homopolymerSafeDnaToAddress(dna);
    const match = whitened.every((b, i) => b === decoded[i]);
    console.log(`  idx=${idx}: whitened=[${Array.from(whitened).map(b=>b.toString(16).padStart(2,'0')).join(',')}] dna="${dna}" decode_ok=${match}`);
  }

  // Test full oligo-level roundtrip
  console.log("\n=== Oligo-level roundtrip ===");
  const dsDir = path.join(__dirname, "..", "datasets");
  const content = fs.readFileSync(path.join(dsDir, "small/uniprot-p00533.fa"), "utf-8");
  const seq = parseFasta(content);
  const data = new Uint8Array(Buffer.from(seq, "utf-8"));
  
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const enc = await encodeFile(data, cfg, { fileName: "test", contentType: "application/octet-stream" });
  
  const primerLen = cfg.primerLength;
  const totalInnerBytes = enc.encoded.metadata.innerRS.n + 2;
  
  let passCount = 0;
  let failCount = 0;
  
  for (let i = 0; i < enc.encoded.oligos.length; i++) {
    const oligo = enc.encoded.oligos[i];
    const innerDna = oligo.sequence.slice(primerLen, oligo.sequence.length - primerLen);
    
    try {
      const result = splitConstrainedDnaToBytesWithErasure(innerDna, 3, 4, totalInnerBytes);
      const decoded = result.data;
      const erasures = result.erasures.filter(e => e).length;
      
      // Check address
      const decodedAddr = decoded.slice(0, 4);
      const expectedAddr = new Uint8Array(4);
      expectedAddr[0] = (oligo.index >> 16) & 0xff;
      expectedAddr[1] = (oligo.index >> 8) & 0xff;
      expectedAddr[2] = oligo.index & 0xff;
      expectedAddr[3] = oligo.seed;
      const expectedWhitened = whitenAddress(expectedAddr);
      const addrMatch = expectedWhitened.every((b, j) => b === decodedAddr[j]);
      
      // Re-encode and compare
      const reEncoded = bytesToSplitConstrainedDna(decoded, 3, 4, 0.4, 0.6);
      const roundtripOk = reEncoded.dna === innerDna;
      
      if (!addrMatch || !roundtripOk) {
        console.log(`  Oligo ${i}: addr_match=${addrMatch} reencode_match=${roundtripOk} erasures=${erasures}`);
        if (!addrMatch) {
          console.log(`    Expected addr: [${Array.from(expectedWhitened).map(b=>b.toString(16).padStart(2,'0')).join(',')}]`);
          console.log(`    Decoded  addr: [${Array.from(decodedAddr).map(b=>b.toString(16).padStart(2,'0')).join(',')}]`);
        }
        failCount++;
      } else {
        passCount++;
      }
    } catch (e: any) {
      console.log;(`  Oligo ${i}: ERROR: ${e.message}`);
      failCount++;
    }
  }
  
  console.log(`\nPass: ${passCount}, Fail: ${failCount}`);
}

main().catch(e => { console.error(e); process.exit(1); });
