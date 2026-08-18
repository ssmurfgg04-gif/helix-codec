/**
 * Diagnostic: encode → decode (no simulation) for SARS-CoV-2
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
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
  console.log("SARS-CoV-2:", data.length, "bytes");

  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const enc = await encodeFile(data, cfg, { fileName: "sars-cov-2", contentType: "application/octet-stream" });
  console.log("Encoded:", enc.stats.oligoCount, "oligos, density:", enc.stats.netDensityBitsPerNt.toFixed(3));
  console.log("Metadata:", JSON.stringify({
    innerRS: enc.encoded.metadata.innerRS,
    outerRS: enc.encoded.metadata.outerRS,
    mappingMode: enc.encoded.metadata.mappingMode,
    parityOligos: enc.encoded.metadata.parityOligos,
    oligoCount: enc.encoded.metadata.oligoCount,
  }, null, 2));

  // Simulate clean reads
  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });
  console.log("Simulated:", sim.reads.length, "reads");

  // Try decode
  try {
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    console.log("Decode result:", {
      dataLen: dec.data?.length,
      hashMatches: dec.hashMatches,
      error: dec.error,
    });
    if (dec.data && dec.data.length === data.length) {
      let mismatches = 0;
      for (let i = 0; i < data.length; i++) {
        if (dec.data[i] !== data[i]) mismatches++;
      }
      console.log("Mismatches:", mismatches, "of", data.length, "bytes");
    } else {
      console.log("Data length mismatch:", dec.data?.length, "vs expected", data.length);
    }
  } catch (e: any) {
    console.log("Decode error:", e.message);
  }
}

main().catch(e => { console.error(e); process.exit(1); });
