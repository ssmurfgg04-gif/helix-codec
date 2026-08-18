/**
 * Detailed decode diagnostic: check clustering and per-oligo decode
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { bytesToSplitConstrainedDna, constrainedDnaToBytes } from "../src/lib/dna/constrained-mapping";
import { dnaToBytes, unwhitenAddress } from "../src/lib/dna/mapping";
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
  console.log("UniProt EGFR:", data.length, "bytes");

  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const enc = await encodeFile(data, cfg, { fileName: "uniprot-egfr", contentType: "application/octet-stream" });
  console.log("Encoded:", enc.stats.oligoCount, "oligos");

  // Manually check roundtrip on first oligo
  const oligo0 = enc.encoded.oligos[0];
  const primerLen = cfg.primerLength;
  const innerDna = oligo0.sequence.slice(primerLen, oligo0.sequence.length - primerLen);
  console.log("\nOligo 0 inner DNA length:", innerDna.length, "nt");
  console.log("Oligo 0 index:", oligo0.index, "gc:", oligo0.gc.toFixed(3), "seed:", oligo0.seed);

  // Try constrained decode of the inner DNA
  try {
    const decoded = constrainedDnaToBytes(innerDna, 3);
    console.log("Constrained decode OK:", decoded.length, "bytes");
    console.log("First 4 bytes (address):", Array.from(decoded.slice(0, 4)).map(b => b.toString(16).padStart(2, '0')).join(' '));
    const unwhitened = unwhitenAddress(decoded.slice(0, 4));
    const index = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
    console.log("Extracted index:", index, "(expected:", oligo0.index, ")");
  } catch (e: any) {
    console.log("Constrained decode FAILED:", e.message);
  }

  // Try split constrained decode
  try {
    const { splitConstrainedDnaToBytesWithErasure } = await import("../src/lib/dna/constrained-mapping");
    const result = splitConstrainedDnaToBytesWithErasure(innerDna, 3, 0, decoded.length || 66);
    console.log("Split constrained decode OK:", result.data.length, "bytes, erasures:", result.erasures);
  } catch (e: any) {
    console.log("Split constrained decode FAILED:", e.message);
  }

  // Full roundtrip with simulation
  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });
  console.log("\nSimulated:", sim.reads.length, "reads");

  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log("Decode: dataLen:", dec.data?.length, "hash:", dec.hashMatches);
  if (dec.data && dec.data.length === data.length) {
    let mm = 0;
    for (let i = 0; i < data.length; i++) if (dec.data[i] !== data[i]) mm++;
    console.log("Mismatches:", mm);
  }
}

main().catch(e => { console.error(e); process.exit(1); });
