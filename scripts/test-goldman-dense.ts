// Test Goldman dense mode encode/decode
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 65536);
  console.log(`Payload: ${payload.length} bytes`);

  for (const mode of ["fast", "dense"] as const) {
    console.log(`\n=== Goldman ${mode} ===`);
    // Goldman fast needs innerNt % 6 == 0; dense needs innerNt % 26 == 0
    // oligoLength=304, primer=20 → innerNt=264 → 264%6=0 ✓, 264%26=10 ✗
    // oligoLength=314, primer=20 → innerNt=274 → 274%26=14 ✗
    // oligoLength=300, primer=20 → innerNt=260 → 260%26=0 ✓
    // But 260%6 = 2 ✗. So we need different oligo lengths per mode.
    const oligoLength = mode === "fast" ? 304 : 300;
    const config = { ...DEFAULT_CONFIG, mappingMode: "goldman" as const, goldmanMode: mode, oligoLength, maxRetries: 1 };
    const enc = await encodeFile(payload, config, { fileName: "test.bin", contentType: "application/octet-stream" });
    console.log(`  Oligos: ${enc.encoded.oligos.length}, Density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);

    const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, config, enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true);
    console.log(`  Decode: ${dec.stats.decodeTimeMs}ms, hash matches: ${dec.hashMatches ? "PASS" : "FAIL"}`);
  }
}
main().catch(e => { console.error(e); process.exit(1); });
