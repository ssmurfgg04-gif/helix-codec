/**
 * Simple encode→decode roundtrip with synthetic data
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";

async function main() {
  // Test with "Hello World"
  const data = new Uint8Array(Buffer.from("Hello World!"));
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  
  console.log("Input:", data.length, "bytes:", new TextDecoder().decode(data));
  
  const enc = await encodeFile(data, cfg, { fileName: "test", contentType: "text/plain" });
  console.log("Encoded:", enc.stats.oligoCount, "oligos, density:", enc.stats.netDensityBitsPerNt.toFixed(3));
  
  // Check constraints
  let gcV=0, hpV=0;
  for (const o of enc.encoded.oligos) {
    if (o.gc < 0.4 || o.gc > 0.6) gcV++;
    if (o.maxHomopolymer > 3) hpV++;
  }
  console.log("Constraints: gcV=", gcV, "hpV=", hpV);
  
  // Clean simulation
  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });
  console.log("Simulated:", sim.reads.length, "reads");
  
  // Decode
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log("Decode: dataLen=", dec.data?.length, "hashOk=", dec.hashMatches, "error=", dec.error);
  
  if (dec.data) {
    const str = new TextDecoder().decode(dec.data);
    console.log("Decoded:", JSON.stringify(str));
    const match = data.length === dec.data.length && data.every((b, i) => b === dec.data[i]);
    console.log("Roundtrip:", match ? "OK" : "FAIL");
  }
}

main().catch(e => { console.error(e); process.exit(1); });
