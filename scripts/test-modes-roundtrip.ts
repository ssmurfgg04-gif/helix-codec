/**
 * Test roundtrip with different mapping modes
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import type { CodecConfig } from "../src/lib/dna/types";

async function testRoundtrip(data: Uint8Array, cfg: CodecConfig, label: string) {
  const enc = await encodeFile(data, cfg, { fileName: label, contentType: "application/octet-stream" });
  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const ok = dec.data && dec.data.length === data.length && data.every((b, i) => b === dec.data![i]);
  console.log(`  ${label}: mapping=${cfg.mappingMode} oligos=${enc.stats.oligoCount} hash=${dec.hashMatches} roundtrip=${ok ? 'OK' : 'FAIL'}`);
  return ok;
}

async function main() {
  const data = new Uint8Array(Buffer.from("Hello World! Test 12345"));
  
  // Test yinyang mode (DEFAULT_CONFIG)
  const yinyangOk = await testRoundtrip(data, DEFAULT_CONFIG, "yinyang");
  
  // Test constrained mode (V51_DEFAULT_CONFIG)  
  const constrainedOk = await testRoundtrip(data, { ...V51_DEFAULT_CONFIG, maxRetries: 10 }, "constrained");
  
  // Test SRT mode
  const srtOk = await testRoundtrip(data, { ...V51_DEFAULT_CONFIG, mappingMode: "srt", maxRetries: 10 }, "srt");
  
  console.log(`\nResults: yinyang=${yinyangOk} constrained=${constrainedOk} srt=${srtOk}`);
}

main().catch(e => { console.error(e); process.exit(1); });
