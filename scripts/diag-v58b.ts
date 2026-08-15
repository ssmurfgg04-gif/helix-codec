import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 8192);
console.log("=== DEFAULT_CONFIG (lowCoverageTrigger=5) - WASM path ===");
const cfg = { ...DEFAULT_CONFIG };
const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
const t0 = Date.now();
const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
console.log("decode ms:", Date.now()-t0, "hash:", dec.hashMatches?"OK":"FAIL", "len:", dec.data?.length, "first 16:", Array.from(dec.data?.slice(0,16) ?? []).map(b=>b.toString(16).padStart(2,'0')).join(' '));
console.log("expected first 16:", Array.from(payload.slice(0,16)).map(b=>b.toString(16).padStart(2,'0')).join(' '));

console.log("\n=== Force JS path (lowCoverageTrigger=0) ===");
const cfg2 = { ...DEFAULT_CONFIG, lowCoverageTrigger: 0 };
const enc2 = await encodeFile(payload, cfg2, { fileName: "test.bin", contentType: "application/octet-stream" });
const sim2 = simulate(enc2.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
const t1 = Date.now();
// Force JS path - use decodeReads directly
const dec2 = await decodeReads(sim2.reads, enc2.encoded.metadata, cfg2, enc2.encoded.forwardPrimer, enc2.encoded.reversePrimer, true);
console.log("decode ms:", Date.now()-t1, "hash:", dec2.hashMatches?"OK":"FAIL", "len:", dec2.data?.length);
console.log("first 16:", Array.from(dec2.data?.slice(0,16) ?? []).map(b=>b.toString(16).padStart(2,'0')).join(' '));

console.log("\n=== Force JS path low coverage trigger = 999 ===");
const cfg3 = { ...DEFAULT_CONFIG, lowCoverageTrigger: 999 };
const enc3 = await encodeFile(payload, cfg3, { fileName: "test.bin", contentType: "application/octet-stream" });
const sim3 = simulate(enc3.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
const t2 = Date.now();
const dec3 = await decodeReadsUltra(sim3.reads, enc3.encoded.metadata, cfg3, enc3.encoded.forwardPrimer, enc3.encoded.reversePrimer);
console.log("decode ms:", Date.now()-t2, "hash:", dec3.hashMatches?"OK":"FAIL", "len:", dec3.data?.length);
console.log("first 16:", Array.from(dec3.data?.slice(0,16) ?? []).map(b=>b.toString(16).padStart(2,'0')).join(' '));
