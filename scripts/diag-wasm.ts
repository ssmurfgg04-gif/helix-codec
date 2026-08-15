// Test WASM full_decode with synthetic input to identify the failure mode
import { encodeFile } from "../src/lib/dna/codec";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG, computeLayoutAuto } from "../src/lib/dna/types";
import { fullDecode } from "../src/lib/dna/wasm-batch-decode";
import * as fs from "fs";

const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 1024); // tiny
const cfg = { ...DEFAULT_CONFIG };
const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
console.log("oligos:", enc.encoded.oligos.length);
console.log("first oligo:", enc.encoded.oligos[0]);
console.log("metadata:", JSON.stringify({
  oligoCount: enc.encoded.metadata.oligoCount,
  mappingMode: enc.encoded.metadata.mappingMode,
  outerRS: enc.encoded.metadata.outerRS,
  fileSize: enc.encoded.metadata.fileSize,
  compression: enc.encoded.metadata.compression,
  interleaveDepth: enc.encoded.metadata.interleaveDepth,
}, null, 2));
console.log("fwd primer:", enc.encoded.forwardPrimer);
console.log("rev primer:", enc.encoded.reversePrimer);

const layout = computeLayoutAuto(cfg);
console.log("layout:", JSON.stringify(layout));
const innerK = layout.addressBytes + layout.payloadBytes;
const innerN = innerK + layout.innerParityBytes;
console.log("innerK:", innerK, "innerN:", innerN, "totalInnerBytes:", layout.totalInnerBytes, "payloadBytes:", layout.payloadBytes);

const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
console.log("reads:", sim.reads.length, "first read seq:", sim.reads[0].sequence.slice(0,40) + "...");

const t0 = Date.now();
const result = fullDecode(
  sim.reads,
  enc.encoded.forwardPrimer,
  enc.encoded.reversePrimer,
  enc.encoded.metadata.oligoCount,
  innerN, innerK, layout.totalInnerBytes,
  enc.encoded.metadata.outerRS.n, enc.encoded.metadata.outerRS.k,
  layout.payloadBytes,
  enc.encoded.metadata.fileSize,
  enc.encoded.metadata.compression === "deflate",
);
console.log("WASM decode ms:", Date.now()-t0, "result len:", result.length);
console.log("result first 16:", Array.from(result.slice(0,16)).map(b=>b.toString(16).padStart(2,'0')).join(' '));
console.log("expected first 16:", Array.from(payload.slice(0,16)).map(b=>b.toString(16).padStart(2,'0')).join(' '));
console.log("result all zero?:", result.every(b=>b===0));
