import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import * as fs from "fs";

const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 8192); // 8KB
console.log("payload bytes:", payload.length);

const cfg = { ...DEFAULT_CONFIG };
const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
console.log("oligos:", enc.encoded.oligos.length, "density:", enc.stats.netDensityBitsPerNt.toFixed(3));

const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
console.log("simulated reads:", sim.reads.length, "avg cluster:", (sim.reads.length/enc.encoded.oligos.length).toFixed(2));

const t0 = Date.now();
const dec = await decodeReadsUltra(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
console.log("decode ms:", Date.now()-t0, "hash:", dec.hashMatches?"OK":"FAIL");
console.log("stats:", JSON.stringify(dec.stats));
