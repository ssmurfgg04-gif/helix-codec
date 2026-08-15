/**
 * v59 Debug path — check which decode path is used.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { computeLayoutAuto } from "../src/lib/dna/types";
import * as fs from "fs";

const TAG = "[v59-path]";

async function main() {
  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin").slice(0, 262144);

  const cfg = ULTIMATE_V55_DENSITY_CONFIG;
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });

  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });

  // Check which path decodeReadsUltra would take
  const layout = computeLayoutAuto(cfg);
  const useArithmetic = (enc.encoded.metadata.mappingMode ?? "direct") === "arithmetic";
  const useConvInner = !!enc.encoded.metadata.useConvolutionalInner;
  const channel = enc.encoded.metadata.channel ?? cfg.channel ?? "illumina";
  const lowCovTrigger = enc.encoded.metadata.lowCoverageTrigger ?? cfg.lowCoverageTrigger ?? 5;
  const avgClusterSize = enc.encoded.metadata.oligoCount > 0 ? sim.reads.length / enc.encoded.metadata.oligoCount : 0;
  const useLowCoveragePath = lowCovTrigger > 0 && avgClusterSize < lowCovTrigger;
  const forceJsPath = useArithmetic || channel === "nanopore";

  console.log(`${TAG} Path decision:`);
  console.log(`${TAG}   mappingMode: ${enc.encoded.metadata.mappingMode}`);
  console.log(`${TAG}   useArithmetic: ${useArithmetic}`);
  console.log(`${TAG}   useConvInner: ${useConvInner}`);
  console.log(`${TAG}   channel: ${channel}`);
  console.log(`${TAG}   lowCovTrigger: ${lowCovTrigger}`);
  console.log(`${TAG}   avgClusterSize: ${avgClusterSize.toFixed(1)}`);
  console.log(`${TAG}   useLowCoveragePath: ${useLowCoveragePath}`);
  console.log(`${TAG}   forceJsPath: ${forceJsPath}`);
  console.log(`${TAG}   → Will use: ${(useLowCoveragePath || useConvInner || forceJsPath) ? "JS path" : "WASM path"}`);

  console.log(`\n${TAG} Layout:`);
  console.log(`${TAG}   innerN: ${layout.addressBytes + layout.payloadBytes + layout.innerParityBytes}`);
  console.log(`${TAG}   innerK: ${layout.addressBytes + layout.payloadBytes}`);
  console.log(`${TAG}   totalInnerBytes: ${layout.totalInnerBytes}`);
  console.log(`${TAG}   payloadBytes: ${layout.payloadBytes}`);

  console.log(`\n${TAG} Metadata:`);
  console.log(`${TAG}   oligoCount: ${enc.encoded.metadata.oligoCount}`);
  console.log(`${TAG}   outerRS: n=${enc.encoded.metadata.outerRS.n}, k=${enc.encoded.metadata.outerRS.k}`);
  console.log(`${TAG}   fileSize: ${enc.encoded.metadata.fileSize}`);
  console.log(`${TAG}   compression: ${enc.encoded.metadata.compression}`);
}
main().catch(e => { console.error(e); process.exit(1); });
