// Profile encoding stages.
import { encodeFile } from "../src/lib/dna/codec";
import { DEFAULT_CONFIG, computeLayout } from "../src/lib/dna/types";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { ReedSolomon216 } from "../src/lib/dna/reedsolomon216";
import * as fs from "fs";

async function main() {
  console.log("=== Profile encoding ===\n");

  const payload = fs.readFileSync("benchmarks/data/erlich_payload/erlich_payload.bin");
  console.log(`Payload: ${payload.length.toLocaleString()} bytes`);

  const testConfig = { ...DEFAULT_CONFIG, oligoLength: 300, primerLength: 20, outerParityRatio: 0.3 };
  const layout = computeLayout(testConfig);
  console.log(`Layout: payloadBytes=${layout.payloadBytes}, innerParityBytes=${layout.innerParityBytes}`);

  // Profile LDPC constructor
  const t0 = Date.now();
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
  console.log(`LDPC constructor: ${Date.now() - t0}ms`);

  // Profile single LDPC encode
  const testData = new Uint8Array(innerK);
  for (let i = 0; i < innerK; i++) testData[i] = i & 0xff;
  const t1 = Date.now();
  for (let i = 0; i < 10000; i++) ldpc.encode(testData);
  console.log(`10K LDPC encodes: ${Date.now() - t1}ms`);

  // Profile full encodeFile
  console.log(`\nEncoding full payload...`);
  const t2 = Date.now();
  const enc = await encodeFile(payload, testConfig, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Full encode: ${Date.now() - t2}ms`);
  console.log(`  Oligos: ${enc.encoded.oligos.length.toLocaleString()}`);
  console.log(`  Screening retries: ${enc.stats.screeningRetries.toLocaleString()}`);
  console.log(`  Avg retries per oligo: ${(enc.stats.screeningRetries / enc.encoded.oligos.length).toFixed(1)}`);
  console.log(`  Compressed: ${enc.stats.compressedSize.toLocaleString()}`);
  console.log(`  Density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
