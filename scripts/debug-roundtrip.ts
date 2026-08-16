/**
 * Debug round-trip failure for default-constrained at 1KB.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { readFile } from "fs/promises";
import { join } from "path";
import { createHash } from "crypto";

const TEST_DATA_DIR = join(import.meta.dirname ?? ".", "..", "test-data");

const cfg = {
  oligoLength: 300,
  primerLength: 12,
  innerCode: "ldpc" as const,
  ldpcDecoder: "auto" as const,
  mappingMode: "constrained" as const,
  innerParityBytes: 4,
  outerParityRatio: 0.1,
  constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
  compress: true,
  maxRetries: 1,
  interleaveDepth: 0,
  channel: "illumina" as const,
  lowCoverageTrigger: 5,
};

async function main() {
  // Load 1KB random
  const buf = await readFile(join(TEST_DATA_DIR, "random_1kb.bin"));
  const data = new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
  const origHash = createHash("sha256").update(data).digest("hex");
  console.log(`Input: ${data.length} bytes, SHA256: ${origHash.slice(0, 32)}...`);

  // Encode
  const enc = await encodeFile(data, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`Encode: ${enc.stats.oligoCount} oligos, density ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt, ${enc.stats.encodeTimeMs}ms`);
  console.log(`  Meta: fileSize=${enc.encoded.metadata.fileSize}, oligoCount=${enc.encoded.metadata.oligoCount}, compression=${enc.encoded.metadata.compression}`);
  console.log(`  Meta: innerRS=(${enc.encoded.metadata.innerRS.n},${enc.encoded.metadata.innerRS.k}), outerRS=(${enc.encoded.metadata.outerRS.n},${enc.encoded.metadata.outerRS.k})`);
  console.log(`  Meta: payloadBytesPerOligo=${enc.encoded.metadata.payloadBytesPerOligo}`);
  console.log(`  Meta: fileHash=${enc.encoded.metadata.fileHash.slice(0, 32)}...`);

  // Simulate clean
  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
  console.log(`\nSimulation: ${sim.reads.length} reads, ${sim.totalErrors} errors`);

  // Decode
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`\nDecode: ${dec.stats.decodeTimeMs}ms`);
  console.log(`  totalReads=${dec.stats.totalReads}, readsUsed=${dec.stats.readsUsed}`);
  console.log(`  clustersFormed=${dec.stats.clustersFormed}, oligosRecovered=${dec.stats.oligosRecovered}`);
  console.log(`  oligosErased=${dec.stats.oligosErased}, oligosFailedInnerRS=${dec.stats.oligosFailedInnerRS}, oligosFailedOuterRS=${dec.stats.oligosFailedOuterRS}`);
  console.log(`  consensusSuccessRate=${dec.stats.consensusSuccessRate.toFixed(3)}`);
  console.log(`  hashMatches=${dec.hashMatches}`);
  console.log(`  data is null: ${dec.data === null}`);
  if (dec.data) {
    const decHash = createHash("sha256").update(dec.data).digest("hex");
    console.log(`  decoded size: ${dec.data.length}, hash: ${decHash.slice(0, 32)}...`);
    console.log(`  hash match: ${decHash === enc.encoded.metadata.fileHash}`);

    // Byte diff
    if (dec.data.length === data.length) {
      let diffs = 0;
      let firstDiff = -1;
      for (let i = 0; i < data.length; i++) {
        if (dec.data[i] !== data[i]) {
          diffs++;
          if (firstDiff === -1) firstDiff = i;
        }
      }
      console.log(`  byte diffs: ${diffs}/${data.length}, first at offset ${firstDiff}`);
    } else {
      console.log(`  SIZE MISMATCH: ${dec.data.length} vs ${data.length}`);
    }
  }

  // Per-oligo strategies
  const strategies = new Map<string, number>();
  for (const p of dec.perOligo) {
    const s = p.strategy ?? "none";
    strategies.set(s, (strategies.get(s) ?? 0) + 1);
  }
  console.log(`  decode strategies:`, Object.fromEntries(strategies));
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
