/**
 * Deeper debug: why does 300nt constrained fail outer RS?
 * Compare 300nt vs 700nt decode for the same data.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { readFile } from "fs/promises";
import { join } from "path";
import { createHash } from "crypto";

const TEST_DATA_DIR = join(import.meta.dirname ?? ".", "..", "test-data");

async function main() {
  const buf = await readFile(join(TEST_DATA_DIR, "random_1kb.bin"));
  const data = new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);

  for (const [label, cfg] of Object.entries({
    "300nt-4parity": {
      oligoLength: 300, primerLength: 12, innerCode: "ldpc" as const,
      ldpcDecoder: "auto" as const, mappingMode: "constrained" as const,
      innerParityBytes: 4, outerParityRatio: 0.1,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: true, maxRetries: 1, interleaveDepth: 0,
      channel: "illumina" as const, lowCoverageTrigger: 999,
    },
    "300nt-8parity": {
      oligoLength: 300, primerLength: 12, innerCode: "ldpc" as const,
      ldpcDecoder: "auto" as const, mappingMode: "constrained" as const,
      innerParityBytes: 8, outerParityRatio: 0.1,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: true, maxRetries: 1, interleaveDepth: 0,
      channel: "illumina" as const, lowCoverageTrigger: 999,
    },
    "300nt-8parity-nocomp": {
      oligoLength: 300, primerLength: 12, innerCode: "ldpc" as const,
      ldpcDecoder: "auto" as const, mappingMode: "constrained" as const,
      innerParityBytes: 8, outerParityRatio: 0.1,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: false, maxRetries: 1, interleaveDepth: 0,
      channel: "illumina" as const, lowCoverageTrigger: 999,
    },
    "700nt-8parity": {
      oligoLength: 700, primerLength: 12, innerCode: "ldpc" as const,
      ldpcDecoder: "auto" as const, mappingMode: "constrained" as const,
      innerParityBytes: 8, outerParityRatio: 0.03,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: true, maxRetries: 1, interleaveDepth: 0,
      channel: "illumina" as const, lowCoverageTrigger: 999,
    },
  })) {
    console.log(`\n── ${label} ──`);
    const enc = await encodeFile(data, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 1 });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);

    console.log(`  Meta: oligos=${enc.encoded.metadata.oligoCount} innerRS=(${enc.encoded.metadata.innerRS.n},${enc.encoded.metadata.innerRS.k}) outerRS=(${enc.encoded.metadata.outerRS.n},${enc.encoded.metadata.outerRS.k})`);
    console.log(`  Decode: recovered=${dec.stats.oligosRecovered} erased=${dec.stats.oligosErased} failInner=${dec.stats.oligosFailedInnerRS} failOuter=${dec.stats.oligosFailedOuterRS}`);
    console.log(`  Hash match: ${dec.hashMatches}`);

    // Show per-oligo: which ones have CRC failure, inner RS failure
    const failOligos = dec.perOligo.filter(p => !p.crcPassed || !p.innerRS.success);
    if (failOligos.length > 0) {
      console.log(`  Failed oligos:`);
      for (const p of failOligos) {
        console.log(`    idx=${p.index} reads=${p.readCount} crc=${p.crcPassed} innerRS=${p.innerRS.success} (corrected=${p.innerRS.corrected}) strategy=${p.strategy} isParity=${p.isParity}`);
      }
    }

    if (dec.data) {
      const decHash = createHash("sha256").update(dec.data).digest("hex");
      console.log(`  Decoded hash: ${decHash.slice(0, 32)}...`);
      console.log(`  Expected hash: ${enc.encoded.metadata.fileHash.slice(0, 32)}...`);
    }
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
