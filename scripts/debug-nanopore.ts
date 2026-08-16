import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";

// Test nanopore without convolutional inner (just LDPC + constrained mapping)
const data = new Uint8Array(1024);
for (let i = 0; i < 1024; i++) data[i] = i & 0xff;

const cfgs = [
  {
    name: "nanopore-no-conv",
    cfg: {
      oligoLength: 300, primerLength: 12, innerCode: "ldpc" as const,
      ldpcDecoder: "auto" as const, mappingMode: "constrained" as const,
      innerParityBytes: 8, outerParityRatio: 0.25,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: true, maxRetries: 1, interleaveDepth: 0,
      channel: "nanopore" as const, lowCoverageTrigger: 5,
      // NO useConvolutionalInner
    },
  },
  {
    name: "nanopore-conv-500nt",
    cfg: {
      oligoLength: 500, primerLength: 12, innerCode: "ldpc" as const,
      ldpcDecoder: "auto" as const, mappingMode: "constrained" as const,
      innerParityBytes: 8, outerParityRatio: 0.25,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: true, maxRetries: 1, interleaveDepth: 0,
      channel: "nanopore" as const, lowCoverageTrigger: 5,
      useConvolutionalInner: true,
    },
  },
];

for (const { name, cfg } of cfgs) {
  console.log(`\n── ${name} ──`);
  try {
    const enc = await encodeFile(data, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
    console.log(`Encoded: ${enc.stats.oligoCount} oligos, ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);

    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: "basic" });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    console.log(`Decode: hash=${dec.hashMatches} rec=${dec.stats.oligosRecovered}/${enc.stats.oligoCount} erased=${dec.stats.oligosErased} failI=${dec.stats.oligosFailedInnerRS} failO=${dec.stats.oligosFailedOuterRS}`);
  } catch (e: any) {
    console.log(`Error: ${e.message}`);
  }
}
