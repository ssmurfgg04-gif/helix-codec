// Test with a larger, more realistic payload.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";

async function main() {
  // Generate a 10KB payload with mixed content
  const text = `The Synthetic DNA Data Storage Codec is a system for encoding digital information into synthetic DNA sequences. DNA has an information density of approximately 455 exabytes per gram, making it an attractive medium for long-term archival storage. Unlike magnetic tape or optical disks, DNA can last for thousands of years if stored in cool, dry conditions.

The codec implements a multi-layer error correction strategy. At the inner level, each oligo is protected by Reed-Solomon codes that can correct substitution errors introduced during synthesis and sequencing. At the outer level, erasure coding across oligos allows the system to recover from complete strand loss (dropouts). A CRC-16 checksum per oligo provides an additional layer of integrity verification.

Biological constraints are enforced through screening: each oligo must have a GC content between 40 and 60 percent, and no homopolymer run longer than 3 nucleotides. Oligos that fail screening are re-encoded with a pseudo-random XOR mask derived from a seed, and the seed is stored in the oligo header for decoder reversal.

The system models realistic sequencing errors from three platforms: Illumina (short-read, substitution-dominant), Oxford Nanopore (long-read, indel-heavy), and Pacific Biosciences (insertion-dominant). Coverage depth can be adjusted from 1x to 50x, with 20x being typical for Illumina DNA storage experiments.

Recovery proceeds through clustering (grouping reads by oligo index), consensus (column-wise plurality vote), inner RS decoding, outer RS erasure decoding, and finally SHA-256 hash verification. The entire pipeline is deterministic given the same random seed.`;

  const payload = new TextEncoder().encode(text.repeat(5));
  console.log(`Test payload: ${payload.length} bytes`);

  const encodeResult = await encodeFile(payload, DEFAULT_CONFIG, {
    fileName: "dna-storage-essay.txt",
    contentType: "text/plain",
  });

  console.log(`Encoded: ${encodeResult.stats.oligoCount} oligos`);
  console.log(`  Compressed: ${encodeResult.stats.compressedSize} bytes`);
  console.log(`  Net density: ${encodeResult.stats.netDensityBitsPerNt.toFixed(3)} bits/nt`);
  console.log(`  Overhead: ${encodeResult.stats.overheadPercent.toFixed(1)}%`);
  console.log(`  Screening retries: ${encodeResult.stats.screeningRetries}`);

  // Test 1: Clean
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 1, substitutionRate: 0, insertionRate: 0, deletionRate: 0 });
    const dec = await decodeReads(sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG, encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer);
    console.log(`\nClean decode: ${dec.hashMatches ? "PASS" : "FAIL"} (${sim.totalReads} reads, ${dec.stats.decodeTimeMs}ms)`);
  }

  // Test 2: Illumina 20x
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_ILLUMINA, seed: 42 });
    const dec = await decodeReads(sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG, encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer);
    console.log(`Illumina 20x: ${dec.hashMatches ? "PASS" : "FAIL"} (${sim.totalReads} reads, ${sim.totalErrors} errors, ${dec.stats.decodeTimeMs}ms)`);
    console.log(`  Recovered: ${dec.stats.oligosRecovered}/${encodeResult.encoded.metadata.oligoCount}, erased: ${dec.stats.oligosErased}, inner RS fails: ${dec.stats.oligosFailedInnerRS}`);
  }

  // Test 3: Nanopore 15x (high indel)
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_NANOPORE, seed: 42 });
    const dec = await decodeReads(sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG, encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer);
    console.log(`Nanopore 15x: ${dec.hashMatches ? "PASS" : "FAIL"} (${sim.totalReads} reads, ${sim.totalErrors} errors, ${dec.stats.decodeTimeMs}ms)`);
    console.log(`  Recovered: ${dec.stats.oligosRecovered}/${encodeResult.encoded.metadata.oligoCount}, erased: ${dec.stats.oligosErased}, inner RS fails: ${dec.stats.oligosFailedInnerRS}`);
  }

  // Test 4: With dropouts
  {
    const sim = simulate(encodeResult.encoded.oligos, { ...PRESET_ILLUMINA, dropoutRate: 0.15, seed: 42 });
    const dec = await decodeReads(sim.reads, encodeResult.encoded.metadata, DEFAULT_CONFIG, encodeResult.encoded.forwardPrimer, encodeResult.encoded.reversePrimer);
    console.log(`Illumina 20x + 15% dropout: ${dec.hashMatches ? "PASS" : "FAIL"} (${sim.totalReads} reads, ${sim.droppedOligos.length} dropped, ${dec.stats.decodeTimeMs}ms)`);
    console.log(`  Recovered: ${dec.stats.oligosRecovered}, erased: ${dec.stats.oligosErased}, outer RS fails: ${dec.stats.oligosFailedOuterRS}`);
  }
}

main().catch(console.error);
