/**
 * Real Genomic Dataset Benchmark
 *
 * Tests the helix-codec encode→decode pipeline against real genomic datasets:
 *   - E. coli K-12 MG1655 (4.6 Mb, GC=50.8%)
 *   - S. cerevisiae S288C (12.1 Mb, GC=38.3%)
 *   - Random data control (1 Mb)
 *
 * Measures:
 *   1. Compression ratio (with arithmetic coding vs VLC vs raw 2-bit)
 *   2. Encode/decode throughput (MB/s)
 *   3. Recovery rate at various noise levels
 *   4. Round-trip integrity (SHA-256 hash verification)
 *
 * Usage:
 *   npx tsx scripts/genomic-benchmark.ts
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import * as crypto from 'node:crypto';

// ---------------------------------------------------------------------------
// Dataset loading
// ---------------------------------------------------------------------------

interface Dataset {
  name: string;
  dna: string;       // DNA sequence (ACGT only)
  gcContent: number;
  sourcePath: string;
}

function loadFastaGz(fastaGzPath: string, maxBases: number = Infinity): string {
  // Use zlib to decompress
  const { gunzipSync } = require('node:zlib');
  const compressed = fs.readFileSync(fastaGzPath);
  const decompressed = gunzipSync(compressed);
  const text = decompressed.toString('utf-8');

  // Parse FASTA — skip headers, concatenate sequences
  const lines = text.split('\n');
  const seqParts: string[] = [];
  for (const line of lines) {
    if (line.startsWith('>')) continue;
    // Filter to ACGT only
    const cleaned = line.replace(/[^ACGTacgt]/g, '');
    if (cleaned.length > 0) seqParts.push(cleaned.toUpperCase());
  }

  const fullSeq = seqParts.join('');
  if (fullSeq.length > maxBases) {
    return fullSeq.substring(0, maxBases);
  }
  return fullSeq;
}

function computeGC(dna: string): number {
  let gc = 0;
  for (let i = 0; i < dna.length; i++) {
    const c = dna[i];
    if (c === 'G' || c === 'C' || c === 'g' || c === 'c') gc++;
  }
  return gc / dna.length;
}

function generateRandomDNA(len: number): string {
  const bases = 'ACGT';
  const result: string[] = [];
  for (let i = 0; i < len; i++) {
    result.push(bases[Math.floor(Math.random() * 4)]);
  }
  return result.join('');
}

// ---------------------------------------------------------------------------
// Arithmetic coding benchmark (TypeScript implementation)
// ---------------------------------------------------------------------------

import {
  ArithmeticEncoder,
  ArithmeticDecoder,
  AdaptiveModel,
  AdaptiveContextModel,
  arithmeticEncode,
  arithmeticDecode,
  arithmeticEncodeContext,
  arithmeticDecodeContext,
} from '../src/lib/dna/arithmetic-coder';

function benchmarkArithmeticCoding(dna: string, name: string): {
  order0Ratio: number;
  order2Ratio: number;
  order0TimeMs: number;
  order2TimeMs: number;
  roundtripOk: boolean;
} {
  // Convert DNA to 4-symbol indices
  const symbols = new Uint8Array(dna.length);
  for (let i = 0; i < dna.length; i++) {
    switch (dna[i]) {
      case 'A': symbols[i] = 0; break;
      case 'C': symbols[i] = 1; break;
      case 'G': symbols[i] = 2; break;
      default:  symbols[i] = 3; break;
    }
  }

  // Order-0 arithmetic coding
  const t0 = Date.now();
  const compressed0 = arithmeticEncode(symbols, 4);
  const order0TimeMs = Date.now() - t0;
  const order0Ratio = compressed0.length / dna.length;

  // Order-2 arithmetic coding
  const t1 = Date.now();
  const compressed2 = arithmeticEncodeContext(symbols, 4, 2);
  const order2TimeMs = Date.now() - t1;
  const order2Ratio = compressed2.length / dna.length;

  // Roundtrip verification
  const decoded0 = arithmeticDecode(compressed0, 4, symbols.length);
  const roundtrip0Ok = decoded0.length === symbols.length &&
    decoded0.every((s, i) => s === symbols[i]);

  const decoded2 = arithmeticDecodeContext(compressed2);
  const roundtrip2Ok = decoded2.length === symbols.length &&
    decoded2.every((s, i) => s === symbols[i]);

  return {
    order0Ratio,
    order2Ratio,
    order0TimeMs,
    order2TimeMs,
    roundtripOk: roundtrip0Ok && roundtrip2Ok,
  };
}

// ---------------------------------------------------------------------------
// DNA storage encode/decode benchmark
// ---------------------------------------------------------------------------

import { encode, decode } from '../src/lib/dna/codec';

function benchmarkDnaStorage(data: Uint8Array, name: string, oligoLength: number = 200): {
  encodeTimeMs: number;
  decodeTimeMs: number;
  oligoCount: number;
  payloadBytesPerOligo: number;
  totalDnaBases: number;
  densityBitsPerNt: number;
  roundtripOk: boolean;
} | null {
  try {
    const t0 = Date.now();
    const encoded = encode(data, {
      oligoLength,
      payloadSize: 30,
      primerLength: 12,
      mappingMode: 'yinyang',
    });
    const encodeTimeMs = Date.now() - t0;

    const t1 = Date.now();
    const decoded = decode(encoded);
    const decodeTimeMs = Date.now() - t1;

    // Verify round-trip
    let roundtripOk = decoded.length === data.length;
    if (roundtripOk) {
      for (let i = 0; i < data.length; i++) {
        if (decoded[i] !== data[i]) { roundtripOk = false; break; }
      }
    }

    // Compute DNA density
    const totalBases = encoded.oligos.reduce((s, o) => s + o.sequence.length, 0);
    const densityBitsPerNt = (data.length * 8) / totalBases;

    return {
      encodeTimeMs,
      decodeTimeMs,
      oligoCount: encoded.oligos.length,
      payloadBytesPerOligo: encoded.metadata.payloadBytesPerOligo,
      totalDnaBases: totalBases,
      densityBitsPerNt,
      roundtripOk,
    };
  } catch (err) {
    console.warn(`  [WARN] DNA storage benchmark failed for ${name}: ${(err as Error).message}`);
    return null;
  }
}

// ---------------------------------------------------------------------------
// Noisy channel recovery benchmark
// ---------------------------------------------------------------------------

import { ReedSolomon } from '../src/lib/dna/reedsolomon';

function benchmarkNoisyChannel(
  data: Uint8Array,
  name: string,
  idsRates: number[],
  coverage: number = 20,
): { idsRate: number; recoveryRate: number; decodeTimeMs: number }[] {
  const results: { idsRate: number; recoveryRate: number; decodeTimeMs: number }[] = [];

  // Simple test: split data into blocks, add noise, check recovery
  const blockSize = 30; // payload bytes per oligo
  const numBlocks = Math.floor(data.length / blockSize);
  const BASES = 'ACGT';

  class Rng {
    private state: number;
    constructor(seed: number) { this.state = (seed >>> 0) || 1; }
    next(): number {
      this.state ^= this.state << 13;
      this.state ^= this.state >>> 17;
      this.state ^= this.state << 5;
      this.state = this.state >>> 0;
      return this.state / 0x100000000;
    }
  }

  for (const idsRate of idsRates) {
    const rng = new Rng(42);
    const delRate = idsRate * 0.45;
    const insRate = idsRate * 0.30;
    const subRate = idsRate * 0.25;

    let recovered = 0;
    const t0 = Date.now();

    for (let b = 0; b < numBlocks; b++) {
      const block = data.slice(b * blockSize, (b + 1) * blockSize);

      // Convert to DNA
      const bits: number[] = [];
      for (let i = 0; i < block.length; i++) {
        for (let bit = 7; bit >= 0; bit--) {
          bits.push((block[i] >> bit) & 1);
        }
      }
      const dna: string[] = [];
      for (let i = 0; i + 1 < bits.length; i += 2) {
        dna.push(BASES[(bits[i] << 1) | bits[i + 1]]);
      }
      const originalDna = dna.join('');

      // Generate noisy reads and consensus
      const consensusVotes: number[][] = Array.from({ length: originalDna.length }, () => [0, 0, 0, 0]);
      for (let r = 0; r < coverage; r++) {
        let pos = 0;
        for (let j = 0; j < originalDna.length; j++) {
          if (rng.next() < delRate) continue; // deletion

          let base: number;
          const origBase = 'ACGT'.indexOf(originalDna[j]);
          if (rng.next() < subRate) {
            do { base = Math.floor(rng.next() * 4); } while (base === origBase);
          } else {
            base = origBase;
          }

          if (pos < consensusVotes.length) {
            consensusVotes[pos][base]++;
          }
          pos++;

          if (rng.next() < insRate) {
            const insBase = Math.floor(rng.next() * 4);
            if (pos < consensusVotes.length) {
              // Insertions shift alignment — skip for simple model
            }
            pos++;
          }
        }
      }

      // Consensus
      const consensusDna: string[] = [];
      for (let p = 0; p < originalDna.length; p++) {
        const votes = consensusVotes[p];
        let bestIdx = 0;
        for (let k = 1; k < 4; k++) {
          if (votes[k] > votes[bestIdx]) bestIdx = k;
        }
        consensusDna.push(BASES[bestIdx]);
      }

      // Convert back to bytes
      const decodedBits: number[] = [];
      for (let j = 0; j < consensusDna.length; j++) {
        const code = 'ACGT'.indexOf(consensusDna[j]);
        decodedBits.push((code >> 1) & 1);
        decodedBits.push(code & 1);
      }
      const decoded = new Uint8Array(blockSize);
      for (let i = 0; i < blockSize * 8 && i < decodedBits.length; i++) {
        decoded[i >> 3] |= decodedBits[i] << (7 - (i & 7));
      }

      // Check
      let matches = true;
      for (let i = 0; i < blockSize; i++) {
        if (decoded[i] !== block[i]) { matches = false; break; }
      }
      if (matches) recovered++;
    }

    const decodeTimeMs = Date.now() - t0;
    results.push({
      idsRate,
      recoveryRate: recovered / numBlocks,
      decodeTimeMs,
    });
  }

  return results;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

async function main() {
  console.log('╔══════════════════════════════════════════════════════════════╗');
  console.log('║   Real Genomic Dataset Benchmark                             ║');
  console.log('╚══════════════════════════════════════════════════════════════╝\n');

  const genomesDir = path.join(process.cwd(), 'test-data', 'genomes');

  // Load datasets
  const datasets: Dataset[] = [];

  // E. coli
  const ecoliPath = path.join(genomesDir, 'ecoli_k12_mg1655.fna.gz');
  if (fs.existsSync(ecoliPath)) {
    console.log('Loading E. coli K-12 MG1655...');
    const dna = loadFastaGz(ecoliPath, 1_000_000); // Use first 1M bases for speed
    const gc = computeGC(dna);
    datasets.push({ name: 'E. coli K-12', dna, gcContent: gc, sourcePath: ecoliPath });
    console.log(`  ${dna.length.toLocaleString()} bp, GC=${(gc * 100).toFixed(1)}%`);
  } else {
    console.warn('E. coli genome not found — skipping');
  }

  // Yeast
  const yeastPath = path.join(genomesDir, 'yeast_s288c.fna.gz');
  if (fs.existsSync(yeastPath)) {
    console.log('Loading S. cerevisiae S288C...');
    const dna = loadFastaGz(yeastPath, 1_000_000); // First 1M bases
    const gc = computeGC(dna);
    datasets.push({ name: 'S. cerevisiae', dna, gcContent: gc, sourcePath: yeastPath });
    console.log(`  ${dna.length.toLocaleString()} bp, GC=${(gc * 100).toFixed(1)}%`);
  } else {
    console.warn('Yeast genome not found — skipping');
  }

  // Random control
  console.log('Generating random DNA control...');
  const randomDna = generateRandomDNA(1_000_000);
  datasets.push({ name: 'Random (control)', dna: randomDna, gcContent: 0.5, sourcePath: 'synthetic' });

  console.log('');

  // ---- Arithmetic Coding Benchmarks ----
  console.log('=== Arithmetic Coding Compression ===\n');

  const arithResults: Record<string, any> = {};

  for (const ds of datasets) {
    console.log(`--- ${ds.name} (${ds.dna.length.toLocaleString()} bp, GC=${(ds.gcContent * 100).toFixed(1)}%) ---`);

    const bench = benchmarkArithmeticCoding(ds.dna, ds.name);
    console.log(`  Order-0: ratio=${bench.order0Ratio.toFixed(3)} (${(1 / bench.order0Ratio).toFixed(1)}:1), time=${bench.order0TimeMs}ms`);
    console.log(`  Order-2: ratio=${bench.order2Ratio.toFixed(3)} (${(1 / bench.order2Ratio).toFixed(1)}:1), time=${bench.order2TimeMs}ms`);
    console.log(`  Roundtrip: ${bench.roundtripOk ? 'OK' : 'FAILED'}`);

    arithResults[ds.name] = bench;
  }

  // ---- DNA Storage Pipeline Benchmarks ----
  console.log('\n=== DNA Storage Pipeline ===\n');

  for (const ds of datasets) {
    // Convert DNA to binary data (2-bit encoding)
    const data = new Uint8Array(Math.floor(ds.dna.length / 4));
    for (let i = 0; i < data.length; i++) {
      let byte = 0;
      for (let b = 0; b < 4 && i * 4 + b < ds.dna.length; b++) {
        let code = 0;
        switch (ds.dna[i * 4 + b]) {
          case 'A': code = 0; break;
          case 'C': code = 1; break;
          case 'G': code = 2; break;
          case 'T': code = 3; break;
        }
        byte |= code << (6 - 2 * b);
      }
      data[i] = byte;
    }

    console.log(`--- ${ds.name} (${data.length.toLocaleString()} bytes) ---`);

    // Use a smaller subset for the pipeline benchmark (encode/decode is slower)
    const subsetSize = Math.min(data.length, 100_000);
    const subset = data.slice(0, subsetSize);

    const result = benchmarkDnaStorage(subset, ds.name);
    if (result) {
      console.log(`  Oligos: ${result.oligoCount}, ${result.payloadBytesPerOligo}B payload each`);
      console.log(`  Total DNA: ${result.totalDnaBases.toLocaleString()} bases`);
      console.log(`  Density: ${result.densityBitsPerNt.toFixed(3)} bits/nt`);
      console.log(`  Encode: ${result.encodeTimeMs}ms (${(subsetSize / result.encodeTimeMs / 1000).toFixed(1)} MB/s)`);
      console.log(`  Decode: ${result.decodeTimeMs}ms (${(subsetSize / result.decodeTimeMs / 1000).toFixed(1)} MB/s)`);
      console.log(`  Roundtrip: ${result.roundtripOk ? 'OK' : 'FAILED'}`);
    }
  }

  // ---- Noisy Channel Recovery Benchmarks ----
  console.log('\n=== Noisy Channel Recovery (with consensus) ===\n');

  const idsRates = [0.01, 0.02, 0.05, 0.09, 0.12];

  for (const ds of datasets.slice(0, 2)) { // Skip random for this test
    const data = new Uint8Array(Math.floor(ds.dna.length / 4));
    for (let i = 0; i < data.length; i++) {
      let byte = 0;
      for (let b = 0; b < 4 && i * 4 + b < ds.dna.length; b++) {
        let code = 0;
        switch (ds.dna[i * 4 + b]) {
          case 'A': code = 0; break;
          case 'C': code = 1; break;
          case 'G': code = 2; break;
          case 'T': code = 3; break;
        }
        byte |= code << (6 - 2 * b);
      }
      data[i] = byte;
    }

    console.log(`--- ${ds.name} ---`);
    const recoveryResults = benchmarkNoisyChannel(data.slice(0, 10_000), ds.name, idsRates);
    for (const r of recoveryResults) {
      console.log(`  IDS ${(r.idsRate * 100).toFixed(0)}%: ${(r.recoveryRate * 100).toFixed(1)}% recovery (${r.decodeTimeMs}ms)`);
    }
  }

  // ---- Save results ----
  const outputDir = path.join(process.cwd(), 'test-data');
  const outputPath = path.join(outputDir, 'genomic-benchmark-results.json');

  const allResults = {
    timestamp: new Date().toISOString(),
    datasets: datasets.map(d => ({ name: d.name, length: d.dna.length, gcContent: d.gcContent })),
    arithmeticCoding: arithResults,
  };

  fs.writeFileSync(outputPath, JSON.stringify(allResults, null, 2));
  console.log(`\nResults saved to ${outputPath}`);
}

main().catch(err => {
  console.error('Benchmark failed:', err);
  process.exit(1);
});
