import { encodeFile } from '../src/lib/dna/codec';
import { decodeReads } from '../src/lib/dna/decode';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';
import type { CodecConfig } from '../src/lib/dna/types';
import { readFile } from 'fs/promises';

async function testRoundtrip(name: string, data: Uint8Array, cfg: CodecConfig) {
  try {
    const enc = await encodeFile(data, cfg, { fileName: name, contentType: 'application/octet-stream' });
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const dataMatch = Buffer.from(dec.data).equals(Buffer.from(data));
    console.log(`${name}: ${dataMatch ? 'PASS' : 'FAIL'} (${data.length}B, ${enc.encoded.oligos.length} oligos, enc=${enc.stats.encodeTimeMs}ms, erased=${dec.stats.oligosErased}/${enc.encoded.oligos.length})`);
    if (!dataMatch) {
      console.log(`  Hash: ${dec.hashMatches}, recovered=${dec.stats.oligosRecovered}, failI=${dec.stats.oligosFailedInnerRS}, failO=${dec.stats.oligosFailedOuterRS}`);
    }
    return dataMatch;
  } catch (e: any) {
    console.log(`${name}: ERROR - ${e.message?.slice(0, 120)}`);
    return false;
  }
}

async function main() {
  const configs: Record<string, CodecConfig> = {
    'rs-direct': {
      oligoLength: 200, primerLength: 20, innerCode: 'rs', mappingMode: 'direct',
      innerParityBytes: 4, outerParityRatio: 0.15,
      constraints: { gcMin: 0.0, gcMax: 1.0, maxHomopolymer: 4 },
      compress: false, maxRetries: 10, channel: 'illumina', lowCoverageTrigger: 5,
    },
    'ldpc-constrained': {
      oligoLength: 300, primerLength: 12, innerCode: 'ldpc', ldpcDecoder: 'auto', mappingMode: 'constrained',
      innerParityBytes: 8, outerParityRatio: 0.15,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: true, maxRetries: 5, channel: 'illumina', lowCoverageTrigger: 5,
    },
  };

  // Generate test data of various sizes
  for (const size of [100, 1000, 10000, 100000]) {
    const data = new Uint8Array(size);
    for (let i = 0; i < size; i++) data[i] = (i * 7 + 13) & 0xFF;
    
    for (const [cfgName, cfg] of Object.entries(configs)) {
      await testRoundtrip(`${size}B-${cfgName}`, data, cfg);
    }
  }

  // Real data test
  try {
    const textData = await readFile('test-data/pride_and_prejudice.txt');
    const data = new Uint8Array(textData.buffer, textData.byteOffset, textData.byteLength);
    for (const [cfgName, cfg] of Object.entries(configs)) {
      await testRoundtrip(`pnp-${cfgName}`, data, cfg);
    }
  } catch (e: any) {
    console.log(`Gutenberg: ${e.message}`);
  }
}

main().catch(console.error);
