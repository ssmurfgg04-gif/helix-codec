import { encodeFile } from '../src/lib/dna/codec';
import { decodeReads } from '../src/lib/dna/decode';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';
import type { CodecConfig } from '../src/lib/dna/types';
import { readFile } from 'fs/promises';

async function rt(name: string, data: Uint8Array, cfg: CodecConfig) {
  try {
    const enc = await encodeFile(data, cfg, { fileName: name, contentType: 'application/octet-stream' });
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const ok = Buffer.from(dec.data).equals(Buffer.from(data));
    console.log(`${ok ? '✓' : '✗'} ${name}: ${data.length}B, ${enc.encoded.oligos.length} oligos, enc=${enc.stats.encodeTimeMs}ms, dec=ms, density=${enc.stats.netDensityBitsPerNt?.toFixed(3)}b/nt, erased=${dec.stats.oligosErased}`);
    return ok;
  } catch (e: any) {
    console.log(`✗ ${name}: ${e.message?.slice(0,80)}`);
    return false;
  }
}

async function main() {
  const data1k = await readFile('test-data/random_1kb.bin');
  const data10k = await readFile('test-data/random_10kb.bin');
  const data100k = await readFile('test-data/random_100kb.bin');
  
  const cfgs: Record<string, CodecConfig> = {
    'ldpc-300': { oligoLength: 300, primerLength: 12, innerCode: 'ldpc', ldpcDecoder: 'auto', mappingMode: 'constrained', innerParityBytes: 8, outerParityRatio: 0.15, constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 }, compress: true, maxRetries: 5, channel: 'illumina', lowCoverageTrigger: 5 },
    'ldpc-700': { oligoLength: 700, primerLength: 12, innerCode: 'ldpc', ldpcDecoder: 'auto', mappingMode: 'constrained', innerParityBytes: 8, outerParityRatio: 0.15, constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 }, compress: true, maxRetries: 5, channel: 'illumina', lowCoverageTrigger: 5 },
    'rs-200': { oligoLength: 200, primerLength: 20, innerCode: 'rs', mappingMode: 'direct', innerParityBytes: 4, outerParityRatio: 0.15, constraints: { gcMin: 0.0, gcMax: 1.0, maxHomopolymer: 4 }, compress: false, maxRetries: 10, channel: 'illumina', lowCoverageTrigger: 5 },
  };
  
  for (const [cn, cfg] of Object.entries(cfgs)) {
    await rt(`1KB-${cn}`, new Uint8Array(data1k.buffer, data1k.byteOffset, data1k.byteLength), cfg);
    await rt(`10KB-${cn}`, new Uint8Array(data10k.buffer, data10k.byteOffset, data10k.byteLength), cfg);
    await rt(`100KB-${cn}`, new Uint8Array(data100k.buffer, data100k.byteOffset, data100k.byteLength), cfg);
  }
}

main().catch(console.error);
