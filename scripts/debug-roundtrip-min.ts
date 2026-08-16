import { encodeFile } from '../src/lib/dna/codec';
import { decodeReads } from '../src/lib/dna/decode';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';
import type { CodecConfig } from '../src/lib/dna/types';

async function main() {
  const data = new Uint8Array(100);
  for (let i = 0; i < 100; i++) data[i] = i;

  // Test 1: RS + direct mapping (simplest)
  const cfg1: CodecConfig = {
    oligoLength: 200,
    primerLength: 20,
    innerCode: 'rs',
    mappingMode: 'direct',
    innerParityBytes: 4,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.0, gcMax: 1.0, maxHomopolymer: 4 },
    compress: false,
    maxRetries: 10,
    channel: 'illumina',
    lowCoverageTrigger: 5,
  };

  console.log('=== Test 1: RS + direct mapping ===');
  try {
    const enc = await encodeFile(data, cfg1, { fileName: 'test.bin', contentType: 'application/octet-stream' });
    console.log(`Encoded: ${enc.encoded.oligos.length} oligos, ${enc.stats.encodeTimeMs}ms`);
    
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
    console.log(`Simulated: ${sim.reads.length} reads`);
    
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg1, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    console.log(`Hash match: ${dec.hashMatches}`);
    const dataMatch = Buffer.from(dec.data).equals(Buffer.from(data));
    console.log(`Data match: ${dataMatch}`);
    if (!dataMatch) {
      console.log(`Decoded length: ${dec.data.length}, expected: ${data.length}`);
      for (let i = 0; i < Math.min(data.length, dec.data.length); i++) {
        if (data[i] !== dec.data[i]) {
          console.log(`First diff at byte ${i}: expected ${data[i]}, got ${dec.data[i]}`);
          break;
        }
      }
    }
  } catch (e: any) {
    console.log(`ERROR: ${e.message}`);
    console.log(e.stack?.slice(0, 500));
  }

  // Test 2: LDPC + constrained mapping
  const cfg2: CodecConfig = {
    oligoLength: 300,
    primerLength: 12,
    innerCode: 'ldpc',
    ldpcDecoder: 'auto',
    mappingMode: 'constrained',
    innerParityBytes: 8,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 5,
    channel: 'illumina',
    lowCoverageTrigger: 5,
  };

  console.log('\n=== Test 2: LDPC + constrained mapping ===');
  try {
    const enc = await encodeFile(data, cfg2, { fileName: 'test.bin', contentType: 'application/octet-stream' });
    console.log(`Encoded: ${enc.encoded.oligos.length} oligos, ${enc.stats.encodeTimeMs}ms`);
    
    const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
    const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg2, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    console.log(`Hash match: ${dec.hashMatches}`);
    const dataMatch = Buffer.from(dec.data).equals(Buffer.from(data));
    console.log(`Data match: ${dataMatch}`);
    if (!dataMatch) {
      console.log(`Decoded length: ${dec.data.length}, expected: ${data.length}`);
      console.log(`Stats: recovered=${dec.stats.oligosRecovered}, erased=${dec.stats.oligosErased}, failInner=${dec.stats.oligosFailedInnerRS}, failOuter=${dec.stats.oligosFailedOuterRS}`);
    }
  } catch (e: any) {
    console.log(`ERROR: ${e.message}`);
    console.log(e.stack?.slice(0, 500));
  }
}

main().catch(console.error);
