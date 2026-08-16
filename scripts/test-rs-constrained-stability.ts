import { encodeFile } from '../src/lib/dna/codec';
import { decodeReads } from '../src/lib/dna/decode';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';
import { randomBytes } from 'crypto';

async function main() {
  // Test RS+constrained 1KB multiple times with 0.15 parity
  for (let run = 0; run < 5; run++) {
    const data = randomBytes(1024);
    const cfg = {
      oligoLength: 200, primerLength: 20, innerCode: 'rs' as const, mappingMode: 'constrained' as const,
      innerParityBytes: 4, outerParityRatio: 0.15,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: false, maxRetries: 1, channel: 'illumina' as const,
    };
    const encoded = await encodeFile(data, cfg, { fileName: 'test.bin', contentType: 'application/octet-stream' });
    const sim = simulate(encoded.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
    const decoded = await decodeReads(sim.reads, encoded.encoded.metadata, cfg, encoded.encoded.forwardPrimer, encoded.encoded.reversePrimer);
    const k = encoded.encoded.metadata.outerRS.k;
    const n = encoded.encoded.metadata.outerRS.n;
    console.log(`Run ${run}: ${decoded.hashMatches ? 'PASS' : 'FAIL'} erased=${decoded.stats.oligosErased} failedInner=${decoded.stats.oligosFailedInnerRS} failedOuter=${decoded.stats.oligosFailedOuterRS} k=${k} n=${n}`);
  }

  // Test with higher parity
  console.log('\n--- With 25% outer parity ---');
  for (let run = 0; run < 5; run++) {
    const data = randomBytes(1024);
    const cfg = {
      oligoLength: 200, primerLength: 20, innerCode: 'rs' as const, mappingMode: 'constrained' as const,
      innerParityBytes: 4, outerParityRatio: 0.25,
      constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
      compress: false, maxRetries: 1, channel: 'illumina' as const,
    };
    const encoded = await encodeFile(data, cfg, { fileName: 'test.bin', contentType: 'application/octet-stream' });
    const sim = simulate(encoded.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
    const decoded = await decodeReads(sim.reads, encoded.encoded.metadata, cfg, encoded.encoded.forwardPrimer, encoded.encoded.reversePrimer);
    console.log(`Run ${run}: ${decoded.hashMatches ? 'PASS' : 'FAIL'} erased=${decoded.stats.oligosErased} failedInner=${decoded.stats.oligosFailedInnerRS} failedOuter=${decoded.stats.oligosFailedOuterRS}`);
  }
}
main().catch(console.error);
