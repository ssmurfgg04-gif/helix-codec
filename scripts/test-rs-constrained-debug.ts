/**
 * Debug RS+constrained 10KB failure.
 */
import { encodeFile } from '../src/lib/dna/codec';
import { decodeReads } from '../src/lib/dna/decode';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';
import { randomBytes } from 'crypto';

async function main() {
  const data = randomBytes(10240);

  const cfg = {
    oligoLength: 200,
    primerLength: 20,
    innerCode: 'rs' as const,
    mappingMode: 'constrained' as const,
    innerParityBytes: 4,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 1,
    channel: 'illumina' as const,
  };

  const encoded = await encodeFile(data, cfg, { fileName: 'test.bin', contentType: 'application/octet-stream' });
  console.log(`Encoded: ${encoded.encoded.oligos.length} oligos, ${encoded.encoded.metadata.outerRS.k} data + ${encoded.encoded.metadata.outerRS.n - encoded.encoded.metadata.outerRS.k} parity`);
  console.log(`payloadBytesPerOligo: ${encoded.encoded.metadata.payloadBytesPerOligo}`);

  const sim = simulate(encoded.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
  
  const decoded = await decodeReads(
    sim.reads,
    encoded.encoded.metadata,
    cfg,
    encoded.encoded.forwardPrimer,
    encoded.encoded.reversePrimer
  );

  console.log(`hash match: ${decoded.hashMatches}`);
  console.log(`stats: recovered=${decoded.stats.oligosRecovered}, erased=${decoded.stats.oligosErased}, failedInner=${decoded.stats.oligosFailedInnerRS}, failedOuter=${decoded.stats.oligosFailedOuterRS}`);

  // Count strategies
  const strategyCounts: Record<string, number> = {};
  for (const p of decoded.perOligo) {
    const s = p.strategy || 'unknown';
    strategyCounts[s] = (strategyCounts[s] || 0) + 1;
  }
  console.log(`strategies: ${JSON.stringify(strategyCounts)}`);

  // Analyze failed oligos
  const failed = decoded.perOligo.filter(p => p.strategy === 'erasure');
  console.log(`\nFailed oligos: ${failed.length}`);
  
  // Check how many reads per failed oligo
  const readCounts = failed.map(f => f.readCount);
  console.log(`Read counts for failed oligos (first 20): ${readCounts.slice(0, 20).join(', ')}`);
  console.log(`Oligos with 0 reads: ${failed.filter(f => f.readCount === 0).length}`);
  console.log(`Oligos with >0 reads: ${failed.filter(f => f.readCount > 0).length}`);

  // Also test with higher parity
  console.log(`\n--- Testing with 25% outer parity ---`);
  const cfg2 = { ...cfg, outerParityRatio: 0.25 };
  const encoded2 = await encodeFile(data, cfg2, { fileName: 'test.bin', contentType: 'application/octet-stream' });
  console.log(`Encoded: ${encoded2.encoded.oligos.length} oligos, ${encoded2.encoded.metadata.outerRS.k} data + ${encoded2.encoded.metadata.outerRS.n - encoded2.encoded.metadata.outerRS.k} parity`);
  
  const sim2 = simulate(encoded2.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
  const decoded2 = await decodeReads(
    sim2.reads,
    encoded2.encoded.metadata,
    cfg2,
    encoded2.encoded.forwardPrimer,
    encoded2.encoded.reversePrimer
  );
  console.log(`hash match: ${decoded2.hashMatches}`);
  console.log(`stats: recovered=${decoded2.stats.oligosRecovered}, erased=${decoded2.stats.oligosErased}, failedInner=${decoded2.stats.oligosFailedInnerRS}, failedOuter=${decoded2.stats.oligosFailedOuterRS}`);
}

main().catch(console.error);
