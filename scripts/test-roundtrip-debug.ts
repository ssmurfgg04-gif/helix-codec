/**
 * Minimal round-trip test to debug encode→decode hash mismatches.
 */
import { encodeFile } from '../src/lib/dna/codec';
import { decodeReads } from '../src/lib/dna/decode';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';

async function testRoundtrip(
  label: string,
  data: Uint8Array,
  cfg: any
) {
  console.log(`\n=== ${label} ===`);
  try {
    const encoded = await encodeFile(data, cfg, { fileName: 'test.bin', contentType: 'application/octet-stream' });
    console.log(`  Encoded: ${encoded.encoded.oligos.length} oligos, ${encoded.encoded.metadata.outerRS.k} data + ${encoded.encoded.metadata.outerRS.n - encoded.encoded.metadata.outerRS.k} parity`);
    console.log(`  payloadBytesPerOligo: ${encoded.encoded.metadata.payloadBytesPerOligo}`);
    console.log(`  mappingMode: ${encoded.encoded.metadata.mappingMode}, innerCode: ${encoded.encoded.metadata.innerCode}`);
    console.log(`  fileHash: ${encoded.encoded.metadata.fileHash}`);

    // Clean simulation (no errors)
    const sim = simulate(encoded.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
    console.log(`  Simulated: ${sim.reads.length} reads`);

    const decoded = await decodeReads(
      sim.reads,
      encoded.encoded.metadata,
      cfg,
      encoded.encoded.forwardPrimer,
      encoded.encoded.reversePrimer
    );

    const dataMatch = decoded.data && Buffer.from(decoded.data).equals(Buffer.from(data));
    console.log(`  hash match: ${decoded.hashMatches}`);
    console.log(`  data match: ${dataMatch}`);
    console.log(`  stats: recovered=${decoded.stats.oligosRecovered}, erased=${decoded.stats.oligosErased}, failedInner=${decoded.stats.oligosFailedInnerRS}, failedOuter=${decoded.stats.oligosFailedOuterRS}`);

    // Count per-oligo strategies
    const strategyCounts: Record<string, number> = {};
    for (const p of decoded.perOligo) {
      const s = p.strategy || 'unknown';
      strategyCounts[s] = (strategyCounts[s] || 0) + 1;
    }
    console.log(`  strategies: ${JSON.stringify(strategyCounts)}`);

    if (!decoded.hashMatches) {
      // Show first few erased/failed oligos
      const failed = decoded.perOligo.filter(p => p.strategy === 'erasure');
      console.log(`  FAILED oligos (first 5): ${failed.slice(0, 5).map(f => f.index).join(', ')}`);
    }

    return decoded.hashMatches;
  } catch (e: any) {
    console.log(`  ERROR: ${e.message}`);
    return false;
  }
}

async function main() {
  // Create test data
  const data = new Uint8Array(100);
  for (let i = 0; i < 100; i++) data[i] = i;

  // Test 1: RS + direct mapping (simplest config)
  await testRoundtrip('RS + direct (200nt/20nt primer)', data, {
    oligoLength: 200,
    primerLength: 20,
    innerCode: 'rs',
    mappingMode: 'direct',
    innerParityBytes: 4,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 5,
    channel: 'illumina',
  });

  // Test 2: RS + constrained mapping
  await testRoundtrip('RS + constrained (200nt/20nt primer)', data, {
    oligoLength: 200,
    primerLength: 20,
    innerCode: 'rs',
    mappingMode: 'constrained',
    innerParityBytes: 4,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 5,
    channel: 'illumina',
  });

  // Test 3: LDPC + constrained (default-like config)
  await testRoundtrip('LDPC + constrained (200nt/20nt primer)', data, {
    oligoLength: 200,
    primerLength: 20,
    innerCode: 'ldpc',
    mappingMode: 'constrained',
    innerParityBytes: 4,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 5,
    channel: 'illumina',
  });

  // Test 4: LDPC + constrained with 12nt primer (300nt oligo)
  await testRoundtrip('LDPC + constrained (300nt/12nt primer)', data, {
    oligoLength: 300,
    primerLength: 12,
    innerCode: 'ldpc',
    mappingMode: 'constrained',
    innerParityBytes: 4,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 5,
    channel: 'illumina',
  });

  // Test 5: RS + direct with maxRetries=1 (mimics DEFAULT_CONFIG issue)
  await testRoundtrip('RS + direct maxRetries=1 (200nt/20nt primer)', data, {
    oligoLength: 200,
    primerLength: 20,
    innerCode: 'rs',
    mappingMode: 'direct',
    innerParityBytes: 4,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 1,
    channel: 'illumina',
  });

  // Test 6: Default config (300nt/12nt/LDPC/yinyang/maxRetries=1)
  await testRoundtrip('Default config (300nt/12nt/LDPC/yinyang)', data, {
    oligoLength: 300,
    primerLength: 12,
    innerCode: 'ldpc',
    ldpcDecoder: 'auto',
    mappingMode: 'yinyang',
    innerParityBytes: 4,
    outerParityRatio: 0.1,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 1,
    channel: 'illumina',
  });

  // Test 7: Nanopore short oligo (150nt/12nt) - should trigger divisibility error
  await testRoundtrip('Nanopore short (150nt/12nt)', data, {
    oligoLength: 150,
    primerLength: 12,
    innerCode: 'ldpc',
    mappingMode: 'constrained',
    innerParityBytes: 4,
    outerParityRatio: 0.3,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 1,
    channel: 'nanopore',
  });

  // Test 8: Larger data (1KB) with RS + direct
  const data1k = new Uint8Array(1024);
  for (let i = 0; i < 1024; i++) data1k[i] = i & 0xff;
  await testRoundtrip('RS + direct 1KB (200nt/20nt primer)', data1k, {
    oligoLength: 200,
    primerLength: 20,
    innerCode: 'rs',
    mappingMode: 'direct',
    innerParityBytes: 4,
    outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false,
    maxRetries: 5,
    channel: 'illumina',
  });
}

main().catch(console.error);
