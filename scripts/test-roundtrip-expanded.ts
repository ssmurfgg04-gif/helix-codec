/**
 * Expanded round-trip tests with larger data, multiple configs, and noisy channels.
 */
import { encodeFile } from '../src/lib/dna/codec';
import { decodeReads } from '../src/lib/dna/decode';
import { simulate, PRESET_CLEAN, PRESET_ILLUMINA } from '../src/lib/dna/simulate';
import { randomBytes } from 'crypto';

let passCount = 0;
let failCount = 0;

async function testRoundtrip(
  label: string,
  data: Uint8Array,
  cfg: any,
  simCfg?: any,
) {
  try {
    const encoded = await encodeFile(data, cfg, { fileName: 'test.bin', contentType: 'application/octet-stream' });
    const sim = simulate(encoded.encoded.oligos, simCfg || { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
    const decoded = await decodeReads(
      sim.reads,
      encoded.encoded.metadata,
      cfg,
      encoded.encoded.forwardPrimer,
      encoded.encoded.reversePrimer
    );

    const dataMatch = decoded.data && Buffer.from(decoded.data).equals(Buffer.from(data));
    const ok = decoded.hashMatches && dataMatch;
    if (ok) {
      passCount++;
      console.log(`  ✓ ${label}`);
    } else {
      failCount++;
      console.log(`  ✗ ${label}: hash=${decoded.hashMatches} data=${dataMatch} erased=${decoded.stats.oligosErased} failedOuter=${decoded.stats.oligosFailedOuterRS}`);
    }
    return ok;
  } catch (e: any) {
    failCount++;
    console.log(`  ✗ ${label}: ERROR ${e.message}`);
    return false;
  }
}

async function main() {
  // Generate test data sizes
  const sizes = [100, 1024, 10 * 1024];  // 100B, 1KB, 10KB

  for (const size of sizes) {
    const data = randomBytes(size);

    // Test all major mapping modes with RS
    for (const mode of ['direct', 'constrained', 'yinyang'] as const) {
      await testRoundtrip(
        `RS+${mode} ${size}B (200nt/20nt)`,
        data,
        {
          oligoLength: 200,
          primerLength: 20,
          innerCode: 'rs',
          mappingMode: mode,
          innerParityBytes: 4,
          outerParityRatio: 0.15,
          constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
          compress: false,
          maxRetries: mode === 'direct' ? 5 : 1,
          channel: 'illumina',
        }
      );
    }

    // Test all major mapping modes with LDPC
    for (const mode of ['direct', 'constrained', 'yinyang'] as const) {
      await testRoundtrip(
        `LDPC+${mode} ${size}B (200nt/20nt)`,
        data,
        {
          oligoLength: 200,
          primerLength: 20,
          innerCode: 'ldpc',
          mappingMode: mode,
          innerParityBytes: 4,
          outerParityRatio: 0.15,
          constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
          compress: false,
          maxRetries: mode === 'direct' ? 5 : 1,
          channel: 'illumina',
        }
      );
    }
  }

  // Test with 300nt/12nt primer (default-like configs)
  for (const size of [100, 1024, 10 * 1024]) {
    const data = randomBytes(size);

    for (const mode of ['constrained', 'yinyang'] as const) {
      await testRoundtrip(
        `LDPC+${mode} ${size}B (300nt/12nt)`,
        data,
        {
          oligoLength: 300,
          primerLength: 12,
          innerCode: 'ldpc',
          mappingMode: mode,
          innerParityBytes: 4,
          outerParityRatio: 0.15,
          constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
          compress: false,
          maxRetries: 1,
          channel: 'illumina',
        }
      );
    }
  }

  // Test with compression enabled
  for (const size of [100, 1024, 10 * 1024]) {
    const data = randomBytes(size);

    await testRoundtrip(
      `LDPC+yinyang+compress ${size}B (300nt/12nt)`,
      data,
      {
        oligoLength: 300,
        primerLength: 12,
        innerCode: 'ldpc',
        mappingMode: 'yinyang',
        innerParityBytes: 4,
        outerParityRatio: 0.15,
        constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
        compress: true,
        maxRetries: 1,
        channel: 'illumina',
      }
    );
  }

  // Test with Illumina noise
  for (const size of [100, 1024]) {
    const data = randomBytes(size);

    await testRoundtrip(
      `LDPC+constrained+noise ${size}B (300nt/12nt)`,
      data,
      {
        oligoLength: 300,
        primerLength: 12,
        innerCode: 'ldpc',
        mappingMode: 'constrained',
        innerParityBytes: 4,
        outerParityRatio: 0.25,
        constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
        compress: false,
        maxRetries: 1,
        channel: 'illumina',
      },
      { ...PRESET_ILLUMINA, coverage: 30, simulator: 'basic' }
    );
  }

  // Test odd oligo lengths (150nt, 12nt primer -> innerNt=126, rounds to 124)
  for (const size of [100, 1024]) {
    const data = randomBytes(size);

    await testRoundtrip(
      `LDPC+constrained ${size}B (150nt/12nt rounded)`,
      data,
      {
        oligoLength: 150,
        primerLength: 12,
        innerCode: 'ldpc',
        mappingMode: 'constrained',
        innerParityBytes: 4,
        outerParityRatio: 0.3,
        constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
        compress: false,
        maxRetries: 1,
        channel: 'illumina',
      }
    );
  }

  console.log(`\n${'='.repeat(60)}`);
  console.log(`Results: ${passCount} passed, ${failCount} failed`);
  console.log(`${'='.repeat(60)}`);
}

main().catch(console.error);
