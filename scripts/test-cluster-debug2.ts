import { encodeFile } from '../src/lib/dna/codec';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';
import { computeLayout } from '../src/lib/dna/types';
import { dnaToBytes, unwhitenAddress } from '../src/lib/dna/mapping';
import { randomBytes } from 'crypto';

async function main() {
  const data = randomBytes(10240);
  const cfg = {
    oligoLength: 200, primerLength: 20, innerCode: 'rs' as const, mappingMode: 'constrained' as const,
    innerParityBytes: 4, outerParityRatio: 0.15,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: false, maxRetries: 1, channel: 'illumina' as const,
  };
  const encoded = await encodeFile(data, cfg, { fileName: 'test.bin', contentType: 'application/octet-stream' });
  const layout = computeLayout(cfg);
  const addressNt = layout.addressBytes * 4;
  const fwd = encoded.encoded.forwardPrimer;
  const rev = encoded.encoded.reversePrimer;

  let nonzeroSeeds = 0;
  for (const o of encoded.encoded.oligos) {
    if (o.seed !== 0) nonzeroSeeds++;
  }
  console.log('Oligos:', encoded.encoded.oligos.length, 'nonzero seeds:', nonzeroSeeds);

  // Check oligos 36 and 37
  for (const idx of [36, 37]) {
    const o = encoded.encoded.oligos[idx];
    const inner = o.sequence.slice(fwd.length, o.sequence.length - rev.length);
    const addr = dnaToBytes(inner.slice(0, addressNt));
    const uw = unwhitenAddress(addr);
    const di = (uw[0] << 16) | (uw[1] << 8) | uw[2];
    console.log('Oligo', idx, 'seed=', o.seed, 'decodedIdx=', di, 'innerLen=', inner.length, 'expectedLen=', layout.totalInnerBytes * 4);
  }

  // Check read sequences
  const sim = simulate(encoded.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
  for (const idx of [36, 37]) {
    const reads = sim.reads.filter(r => r.oligoIndex === idx);
    console.log('Reads for oligo', idx, ':', reads.length);
    if (reads.length > 0) {
      const r0 = reads[0];
      const inner = r0.sequence.slice(fwd.length, r0.sequence.length - rev.length);
      const addr = dnaToBytes(inner.slice(0, addressNt));
      const uw = unwhitenAddress(addr);
      const di = (uw[0] << 16) | (uw[1] << 8) | uw[2];
      console.log('  Read0: oligoIndex=', r0.oligoIndex, 'decodedIdx=', di);
    }
  }
}
main().catch(console.error);
