/**
 * Debug clustering for RS+constrained 10KB.
 */
import { encodeFile } from '../src/lib/dna/codec';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';
import { computeLayout } from '../src/lib/dna/types';
import { dnaToBytes, unwhitenAddress } from '../src/lib/dna/mapping';
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
  console.log(`Encoded: ${encoded.encoded.oligos.length} oligos`);

  const layout = computeLayout(cfg);
  const addressNt = layout.addressBytes * 4;

  const sim = simulate(encoded.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
  console.log(`Total reads: ${sim.reads.length}`);

  // Manually check clustering for first few reads
  const fwd = encoded.encoded.forwardPrimer;
  const rev = encoded.encoded.reversePrimer;

  // Check a few oligos
  for (let i = 0; i < Math.min(5, encoded.encoded.oligos.length); i++) {
    const oligo = encoded.encoded.oligos[i];
    const fullSeq = oligo.sequence;
    const innerDna = fullSeq.slice(fwd.length, fullSeq.length - rev.length);
    
    // Extract address using direct mapping (for constrained, address is direct)
    const addressDna = innerDna.slice(0, addressNt);
    const addressBytes = dnaToBytes(addressDna);
    const unwhitened = unwhitenAddress(addressBytes);
    const index = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
    
    console.log(`Oligo ${i}: address=${Array.from(addressBytes).map(b => b.toString(16).padStart(2,"0")).join("")} -> index=${index}, seed=${unwhitened[3]}, expected=${i}`);
    if (index !== i) {
      console.log(`  MISMATCH! Address decode gives ${index} but expected ${i}`);
    }
  }

  // Check how reads are distributed
  const indexCounts = new Map<number, number>();
  for (const read of sim.reads) {
    const inner = read.sequence.slice(fwd.length, read.sequence.length - rev.length);
    if (inner.length < addressNt) continue;
    const addressDna = inner.slice(0, addressNt);
    try {
      const addressBytes = dnaToBytes(addressDna);
      const unwhitened = unwhitenAddress(addressBytes);
      const index = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
      indexCounts.set(index, (indexCounts.get(index) || 0) + 1);
    } catch {
      console.log(`  Failed to decode address from read`);
    }
  }

  // Find indices with 0 reads or too many reads
  const zeroReadIndices: number[] = [];
  const tooManyReadIndices: number[] = [];
  for (let i = 0; i < encoded.encoded.oligos.length; i++) {
    const count = indexCounts.get(i) || 0;
    if (count === 0) zeroReadIndices.push(i);
    if (count > 15) tooManyReadIndices.push(i);
  }
  console.log(`\nIndices with 0 reads: ${zeroReadIndices.length} (first 10: ${zeroReadIndices.slice(0, 10).join(', ')})`);
  console.log(`Indices with >15 reads: ${tooManyReadIndices.length}`);
  if (tooManyReadIndices.length > 0) {
    console.log(`  Details: ${tooManyReadIndices.slice(0, 10).map(i => `${i}:${indexCounts.get(i)}`).join(', ')}`);
  }

  // Check total reads assigned
  let totalAssigned = 0;
  for (const [_, count] of indexCounts) totalAssigned += count;
  console.log(`Total reads assigned: ${totalAssigned} (expected ${sim.reads.length})`);
}

main().catch(console.error);
