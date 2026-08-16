/**
 * Targeted debug test for yinyang round-trip.
 */
import { encodeFile } from '../src/lib/dna/codec';
import { decodeReads } from '../src/lib/dna/decode';
import { simulate, PRESET_CLEAN } from '../src/lib/dna/simulate';
import { computeLayout, computeLayoutAuto } from '../src/lib/dna/types';
import { yinyangEncode, yinyangDecode } from '../src/lib/dna/yinyang';
import { dnaToBytes, bytesToDna, unwhitenAddress } from '../src/lib/dna/mapping';

async function main() {
  const data = new Uint8Array(100);
  for (let i = 0; i < 100; i++) data[i] = i;

  const cfg = {
    oligoLength: 300,
    primerLength: 12,
    innerCode: 'ldpc' as const,
    ldpcDecoder: 'auto' as const,
    mappingMode: 'yinyang' as const,
    innerParityBytes: 4,
    outerParityRatio: 0.1,
    constraints: { gcMin: 0.4, gcMax: 0.6, maxHomopolymer: 3 },
    compress: true,
    maxRetries: 1,
    channel: 'illumina' as const,
  };

  // Check layout
  const layout = computeLayout(cfg);
  console.log('Layout:', JSON.stringify(layout, null, 2));
  const totalUsed = layout.addressBytes + layout.payloadBytes + layout.innerParityBytes + layout.crcBytes;
  console.log(`Total used bytes: ${totalUsed}, totalInnerBytes: ${layout.totalInnerBytes}`);
  console.log(`Gap: ${layout.totalInnerBytes - totalUsed} bytes = ${(layout.totalInnerBytes - totalUsed) * 4} nt`);
  console.log(`Expected DNA len: ${layout.totalInnerBytes * 4}`);
  console.log(`Actual used DNA len: ${totalUsed * 4}`);

  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  console.log(`innerK: ${innerK}, innerN: ${innerN}`);
  console.log(`CRC starts at byte: ${innerN}, CRC bytes: ${layout.crcBytes}`);
  console.log(`CRC end at byte: ${innerN + layout.crcBytes}`);
  console.log(`Padding bytes: ${layout.totalInnerBytes - innerN - layout.crcBytes}`);

  // Encode
  const encoded = await encodeFile(data, cfg, { fileName: 'test.bin', contentType: 'application/octet-stream' });
  console.log(`\nEncoded: ${encoded.encoded.oligos.length} oligos`);
  console.log(`Forward primer: ${encoded.encoded.forwardPrimer} (len=${encoded.encoded.forwardPrimer.length})`);
  console.log(`Reverse primer: ${encoded.encoded.reversePrimer} (len=${encoded.encoded.reversePrimer.length})`);

  // Check first oligo
  const oligo0 = encoded.encoded.oligos[0];
  console.log(`\nOligo 0:`);
  console.log(`  Full sequence length: ${oligo0.sequence.length}`);
  console.log(`  Full sequence: ${oligo0.sequence}`);

  // Strip primers manually
  const fwd = encoded.encoded.forwardPrimer;
  const rev = encoded.encoded.reversePrimer;
  const innerDna = oligo0.sequence.slice(fwd.length, oligo0.sequence.length - rev.length);
  console.log(`  Inner DNA length: ${innerDna.length}`);
  console.log(`  Inner DNA: ${innerDna}`);

  // Try yinyang decode
  const expectedDnaLen = layout.totalInnerBytes * 4;
  console.log(`  Expected inner DNA length: ${expectedDnaLen}`);
  console.log(`  Match: ${innerDna.length === expectedDnaLen}`);

  try {
    const decoded = yinyangDecode(innerDna, layout.totalInnerBytes);
    console.log(`  Decoded bytes length: ${decoded.length}`);
    console.log(`  First 4 bytes (address): ${Array.from(decoded.slice(0, 4)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    
    // Unwhiten address
    const unwhitened = unwhitenAddress(decoded.slice(0, 4));
    const idx = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
    console.log(`  Decoded index: ${idx}, seed: ${unwhitened[3]}`);
    
    // Check rsCodeword
    const rsCodeword = decoded.slice(0, innerN);
    console.log(`  RS codeword length: ${rsCodeword.length}`);
    console.log(`  CRC bytes: ${Array.from(decoded.slice(innerN, innerN + 2)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    if (decoded.length > innerN + 2) {
      console.log(`  Extra bytes after CRC: ${decoded.length - innerN - 2}`);
      console.log(`  Extra bytes: ${Array.from(decoded.slice(innerN + 2)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    }

    // Try LDPC decode
    const { getCachedLDPCInner } = await import('../src/lib/dna/ldpc-codec');
    const ldpc = getCachedLDPCInner(innerN, innerK);
    const ldpcResult = ldpc.decode(rsCodeword);
    console.log(`  LDPC decode: corrected=${ldpcResult.corrected}`);
    console.log(`  LDPC decoded data length: ${ldpcResult.data.length}`);
    
    // Verify CRC
    const { crc16Bytes } = await import('../src/lib/dna/crc16');
    const reEncoded = ldpc.encode(ldpcResult.data);
    const expectedCrc = crc16Bytes(reEncoded);
    const actualCrc = decoded.slice(innerN, innerN + 2);
    console.log(`  Re-encoded CRC: ${Array.from(expectedCrc).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    console.log(`  Actual CRC: ${Array.from(actualCrc).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
    console.log(`  CRC match: ${expectedCrc[0] === actualCrc[0] && expectedCrc[1] === actualCrc[1]}`);
  } catch (e: any) {
    console.log(`  Decode error: ${e.message}`);
  }

  // Now try the full decode pipeline with HELIX_DEBUG=1
  const sim = simulate(encoded.encoded.oligos, { ...PRESET_CLEAN, coverage: 10, simulator: 'basic' });
  console.log(`\nSimulated ${sim.reads.length} reads`);
  
  // Check first read
  console.log(`First read: index=${sim.reads[0].oligoIndex}, seq_len=${sim.reads[0].sequence.length}, seq=${sim.reads[0].sequence}`);

  const decoded = await decodeReads(
    sim.reads,
    encoded.encoded.metadata,
    cfg,
    encoded.encoded.forwardPrimer,
    encoded.encoded.reversePrimer
  );

  console.log(`\nFull decode result:`);
  console.log(`  hash match: ${decoded.hashMatches}`);
  console.log(`  data match: ${decoded.data && Buffer.from(decoded.data).equals(Buffer.from(data))}`);
  console.log(`  stats: ${JSON.stringify(decoded.stats)}`);
  
  // Show per-oligo details
  for (const p of decoded.perOligo) {
    console.log(`  Oligo ${p.index}: reads=${p.readCount}, crc=${p.crcPassed}, innerRS=${JSON.stringify(p.innerRS)}, strategy=${p.strategy}`);
  }
}

main().catch(console.error);
