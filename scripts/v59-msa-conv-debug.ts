/**
 * v59 MSA + Conv debug — trace why MSA + conv Viterbi fails.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { simulate, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { ULTIMATE_NANOPORE_V52_CONFIG } from "../src/lib/dna/presets";
import { computeLayoutAuto } from "../src/lib/dna/types";
import { ConvolutionalInnerCode } from "../src/lib/dna/convolutional";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { dnaToBytes, unwhitenAddress, whitenAddress, bytesToDna } from "../src/lib/dna/mapping";
import { buildReferenceKmerIndex, matchReadToReference, extractKmers } from "../src/lib/dna/kmer";
import { progressiveMSA, msaConsensus, DEFAULT_MSA_CONFIG } from "../src/lib/dna/progressive-msa";
import { crc16Bytes } from "../src/lib/dna/crc16";

const TAG = "[v59-msa-conv]";

async function main() {
  const payload = new Uint8Array(16 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;

  const cfg = ULTIMATE_NANOPORE_V52_CONFIG;
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`${TAG} Encoded ${enc.encoded.oligos.length} oligos`);

  const sim = simulate(enc.encoded.oligos, { ...PRESET_NANOPORE, coverage: 10, seed: 42 });

  const layout = computeLayoutAuto(cfg);
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  console.log(`${TAG} Layout: innerN=${innerN}, innerK=${innerK}, convEncoded=${layout.convEncodedBytes}, totalInner=${layout.totalInnerBytes}`);

  const fwd = enc.encoded.forwardPrimer;
  const rev = enc.encoded.reversePrimer;
  const expectedDnaLen = layout.totalInnerBytes * 4;

  // Build reference addresses
  const refs: string[] = [];
  for (let idx = 0; idx < enc.encoded.oligos.length; idx++) {
    const raw = new Uint8Array(4);
    raw[0] = (idx >> 16) & 0xff;
    raw[1] = (idx >> 8) & 0xff;
    raw[2] = idx & 0xff;
    raw[3] = 0;
    refs.push(bytesToDna(whitenAddress(raw)));
  }
  const kmerIdx = buildReferenceKmerIndex(refs, 4);

  // Cluster reads using k-mer matching
  const clusters: Map<number, any[]> = new Map();
  for (const read of sim.reads) {
    const seq = read.sequence;
    // K-mer find forward primer
    const k = 5;
    const fwdKmers = new Set(extractKmers(fwd, k));
    let fwdBestPos = 0;
    let fwdBestOverlap = 0;
    for (let p = 0; p <= 20; p++) {
      if (p + fwd.length > seq.length) break;
      const cand = new Set(extractKmers(seq.slice(p, p + fwd.length), k));
      let overlap = 0;
      for (const km of fwdKmers) if (cand.has(km)) overlap++;
      if (overlap > fwdBestOverlap) { fwdBestOverlap = overlap; fwdBestPos = p; if (overlap >= fwdKmers.size - 1) break; }
    }
    if (fwdBestOverlap < 3) continue;

    const innerStart = fwdBestPos + fwd.length;
    const inner = seq.slice(innerStart, innerStart + expectedDnaLen);
    if (inner.length < 16) continue;

    const addressDna = inner.slice(0, 16);
    let assignedIdx = -1;
    try {
      const addrBytes = dnaToBytes(addressDna);
      const unwh = unwhitenAddress(addrBytes);
      const idx = (unwh[0] << 16) | (unwh[1] << 8) | unwh[2];
      if (idx >= 0 && idx < refs.length) {
        const refDna = bytesToDna(whitenAddress(unwh));
        let dist = 0;
        for (let j = 0; j < 16; j++) if (addressDna[j] !== refDna[j]) dist++;
        if (dist <= 2) assignedIdx = idx;
      }
    } catch {}
    if (assignedIdx === -1) {
      const result = matchReadToReference(addressDna, kmerIdx, 4, 2);
      if (result.bestIdx >= 0) assignedIdx = result.bestIdx;
    }
    if (assignedIdx === -1) continue;

    if (!clusters.has(assignedIdx)) clusters.set(assignedIdx, []);
    clusters.get(assignedIdx)!.push({ sequence: inner });
  }

  // Pick oligo 2 (which has 13 reads but fails)
  const oligoIdx = 2;
  const clusterReads = clusters.get(oligoIdx) ?? [];
  console.log(`\n${TAG} Oligo ${oligoIdx}: ${clusterReads.length} reads`);

  if (clusterReads.length < 3) {
    console.log(`${TAG} Not enough reads`);
    return;
  }

  // Run MSA
  const trimmedSeqs = clusterReads.slice(0, 15).map(r => {
    let s = r.sequence;
    const maxLen = Math.floor(expectedDnaLen * 1.2);
    const minLen = Math.floor(expectedDnaLen * 0.8);
    if (s.length > maxLen) s = s.slice(0, maxLen);
    else if (s.length < minLen) s = s + "A".repeat(minLen - s.length);
    return s;
  });

  console.log(`${TAG} Read lengths: ${trimmedSeqs.map(s => s.length).join(", ")}`);

  const aligned = progressiveMSA(trimmedSeqs, DEFAULT_MSA_CONFIG);
  const consensusResult = msaConsensus(aligned);
  console.log(`${TAG} MSA consensus length: ${consensusResult.sequence.length} (expected ${expectedDnaLen})`);

  let consensus = consensusResult.sequence;
  if (consensus.length > expectedDnaLen) consensus = consensus.slice(0, expectedDnaLen);
  else if (consensus.length < expectedDnaLen) consensus = consensus + "A".repeat(expectedDnaLen - consensus.length);

  // Convert to bytes
  const innerBlock = dnaToBytes(consensus);
  console.log(`${TAG} Inner block bytes: ${innerBlock.length}`);

  // Extract address
  const addrBytes = innerBlock.slice(0, 4);
  const unwh = unwhitenAddress(addrBytes);
  const decIdx = (unwh[0] << 16) | (unwh[1] << 8) | unwh[2];
  console.log(`${TAG} Decoded address index: ${decIdx} (expected ${oligoIdx})`);

  // Conv decode
  const convInner = new ConvolutionalInnerCode(innerN);
  const convBytes = innerBlock.slice(4, 4 + layout.convEncodedBytes);
  console.log(`${TAG} Conv bytes: ${convBytes.length} (expected ${convInner.outputBytes})`);

  try {
    const rsCodewordDecoded = convInner.decode(convBytes);
    console.log(`${TAG} Conv decoded: ${rsCodewordDecoded.length} bytes (expected ${innerN})`);

    // LDPC decode
    const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
    const r = ldpc.decode(rsCodewordDecoded);
    console.log(`${TAG} LDPC decoded: corrected=${r.corrected}`);

    // Verify address
    const reEncoded = ldpc.encode(r.data);
    const addr2 = unwhitenAddress(reEncoded.slice(0, 4));
    const decIdx2 = (addr2[0] << 16) | (addr2[1] << 8) | addr2[2];
    console.log(`${TAG} LDPC address index: ${decIdx2} (expected ${oligoIdx})`);

    // CRC check
    const crc = crc16Bytes(reEncoded);
    const crcFromRead = innerBlock.slice(innerBlock.length - 2);
    console.log(`${TAG} CRC: expected=[${crc[0]},${crc[1]}] got=[${crcFromRead[0]},${crcFromRead[1]}] match=${crc[0] === crcFromRead[0] && crc[1] === crcFromRead[1]}`);
  } catch (e) {
    console.log(`${TAG} Conv/LDPC failed: ${(e as Error).message}`);
  }
}
main().catch(e => { console.error(e); process.exit(1); });
