/**
 * v59 Misclustering check — are k-mer matches assigning reads to the right oligo?
 */
import { encodeFile } from "../src/lib/dna/codec";
import { simulate, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { ULTIMATE_NANOPORE_V52_CONFIG } from "../src/lib/dna/presets";
import { whitenAddress, bytesToDna, dnaToBytes, unwhitenAddress } from "../src/lib/dna/mapping";
import { buildReferenceKmerIndex, matchReadToReference, extractKmers } from "../src/lib/dna/kmer";

const TAG = "[v59-misclust]";

async function main() {
  const payload = new Uint8Array(16 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;

  const cfg = ULTIMATE_NANOPORE_V52_CONFIG;
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const sim = simulate(enc.encoded.oligos, { ...PRESET_NANOPORE, coverage: 10, seed: 42 });

  const fwd = enc.encoded.forwardPrimer;
  const expectedDnaLen = 260;

  const refs: string[] = [];
  for (let idx = 0; idx < enc.encoded.oligos.length; idx++) {
    const raw = new Uint8Array(4);
    raw[0] = (idx >> 16) & 0xff;
    raw[1] = (idx >> 8) & 0xff;
    raw[2] = idx & 0xff;
    raw[3] = 0;
    refs.push(bytesToDna(whitenAddress(raw)));
  }
  const kmerIdx = buildReferenceKmerIndex(refs, 5);

  let correct = 0, wrong = 0, unassigned = 0;
  const wrongPairs: { true: number; assigned: number }[] = [];

  for (const read of sim.reads) {
    const seq = read.sequence;
    const trueIdx = read.oligoIndex;

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
    if (fwdBestOverlap < 3) { unassigned++; continue; }

    const innerStart = fwdBestPos + fwd.length;
    const inner = seq.slice(innerStart, innerStart + expectedDnaLen);
    if (inner.length < 16) { unassigned++; continue; }

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
      const result = matchReadToReference(addressDna, kmerIdx, 5, 3);
      if (result.bestIdx >= 0) assignedIdx = result.bestIdx;
    }

    if (assignedIdx === -1) {
      unassigned++;
    } else if (assignedIdx === trueIdx) {
      correct++;
    } else {
      wrong++;
      wrongPairs.push({ true: trueIdx, assigned: assignedIdx });
    }
  }

  console.log(`${TAG} Clustering accuracy:`);
  console.log(`${TAG}   Correct: ${correct}/${sim.reads.length} (${(correct / sim.reads.length * 100).toFixed(1)}%)`);
  console.log(`${TAG}   Wrong: ${wrong}/${sim.reads.length} (${(wrong / sim.reads.length * 100).toFixed(1)}%)`);
  console.log(`${TAG}   Unassigned: ${unassigned}/${sim.reads.length} (${(unassigned / sim.reads.length * 100).toFixed(1)}%)`);

  if (wrongPairs.length > 0) {
    console.log(`\n${TAG} Wrong assignments (first 10):`);
    for (const p of wrongPairs.slice(0, 10)) {
      console.log(`${TAG}   true=${p.true} → assigned=${p.assigned}`);
    }
  }
}
main().catch(e => { console.error(e); process.exit(1); });
