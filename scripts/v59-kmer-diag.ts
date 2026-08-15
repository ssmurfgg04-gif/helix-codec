/**
 * v59 K-mer Clustering Diagnostic — understand why reads are discarded.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { simulate, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { ULTIMATE_NANOPORE_V52_CONFIG } from "../src/lib/dna/presets";
import { whitenAddress, bytesToDna, dnaToBytes, unwhitenAddress } from "../src/lib/dna/mapping";
import { buildReferenceKmerIndex, matchReadToReference, kmerToBits, extractKmers } from "../src/lib/dna/kmer";

const TAG = "[v59-diag]";

async function main() {
  const payload = new Uint8Array(16 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = (i * 31 + 17) & 0xff;

  const cfg = ULTIMATE_NANOPORE_V52_CONFIG;
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  console.log(`${TAG} Encoded ${enc.encoded.oligos.length} oligos`);

  const sim = simulate(enc.encoded.oligos, { ...PRESET_NANOPORE, coverage: 10, seed: 42 });
  console.log(`${TAG} Simulated ${sim.reads.length} reads`);

  const fwd = enc.encoded.forwardPrimer;
  const rev = enc.encoded.reversePrimer;
  console.log(`${TAG} Fwd primer: ${fwd} (${fwd.length}nt)`);
  console.log(`${TAG} Rev primer: ${rev} (${rev.length}nt)`);

  // Build reference addresses
  const oligoCount = enc.encoded.oligos.length;
  const refs: string[] = [];
  for (let idx = 0; idx < oligoCount; idx++) {
    const raw = new Uint8Array(4);
    raw[0] = (idx >> 16) & 0xff;
    raw[1] = (idx >> 8) & 0xff;
    raw[2] = idx & 0xff;
    raw[3] = 0;
    refs.push(bytesToDna(whitenAddress(raw)));
  }
  const kmerIdx = buildReferenceKmerIndex(refs, 4);
  console.log(`${TAG} Built k-mer index: ${kmerIdx.size} unique k-mers from ${refs.length} refs`);

  // Stats
  let noFwdPrimer = 0;
  let noRevPrimer = 0;
  let tooShort = 0;
  let exactMatch = 0;
  let kmerMatch = 0;
  let noMatch = 0;
  let kmerLowOverlap = 0;

  // Distribution of overlap scores
  const overlapHist = new Map<number, number>();

  for (let i = 0; i < Math.min(sim.reads.length, 300); i++) {
    const read = sim.reads[i];
    const seq = read.sequence;
    const trueIdx = read.oligoIndex;

    // Forward primer search (tolerant)
    const fwdLen = fwd.length;
    let fwdBestDist = Infinity;
    let fwdBestPos = 0;
    for (let p = 0; p <= 15; p++) {
      if (p + fwdLen > seq.length) break;
      const candidate = seq.slice(p, p + fwdLen);
      let d = 0;
      for (let j = 0; j < fwdLen; j++) {
        if (candidate[j] !== fwd[j]) { d++; if (d > 6) break; }
      }
      if (d < fwdBestDist) { fwdBestDist = d; fwdBestPos = p; if (d <= 2) break; }
    }
    if (fwdBestDist > 5) { noFwdPrimer++; continue; }

    // Get inner DNA (use expected length)
    const innerStart = fwdBestPos + fwdLen;
    const expectedInnerLen = 260; // 300nt - 2*20nt primers
    const inner = seq.slice(innerStart, innerStart + expectedInnerLen);
    if (inner.length < 16) { tooShort++; continue; }

    // Try exact address match
    const addressDna = inner.slice(0, 16);
    let assignedIdx = -1;
    try {
      const addrBytes = dnaToBytes(addressDna);
      const unwhitened = unwhitenAddress(addrBytes);
      const idx = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
      if (idx >= 0 && idx < oligoCount) {
        const refDna = bytesToDna(whitenAddress(unwhitened));
        let dist = 0;
        for (let j = 0; j < 16; j++) if (addressDna[j] !== refDna[j]) dist++;
        if (dist <= 2) {
          assignedIdx = idx;
          exactMatch++;
        }
      }
    } catch {}

    if (assignedIdx === -1) {
      // Try k-mer matching
      const result = matchReadToReference(addressDna, kmerIdx, 4, 2);
      overlapHist.set(result.bestOverlap, (overlapHist.get(result.bestOverlap) ?? 0) + 1);
      if (result.bestIdx >= 0) {
        assignedIdx = result.bestIdx;
        kmerMatch++;
        if (assignedIdx !== trueIdx) {
          // Miscluster!
          // console.log(`${TAG} MISCLUSTER read ${i}: true=${trueIdx} assigned=${assignedIdx} overlap=${result.bestOverlap}`);
        }
      } else {
        noMatch++;
        if (result.bestOverlap > 0) kmerLowOverlap++;
      }
    }
  }

  console.log(`\n${TAG} Clustering stats (300 reads):`);
  console.log(`${TAG}   Exact match: ${exactMatch}`);
  console.log(`${TAG}   K-mer match: ${kmerMatch}`);
  console.log(`${TAG}   No match: ${noMatch} (low overlap: ${kmerLowOverlap})`);
  console.log(`${TAG}   No fwd primer: ${noFwdPrimer}`);
  console.log(`${TAG}   Too short: ${tooShort}`);
  console.log(`${TAG}   Total clustered: ${exactMatch + kmerMatch}/${300}`);

  console.log(`\n${TAG} Overlap histogram (k-mer matches):`);
  for (const [overlap, count] of Array.from(overlapHist.entries()).sort((a, b) => a[0] - b[0])) {
    console.log(`${TAG}   overlap=${overlap}: ${count} reads`);
  }

  // Also check: what's the address length distribution?
  console.log(`\n${TAG} Read length stats:`);
  const lengths = sim.reads.map(r => r.sequence.length).sort((a, b) => a - b);
  console.log(`${TAG}   min: ${lengths[0]}, max: ${lengths[lengths.length - 1]}, median: ${lengths[Math.floor(lengths.length / 2)]}`);
  console.log(`${TAG}   expected: ~300nt (with indels: 270-330)`);
}
main().catch(e => { console.error(e); process.exit(1); });
