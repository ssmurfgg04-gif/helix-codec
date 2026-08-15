// Direct OSD-2 test: inject errors that hard-decision can't fix, verify OSD-2 recovers.
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";

function randomBytes(n: number, seed: number): Uint8Array {
  let s = seed >>> 0;
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    s ^= s << 13; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
    out[i] = s & 0xff;
  }
  return out;
}

function bytesToBits(data: Uint8Array): Uint8Array {
  const bits = new Uint8Array(data.length * 8);
  for (let i = 0; i < data.length; i++) {
    for (let bit = 0; bit < 8; bit++) {
      bits[i * 8 + bit] = (data[i] >> (7 - bit)) & 1;
    }
  }
  return bits;
}

function bitsToBytes(bits: Uint8Array): Uint8Array {
  const out = new Uint8Array(bits.length / 8);
  for (let i = 0; i < out.length; i++) {
    let b = 0;
    for (let bit = 0; bit < 8; bit++) b |= bits[i * 8 + bit] << (7 - bit);
    out[i] = b;
  }
  return out;
}

function injectBitErrors(bits: Uint8Array, numErrors: number, seed: number): { bits: Uint8Array; positions: number[] } {
  let s = seed >>> 0;
  const positions = new Set<number>();
  while (positions.size < numErrors) {
    s ^= s << 13; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
    positions.add(s % bits.length);
  }
  const out = bits.slice();
  for (const p of positions) out[p] ^= 1;
  return { bits: out, positions: Array.from(positions) };
}

async function main() {
  console.log("=== Direct OSD-2 Test ===\n");

  // LDPC(304 bits, 272 bits, 32 bits parity) — same as 200nt oligo
  const k = 34, n = 38;
  const ldpc = new LDPCInnerCode({ n, k });

  console.log("Testing OSD-2 recovery with REALISTIC Q-scores (low for error positions):\n");
  console.log("(In real sequencing, substituted bases get Q5-Q15, correct bases get Q30-Q40)\n");

  for (const numErrors of [3, 4, 5, 6]) {
    let hardPass = 0, osd2Pass = 0, total = 0;
    for (let trial = 0; trial < 100; trial++) {
      const info = randomBytes(k, trial);
      const codeword = ldpc.encode(info);
      const bits = bytesToBits(codeword);
      const { bits: noisyBits, positions } = injectBitErrors(bits, numErrors, trial + numErrors * 1000);
      const noisyBytes = bitsToBytes(noisyBits);

      // Build realistic Q-scores: low (Q10) for error positions, high (Q35) for correct
      // In direct mapping: bit i → base i/2. Error at bit position p → base p/2.
      const qScores = new Uint8Array(n * 4);
      for (let i = 0; i < qScores.length; i++) qScores[i] = 35;
      for (const p of positions) {
        const baseIdx = Math.floor(p / 2);
        if (baseIdx < qScores.length) qScores[baseIdx] = 10; // low Q for error
      }

      total++;
      try {
        const r = ldpc.decode(noisyBytes);
        const ok = r.data.every((b, i) => b === info[i]);
        if (ok) hardPass++;
      } catch {
        try {
          const r = ldpc.decodeWithSoftInfo(noisyBytes, qScores, false);
          const ok = r.data.every((b, i) => b === info[i]);
          if (ok) osd2Pass++;
        } catch {}
      }
    }
    console.log(`${numErrors}-bit errors: hard=${hardPass}/100, OSD-2=${osd2Pass}/100, total=${hardPass + osd2Pass}/100`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
