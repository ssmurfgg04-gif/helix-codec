// Test LDPC Belief-Propagation decoder vs hard-decision at various error rates.
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
  console.log("=== LDPC Belief-Propagation Decoder Test ===\n");

  const k = 34, n = 38;
  const ldpc = new LDPCInnerCode({ n, k });
  console.log(`LDPC(${n * 8} bits, ${k * 8} bits, ${(n - k) * 8} bits parity)`);
  console.log(`Comparing: hard-decision syndrome lookup vs belief-propagation\n`);

  for (const numErrors of [3, 4, 5, 6, 8, 10]) {
    let hardPass = 0, bpPass = 0, total = 100;
    let bpTime = 0, hardTime = 0;

    for (let trial = 0; trial < total; trial++) {
      const info = randomBytes(k, trial);
      const codeword = ldpc.encode(info);
      const bits = bytesToBits(codeword);
      const { bits: noisyBits, positions } = injectBitErrors(bits, numErrors, trial + numErrors * 1000);
      const noisyBytes = bitsToBytes(noisyBits);

      // Build realistic Q-scores: low (Q10) for error positions, high (Q35) for correct
      const qScores = new Uint8Array(n * 4);
      for (let i = 0; i < qScores.length; i++) qScores[i] = 35;
      for (const p of positions) {
        const baseIdx = Math.floor(p / 2);
        if (baseIdx < qScores.length) qScores[baseIdx] = 10;
      }

      // Hard-decision
      const t0 = Date.now();
      try {
        const r = ldpc.decode(noisyBytes);
        if (r.data.every((b, i) => b === info[i])) hardPass++;
      } catch {}
      hardTime += Date.now() - t0;

      // Belief-propagation
      const t1 = Date.now();
      try {
        const r = ldpc.decodeBeliefPropagation(noisyBytes, qScores, false, 20);
        if (r.data.every((b, i) => b === info[i])) bpPass++;
      } catch {}
      bpTime += Date.now() - t1;
    }

    console.log(`${numErrors} errors: hard=${hardPass}/${total} (${hardTime}ms), BP=${bpPass}/${total} (${bpTime}ms)`);
  }
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
