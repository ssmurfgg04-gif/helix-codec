/**
 * K=9 Penalty Tuning — sweep insertionPenalty/deletionPenalty/maxDrift to find
 * the configuration with the lowest BER on a noisy Nanopore-like channel.
 *
 * Test matrix:
 *   - Ins/Del penalties: [0.5, 1.0, 1.5, 2.0, 2.5]
 *   - Max drift: [5, 10, 15, 20, 25]
 *   - For each combo, run K=9 encode + indel decode on 32-byte payloads with
 *     simulated 1-3 random indels (insertions and deletions)
 *   - Track: decode success rate, mean BER, mean decode time
 *
 * Output: ranked config list + recommended config saved to disk.
 */
import * as fs from 'node:fs';
import * as path from 'node:path';

import {
  enableNativeViterbi,
  nativeConvK9Encode,
  nativeViterbiK9Decode,
  nativeViterbiK9DecodeStandard,
} from '../src/lib/dna/native/viterbi-napi';

interface TestResult {
  ins: number;
  del: number;
  drift: number;
  success: number;       // 0..1 fraction
  meanBer: number;       // mean bit error rate
  meanMs: number;        // mean decode time ms
  samples: number;
}

function hex(b: Uint8Array): string {
  return Buffer.from(b).toString('hex');
}

/** Simulate `nIns` insertions + `nDel` deletions at random bit positions. */
function applyIndels(bits: number[], nIns: number, nDel: number, rng: () => number): number[] {
  const out = bits.slice();
  for (let i = 0; i < nIns; i++) {
    const pos = Math.floor(rng() * (out.length + 1));
    out.splice(pos, 0, Math.floor(rng() * 2));
  }
  for (let i = 0; i < nDel; i++) {
    if (out.length === 0) break;
    const pos = Math.floor(rng() * out.length);
    out.splice(pos, 1);
  }
  return out;
}

function bytesToBits(b: Uint8Array): number[] {
  const out: number[] = [];
  for (let i = 0; i < b.length; i++) {
    for (let k = 7; k >= 0; k--) out.push((b[i] >> k) & 1);
  }
  return out;
}

function bitsToBytes(bits: number[]): Uint8Array {
  const n = Math.floor(bits.length / 8);
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    let v = 0;
    for (let k = 7; k >= 0; k--) v |= (bits[i * 8 + (7 - k)] & 1) << k;
    out[i] = v;
  }
  return out;
}

// Mulberry32 PRNG for reproducibility
function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

async function main() {
  console.log('=== K=9 Penalty Tuning ===\n');

  const enabled = await enableNativeViterbi();
  if (!enabled) {
    console.error('Native Viterbi addon not loaded — aborting.');
    process.exit(1);
  }

  // Test payloads — 32 bytes each, multiple samples
  const N_SAMPLES = 5;
  const payloads: Uint8Array[] = [];
  const rng0 = mulberry32(42);
  for (let i = 0; i < N_SAMPLES; i++) {
    const p = new Uint8Array(32);
    for (let j = 0; j < 32; j++) p[j] = Math.floor(rng0() * 256);
    payloads.push(p);
  }

  // Pre-encode all payloads once
  console.log(`Pre-encoding ${N_SAMPLES} payloads (32 bytes each)...`);
  const encoded = payloads.map(p => nativeConvK9Encode(p));
  console.log(`Encoded ${encoded.length} payloads, each ${encoded[0].length} bytes\n`);

  // Test matrix — reduced for fast iteration
  const penalties = [1.0, 1.5, 2.0];
  const drifts = [10, 15, 20];

  // Indel scenarios: (insCount, delCount) — represent realistic ONT R10.4.1
  const scenarios = [
    { name: 'clean', ins: 0, del: 0 },
    { name: 'light-1i-1d', ins: 1, del: 1 },
    { name: 'med-3i-2d', ins: 3, del: 2 },
  ];

  const results: TestResult[] = [];

  console.log('Sweeping configs...\n');
  console.log('Config            | clean    | light    | med      | heavy    ');
  console.log('------------------|----------|----------|----------|----------');

  for (const ins of penalties) {
    for (const del of penalties) {
      for (const drift of drifts) {
        const row: string[] = [];
        let totalSuccess = 0;
        let totalSamples = 0;
        let totalBer = 0;
        let totalMs = 0;
        let measuredSamples = 0;

        for (const sc of scenarios) {
          let ok = 0;
          let berSum = 0;
          let timeSum = 0;
          let measured = 0;
          for (let i = 0; i < N_SAMPLES; i++) {
            // Apply indels at the BIT level
            const bits = bytesToBits(encoded[i]);
            const corrupted = applyIndels(bits, sc.ins, sc.del, mulberry32(1000 + i * 17 + sc.ins * 31 + sc.del * 7));
            const corruptedBytes = bitsToBytes(corrupted);

            const t0 = Date.now();
            let decoded: Uint8Array;
            try {
              decoded = nativeViterbiK9Decode(corruptedBytes, {
                maxDrift: drift,
                insertionPenalty: ins,
                deletionPenalty: del,
                numInfoBits: 32 * 8,
              });
            } catch {
              timeSum += Date.now() - t0;
              continue;
            }
            const dt = Date.now() - t0;
            timeSum += dt;

            // Compute BER vs original payload
            const trimmed = decoded.length > 32 ? decoded.slice(0, 32) : decoded;
            let errs = 0;
            for (let j = 0; j < 32; j++) {
              if (trimmed[j] !== payloads[i][j]) errs++;
            }
            const ber = errs / 32;
            berSum += ber;
            if (ber === 0) ok++;
            if (sc.name !== 'clean') {
              totalBer += ber;
              totalMs += dt;
              measured++;
            }
          }
          const succRate = ok / N_SAMPLES;
          totalSuccess += ok;
          totalSamples += N_SAMPLES;
          row.push(`${(succRate * 100).toFixed(0).padStart(3)}%`);
          if (sc.name !== 'clean') {
            measuredSamples += measured;
          }
        }

        // Aggregate stats for noisy scenarios only
        const overallSuccess = totalSuccess / totalSamples;
        const meanBer = measuredSamples > 0 ? totalBer / measuredSamples : 0;
        const meanMs = measuredSamples > 0 ? totalMs / measuredSamples : 0;

        results.push({
          ins, del, drift,
          success: overallSuccess,
          meanBer,
          meanMs,
          samples: measuredSamples,
        });

        // Print row only for selected combos to keep table compact
        if ((ins === 1.5 || ins === 2.0) && (del === 1.5 || del === 2.0)) {
          console.log(`ins=${ins} del=${del} d=${String(drift).padStart(2)} | ${row.join(' | ')}`);
        }
      }
    }
  }

  // Print top-10 by lowest mean BER
  console.log('\n=== Top 10 configs by lowest mean BER (across noisy scenarios) ===');
  const sorted = results.slice().sort((a, b) => a.meanBer - b.meanBer);
  console.log('Rank | ins  | del  | drift | succ%  | meanBER | meanMs');
  console.log('-----|------|------|-------|--------|---------|-------');
  for (let i = 0; i < Math.min(10, sorted.length); i++) {
    const r = sorted[i];
    console.log(
      `${String(i + 1).padStart(4)} | ${r.ins.toFixed(1)} | ${r.del.toFixed(1)} | ${String(r.drift).padStart(5)} | ${(r.success * 100).toFixed(1).padStart(5)}% | ${r.meanBer.toFixed(4)}  | ${r.meanMs.toFixed(1)}`,
    );
  }

  // Save results
  const outPath = '/home/z/my-project/datasets/k9-penalty-tuning.json';
  fs.writeFileSync(outPath, JSON.stringify({
    generatedAt: new Date().toISOString(),
    samplesPerScenario: N_SAMPLES,
    scenarios: scenarios.map(s => s.name),
    results: sorted,
  }, null, 2));
  console.log(`\nFull results saved to ${outPath}`);

  // Print recommended config
  const best = sorted[0];
  console.log(`\n=== Recommended config ===`);
  console.log(`  insertionPenalty: ${best.ins}`);
  console.log(`  deletionPenalty:  ${best.del}`);
  console.log(`  maxDrift:         ${best.drift}`);
  console.log(`  success rate:     ${(best.success * 100).toFixed(1)}%`);
  console.log(`  mean BER:         ${best.meanBer.toFixed(4)}`);
  console.log(`  mean decode ms:   ${best.meanMs.toFixed(2)}`);
}

main().catch(e => { console.error(e); process.exit(1); });
