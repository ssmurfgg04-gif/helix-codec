/**
 * v63 Nanopore K=9 IDS Recovery Test
 *
 * Tests the K=9 convolutional code (NASA standard, d_free=24) + indel-tolerant
 * Viterbi decoder at various IDS levels. Target: 90% recovery at 9% IDS.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, MutationConfig } from "../src/lib/dna/simulate";
import { ULTIMATE_V61_NANOPORE_CONFIG } from "../src/lib/dna/presets";

async function main() {
  const cfg = ULTIMATE_V61_NANOPORE_CONFIG;
  const payloadSize = 256; // 256B — K=9 indel Viterbi is ~1.4s/decode, keep tiny
  const payload = new Uint8Array(payloadSize);
  for (let i = 0; i < payloadSize; i++) payload[i] = Math.floor(Math.random() * 256);

  console.log("=== v63 Nanopore K=9 IDS Recovery Test ===");
  console.log(`config: 300nt oligos, K=9 conv (d_free=24), 8B LDPC, 25% outer RS`);
  console.log(`payload: 256B random data, 10× coverage\n`);

  const t0 = performance.now();
  const enc = await encodeFile(payload, cfg, { fileName: "b.bin", contentType: "application/octet-stream" });
  const encMs = performance.now() - t0;
  const { oligos, metadata } = enc.encoded;
  console.log(`encoded: ${oligos.length} oligos in ${encMs.toFixed(0)}ms\n`);

  const tests: Array<[string, number, number, number]> = [
    ["0% IDS (clean)",       0.00, 0.00, 0.00],
    ["1% IDS (minimal)",     0.005, 0.003, 0.002],
    ["3% IDS (mild)",        0.01, 0.01, 0.01],
    ["6% IDS (moderate)",    0.02, 0.02, 0.02],
    ["9% IDS (PRESET_NANOPORE)", 0.02, 0.03, 0.04],
    ["12% IDS (heavy)",      0.03, 0.04, 0.05],
    ["15% IDS (extreme)",    0.04, 0.05, 0.06],
  ];

  console.log("IDS Level                    | hash   | data   | decode  | recovered");
  console.log("-----------------------------|--------|--------|---------|----------");

  for (const [label, sub, ins, del] of tests) {
    const simCfg: MutationConfig = {
      substitutionRate: sub,
      insertionRate: ins,
      deletionRate: del,
      coverage: 10,
      dropoutRate: 0,
      seed: 42,
    };
    const { reads } = simulate(oligos, simCfg);

    const t1 = performance.now();
    let result;
    try {
      result = await decodeReads(reads, metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    } catch (e: any) {
      console.log(`${label.padEnd(29)} | ERROR  | -      | -       | ${e.message.slice(0, 40)}`);
      continue;
    }
    const decMs = performance.now() - t1;

    const hashMatch = result.hashMatches;
    const dataMatch = result.data
      ? Buffer.compare(Buffer.from(payload), Buffer.from(result.data)) === 0
      : false;

    console.log(
      `${label.padEnd(29)} | ${hashMatch ? "OK ✅ " : "FAIL ❌"} | ${dataMatch ? "OK ✅ " : "FAIL ❌"} | ${decMs.toFixed(0).padStart(5)}ms | ${result.stats.oligosRecovered}/${oligos.length}`,
    );
  }
}

main().catch(e => { console.error(e); process.exit(1); });
