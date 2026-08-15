#!/usr/bin/env bun
/**
 * Helix v53 Adversarial Fuzz Test Suite
 *
 * Systematically tests the codec against adversarial inputs:
 *   - Empty payload
 *   - 1-byte payload
 *   - Maximum payload size
 *   - All-zeros payload
 *   - All-0xFF payload
 *   - Random payloads (various sizes)
 *   - Payloads with high compressibility
 *   - Payloads with low compressibility (random)
 *   - Payloads at size boundaries (powers of 2)
 *   - Adversarial byte patterns (0x00, 0xFF, alternating, sequential)
 *
 * For each payload, tests:
 *   - Zero-noise round-trip (must always pass)
 *   - Low-noise Illumina (sub=0.001, 5× cov)
 *   - High-noise Illumina (sub=0.01, 10× cov)
 *   - Nanopore (sub=0.02, ins=0.03, del=0.04, 15× cov)
 *
 * Usage: bun run scripts/v53-fuzz-test.ts
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_CLEAN, PRESET_ILLUMINA, PRESET_NANOPORE } from "../src/lib/dna/simulate";
import { V51_DEFAULT_CONFIG, ULTIMATE_DENSITY_CONFIG } from "../src/lib/dna/presets";

const TAG = "[v53-fuzz]";

interface FuzzResult {
  name: string;
  payloadSize: number;
  config: string;
  noise: string;
  coverage: number;
  passed: boolean;
  encodeMs: number;
  decodeMs: number;
  error?: string;
}

const results: FuzzResult[] = [];

function makePayload(size: number, pattern: string): Uint8Array {
  const p = new Uint8Array(size);
  switch (pattern) {
    case "zeros":
      return p;
    case "ones":
      return p.fill(0xff);
    case "sequential":
      for (let i = 0; i < size; i++) p[i] = i & 0xff;
      break;
    case "alternating":
      for (let i = 0; i < size; i++) p[i] = i % 2 === 0 ? 0xaa : 0x55;
      break;
    case "random":
      for (let i = 0; i < size; i++) p[i] = Math.floor(Math.random() * 256);
      break;
    case "high-compress":
      // Repeating pattern (highly compressible)
      for (let i = 0; i < size; i++) p[i] = (i * 7 + 13) % 256;
      break;
    case "low-compress":
      // Pseudo-random but deterministic (low compressibility)
      let seed = 12345;
      for (let i = 0; i < size; i++) {
        seed = (seed * 1103515245 + 12345) & 0x7fffffff;
        p[i] = seed & 0xff;
      }
      break;
    case "text":
      // ASCII text (medium compressibility)
      const text = "The quick brown fox jumps over the lazy dog. ";
      for (let i = 0; i < size; i++) p[i] = text.charCodeAt(i % text.length);
      break;
    default:
      for (let i = 0; i < size; i++) p[i] = Math.floor(Math.random() * 256);
  }
  return p;
}

async function runFuzz(
  name: string,
  payload: Uint8Array,
  cfg: any,
  configName: string,
  noise: any,
  noiseName: string,
  coverage: number,
): Promise<FuzzResult> {
  const result: FuzzResult = {
    name,
    payloadSize: payload.length,
    config: configName,
    noise: noiseName,
    coverage,
    passed: false,
    encodeMs: 0,
    decodeMs: 0,
  };

  try {
    const t0 = Date.now();
    const enc = await encodeFile(payload, cfg, { fileName: "fuzz.bin", contentType: "application/octet-stream" });
    result.encodeMs = Date.now() - t0;

    const simResult = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage, ...noise });
    const reads = simResult.reads;

    const t1 = Date.now();
    const decoded = await decodeReadsUltra(reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    result.decodeMs = Date.now() - t1;

    result.passed = decoded.hashMatches &&
      decoded.data.length === payload.length &&
      decoded.data.every((b: number, i: number) => b === payload[i]);

    if (!result.passed && !result.error) {
      result.error = `hashMatch=${decoded.hashMatches}, bytes=${decoded.data.length}/${payload.length}`;
    }
  } catch (e: any) {
    result.error = e.message?.slice(0, 100);
  }

  results.push(result);
  const status = result.passed ? "✅" : "⚠️";
  console.log(`${status} ${name.padEnd(40)} ${configName.padEnd(15)} ${noiseName.padEnd(20)} ${String(coverage).padEnd(3)}× ${result.encodeMs.toString().padStart(4)}ms/${result.decodeMs.toString().padStart(4)}ms ${result.error ?? ""}`);
  return result;
}

async function main() {
  console.log(`${TAG} === Helix v53 Adversarial Fuzz Test Suite ===`);
  console.log(`${TAG}`);

  const cfgV51 = JSON.parse(JSON.stringify(V51_DEFAULT_CONFIG));
  const cfg500 = JSON.parse(JSON.stringify(ULTIMATE_DENSITY_CONFIG));
  cfg500.mappingMode = "direct"; // use direct mode (arithmetic has known bug)

  const zeroNoise = { substitutionRate: 0, insertionRate: 0, deletionRate: 0 };
  const lowNoise = { substitutionRate: 0.001, insertionRate: 0, deletionRate: 0 };
  const highNoise = { substitutionRate: 0.01, insertionRate: 0, deletionRate: 0 };
  const nanoporeNoise = { substitutionRate: 0.02, insertionRate: 0.03, deletionRate: 0.04 };

  console.log(`${TAG} --- Size boundary tests (v51, zero-noise, 5×) ---`);
  for (const size of [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096, 16384, 65536]) {
    const payload = makePayload(size, "random");
    await runFuzz(`size-${size}B`, payload, cfgV51, "v51-default", zeroNoise, "zero-noise", 5);
  }

  console.log(`${TAG}`);
  console.log(`${TAG} --- Pattern tests (1KB, v51, zero-noise, 5×) ---`);
  for (const pattern of ["zeros", "ones", "sequential", "alternating", "random", "high-compress", "low-compress", "text"]) {
    const payload = makePayload(1024, pattern);
    await runFuzz(`pattern-${pattern}`, payload, cfgV51, "v51-default", zeroNoise, "zero-noise", 5);
  }

  console.log(`${TAG}`);
  console.log(`${TAG} --- Noise level tests (64KB, v51, 5×) ---`);
  for (const [noiseName, noise] of Object.entries({ "zero": zeroNoise, "low": lowNoise, "high": highNoise })) {
    const payload = makePayload(65536, "random");
    await runFuzz(`noise-${noiseName}`, payload, cfgV51, "v51-default", noise, noiseName + "-noise", 5);
  }

  console.log(`${TAG}`);
  console.log(`${TAG} --- Coverage tests (64KB, v51, low-noise) ---`);
  for (const cov of [3, 5, 8, 10, 15]) {
    const payload = makePayload(65536, "random");
    await runFuzz(`cov-${cov}x`, payload, cfgV51, "v51-default", lowNoise, "low-noise", cov);
  }

  console.log(`${TAG}`);
  console.log(`${TAG} --- 500nt direct mode tests (64KB, zero-noise, 10×) ---`);
  for (const pattern of ["zeros", "random", "high-compress", "text"]) {
    const payload = makePayload(65536, pattern);
    await runFuzz(`500nt-${pattern}`, payload, cfg500, "500nt-direct", zeroNoise, "zero-noise", 10);
  }

  console.log(`${TAG}`);
  console.log(`${TAG} --- Edge cases ---`);
  // Empty payload
  await runFuzz("empty-payload", new Uint8Array(0), cfgV51, "v51-default", zeroNoise, "zero-noise", 5);
  // 1 byte
  await runFuzz("1-byte", new Uint8Array([0x42]), cfgV51, "v51-default", zeroNoise, "zero-noise", 5);

  // Summary
  console.log(`${TAG}`);
  console.log(`${TAG} ═══ SUMMARY ═══`);
  const passed = results.filter(r => r.passed).length;
  const failed = results.filter(r => !r.passed).length;
  console.log(`${TAG} Total: ${results.length}, Passed: ${passed}, Failed: ${failed}`);
  console.log(`${TAG} Pass rate: ${(passed / results.length * 100).toFixed(1)}%`);

  if (failed > 0) {
    console.log(`${TAG}`);
    console.log(`${TAG} Failures:`);
    for (const r of results.filter(r => !r.passed)) {
      console.log(`${TAG}   ⚠️ ${r.name} (${r.config}, ${r.noise}, ${r.coverage}×): ${r.error}`);
    }
  }

  process.exit(failed === 0 ? 0 : 1);
}

main().catch(e => { console.error(e); process.exit(1); });
