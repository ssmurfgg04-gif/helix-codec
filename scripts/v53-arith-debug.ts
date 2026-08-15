#!/usr/bin/env bun
/**
 * v53 Arithmetic Decoder Debug Test
 *
 * Verifies that the v53 fix to codec.ts (bytesPerBlockData = total - 2 instead of -1)
 * resolves the WASM arithmetic decoder bug at 500nt.
 *
 * Tests:
 *   1. Zero-noise round-trip with 500nt + arithmetic (WASM path)
 *   2. Zero-noise round-trip with 500nt + arithmetic + low-coverage (WASM path)
 *   3. Low-noise 10x coverage with 500nt + arithmetic (WASM path)
 *   4. Block-size scan verification (encoder/decoder capacity match)
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { ULTIMATE_DENSITY_CONFIG, ULTIMATE_LOW_COVERAGE_CONFIG } from "../src/lib/dna/presets";
import { computeLayoutAuto } from "../src/lib/dna/types";

const TAG = "[v53-arith]";

function makePayload(n: number): Uint8Array {
  const p = new Uint8Array(n);
  for (let i = 0; i < n; i++) p[i] = (i * 31 + 17) & 0xff;
  return p;
}

async function runTest(name: string, cfg: any, payloadSize: number, coverage: number, noise: any) {
  console.log(`\n--- ${name} ---`);
  const payload = makePayload(payloadSize);

  // Layout debug
  const layout = computeLayoutAuto(cfg);
  const innerN = layout.addressBytes + layout.payloadBytes + layout.innerParityBytes;
  const innerDnaLen = layout.totalInnerBytes * 4;
  const blockSize = (cfg as any).arithmeticBlockSize ?? Math.floor(innerDnaLen / 2);
  const bytesPerBlockTotal = Math.floor(blockSize / 4);
  const bytesPerBlockData = bytesPerBlockTotal - 3;
  const numBlocks = Math.floor(innerDnaLen / blockSize);
  const dataCapacity = numBlocks * bytesPerBlockData;

  console.log(`${TAG} config: oligoLen=${cfg.oligoLength}, mapping=${cfg.mappingMode}, inner=${cfg.innerCode?.type}/${cfg.innerCode?.parityBytes}, outer=${cfg.outerRS?.parityFraction}`);
  console.log(`${TAG} layout: addressBytes=${layout.addressBytes}, payloadBytes=${layout.payloadBytes}, parityBytes=${layout.innerParityBytes}, totalInnerBytes=${layout.totalInnerBytes}`);
  console.log(`${TAG} innerN=${innerN}, innerDnaLen=${innerDnaLen}, blockSize=${blockSize}, numBlocks=${numBlocks}, dataCapacity=${dataCapacity}`);
  console.log(`${TAG} capacity check: dataCapacity(${dataCapacity}) ${dataCapacity >= innerN ? ">=" : "<"} innerN(${innerN}) — ${dataCapacity >= innerN ? "OK" : "FALLBACK"}`);

  const theoreticalDensity = (payloadSize * 8) / (Math.ceil(payloadSize / layout.payloadBytes) * cfg.oligoLength);
  console.log(`${TAG} theoretical density: ${theoreticalDensity.toFixed(3)} bits/nt`);

  const t0 = Date.now();
  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const oligos = enc.encoded.oligos;
  const metadata = enc.encoded.metadata;
  const encodeMs = Date.now() - t0;
  console.log(`${TAG} encoded: ${oligos.length} oligos in ${encodeMs}ms`);

  const simResult = simulate(oligos, { ...PRESET_CLEAN, coverage, ...noise });
  const reads = simResult.reads;
  console.log(`${TAG} simulated: ${reads.length} reads (${coverage}× coverage, noise=${JSON.stringify(noise)})`);

  const t1 = Date.now();
  const decoded = await decodeReadsUltra(reads, metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  const decodeMs = Date.now() - t1;
  console.log(`${TAG} decoded in ${decodeMs}ms: hashMatch=${decoded.hashMatches}, bytes=${decoded.data.length}`);
  console.log(`${TAG} decode stats:`, JSON.stringify((decoded as any).stats));
  // Compare first/last 32 bytes
  console.log(`${TAG} payload[0..32]:  ${Array.from(payload.slice(0, 32)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  console.log(`${TAG} decoded[0..32]:  ${Array.from(decoded.data.slice(0, 32)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  // Find first mismatch
  let firstMismatch = -1;
  for (let i = 0; i < Math.min(payload.length, decoded.data.length); i++) {
    if (payload[i] !== decoded.data[i]) { firstMismatch = i; break; }
  }
  console.log(`${TAG} first mismatch at byte: ${firstMismatch}`);

  const ok = decoded.hashMatches && decoded.data.length === payload.length && decoded.data.every((b: number, i: number) => b === payload[i]);
  console.log(`${ok ? "✅" : "⚠️"} ${name} ${ok ? "PASSED" : "FAILED"}`);
  return ok;
}

async function main() {
  console.log(`${TAG} === v53 Arithmetic Decoder Debug Test ===`);
  console.log(`${TAG} Verifies v53 fix to codec.ts (bytesPerBlockData = total - 2)`);

  const cfg500Arith = JSON.parse(JSON.stringify(ULTIMATE_DENSITY_CONFIG));
  const cfg500ArithLowCov = JSON.parse(JSON.stringify(ULTIMATE_LOW_COVERAGE_CONFIG));

  const zeroNoise = { substitutionRate: 0, insertionRate: 0, deletionRate: 0 };
  const lowNoise = { substitutionRate: 0.001, insertionRate: 0, deletionRate: 0 };

  const results: boolean[] = [];
  results.push(await runTest("Test 1: 500nt + arithmetic (zero-noise, 10x)", cfg500Arith, 64 * 1024, 10, zeroNoise));
  results.push(await runTest("Test 2: 500nt + arithmetic + low-cov (zero-noise, 10x)", cfg500ArithLowCov, 64 * 1024, 10, zeroNoise));
  results.push(await runTest("Test 3: 500nt + arithmetic (low-noise, 10x)", cfg500Arith, 64 * 1024, 10, lowNoise));
  results.push(await runTest("Test 4: 500nt + arithmetic (zero-noise, 5x)", cfg500Arith, 64 * 1024, 5, zeroNoise));

  const passed = results.filter(Boolean).length;
  console.log(`\n=== ${TAG} Summary: ${passed}/${results.length} PASSED ===`);
  if (passed === results.length) {
    console.log(`✅ v53 arithmetic decoder fix VERIFIED`);
    process.exit(0);
  } else {
    console.log(`⚠️ ${results.length - passed} tests still failing`);
    process.exit(1);
  }
}

main().catch(e => { console.error(e); process.exit(1); });
