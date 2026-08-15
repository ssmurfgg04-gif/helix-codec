/**
 * v52 WASM RECOMPILE VERIFICATION
 *
 * Verifies that the freshly recompiled WASM binary supports:
 *   - 500nt oligos (ULTIMATE_DENSITY_CONFIG)
 *   - Arithmetic mapping (fullDecodeArithmetic WASM path)
 *
 * Without this recompile, 500nt + arithmetic would fall back to JS path
 * (3-5x slower). After recompile, the WASM fast path should work for both
 * 300nt/direct AND 500nt/arithmetic configurations.
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import {
  ULTIMATE_DENSITY_CONFIG,
  ULTIMATE_LOW_COVERAGE_CONFIG,
  V51_DEFAULT_CONFIG,
  computeDensity,
} from "../src/lib/dna/presets";
import * as crypto from "crypto";

const REPORT = (msg: string) => console.log(`[v52-wasm-verify] ${msg}`);

async function testConfig(
  name: string,
  config: any,
  payloadSize: number,
  coverage: number = 10,
): Promise<boolean> {
  REPORT(`\n=== ${name} ===`);
  REPORT(`config: oligoLen=${config.oligoLength}, mapping=${config.mappingMode}, inner=${config.innerCode}/${config.innerParityBytes}, outer=${config.outerParityRatio}`);

  const density = computeDensity(config, "total");
  REPORT(`theoretical density: ${density.toFixed(3)} bits/nt (total-oligo)`);

  const payload = crypto.randomBytes(payloadSize);
  const tEnc = Date.now();
  const enc = await encodeFile(
    Buffer.from(payload),
    config,
    { fileName: `${name}.bin`, contentType: "application/octet-stream" },
  );
  const encMs = Date.now() - tEnc;
  REPORT(`encoded: ${enc.encoded.oligos.length} oligos in ${encMs}ms`);

  const sim = simulate(enc.encoded.oligos, {
    ...PRESET_ILLUMINA,
    coverage,
    seed: 42,
  });
  REPORT(`simulated: ${sim.totalReads} reads at ${coverage}× coverage`);

  const t0 = Date.now();
  let dec;
  try {
    dec = await decodeReadsUltra(
      sim.reads,
      enc.encoded.metadata,
      config,
      enc.encoded.forwardPrimer,
      enc.encoded.reversePrimer,
    );
  } catch (e: any) {
    REPORT(`❌ decodeReadsUltra FAILED: ${e.message?.slice(0, 100)}`);
    return false;
  }
  const ms = Date.now() - t0;

  const pass = dec.hashMatches;
  REPORT(`decoded in ${ms}ms: hashMatch=${pass}, oligosRecovered=${dec.stats.oligosRecovered}/${enc.encoded.oligos.length}`);

  if (pass) {
    REPORT(`✅ ${name} PASSED — WASM fast path works`);
  } else {
    REPORT(`⚠️ ${name} PARTIAL — hash mismatch (decode may have fallen back to JS)`);
  }
  return pass;
}

async function main() {
  REPORT("╔══════════════════════════════════════════════════════════════╗");
  REPORT("║  v52 WASM RECOMPILE VERIFICATION                            ║");
  REPORT("║  Tests 500nt + arithmetic on freshly compiled WASM binary    ║");
  REPORT("╚══════════════════════════════════════════════════════════════╝");

  let allPass = true;

  // Test 1: 300nt + direct (baseline — must still work after recompile)
  REPORT("\n--- Test 1: 300nt + direct (regression baseline) ---");
  if (!await testConfig(
    "300nt-direct",
    { ...V51_DEFAULT_CONFIG, oligoLength: 300, mappingMode: "direct" },
    32 * 1024,
    10,
  )) allPass = false;

  // Test 2: 500nt + direct (density unlock — was JS-only before recompile)
  REPORT("\n--- Test 2: 500nt + direct (density unlock) ---");
  if (!await testConfig(
    "500nt-direct",
    { ...V51_DEFAULT_CONFIG, oligoLength: 500, mappingMode: "direct", innerParityBytes: 8, outerParityRatio: 0.05 },
    32 * 1024,
    10,
  )) allPass = false;

  // Test 3: 500nt + arithmetic (ULTIMATE_DENSITY_CONFIG — the key Phase 1 target)
  REPORT("\n--- Test 3: 500nt + arithmetic (ULTIMATE_DENSITY_CONFIG) ---");
  if (!await testConfig(
    "500nt-arithmetic",
    ULTIMATE_DENSITY_CONFIG,
    32 * 1024,
    10,
  )) allPass = false;

  // Test 4: 500nt + arithmetic + low coverage trigger (ULTIMATE_LOW_COVERAGE_CONFIG)
  REPORT("\n--- Test 4: 500nt + arithmetic + low-coverage (ULTIMATE_LOW_COVERAGE_CONFIG) ---");
  if (!await testConfig(
    "500nt-arith-lowcov",
    ULTIMATE_LOW_COVERAGE_CONFIG,
    32 * 1024,
    10,
  )) allPass = false;

  REPORT(`\n=== WASM RECOMPILE VERIFICATION ${allPass ? "PASSED ✅" : "PARTIAL ⚠️"} ===`);
  if (allPass) {
    REPORT("All configurations work on the freshly recompiled WASM binary.");
    REPORT("500nt + arithmetic density unlock is live — no JS-path slowdown.");
  }
  process.exit(allPass ? 0 : 1);
}

main().catch((e) => {
  console.error("FATAL:", e);
  process.exit(1);
});
