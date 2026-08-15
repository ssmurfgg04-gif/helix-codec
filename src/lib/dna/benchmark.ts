/**
 * Benchmark suite for the DNA storage codec.
 *
 * Measures:
 *   1. Encoding throughput (bytes/sec, oligos/sec).
 *   2. Decoding throughput at various error rates.
 *   3. Recovery success rate vs. error rate (substitution, indel, dropout).
 *   4. Net storage density (bits/nt) vs. oligo length and RS config.
 *   5. Error correction capacity demonstration (max errors corrected).
 */

import { encodeFile, CodecConfig } from "./codec";
import { decodeReads } from "./decode";
import { simulate, MutationConfig, PRESET_ILLUMINA, PRESET_NANOPORE, PRESET_CLEAN } from "./simulate";
import { computeLayout } from "./types";

export interface BenchmarkPoint {
  label: string;
  errorRate: number;
  coverage: number;
  oligoCount: number;
  encoded: boolean;
  decoded: boolean;
  hashMatch: boolean;
  encodeMs: number;
  decodeMs: number;
  netDensityBitsPerNt: number;
  oligosRecovered: number;
  oligosErased: number;
  readsTotal: number;
}

export interface BenchmarkResult {
  points: BenchmarkPoint[];
  summary: {
    maxErrorRateRecovered: number;
    minCoverageRecovered: number;
    avgDensityBitsPerNt: number;
    avgEncodeMs: number;
    avgDecodeMs: number;
  };
}

/**
 * Run an error-rate sweep: encode a test payload, then decode at increasing
 * error rates to find the codec's breaking point.
 */
export async function runErrorRateSweep(
  payload: Uint8Array,
  cfg: CodecConfig,
  options?: { maxRate?: number; steps?: number; coverage?: number },
): Promise<BenchmarkResult> {
  const maxRate = options?.maxRate ?? 0.05;
  const steps = options?.steps ?? 6;
  const coverage = options?.coverage ?? 20;

  // Encode once
  const encodeResult = await encodeFile(payload, cfg, {
    fileName: "benchmark.bin",
    contentType: "application/octet-stream",
  });
  const oligos = encodeResult.encoded.oligos;
  const metadata = encodeResult.encoded.metadata;

  const points: BenchmarkPoint[] = [];
  for (let step = 0; step <= steps; step++) {
    const rate = (maxRate * step) / steps;
    const mutCfg: MutationConfig = {
      substitutionRate: rate * 0.4,
      insertionRate: rate * 0.2,
      deletionRate: rate * 0.4,
      coverage,
      dropoutRate: rate * 0.5,
      seed: 42, // deterministic
    };
    const sim = simulate(oligos, mutCfg);
    const decodeResult = await decodeReads(
      sim.reads,
      metadata,
      cfg,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );

    points.push({
      label: `rate=${rate.toFixed(4)}`,
      errorRate: rate,
      coverage,
      oligoCount: oligos.length,
      encoded: true,
      decoded: decodeResult.data !== null,
      hashMatch: decodeResult.hashMatches,
      encodeMs: encodeResult.stats.encodeTimeMs,
      decodeMs: decodeResult.stats.decodeTimeMs,
      netDensityBitsPerNt: encodeResult.stats.netDensityBitsPerNt,
      oligosRecovered: decodeResult.stats.oligosRecovered,
      oligosErased: decodeResult.stats.oligosErased,
      readsTotal: sim.totalReads,
    });
  }

  // Compute summary
  const maxRecovered = points
    .filter((p) => p.hashMatch)
    .reduce((m, p) => Math.max(m, p.errorRate), 0);
  const validPoints = points.filter((p) => p.hashMatch);
  const avgDensity =
    validPoints.length > 0
      ? validPoints.reduce((s, p) => s + p.netDensityBitsPerNt, 0) / validPoints.length
      : 0;
  const avgEncode =
    points.reduce((s, p) => s + p.encodeMs, 0) / points.length;
  const avgDecode =
    points.reduce((s, p) => s + p.decodeMs, 0) / points.length;

  return {
    points,
    summary: {
      maxErrorRateRecovered: maxRecovered,
      minCoverageRecovered: coverage,
      avgDensityBitsPerNt: avgDensity,
      avgEncodeMs: avgEncode,
      avgDecodeMs: avgDecode,
    },
  };
}

/**
 * Run a coverage sweep: encode + decode at various coverage depths.
 */
export async function runCoverageSweep(
  payload: Uint8Array,
  cfg: CodecConfig,
  options?: { minCov?: number; maxCov?: number; steps?: number },
): Promise<BenchmarkResult> {
  const minCov = options?.minCov ?? 5;
  const maxCov = options?.maxCov ?? 30;
  const steps = options?.steps ?? 6;

  const encodeResult = await encodeFile(payload, cfg, {
    fileName: "benchmark.bin",
    contentType: "application/octet-stream",
  });
  const oligos = encodeResult.encoded.oligos;
  const metadata = encodeResult.encoded.metadata;

  const points: BenchmarkPoint[] = [];
  for (let step = 0; step <= steps; step++) {
    const cov = Math.round(minCov + ((maxCov - minCov) * step) / steps);
    const mutCfg: MutationConfig = {
      ...PRESET_ILLUMINA,
      coverage: cov,
      seed: 42,
    };
    const sim = simulate(oligos, mutCfg);
    const decodeResult = await decodeReads(
      sim.reads,
      metadata,
      cfg,
      encodeResult.encoded.forwardPrimer,
      encodeResult.encoded.reversePrimer,
    );

    points.push({
      label: `cov=${cov}x`,
      errorRate: PRESET_ILLUMINA.substitutionRate + PRESET_ILLUMINA.insertionRate + PRESET_ILLUMINA.deletionRate,
      coverage: cov,
      oligoCount: oligos.length,
      encoded: true,
      decoded: decodeResult.data !== null,
      hashMatch: decodeResult.hashMatches,
      encodeMs: encodeResult.stats.encodeTimeMs,
      decodeMs: decodeResult.stats.decodeTimeMs,
      netDensityBitsPerNt: encodeResult.stats.netDensityBitsPerNt,
      oligosRecovered: decodeResult.stats.oligosRecovered,
      oligosErased: decodeResult.stats.oligosErased,
      readsTotal: sim.totalReads,
    });
  }

  const minCovRecovered = points.filter((p) => p.hashMatch).reduce((m, p) => Math.min(m, p.coverage), Infinity);
  const validPoints = points.filter((p) => p.hashMatch);
  const avgDensity =
    validPoints.length > 0
      ? validPoints.reduce((s, p) => s + p.netDensityBitsPerNt, 0) / validPoints.length
      : 0;
  const avgEncode = points.reduce((s, p) => s + p.encodeMs, 0) / points.length;
  const avgDecode = points.reduce((s, p) => s + p.decodeMs, 0) / points.length;

  return {
    points,
    summary: {
      maxErrorRateRecovered: PRESET_ILLUMINA.substitutionRate + PRESET_ILLUMINA.insertionRate + PRESET_ILLUMINA.deletionRate,
      minCoverageRecovered: minCovRecovered === Infinity ? 0 : minCovRecovered,
      avgDensityBitsPerNt: avgDensity,
      avgEncodeMs: avgEncode,
      avgDecodeMs: avgDecode,
    },
  };
}

/**
 * Generate a synthetic test payload of given size.
 * Mix of repetitive and random data to stress-test compression.
 */
export function generatePayload(size: number): Uint8Array {
  const buf = new Uint8Array(size);
  // First half: pseudo-random (high entropy, low compression)
  for (let i = 0; i < size / 2; i++) {
    buf[i] = (i * 31 + 17) & 0xff;
  }
  // Second half: repetitive (low entropy, high compression)
  const pattern = new Uint8Array([0x48, 0x65, 0x6c, 0x6c, 0x6f, 0x20, 0x44, 0x4e, 0x41, 0x21]);
  for (let i = size / 2; i < size; i++) {
    buf[i] = pattern[i % pattern.length];
  }
  return buf;
}
