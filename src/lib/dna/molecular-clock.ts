/**
 * Molecular Clock — Generational Drift Tracker
 *
 * Models the depurination kinetics of DNA over time to predict how error
 * rates grow as an archive ages. This allows the decoder to adjust its
 * priors based on the archive's age — a novel "age-aware decoding" feature.
 *
 * Depurination is the primary DNA degradation mechanism: purine bases (A, G)
 * are lost from the sugar-phosphate backbone, leading to strand breakage and
 * errors. The rate depends on temperature, pH, and humidity.
 *
 * Model (from Allentoft 2012, Nature):
 *   depurination_rate(t) = k * t
 *   where k = 5.21e-6 per nt per year (at 25°C, dry)
 *
 * Mahoraga projects 282 years at 17.1 EB/g (25°C dry). At cooler temperatures,
 * the rate drops dramatically (Arrhenius: ~10× slower per 10°C decrease).
 *
 * Reference:
 *   - Allentoft et al. (2012). "The half-life of DNA in bone."
 *     Proc R Soc B 280:20120990.
 *   - Banal et al. (2026). arXiv:2604.20810. (Mahoraga longevity projection)
 *   - Arrhenius (1889). "Über die Reaktionsgeschwindigkeit..."
 */

export interface ArchiveAge {
  /** Age in years. */
  years: number;
  /** Storage temperature in °C. */
  temperatureC: number;
  /** Whether the archive was stored dry (vs. aqueous). */
  dry: boolean;
}

export interface DriftPrediction {
  /** Expected substitution rate. */
  substitutionRate: number;
  /** Expected deletion rate (dominated by depurination). */
  deletionRate: number;
  /** Expected insertion rate. */
  insertionRate: number;
  /** Expected total error rate. */
  totalErrorRate: number;
  /** Expected strand breakage probability. */
  strandBreakageProb: number;
  /** Recommended coverage for recovery. */
  recommendedCoverage: number;
  /** Recommended outer parity ratio. */
  recommendedParityRatio: number;
  /** Whether recovery is feasible. */
  recoveryFeasible: boolean;
}

// Base depurination rate at 25°C (per nt per year)
const BASE_DEPURINATION_RATE = 5.21e-6;

/**
 * Calculate the Arrhenius temperature correction factor.
 * DNA degradation rate ~10× slower per 10°C decrease.
 */
function temperatureFactor(tempC: number): number {
  const refTemp = 25; // reference temperature
  const tempDiff = tempC - refTemp;
  return Math.pow(10, tempDiff / 10);
}

/**
 * Dry storage reduces depurination by ~100× (water is needed for hydrolysis).
 */
function dryFactor(dry: boolean): number {
  return dry ? 0.01 : 1.0;
}

/**
 * Predict error rates based on archive age and storage conditions.
 *
 * @param age Archive age and storage conditions
 * @returns Predicted error rates and recommendations
 */
export function predictDrift(age: ArchiveAge): DriftPrediction {
  const tempFactor = temperatureFactor(age.temperatureC);
  const dry = dryFactor(age.dry);

  // Depurination rate (per nt per year) at actual conditions
  const depurRate = BASE_DEPURINATION_RATE * tempFactor * dry;

  // After `years` years, expected number of depurinations per nt
  const depurPerNt = depurRate * age.years;

  // Depurination leads to:
  //   - Deletions (purine base lost) — most common
  //   - Substitutions (misincorporation during repair/replication)
  //   - Insertions (rare)
  const deletionRate = Math.min(0.9, depurPerNt * 0.7); // 70% of depurinations become deletions
  const substitutionRate = Math.min(0.3, depurPerNt * 0.25); // 25% become substitutions
  const insertionRate = Math.min(0.1, depurPerNt * 0.05); // 5% become insertions
  const totalErrorRate = deletionRate + substitutionRate + insertionRate;

  // Strand breakage: probability that at least one break occurs in a 200nt strand
  // P(break) = 1 - (1 - depurRate)^(200 * years)
  const strandBreakageProb = 1 - Math.pow(1 - depurRate, 200 * age.years);

  // Recommend coverage and parity based on error rate
  // Higher errors → need more coverage and more parity
  let recommendedCoverage = 20;
  if (totalErrorRate > 0.01) recommendedCoverage = 30;
  if (totalErrorRate > 0.05) recommendedCoverage = 50;
  if (totalErrorRate > 0.10) recommendedCoverage = 100;

  let recommendedParityRatio = 0.2;
  if (totalErrorRate > 0.01) recommendedParityRatio = 0.3;
  if (totalErrorRate > 0.05) recommendedParityRatio = 0.5;
  if (totalErrorRate > 0.10) recommendedParityRatio = 0.8;

  // Recovery is feasible if total error rate < ~15% (codec's max tolerance)
  const recoveryFeasible = totalErrorRate < 0.15;

  return {
    substitutionRate,
    deletionRate,
    insertionRate,
    totalErrorRate,
    strandBreakageProb,
    recommendedCoverage,
    recommendedParityRatio,
    recoveryFeasible,
  };
}

/**
 * Estimate the maximum archival lifetime for a given codec configuration.
 *
 * @param maxErrorRate Maximum tolerable error rate (e.g., 0.10 = 10%)
 * @param tempC Storage temperature in °C
 * @param dry Whether storage is dry
 * @returns Estimated lifetime in years
 */
export function estimateLifetime(
  maxErrorRate: number,
  tempC: number = 25,
  dry: boolean = true,
): number {
  const tempFactor = temperatureFactor(tempC);
  const dryF = dryFactor(dry);
  const depurRate = BASE_DEPURINATION_RATE * tempFactor * dryF;

  // Solve: depurRate * years * 1.0 (total→error conversion) = maxErrorRate
  return maxErrorRate / depurRate;
}

/**
 * Generate a longevity report for an archive.
 */
export function longevityReport(
  age: ArchiveAge,
): {
  prediction: DriftPrediction;
  lifetime10pct: number; // years until 10% error rate
  lifetime15pct: number; // years until 15% (codec limit)
  summary: string;
} {
  const prediction = predictDrift(age);
  const lifetime10 = estimateLifetime(0.10, age.temperatureC, age.dry);
  const lifetime15 = estimateLifetime(0.15, age.temperatureC, age.dry);

  let summary: string;
  if (prediction.recoveryFeasible) {
    summary = `Archive age ${age.years}y at ${age.temperatureC}°C${age.dry ? " (dry)" : ""}: ${(prediction.totalErrorRate * 100).toFixed(2)}% total error. Recovery feasible. Can survive ${(lifetime15 - age.years).toFixed(0)} more years.`;
  } else {
    summary = `Archive age ${age.years}y at ${age.temperatureC}°C${age.dry ? " (dry)" : ""}: ${(prediction.totalErrorRate * 100).toFixed(2)}% total error. Recovery may fail. Exceeded codec tolerance.`;
  }

  return {
    prediction,
    lifetime10pct: lifetime10,
    lifetime15pct: lifetime15,
    summary,
  };
}
