/**
 * SantaLucia Nearest-Neighbor Thermodynamics + Primer Design
 *
 * Implements the SantaLucia 1998 unified nearest-neighbor model for DNA
 * melting temperature (Tm) prediction, plus hairpin and self-dimer ΔG
 * calculation. This is the core algorithm used by Primer3.
 *
 * Tm = ΔH / (ΔS + R * ln(C_T / x)) - 273.15
 * where:
 *   ΔH = sum of nearest-neighbor enthalpies (kcal/mol)
 *   ΔS = sum of nearest-neighbor entropies (cal/(mol·K))
 *   R = 1.987 cal/(mol·K) (gas constant)
 *   C_T = total strand concentration (M)
 *   x = 4 for non-self-complementary, 1 for self-complementary
 *
 * Nearest-neighbor parameters (SantaLucia 1998, Table 1):
 *   Each NN pair has ΔH (kcal/mol) and ΔS (cal/(mol·K))
 *
 * Reference:
 *   - SantaLucia (1998). "A unified view of polymer, dumbbell, and
 *     oligonucleotide DNA nearest-neighbor thermodynamics."
 *     PNAS 95:1460-1465.
 *   - Owczarzy et al. (2004). "Predicting sequence-dependent and
 *     concentration-dependent DNA melting temperatures." Biochemistry 43:12.
 *   - Primer3: github.com/primer3-org/primer3
 */

// Nearest-neighbor ΔH (kcal/mol) and ΔS (cal/(mol·K)) from SantaLucia 1998
const NN_PARAMS: Record<string, { dH: number; dS: number }> = {
  // 5'-3' / 3'-5' pairs
  "AA/TT": { dH: -7.9, dS: -22.2 },
  "AT/AT": { dH: -7.2, dS: -20.4 },
  "TA/TA": { dH: -7.2, dS: -21.3 },
  "CA/GT": { dH: -8.5, dS: -22.7 },
  "GT/CA": { dH: -8.4, dS: -22.4 },
  "CT/GA": { dH: -7.8, dS: -21.0 },
  "GA/CT": { dH: -8.2, dS: -22.2 },
  "CG/CG": { dH: -10.6, dS: -27.2 },
  "GC/GC": { dH: -9.8, dS: -24.4 },
  "GG/CC": { dH: -8.0, dS: -19.9 },
};

// Initiation parameters
const INIT_GC = { dH: 0.1, dS: -2.8 }; // initiation with terminal G·C
const INIT_AT = { dH: 2.3, dS: 4.1 }; // initiation with terminal A·T

const R = 1.987; // gas constant cal/(mol·K)
const DEFAULT_CT = 250e-9; // 250 nM total strand concentration
const DEFAULT_NA = 50e-3; // 50 mM Na+

/**
 * Get the nearest-neighbor key for a dinucleotide step.
 * Returns the canonical form (handles both strands).
 */
function getNNKey(dimer: string): string {
  // dimer is 2 bases, e.g. "AC"
  // The NN parameter key is "AC/GT" but we need to find the canonical form
  const complement: Record<string, string> = { A: "T", T: "A", C: "G", G: "C" };
  const top = dimer;
  const bottom = dimer
    .split("")
    .map((b) => complement[b] ?? "N")
    .reverse()
    .join("");
  const key1 = `${top}/${bottom}`;
  // Also check reverse complement
  const rcTop = top
    .split("")
    .map((b) => complement[b] ?? "N")
    .reverse()
    .join("");
  const rcBottom = rcTop
    .split("")
    .map((b) => complement[b] ?? "N")
    .reverse()
    .join("");
  const key2 = `${rcTop}/${rcBottom}`;

  if (NN_PARAMS[key1]) return key1;
  if (NN_PARAMS[key2]) return key2;
  return key1; // default
}

/**
 * Calculate melting temperature (Tm) using SantaLucia 1998 NN model.
 *
 * @param seq DNA sequence (5'→3')
 * @param naConc Sodium concentration in M (default 50 mM)
 * @param strandConc Total strand concentration in M (default 250 nM)
 * @returns Tm in °C
 */
export function calculateTm(
  seq: string,
  naConc: number = DEFAULT_NA,
  strandConc: number = DEFAULT_CT,
): number {
  if (seq.length < 2) return 0;

  let dH = 0; // kcal/mol
  let dS = 0; // cal/(mol·K)

  // Sum nearest-neighbor parameters
  for (let i = 0; i < seq.length - 1; i++) {
    const dimer = seq.slice(i, i + 2);
    const key = getNNKey(dimer);
    const params = NN_PARAMS[key];
    if (params) {
      dH += params.dH;
      dS += params.dS;
    }
  }

  // Initiation terms (terminal bases)
  const firstBase = seq[0];
  const lastBase = seq[seq.length - 1];
  if (firstBase === "G" || firstBase === "C") {
    dH += INIT_GC.dH;
    dS += INIT_GC.dS;
  } else {
    dH += INIT_AT.dH;
    dS += INIT_AT.dS;
  }
  if (lastBase === "G" || lastBase === "C") {
    dH += INIT_GC.dH;
    dS += INIT_GC.dS;
  } else {
    dH += INIT_AT.dH;
    dS += INIT_AT.dS;
  }

  // Self-complementary check
  const isSelfComp = seq === reverseComplement(seq);
  const x = isSelfComp ? 1 : 4;

  // Tm = ΔH / (ΔS + R * ln(C_T / x)) - 273.15
  // ΔH in kcal/mol, ΔS in cal/(mol·K) — convert ΔH to cal
  const dH_cal = dH * 1000;
  const tm = dH_cal / (dS + R * Math.log(strandConc / x)) - 273.15;

  // Salt correction (Owczarzy 2004)
  const saltCorrectedTm = tm + 16.6 * Math.log10(naConc);

  return saltCorrectedTm;
}

function reverseComplement(seq: string): string {
  const comp: Record<string, string> = { A: "T", T: "A", C: "G", G: "C" };
  return seq
    .split("")
    .reverse()
    .map((b) => comp[b] ?? "N")
    .join("");
}

/**
 * Calculate GC content of a sequence.
 */
export function gcContent(seq: string): number {
  if (seq.length === 0) return 0;
  let gc = 0;
  for (const b of seq) if (b === "G" || b === "C") gc++;
  return gc / seq.length;
}

/**
 * Detect hairpins in a primer sequence.
 * Returns the most stable hairpin's ΔG (kcal/mol), or 0 if none.
 *
 * A hairpin requires a stem ≥ 3 bp and loop ≥ 3 nt.
 */
export function hairpinDG(seq: string): number {
  let minDG = 0;
  const minStem = 3;
  const minLoop = 3;

  for (let i = 0; i < seq.length - 2 * minStem - minLoop; i++) {
    for (let stemLen = minStem; stemLen <= 8; stemLen++) {
      for (let loopLen = minLoop; loopLen <= 12; loopLen++) {
        const arm1Start = i;
        const arm2Start = i + stemLen + loopLen;
        if (arm2Start + stemLen > seq.length) break;

        const arm1 = seq.slice(arm1Start, arm1Start + stemLen);
        const arm2 = seq.slice(arm2Start, arm2Start + stemLen);

        if (arm1 === reverseComplement(arm2)) {
          // Found a hairpin — estimate ΔG from stem stability
          let dg = 0;
          for (let k = 0; k < stemLen - 1; k++) {
            const key = getNNKey(arm1.slice(k, k + 2));
            const params = NN_PARAMS[key];
            if (params) {
              // ΔG = ΔH - T*ΔS at 37°C (310.15 K)
              dg += params.dH - 310.15 * (params.dS / 1000);
            }
          }
          // Loop penalty (approximate: +3 kcal/mol for loops 3-12 nt)
          dg += 3.0 + Math.abs(loopLen - 5) * 0.3;
          if (dg < minDG) minDG = dg;
        }
      }
    }
  }

  return minDG;
}

/**
 * Check for self-dimers (primer self-complementarity at 3' end).
 * Returns the ΔG of the most stable self-dimer, or 0 if none.
 */
export function selfDimerDG(seq: string): number {
  let minDG = 0;
  // Check all alignments of seq with its reverse complement
  const rc = reverseComplement(seq);
  for (let offset = -seq.length + 4; offset < seq.length - 4; offset++) {
    let matches = 0;
    let dg = 0;
    for (let i = 0; i < seq.length; i++) {
      const j = i + offset;
      if (j < 0 || j >= rc.length) continue;
      if (seq[i] === rc[j]) {
        matches++;
        // Add NN ΔG for matched pair
        if (i > 0 && i - 1 + offset >= 0 && i - 1 + offset < rc.length) {
          if (seq[i - 1] === rc[j - 1]) {
            const key = getNNKey(seq.slice(i - 1, i + 1));
            const params = NN_PARAMS[key];
            if (params) {
              dg += params.dH - 310.15 * (params.dS / 1000);
            }
          }
        }
      }
    }
    if (matches >= 4 && dg < minDG) minDG = dg;
  }
  return minDG;
}

export interface PrimerScore {
  /** Overall quality score (0-100, higher = better). */
  score: number;
  /** Melting temperature in °C. */
  tm: number;
  /** GC content (0-1). */
  gc: number;
  /** Hairpin ΔG (kcal/mol, more negative = worse). */
  hairpinDG: number;
  /** Self-dimer ΔG (kcal/mol, more negative = worse). */
  selfDimerDG: number;
  /** Length in nucleotides. */
  length: number;
  /** Issues found. */
  issues: string[];
}

/**
 * Score a primer sequence for quality.
 * Returns a score (0-100) and detailed metrics.
 *
 * Ideal primer (SantaLucia / Primer3 guidelines):
 *   - Tm: 55-62°C (optimal 58°C)
 *   - GC: 40-60% (optimal 50%)
 *   - Length: 18-22 nt (optimal 20 nt)
 *   - No hairpins (ΔG > -2 kcal/mol)
 *   - No self-dimers (ΔG > -2 kcal/mol)
 *   - Ends in G or C (GC clamp)
 */
export function scorePrimer(seq: string): PrimerScore {
  const tm = calculateTm(seq);
  const gc = gcContent(seq);
  const hpDG = hairpinDG(seq);
  const sdDG = selfDimerDG(seq);
  const issues: string[] = [];
  let score = 100;

  // Tm scoring (optimal 58°C)
  const tmDiff = Math.abs(tm - 58);
  if (tmDiff > 10) {
    score -= 30;
    issues.push(`Tm ${tm.toFixed(1)}°C far from optimal 58°C`);
  } else if (tmDiff > 5) {
    score -= 15;
    issues.push(`Tm ${tm.toFixed(1)}°C suboptimal`);
  }

  // GC content
  if (gc < 0.4 || gc > 0.6) {
    score -= 20;
    issues.push(`GC ${(gc * 100).toFixed(0)}% outside 40-60%`);
  } else if (gc < 0.45 || gc > 0.55) {
    score -= 10;
    issues.push(`GC ${(gc * 100).toFixed(0)}% suboptimal`);
  }

  // Length
  if (seq.length < 18 || seq.length > 25) {
    score -= 15;
    issues.push(`Length ${seq.length} outside 18-25 nt`);
  }

  // Hairpin
  if (hpDG < -2) {
    score -= 25;
    issues.push(`Hairpin ΔG ${hpDG.toFixed(1)} kcal/mol (too stable)`);
  }

  // Self-dimer
  if (sdDG < -2) {
    score -= 20;
    issues.push(`Self-dimer ΔG ${sdDG.toFixed(1)} kcal/mol (too stable)`);
  }

  // GC clamp (3' end should be G or C)
  const lastBase = seq[seq.length - 1];
  if (lastBase !== "G" && lastBase !== "C") {
    score -= 10;
    issues.push("No GC clamp at 3' end");
  }

  // Homopolymer check
  let maxRun = 1;
  let currentRun = 1;
  for (let i = 1; i < seq.length; i++) {
    if (seq[i] === seq[i - 1]) {
      currentRun++;
      if (currentRun > maxRun) maxRun = currentRun;
    } else {
      currentRun = 1;
    }
  }
  if (maxRun >= 4) {
    score -= 15;
    issues.push(`Homopolymer run of ${maxRun} nt`);
  }

  return {
    score: Math.max(0, score),
    tm,
    gc,
    hairpinDG: hpDG,
    selfDimerDG: sdDG,
    length: seq.length,
    issues,
  };
}
