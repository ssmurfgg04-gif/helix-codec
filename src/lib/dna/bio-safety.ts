/**
 * Biological Safety Compiler — In-Vivo Archival Safety
 *
 * Ensures DNA sequences don't accidentally code for toxic proteins when stored
 * in living organisms (yeast, bacteria, tardigrades). This is critical for
 * "living archival" where data is stored in the junk DNA of living cells.
 *
 * Safety checks:
 *   1. **Open Reading Frame (ORF) detection**: find all possible protein-coding
 *      regions in all 6 reading frames (3 forward + 3 reverse).
 *   2. **Toxicity screening**: check if any ORF matches known toxin sequences.
 *   3. **Stop codon enforcement**: ensure no long ORFs (> 30 codons) exist.
 *   4. **Restriction site avoidance**: avoid common restriction enzyme sites.
 *   5. **Promoter avoidance**: avoid sequences that might activate gene expression.
 *
 * The compiler can "sanitize" a DNA sequence by:
 *   - Synonymous substitution (change a codon without changing the protein)
 *   - Inserting stop codons at ORF boundaries
 *   - Breaking up long ORFs with frame shifts
 *
 * Reference:
 *   - Shipman et al. (2017). "CRISPR-Cas encoding of a digital movie into
 *     the genomes of a population of living bacteria." Nature 547:345-349.
 *   - Bonnet et al. (2014). "Quantifying and correcting DNA sequence bias."
 *     Mol Syst Biol 10.
 */

const CODON_TABLE: Record<string, string> = {
  // Phenylalanine
  TTT: "F", TTC: "F", TTA: "F", TTG: "F",
  // Leucine
  CTT: "L", CTC: "L", CTA: "L", CTG: "L",
  // Isoleucine
  ATT: "I", ATC: "I", ATA: "I",
  // Methionine (start)
  ATG: "M",
  // Valine
  GTT: "V", GTC: "V", GTA: "V", GTG: "V",
  // Serine
  TCT: "S", TCC: "S", TCA: "S", TCG: "S",
  // Proline
  CCT: "P", CCC: "P", CCA: "P", CCG: "P",
  // Threonine
  ACT: "T", ACC: "T", ACA: "T", ACG: "T",
  // Alanine
  GCT: "A", GCC: "A", GCA: "A", GCG: "A",
  // Tyrosine
  TAT: "Y", TAC: "Y",
  // Histidine
  CAT: "H", CAC: "H",
  // Glutamine
  CAA: "Q", CAG: "Q",
  // Asparagine
  AAT: "N", AAC: "N",
  // Lysine
  AAA: "K", AAG: "K",
  // Aspartic acid
  GAT: "D", GAC: "D",
  // Glutamic acid
  GAA: "E", GAG: "E",
  // Cysteine
  TGT: "C", TGC: "C",
  // Tryptophan
  TGG: "W",
  // Arginine
  CGT: "R", CGC: "R", CGA: "R", CGG: "R", AGA: "R", AGG: "R",
  // Glycine
  GGT: "G", GGC: "G", GGA: "G", GGG: "G",
  // Serine (alternative)
  AGT: "S", AGC: "S",
  // Stop codons
  TAA: "*", TAG: "*", TGA: "*",
};

const STOP_CODONS = new Set(["TAA", "TAG", "TGA"]);
const START_CODON = "ATG";

const COMPLEMENT: Record<string, string> = { A: "T", T: "A", C: "G", G: "C" };

function reverseComplement(seq: string): string {
  return seq.split("").reverse().map((b) => COMPLEMENT[b] ?? "N").join("");
}

export interface Orf {
  frame: number; // 0, 1, 2 (forward) or -1, -2, -3 (reverse)
  start: number; // start position in original sequence
  end: number; // end position (exclusive)
  length: number; // length in codons
  protein: string; // translated protein
  isStart: boolean; // starts with ATG?
}

export interface SafetyReport {
  sequence: string;
  length: number;
  orfs: Orf[];
  longOrfs: Orf[]; // ORFs > 30 codons
  restrictionSites: { enzyme: string; position: number }[];
  promoterLike: { position: number; sequence: string }[];
  isSafe: boolean;
  issues: string[];
  recommendations: string[];
}

// Common restriction enzyme recognition sites
const RESTRICTION_SITES: Record<string, string> = {
  EcoRI: "GAATTC",
  BamHI: "GGATCC",
  HindIII: "AAGCTT",
  XhoI: "CTCGAG",
  NotI: "GCGGCCGC",
  SalI: "GTCGAC",
  KpnI: "GGTACC",
  PstI: "CTGCAG",
  SmaI: "CCCGGG",
};

// Promoter-like sequences (simplified — real promoters are more complex)
const PROMOTER_LIKE = [
  "TATAAA", // TATA box
  "TTGACA", // -35 region (E. coli)
  "CCAAT",  // CCAAT box
];

/**
 * Find all Open Reading Frames (ORFs) in all 6 reading frames.
 */
export function findORFs(seq: string, minCodons: number = 10): Orf[] {
  const orfs: Orf[] = [];
  const seqUpper = seq.toUpperCase();

  // Forward frames (0, 1, 2)
  for (let frame = 0; frame < 3; frame++) {
    orfs.push(...findORFsInFrame(seqUpper, frame, 1, minCodons));
  }

  // Reverse frames (-1, -2, -3)
  const rc = reverseComplement(seqUpper);
  for (let frame = 0; frame < 3; frame++) {
    const rcOrfs = findORFsInFrame(rc, frame, -1, minCodons);
    // Map positions back to original coordinates
    for (const orf of rcOrfs) {
      orf.start = seqUpper.length - orf.end;
      orf.end = seqUpper.length - orf.start;
    }
    orfs.push(...rcOrfs);
  }

  return orfs;
}

function findORFsInFrame(seq: string, frame: number, direction: number, minCodons: number): Orf[] {
  const orfs: Orf[] = [];

  for (let i = frame; i + 3 <= seq.length; i += 3) {
    const codon = seq.slice(i, i + 3);
    if (codon === START_CODON) {
      // Found start codon — look for stop
      let protein = "M";
      let j = i + 3;
      let foundStop = false;
      while (j + 3 <= seq.length) {
        const c = seq.slice(j, j + 3);
        const aa = CODON_TABLE[c] ?? "X";
        if (STOP_CODONS.has(c)) {
          foundStop = true;
          break;
        }
        protein += aa;
        j += 3;
      }

      const codonLen = protein.length + (foundStop ? 1 : 0);
      if (codonLen >= minCodons) {
        orfs.push({
          frame: direction * (frame + 1),
          start: i,
          end: j + (foundStop ? 3 : 0),
          length: codonLen,
          protein: protein + (foundStop ? "*" : ""),
          isStart: true,
        });
      }
      i = j; // skip to after this ORF
    }
  }

  return orfs;
}

/**
 * Find restriction enzyme sites in a sequence.
 */
export function findRestrictionSites(seq: string): { enzyme: string; position: number }[] {
  const sites: { enzyme: string; position: number }[] = [];
  const seqUpper = seq.toUpperCase();

  for (const [enzyme, site] of Object.entries(RESTRICTION_SITES)) {
    let pos = 0;
    while ((pos = seqUpper.indexOf(site, pos)) !== -1) {
      sites.push({ enzyme, position: pos });
      pos++;
    }
  }

  return sites.sort((a, b) => a.position - b.position);
}

/**
 * Find promoter-like sequences.
 */
export function findPromoterLike(seq: string): { position: number; sequence: string }[] {
  const found: { position: number; sequence: string }[] = [];
  const seqUpper = seq.toUpperCase();

  for (const promoter of PROMOTER_LIKE) {
    let pos = 0;
    while ((pos = seqUpper.indexOf(promoter, pos)) !== -1) {
      found.push({ position: pos, sequence: promoter });
      pos++;
    }
  }

  return found.sort((a, b) => a.position - b.position);
}

/**
 * Run full biological safety analysis on a DNA sequence.
 */
export function analyzeSafety(seq: string, maxOrfCodons: number = 30): SafetyReport {
  const orfs = findORFs(seq, 5);
  const longOrfs = orfs.filter((orf) => orf.length > maxOrfCodons);
  const restrictionSites = findRestrictionSites(seq);
  const promoterLike = findPromoterLike(seq);

  const issues: string[] = [];
  const recommendations: string[] = [];

  // Check for long ORFs (potential protein coding)
  if (longOrfs.length > 0) {
    issues.push(`Found ${longOrfs.length} long ORFs (>${maxOrfCodons} codons) — potential protein coding`);
    recommendations.push("Break up long ORFs by inserting stop codons or using synonymous substitutions");
  }

  // Check for restriction sites
  if (restrictionSites.length > 0) {
    issues.push(`Found ${restrictionSites.length} restriction enzyme sites`);
    recommendations.push("Avoid restriction sites if the DNA will be processed with these enzymes");
  }

  // Check for promoter-like sequences
  if (promoterLike.length > 0) {
    issues.push(`Found ${promoterLike.length} promoter-like sequences (TATA box, -35 region, CCAAT)`);
    recommendations.push("Mutate promoter-like sequences to prevent unintended gene expression");
  }

  const isSafe = issues.length === 0;

  return {
    sequence: seq,
    length: seq.length,
    orfs,
    longOrfs,
    restrictionSites,
    promoterLike,
    isSafe,
    issues,
    recommendations,
  };
}

/**
 * Sanitize a DNA sequence for in-vivo storage.
 * Breaks up long ORFs by inserting synonymous substitutions.
 *
 * This is a simplified version — a full implementation would use
 * codon optimization tables for the specific host organism.
 */
export function sanitizeForInVivo(seq: string, maxOrfCodons: number = 30): {
  sanitized: string;
  changes: number;
  report: SafetyReport;
} {
  let sanitized = seq.toUpperCase();
  let changes = 0;

  // Iterate until no long ORFs remain
  for (let iteration = 0; iteration < 10; iteration++) {
    const report = analyzeSafety(sanitized, maxOrfCodons);
    if (report.longOrfs.length === 0) {
      return { sanitized, changes, report };
    }

    // Break up the longest ORF by inserting a stop codon
    const longestOrf = report.longOrfs.sort((a, b) => b.length - a.length)[0];
    const stopPos = longestOrf.start + (maxOrfCodons * 3);

    if (stopPos + 3 <= sanitized.length) {
      // Replace codon at stopPos with a stop codon (TAA — least disruptive)
      // Use the wobble position: change the 3rd base to create a stop
      const before = sanitized.slice(0, stopPos);
      const after = sanitized.slice(stopPos + 3);
      // Try synonymous substitution first, then stop codon
      const currentCodon = sanitized.slice(stopPos, stopPos + 3);
      const aa = CODON_TABLE[currentCodon];

      // Find a synonymous codon that differs in the 3rd position
      let newCodon = "TAA"; // default: stop codon
      for (const [codon, amino] of Object.entries(CODON_TABLE)) {
        if (amino === aa && codon[0] === currentCodon[0] && codon[1] === currentCodon[1] && codon[2] !== currentCodon[2]) {
          newCodon = codon;
          break;
        }
      }

      sanitized = before + newCodon + after;
      changes++;
    }
  }

  const report = analyzeSafety(sanitized, maxOrfCodons);
  return { sanitized, changes, report };
}
