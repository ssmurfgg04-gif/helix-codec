/**
 * 
 * ⚠️ EXPERIMENTAL / ROADMAP — This module is a scaffolding prototype.
 * It is NOT wired into the production encode/decode pipeline.
 * See CHANGELOG.md for v59 status and remaining work.
 *
 * Hardware API Stubs — Wet-Lab Bridge
 *
 * API client stubs for DNA synthesis and sequencing hardware:
 *   1. IDT (Integrated DNA Technologies) — oligo ordering API
 *   2. Twist Bioscience — silicon-based DNA synthesis API
 *   3. Oxford Nanopore MinKNOW — sequencing software output format
 *
 * These are STUB implementations that generate the correct file formats
 * and calculate costs. Real API integration requires vendor accounts.
 *
 * Commands:
 *   helix synthesize dna://pool/key --vendor idt
 *   helix synthesize dna://pool/key --vendor twist
 *   helix export-minknow dna://pool/key
 */

export interface OligoOrder {
  vendor: "idt" | "twist" | "elegen";
  oligos: { name: string; sequence: string; scale: string; purification: string }[];
  totalLength: number;
  totalOligos: number;
  estimatedCostUSD: number;
  estimatedTurnaroundDays: number;
  orderFormat: "csv" | "fasta" | "xlsx";
  orderContent: string;
}

export interface VendorPricing {
  costPerBase: number;
  minimumOrder: number;
  setupFee: number;
  turnaroundDays: number;
  maxOligoLength: number;
  maxOligosPerOrder: number;
}

const VENDOR_PRICING: Record<string, VendorPricing> = {
  idt: {
    costPerBase: 0.10, // $0.10 per base (standard IDT pricing)
    minimumOrder: 50,
    setupFee: 25,
    turnaroundDays: 3,
    maxOligoLength: 230,
    maxOligosPerOrder: 50000,
  },
  twist: {
    costPerBase: 0.07, // Twist is cheaper per base
    minimumOrder: 100,
    setupFee: 50,
    turnaroundDays: 7,
    maxOligoLength: 300,
    maxOligosPerOrder: 100000,
  },
  elegen: {
    costPerBase: 0.05, // Elegen is cheapest (newer technology)
    minimumOrder: 200,
    setupFee: 75,
    turnaroundDays: 10,
    maxOligoLength: 300,
    maxOligosPerOrder: 200000,
  },
};

/**
 * Generate an oligo order for a DNA synthesis vendor.
 *
 * @param oligos Array of oligo sequences to synthesize
 * @param vendor Target vendor ("idt", "twist", "elegen")
 * @param scale Synthesis scale ("25nm", "100nm", "250nm", "1um")
 * @param purification Purification method ("STD", "PAGE", "HPLC")
 */
export function createSynthesisOrder(
  oligos: { name: string; sequence: string }[],
  vendor: "idt" | "twist" | "elegen" = "idt",
  scale: string = "25nm",
  purification: string = "STD",
): OligoOrder {
  const pricing = VENDOR_PRICING[vendor];

  // Validate oligo lengths
  for (const oligo of oligos) {
    if (oligo.sequence.length > pricing.maxOligoLength) {
      throw new Error(
        `Oligo ${oligo.name} is ${oligo.sequence.length}nt, exceeds ${vendor} max of ${pricing.maxOligoLength}nt`,
      );
    }
  }

  if (oligos.length > pricing.maxOligosPerOrder) {
    throw new Error(
      `Order has ${oligos.length} oligos, exceeds ${vendor} max of ${pricing.maxOligosPerOrder}`,
    );
  }

  // Calculate cost
  const totalBases = oligos.reduce((sum, o) => sum + o.sequence.length, 0);
  const baseCost = totalBases * pricing.costPerBase;
  const totalCost = Math.max(baseCost + pricing.setupFee, pricing.minimumOrder);

  // Generate order file
  const orderOligos = oligos.map(o => ({
    name: o.name,
    sequence: o.sequence,
    scale,
    purification,
  }));

  let orderContent: string;
  let orderFormat: "csv" | "fasta" | "xlsx";

  switch (vendor) {
    case "idt":
      orderFormat = "csv";
      orderContent = generateIDTCsv(orderOligos);
      break;
    case "twist":
      orderFormat = "fasta";
      orderContent = generateTwistFasta(orderOligos);
      break;
    case "elegen":
      orderFormat = "csv";
      orderContent = generateElegenCsv(orderOligos);
      break;
  }

  return {
    vendor,
    oligos: orderOligos,
    totalLength: totalBases,
    totalOligos: oligos.length,
    estimatedCostUSD: totalCost,
    estimatedTurnaroundDays: pricing.turnaroundDays,
    orderFormat,
    orderContent,
  };
}

/**
 * Generate IDT-format CSV order file.
 * IDT accepts CSV with columns: Name, Sequence, Scale, Purification
 */
function generateIDTCsv(oligos: { name: string; sequence: string; scale: string; purification: string }[]): string {
  const lines = ["Name,Sequence,Scale,Purification"];
  for (const o of oligos) {
    lines.push(`${o.name},${o.sequence},${o.scale},${o.purification}`);
  }
  return lines.join("\n");
}

/**
 * Generate Twist Bioscience FASTA order file.
 * Twist accepts FASTA format with oligo names.
 */
function generateTwistFasta(oligos: { name: string; sequence: string; scale: string; purification: string }[]): string {
  const lines: string[] = [];
  for (const o of oligos) {
    lines.push(`>${o.name}`);
    lines.push(o.sequence);
  }
  return lines.join("\n");
}

/**
 * Generate Elegen CSV order file.
 */
function generateElegenCsv(oligos: { name: string; sequence: string; scale: string; purification: string }[]): string {
  const lines = ["Oligo Name,Sequence,Length,Scale,Purification"];
  for (const o of oligos) {
    lines.push(`${o.name},${o.sequence},${o.sequence.length},${o.scale},${o.purification}`);
  }
  return lines.join("\n");
}

/**
 * MinKNOW output formatter — generates a file compatible with Oxford Nanopore's
 * MinKNOW sequencing software.
 *
 * MinKNOW expects FASTQ files with specific read naming conventions.
 */
export function exportMinKnowFormat(
  reads: { id: string; sequence: string; quality: Uint8Array }[],
): string {
  const lines: string[] = [];
  const PHRED_OFFSET = 33;

  for (const read of reads) {
    // MinKNOW read ID format: runid=readid read=number ch=channel start_time=time
    const minknowId = `runid=helix_${Date.now()} read=${read.id} ch=1 start_time=2026-08-09T00:00:00Z`;
    lines.push(`@${minknowId}`);
    lines.push(read.sequence);
    lines.push("+");
    lines.push(
      Array.from(read.quality)
        .map((q) => String.fromCharCode(q + PHRED_OFFSET))
        .join(""),
    );
  }

  return lines.join("\n");
}

/**
 * Calculate the cost of synthesizing a DNA archive.
 *
 * @param totalOligos Number of oligos in the archive
 * @param oligoLength Length of each oligo in nucleotides
 * @param vendor Synthesis vendor
 * @param copies Number of physical copies (for redundancy)
 */
export function calculateSynthesisCost(
  totalOligos: number,
  oligoLength: number,
  vendor: "idt" | "twist" | "elegen" = "idt",
  copies: number = 1,
): {
  perCopyCost: number;
  totalCost: number;
  costPerKB: number;
  costPerGB: number;
  vendor: string;
  turnaroundDays: number;
} {
  const pricing = VENDOR_PRICING[vendor];
  const totalBases = totalOligos * oligoLength * copies;
  const baseCost = totalBases * pricing.costPerBase;
  const totalCost = Math.max(baseCost + pricing.setupFee, pricing.minimumOrder);
  const dataBytes = totalOligos * 26; // ~26 bytes payload per oligo
  const costPerKB = totalCost / (dataBytes / 1024);
  const costPerGB = costPerKB * 1024 * 1024;

  return {
    perCopyCost: totalCost / copies,
    totalCost,
    costPerKB,
    costPerGB,
    vendor,
    turnaroundDays: pricing.turnaroundDays,
  };
}

/**
 * Format a synthesis order for CLI display.
 */
export function formatOrderSummary(order: OligoOrder): string {
  return [
    `=== Synthesis Order: ${order.vendor.toUpperCase()} ===`,
    `Oligos:       ${order.totalOligos}`,
    `Total length: ${order.totalLength.toLocaleString()} nt`,
    `Est. cost:    $${order.estimatedCostUSD.toFixed(2)}`,
    `Turnaround:   ${order.estimatedTurnaroundDays} days`,
    `Format:       ${order.orderFormat.toUpperCase()}`,
    ``,
    `Order content preview:`,
    order.orderContent.split("\n").slice(0, 5).join("\n"),
    `... (${order.orderContent.split("\n").length} lines total)`,
  ].join("\n");
}
