// Test v10.0: Ghost benchmark + Read Until + Hardware API + Chamaeleo bridge.
import { runGhostBenchmark, runComprehensiveBenchmark } from "../src/lib/dna/ghost-benchmark";
import { simulateReadUntil, generateFileHeaders, generateDnaPool, DEFAULT_READUNTIL_CONFIG } from "../src/lib/dna/read-until";
import { createSynthesisOrder, exportMinKnowFormat, calculateSynthesisCost, formatOrderSummary, VENDOR_PRICING } from "../src/lib/dna/hardware-api";

function assert(cond: boolean, msg: string) {
  if (!cond) { console.error("FAIL:", msg); process.exit(1); }
  else console.log("PASS:", msg);
}

async function main() {
  console.log("=== v10.0 Module Tests ===\n");

  // --- Ghost Benchmark ---
  console.log("--- Ghost Benchmark: Erlich 2017 Replay ---\n");
  {
    const result = await runGhostBenchmark(2048);
    console.log(`\n${result.summary}\n`);
    assert(result.results.length > 0, "Ghost benchmark produced results");
    assert(result.markdownTable.length > 0, "Markdown table generated");

    // Save the full benchmark report
    const fullReport = await runComprehensiveBenchmark();
    console.log("\n--- Full Benchmark Report ---\n");
    console.log(fullReport.slice(0, 500) + "...\n");
  }

  // --- Read Until Simulator ---
  console.log("\n--- Read Until / Adaptive Sampling ---\n");
  {
    // Generate a DNA pool with 3 files
    const { strands, fileMap } = generateDnaPool([
      { name: "README.txt", numOligos: 100, oligoLength: 200 },
      { name: "data.bin", numOligos: 200, oligoLength: 200 },
      { name: "image.png", numOligos: 150, oligoLength: 200 },
    ]);

    console.log(`  Pool: ${strands.length} strands, 3 files`);
    assert(strands.length === 450, "Generated 450 strands");

    // Generate file headers
    const files = [
      { name: "README.txt", oligos: [{ sequence: strands[0] }] },
      { name: "data.bin", oligos: [{ sequence: strands[100] }] },
      { name: "image.png", oligos: [{ sequence: strands[300] }] },
    ];
    const headers = generateFileHeaders(files);
    assert(headers.length === 3, "Generated 3 file headers");

    // Simulate Read Until targeting "data.bin"
    const result = simulateReadUntil(strands, fileMap, headers, {
      ...DEFAULT_READUNTIL_CONFIG,
      targetFile: "data.bin",
      totalStrands: strands.length,
      decisionBases: 100, // smaller for faster simulation
      matchThreshold: 100,
    });

    console.log(`  ${result.summary}`);
    console.log(`  Kept: ${result.strandsKept}, Ejected: ${result.strandsEjected}`);
    console.log(`  Bases saved: ${result.basesSaved.toLocaleString()}`);
    console.log(`  Time saved: ${result.timeSavedSeconds.toFixed(1)}s`);
    console.log(`  Cost saved: $${result.costSavedUSD.toFixed(4)}`);

    assert(result.strandsEjected > 0, "Some strands ejected");
    assert(result.basesSaved > 0, "Bases saved > 0");
  }

  // --- Hardware API Stubs ---
  console.log("\n--- Hardware API (Synthesis + MinKNOW) ---\n");
  {
    // Create synthesis order for IDT
    const oligos = [
      { name: "oligo_001", sequence: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" },
      { name: "oligo_002", sequence: "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA" },
      { name: "oligo_003", sequence: "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCC" },
    ];

    const idtOrder = createSynthesisOrder(oligos, "idt", "25nm", "STD");
    console.log(`  IDT order: ${idtOrder.totalOligos} oligos, ${idtOrder.totalLength}nt, $${idtOrder.estimatedCostUSD.toFixed(2)}`);
    assert(idtOrder.vendor === "idt", "IDT vendor set");
    assert(idtOrder.orderFormat === "csv", "IDT uses CSV format");
    assert(idtOrder.estimatedCostUSD > 0, "Cost calculated");

    // Twist order
    const twistOrder = createSynthesisOrder(oligos, "twist", "100nm", "PAGE");
    console.log(`  Twist order: $${twistOrder.estimatedCostUSD.toFixed(2)}, ${twistOrder.estimatedTurnaroundDays} days`);
    assert(twistOrder.orderFormat === "fasta", "Twist uses FASTA format");

    // Cost calculation
    const cost = calculateSynthesisCost(10000, 200, "twist", 4);
    console.log(`  10K oligos × 200nt × 4 copies (Twist): $${cost.totalCost.toFixed(2)}`);
    console.log(`  Cost per GB: $${cost.costPerGB.toFixed(2)}`);
    assert(cost.totalCost > 0, "Total cost > 0");
    assert(cost.costPerGB > 0, "Cost per GB > 0");

    // MinKNOW export
    const reads = [
      { id: "read_001", sequence: "ACGTACGT", quality: new Uint8Array([30, 30, 30, 30, 30, 30, 30, 30]) },
      { id: "read_002", sequence: "TGCATGCA", quality: new Uint8Array([25, 30, 35, 30, 25, 30, 35, 30]) },
    ];
    const minknowOutput = exportMinKnowFormat(reads);
    assert(minknowOutput.includes("@runid="), "MinKNOW format has @runid header");
    assert(minknowOutput.includes("ACGTACGT"), "MinKNOW format has sequence");
    console.log(`  MinKNOW export: ${minknowOutput.split("\n").length} lines`);

    // Format order summary
    console.log(`\n${formatOrderSummary(idtOrder)}`);
  }

  // --- Chamaeleo Bridge ---
  console.log("\n--- Chamaeleo Bridge (Python wrapper) ---\n");
  {
    const fs = await import("fs");
    const pyWrapper = fs.readFileSync("python/helix_codec.py", "utf-8");
    assert(pyWrapper.includes("class HelixCodec"), "Python wrapper has HelixCodec class");
    assert(pyWrapper.includes("class HelixChamaeleoMethod"), "Python wrapper has Chamaeleo method");
    assert(pyWrapper.includes("def transcode"), "Chamaeleo transcode method present");
    assert(pyWrapper.includes("def retrieve"), "Chamaeleo retrieve method present");
    assert(pyWrapper.includes("wasmtime"), "Uses wasmtime for WASM execution");
    console.log("  Python wrapper structure verified");
    console.log("  Classes: HelixCodec, HelixChamaeleoMethod");
    console.log("  Methods: encode, decode, bytes_to_dna, benchmark, transcode, retrieve");
  }

  console.log("\n=== All v10.0 module tests passed ===");
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
