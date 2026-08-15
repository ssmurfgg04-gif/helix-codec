// Test v7.0 modules.
import { WorkerPool, parallelBatch } from "../src/lib/dna/worker-pool";
import { progressiveMSA, msaConsensus, DEFAULT_MSA_CONFIG } from "../src/lib/dna/progressive-msa";
import { LeopardRS } from "../src/lib/dna/leopard-rs";
import { NeuralDecoder, createNeuralPolarDecoder } from "../src/lib/dna/neural-decoder";
import { generateSquiggle, dtwDistance, matchSquiggleToReference } from "../src/lib/dna/squiggle";
import { transformerConsensus, multiReadTransformerConsensus } from "../src/lib/dna/transformer-consensus";
import { generatePQKeyPair, pqSign, pqVerify, signArchivePQ, verifyArchivePQ, getSecurityLevel } from "../src/lib/dna/post-quantum";
import { designGuideRNA, scanForMatches, searchArchive, buildSearchIndex, indexedSearch } from "../src/lib/dna/crispr-search";

function assert(cond: boolean, msg: string) {
  if (!cond) { console.error("FAIL:", msg); process.exit(1); }
  else console.log("PASS:", msg);
}

async function main() {
  console.log("=== v7.0 Module Tests ===\n");

  // --- Worker Pool ---
  console.log("--- Web Worker Pool ---\n");
  {
    const pool = new WorkerPool({ size: 4 });
    const items = [1, 2, 3, 4, 5, 6, 7, 8];
    const results = await pool.map(items, (x) => x * 2);
    assert(JSON.stringify(results) === "[2,4,6,8,10,12,14,16]", "Worker pool map works");
    console.log(`  Pool stats: ${JSON.stringify(pool.getStats())}`);

    // Parallel batch
    const batch = await parallelBatch(items, (x) => x * x);
    assert(JSON.stringify(batch) === "[1,4,9,16,25,36,49,64]", "Parallel batch works");
  }

  // --- Progressive MSA ---
  console.log("\n--- Progressive MSA ---\n");
  {
    const seqs = ["ACGTACGTACGT", "ACGTACGTACGA", "ACGTACGTTCGT"];
    const aligned = progressiveMSA(seqs);
    assert(aligned.length === 3, "MSA produces 3 aligned sequences");
    assert(aligned.every((s) => s.length === aligned[0].length), "All aligned sequences same length");

    const consensus = msaConsensus(aligned);
    assert(consensus.sequence.length > 0, "MSA consensus produced");
    console.log(`  Consensus: ${consensus.sequence} (len ${consensus.sequence.length})`);
  }

  // --- Leopard-RS ---
  console.log("\n--- Leopard-RS (NTT-based) ---\n");
  {
    const rs = new LeopardRS({ n: 20, k: 16 });
    const data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
    const encoded = rs.encode(data);
    assert(encoded.length === 20, "Leopard-RS encodes to correct length");

    // Benchmark
    const bench = LeopardRS.benchmark(256);
    console.log(`  NTT vs naive: ${bench.ntt}ms vs ${bench.naive}ms (speedup: ${bench.speedup.toFixed(2)}x)`);
    assert(bench.ntt >= 0 && bench.naive >= 0, "NTT benchmark produces valid results");
  }

  // --- Neural Decoder ---
  console.log("\n--- Neural Decoder Interface ---\n");
  {
    const decoder = createNeuralPolarDecoder();
    // Should fall back since no model is loaded
    const input = new Float32Array(252).fill(1);
    const result = await decoder.decode(input, (llr) => {
      // Fallback: hard decision
      const out = new Float32Array(208);
      for (let i = 0; i < 208; i++) out[i] = llr[i] < 0 ? 1 : 0;
      return out;
    });
    assert(!result.usedModel, "Neural decoder falls back when no model");
    assert(result.output.length === 208, "Fallback produces correct output length");
    console.log(`  Used model: ${result.usedModel}, inference: ${result.inferenceMs}ms`);
  }

  // --- Squiggle ---
  console.log("\n--- Nanopore Squiggle-Native Decoding ---\n");
  {
    const seq = "ACGTACGTACGTACGTACGT";
    const squiggle = generateSquiggle(seq, undefined, 0.5);
    assert(squiggle.length > 0, "Squiggle generated");
    console.log(`  Generated ${squiggle.length} samples for ${seq.length}nt sequence`);

    // DTW distance
    const sig1 = new Float32Array([1, 2, 3, 4, 5]);
    const sig2 = new Float32Array([1, 2, 3, 4, 5]);
    const dist = dtwDistance(sig1, sig2);
    assert(dist === 0, "DTW distance of identical signals is 0");

    const sig3 = new Float32Array([5, 4, 3, 2, 1]);
    const dist2 = dtwDistance(sig1, sig3);
    assert(dist2 > 0, "DTW distance of different signals > 0");

    // Match squiggle to reference
    const matches = matchSquiggleToReference(squiggle, seq);
    assert(matches.length > 0, "Squiggle matching produces results");
    console.log(`  Found ${matches.length} matches`);
  }

  // --- Transformer Consensus ---
  console.log("\n--- Transformer-Based Consensus ---\n");
  {
    const read = "ACGTACGTACGTACGTACGT";
    const quality = new Uint8Array(read.length).fill(30);
    const result = transformerConsensus(read, quality);
    assert(result.sequence.length === read.length, "Transformer produces same-length output");
    assert(result.confidence.length === read.length, "Confidence array correct length");

    // Multi-read
    const reads = ["ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGA", "ACGTACGTACGTACGTACGT"];
    const multi = multiReadTransformerConsensus(reads);
    assert(multi.sequence.length > 0, "Multi-read consensus produces output");
    console.log(`  Multi-read consensus: ${multi.sequence.slice(0, 20)}...`);
  }

  // --- Post-Quantum Signatures ---
  console.log("\n--- Post-Quantum Signatures (ML-DSA) ---\n");
  {
    const keypair = generatePQKeyPair("ml-dsa-65");
    console.log(`  Algorithm: ${keypair.algorithm}, pubkey size: ${keypair.publicKey.length}`);

    const data = { archive: "test", version: 7 };
    const signature = pqSign(data, keypair.privateKey, keypair.algorithm);
    assert(signature.length > 0, "PQ signature produced");

    const valid = pqVerify(data, signature, keypair.publicKey, keypair.algorithm);
    assert(valid, "PQ signature verifies");

    // Tamper detection
    const tampered = { ...data, archive: "hacked" };
    const tamperedValid = pqVerify(tampered, signature, keypair.publicKey, keypair.algorithm);
    assert(!tamperedValid, "Tampered data fails verification");

    // Security level
    const sec = getSecurityLevel(keypair.algorithm);
    console.log(`  Security: ${sec.classicalBits}-bit classical, ${sec.quantumBits}-bit quantum, quantumSafe=${sec.quantumSafe}`);

    // Sign archive
    const signed = signArchivePQ(data, keypair);
    assert(verifyArchivePQ(signed), "Signed archive verifies");
  }

  // --- CRISPR Search ---
  console.log("\n--- CRISPR Search Simulation ---\n");
  {
    const oligos = [
      { index: 0, sequence: "ACGTACGTAGGACGTACGTAGGACGTACGTAGG" },
      { index: 1, sequence: "TTTTACGTACGTAGGCCCCACGTACGTAGG" },
      { index: 2, sequence: "GGGGACGTACGTAGGAAAACGTACGTAGG" },
    ];

    const gRNA = designGuideRNA("ACGTACGT");
    assert(gRNA.keyword === "ACGTACGT", "Guide RNA targets correct keyword");
    assert(gRNA.pam === "NGG", "PAM is NGG");
    console.log(`  gRNA: ${gRNA.sequence} (target: ${gRNA.keyword})`);

    const result = searchArchive(oligos, "ACGTACGT");
    console.log(`  Found ${result.totalMatches} matches in ${result.matchedOligos.length} oligos`);
    assert(result.totalMatches > 0, "CRISPR search finds matches");
    assert(result.matchedOligos.length === 3, "All oligos have matches");

    // Build index
    const index = buildSearchIndex(oligos, 8);
    assert(index.size > 0, "Search index built");

    // Indexed search
    const fastMatches = indexedSearch(index, "ACGTACGT", 8);
    assert(fastMatches.length > 0, "Indexed search finds matches");
    console.log(`  Indexed search: ${fastMatches.length} matches`);
  }

  console.log("\n=== All v7.0 module tests passed ===");
}

main().catch((e) => { console.error("FATAL:", e); process.exit(1); });
