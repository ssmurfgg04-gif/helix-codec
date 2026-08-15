/**
 * v60: Verify profileHmm3 produces correct posteriors.
 */
import { forwardBackward3, DEFAULT_HMM3_PARAMS } from "../src/lib/dna/profileHmm3";

// Test 1: Perfect match — posteriors should be peaked on the correct base
{
  const ref = "ACGTACGTACGTACGTACGT"; // 20nt
  const read = ref; // perfect match
  const quality = new Uint8Array(20).fill(30); // Q30

  const result = forwardBackward3(read, ref, quality, DEFAULT_HMM3_PARAMS, 10);
  
  console.log("Test 1: Perfect match");
  console.log(`  logLikelihood: ${result.logLikelihood.toFixed(2)}`);
  console.log(`  path length: ${result.path.length}`);
  
  // Check posteriors: each should be peaked on the correct base
  let allCorrect = true;
  for (let j = 0; j < ref.length; j++) {
    const pA = result.matchPosteriors[j * 4 + 0];
    const pC = result.matchPosteriors[j * 4 + 1];
    const pG = result.matchPosteriors[j * 4 + 2];
    const pT = result.matchPosteriors[j * 4 + 3];
    const correctBase = ref[j];
    const correctIdx = "ACGT".indexOf(correctBase);
    const maxIdx = [pA, pC, pG, pT].indexOf(Math.max(pA, pC, pG, pT));
    const correct = maxIdx === correctIdx;
    if (!correct) allCorrect = false;
    if (j < 4) {
      console.log(`  pos ${j}: ref=${correctBase} post=[${pA.toFixed(3)}, ${pC.toFixed(3)}, ${pG.toFixed(3)}, ${pT.toFixed(3)}] ${correct ? "OK" : "WRONG"}`);
    }
  }
  console.log(`  All positions correct: ${allCorrect ? "YES" : "NO"}`);
  console.log();
}

// Test 2: One substitution — posteriors should still be peaked (Q30 overcomes 1 sub)
{
  const ref = "ACGTACGTACGTACGTACGT";
  const read = "ACGTACGAACGTACGTACGT"; // 1 sub at position 7 (T→A)
  const quality = new Uint8Array(20).fill(30);

  const result = forwardBackward3(read, ref, quality, DEFAULT_HMM3_PARAMS, 10);
  
  console.log("Test 2: 1 substitution at position 7");
  console.log(`  logLikelihood: ${result.logLikelihood.toFixed(2)}`);
  
  // Position 7: ref=T, read=A. With Q30, the posterior should still favor T (ref)
  // because P(obs=A | true=T, Q30) = 0.001/3 ≈ 0.0003, while P(obs=A | true=A, Q30) = 0.999.
  // But the HMM combines the emission with the transition. If we trust the read (Q30),
  // the posterior favors A. If we trust the ref, it favors T.
  // The HMM posterior P(true base | read, ref) should favor A (the read's call).
  const pA = result.matchPosteriors[7 * 4 + 0];
  const pC = result.matchPosteriors[7 * 4 + 1];
  const pG = result.matchPosteriors[7 * 4 + 2];
  const pT = result.matchPosteriors[7 * 4 + 3];
  console.log(`  pos 7: ref=T, read=A, post=[A=${pA.toFixed(3)}, C=${pC.toFixed(3)}, G=${pG.toFixed(3)}, T=${pT.toFixed(3)}]`);
  console.log();
}

// Test 3: One insertion — read has an extra base
{
  const ref = "ACGTACGTACGTACGTACGT"; // 20nt
  const read = "ACGTACGXTACGTACGTACGT".replace("X", "G"); // 21nt, extra G at position 7
  const quality = new Uint8Array(21).fill(30);

  const result = forwardBackward3(read, ref, quality, DEFAULT_HMM3_PARAMS, 10);
  
  console.log("Test 3: 1 insertion (extra base at position 7)");
  console.log(`  logLikelihood: ${result.logLikelihood.toFixed(2)}`);
  console.log(`  path length: ${result.path.length}`);
  
  // Count M, I, D states in path
  let mCount = 0, iCount = 0, dCount = 0;
  for (const step of result.path) {
    if (step.state === "M") mCount++;
    else if (step.state === "I") iCount++;
    else dCount++;
  }
  console.log(`  Path: ${mCount} M, ${iCount} I, ${dCount} D`);
  console.log();
}

// Test 4: One deletion — read is missing a base
{
  const ref = "ACGTACGTACGTACGTACGT"; // 20nt
  const read = "ACGTACGTACGTACGTACGT".slice(0, 7) + "ACGTACGTACGTACGT".slice(0, 12); // remove 1 base → 19nt
  const quality = new Uint8Array(19).fill(30);

  const result = forwardBackward3(read, ref, quality, DEFAULT_HMM3_PARAMS, 10);
  
  console.log("Test 4: 1 deletion (missing base at position 7)");
  console.log(`  logLikelihood: ${result.logLikelihood.toFixed(2)}`);
  console.log(`  path length: ${result.path.length}`);
  
  let mCount = 0, iCount = 0, dCount = 0;
  for (const step of result.path) {
    if (step.state === "M") mCount++;
    else if (step.state === "I") iCount++;
    else dCount++;
  }
  console.log(`  Path: ${mCount} M, ${iCount} I, ${dCount} D`);
  console.log();
}
