// Debug script to investigate homopolymer issue.
import { encodeFile } from "../src/lib/dna/codec";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";
import { bytesToDna, gcContent, maxHomopolymerRun, countHomopolymers, whitenAddress } from "../src/lib/dna/mapping";

async function main() {
  const testPayload = new TextEncoder().encode(
    "Hello, DNA! This is a test of synthetic DNA data storage. ".repeat(20),
  );
  console.log(`Test payload: ${testPayload.length} bytes`);

  const encodeResult = await encodeFile(testPayload, DEFAULT_CONFIG, {
    fileName: "hello.txt",
    contentType: "text/plain",
  });

  for (const oligo of encodeResult.encoded.oligos) {
    const inner = oligo.sequence.slice(
      DEFAULT_CONFIG.primerLength,
      oligo.sequence.length - DEFAULT_CONFIG.primerLength,
    );
    console.log(`\nOligo ${oligo.index} (seed=${oligo.seed}, GC=${oligo.gc.toFixed(2)}, maxHp=${oligo.maxHomopolymer}):`);
    console.log(`  Inner DNA (${inner.length} nt): ${inner}`);
    
    // Find homopolymer runs
    let runStart = 0;
    let runLen = 1;
    for (let i = 1; i <= inner.length; i++) {
      if (i < inner.length && inner[i] === inner[i - 1]) {
        runLen++;
      } else {
        if (runLen >= 3) {
          console.log(`  Homopolymer run: ${inner[runStart]} x ${runLen} at pos ${runStart}-${runStart + runLen - 1}`);
        }
        runStart = i;
        runLen = 1;
      }
    }
  }
}

main().catch(console.error);
