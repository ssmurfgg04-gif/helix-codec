// Test the aligned consensus specifically at moderate indel rates.
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate } from "../src/lib/dna/simulate";
import { DEFAULT_CONFIG } from "../src/lib/dna/types";

async function main() {
  const text = `DNA storage is a promising technology for long-term archival. `;
  const payload = new TextEncoder().encode(text.repeat(20));
  console.log(`Payload: ${payload.length} bytes`);

  const enc = await encodeFile(payload, DEFAULT_CONFIG, {
    fileName: "align-test.txt",
    contentType: "text/plain",
  });
  console.log(`Encoded: ${enc.stats.oligoCount} oligos\n`);

  // Test at various indel rates
  const tests = [
    { name: "0.5% indel (2x)", sub: 0, ins: 0.003, del: 0.002, cov: 10 },
    { name: "1% indel (10x)", sub: 0, ins: 0.005, del: 0.005, cov: 10 },
    { name: "2% indel (15x)", sub: 0, ins: 0.01, del: 0.01, cov: 15 },
    { name: "3% indel (20x)", sub: 0, ins: 0.015, del: 0.015, cov: 20 },
    { name: "5% indel (25x)", sub: 0, ins: 0.025, del: 0.025, cov: 25 },
    { name: "1% sub + 1% indel (15x)", sub: 0.01, ins: 0.005, del: 0.005, cov: 15 },
    { name: "2% sub + 2% indel (20x)", sub: 0.02, ins: 0.01, del: 0.01, cov: 20 },
  ];

  console.log("Profile                          | Standard  | Aligned (soft-info)");
  console.log("---------------------------------|-----------|-------------------");

  for (const t of tests) {
    const sim = simulate(enc.encoded.oligos, {
      substitutionRate: t.sub,
      insertionRate: t.ins,
      deletionRate: t.del,
      coverage: t.cov,
      dropoutRate: 0,
      seed: 42,
    });

    const decStd = await decodeReads(
      sim.reads, enc.encoded.metadata, DEFAULT_CONFIG,
      enc.encoded.forwardPrimer, enc.encoded.reversePrimer, false,
    );

    const decAlign = await decodeReads(
      sim.reads, enc.encoded.metadata, DEFAULT_CONFIG,
      enc.encoded.forwardPrimer, enc.encoded.reversePrimer, true,
    );

    console.log(
      `${t.name.padEnd(33)}| ${decStd.hashMatches ? "PASS" : "FAIL"} (${decStd.stats.oligosRecovered}/${enc.stats.oligoCount})  | ${decAlign.hashMatches ? "PASS" : "FAIL"} (${decAlign.stats.oligosRecovered}/${enc.stats.oligoCount})`,
    );
  }
}

main().catch(console.error);
