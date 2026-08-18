/**
 * Quick smoke test: encode SARS-CoV-2 with v51-default, verify roundtrip.
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

function parseFasta(content: string): string {
  return content.split("\n")
    .filter(l => !l.startsWith(">") && l.trim())
    .join("");
}

async function main() {
  const fastaPath = path.join(__dirname, "..", "datasets", "small", "sars-cov-2.fa");
  const content = fs.readFileSync(fastaPath, "utf-8");
  const seq = parseFasta(content);
  const data = new Uint8Array(Buffer.from(seq, "utf-8"));
  console.log(`SARS-CoV-2: ${data.length} bytes (${seq.length} nt)`);

  const cfg = V51_DEFAULT_CONFIG;
  console.log(`Encoding with v51-default (300nt, constrained, LDPC 4B, 10% RS)...`);

  const t0 = Date.now();
  const encodeResult = await encodeFile(data, cfg, {
    fileName: "sars-cov-2",
    contentType: "application/octet-stream",
  });
  const encodeMs = Date.now() - t0;
  const { encoded, stats } = { encoded: encodeResult.encoded, stats: encodeResult.stats };

  console.log(`Encoded: ${stats.oligoCount} oligos, density=${stats.netDensityBitsPerNt.toFixed(3)} b/nt, ${encodeMs}ms`);
  console.log(`Screening retries: ${stats.screeningRetries}`);

  // Check constraints
  let violations = 0;
  for (const oligo of encoded.oligos) {
    if (oligo.gc < cfg.constraints.gcMin || oligo.gc > cfg.constraints.gcMax) {
      console.log(`  Oligo ${oligo.index}: GC=${oligo.gc.toFixed(3)} OUT OF RANGE`);
      violations++;
    }
    if (oligo.maxHomopolymer > cfg.constraints.maxHomopolymer) {
      console.log(`  Oligo ${oligo.index}: maxHp=${oligo.maxHomopolymer} OUT OF RANGE`);
      violations++;
    }
  }
  console.log(`Constraint violations: ${violations}`);

  // Simulate clean reads
  console.log(`Simulating 30x clean reads...`);
  const simResult = simulate(encoded.oligos, {
    ...PRESET_CLEAN,
    coverage: 30,
    simulator: "basic",
  });
  console.log(`Generated ${simResult.totalReads} reads`);

  // Decode
  console.log(`Decoding...`);
  const decStart = Date.now();
  const decResult = await decodeReads(
    simResult.reads,
    encoded.metadata,
    cfg,
    encoded.forwardPrimer,
    encoded.reversePrimer,
  );
  const decMs = Date.now() - decStart;
  console.log(`Decoded in ${decMs}ms`);

  // Verify
  const decoded = decResult.data;
  if (!decoded) {
    console.log(`FAIL: decode returned null`);
    process.exit(1);
  }

  let match = decoded.length === data.length;
  if (match) {
    for (let i = 0; i < data.length; i++) {
      if (decoded[i] !== data[i]) { match = false; break; }
    }
  }
  console.log(`Roundtrip: ${match ? "PASS" : "FAIL"} (${decoded.length} vs ${data.length} bytes)`);
  console.log(`Hash match: ${decResult.hashMatches ? "PASS" : "FAIL"}`);

  if (!match || !decResult.hashMatches || violations > 0) {
    process.exit(1);
  }
}

main().catch(e => { console.error(e); process.exit(1); });
