import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_ILLUMINA } from "../src/lib/dna/simulate";
import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { dnaToBytes } from "../src/lib/dna/mapping";
import { computeLayout } from "../src/lib/dna/types";

async function main() {
  const payload = new Uint8Array(256 * 1024);
  for (let i = 0; i < payload.length; i++) payload[i] = ((i * 31 + 17) ^ (i >> 8)) & 0xff;
  
  const enc = await encodeFile(payload, ULTIMATE_V55_DENSITY_CONFIG, { fileName: "t.bin", contentType: "application/octet-stream" });
  const cfg = ULTIMATE_V55_DENSITY_CONFIG;
  const layout = computeLayout(cfg);
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
  
  const sim = simulate(enc.encoded.oligos, { ...PRESET_ILLUMINA, coverage: 10, seed: 42 });
  
  // Time just the LDPC decode for all reads
  const fwd = enc.encoded.forwardPrimer;
  const rev = enc.encoded.reversePrimer;
  const expectedDnaLen = layout.totalInnerBytes * 4;
  
  let decodeCount = 0;
  let slowPathCount = 0;
  const t0 = Date.now();
  for (const read of sim.reads) {
    if (read.sequence.length < fwd.length + rev.length) continue;
    const inner = read.sequence.slice(fwd.length, read.sequence.length - rev.length);
    let dna = inner.length === expectedDnaLen ? inner : 
              inner.length > expectedDnaLen ? inner.slice(0, expectedDnaLen) :
              inner + "A".repeat(expectedDnaLen - inner.length);
    let innerBlock;
    try { innerBlock = dnaToBytes(dna); } catch { continue; }
    const rsCodeword = innerBlock.slice(0, innerN);
    try {
      const result = ldpc.decode(rsCodeword);
      decodeCount++;
      if (result.corrected > 0) slowPathCount++;
    } catch { /* decode failed */ }
  }
  const ldpcMs = Date.now() - t0;
  console.log(`LDPC decode: ${ldpcMs}ms, ${decodeCount} decodes, ${slowPathCount} slow path`);
  
  // Now time the full decode
  const t1 = Date.now();
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, fwd, rev);
  const fullMs = Date.now() - t1;
  console.log(`Full decode: ${fullMs}ms, hash=${dec.hashMatches}`);
  console.log(`Overhead (full - LDPC): ${fullMs - ldpcMs}ms`);
}
main().catch(e => { console.error(e); process.exit(1); });
