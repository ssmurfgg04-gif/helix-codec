/**
 * v63-hd quick test — 64KB only, verify density beats SOTA 1.815
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, MutationConfig } from "../src/lib/dna/simulate";
import { ULTIMATE_V63_HD_CONFIG, ULTIMATE_V63_MAXDENSITY_CONFIG, computeDensity } from "../src/lib/dna/presets";

async function main() {
  const configs = [
    { name: "v63-hd (1100nt)", cfg: ULTIMATE_V63_HD_CONFIG, payloadSize: 256 * 1024 },
    { name: "v63-maxdensity (1500nt)", cfg: ULTIMATE_V63_MAXDENSITY_CONFIG, payloadSize: 256 * 1024 },
  ];

  for (const { name, cfg, payloadSize } of configs) {
    const payload = new Uint8Array(payloadSize);
    for (let i = 0; i < payloadSize; i++) payload[i] = Math.floor(Math.random() * 256);

    const theoreticalDensity = computeDensity(cfg);
    console.log(`\n=== ${name} ===`);
    console.log(`  theoretical density: ${theoreticalDensity.toFixed(3)} b/nt`);
    console.log(`  SOTA (Yi Ding 2024): 1.815 b/nt`);

    const t0 = performance.now();
    const enc = await encodeFile(payload, cfg, { fileName: "b.bin", contentType: "application/octet-stream" });
    const encMs = performance.now() - t0;
    const { oligos, metadata } = enc.encoded;

    const totalNt = oligos.reduce((s, o) => s + o.sequence.length, 0);
    const realizedDensity = (payloadSize * 8) / totalNt;
    console.log(`  encoded: ${oligos.length} oligos in ${encMs.toFixed(0)}ms`);
    console.log(`  total nt: ${totalNt.toLocaleString()}`);
    console.log(`  reported density: ${enc.stats.netDensityBitsPerNt.toFixed(3)} b/nt`);
    console.log(`  realized density: ${realizedDensity.toFixed(3)} b/nt`);
    console.log(`  gap vs SOTA: ${((realizedDensity - 1.815) / 1.815 * 100).toFixed(1)}%`);

    const simCfg: MutationConfig = {
      substitutionRate: 0.001,
      insertionRate: 0,
      deletionRate: 0,
      coverage: 10,
      dropoutRate: 0,
      seed: 42,
    };
    const { reads } = simulate(oligos, simCfg);

    const t1 = performance.now();
    const dec = await decodeReads(reads, metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
    const decMs = performance.now() - t1;

    console.log(`  decode: ${decMs.toFixed(0)}ms = ${(payloadSize / 1024 / 1024 / (decMs / 1000)).toFixed(2)} MB/s`);
    console.log(`  hash: ${dec.hashMatches ? "OK ✅" : "FAIL ❌"}`);
    console.log(`  recovered: ${dec.stats.oligosRecovered}/${oligos.length} oligos, erased: ${dec.stats.oligosErased}`);

    if (dec.data) {
      const match = Buffer.compare(Buffer.from(payload), Buffer.from(dec.data)) === 0;
      console.log(`  data: ${match ? "OK ✅" : "FAIL ❌"}`);
    }

    if (dec.hashMatches && realizedDensity > 1.815) {
      console.log(`  🎉 ${name} BEATS SOTA: ${realizedDensity.toFixed(3)} > 1.815 b/nt with hash OK`);
    } else if (dec.hashMatches) {
      console.log(`  ✅ ${name} hash OK, density ${realizedDensity.toFixed(3)} b/nt`);
    } else {
      console.log(`  ❌ ${name} hash FAIL`);
    }
  }
}

main().catch(e => { console.error(e); process.exit(1); });
