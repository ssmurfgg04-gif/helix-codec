#!/usr/bin/env bun
/**
 * v53 Single-Oligo Arithmetic Test
 *
 * Encode 1 oligo with arithmetic mode at 500nt, decode it without any noise,
 * and verify byte-for-byte match. This isolates the encoder/decoder mismatch.
 */

import { encodeFile } from "../src/lib/dna/codec";
import { decodeReadsUltra } from "../src/lib/dna/ultra-decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { ULTIMATE_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { computeLayoutAuto } from "../src/lib/dna/types";
import * as fs from "fs";

const TAG = "[v53-single]";

async function main() {
  // Small payload that fits in 1 oligo (98 bytes max with current layout)
  const payloadSize = 50;
  const payload = new Uint8Array(payloadSize);
  for (let i = 0; i < payloadSize; i++) payload[i] = (i * 37 + 11) & 0xff;

  const cfg = JSON.parse(JSON.stringify(ULTIMATE_DENSITY_CONFIG));
  const layout = computeLayoutAuto(cfg);
  console.log(`${TAG} layout: addressBytes=${layout.addressBytes}, payloadBytes=${layout.payloadBytes}, parityBytes=${layout.innerParityBytes}, totalInnerBytes=${layout.totalInnerBytes}`);

  const enc = await encodeFile(payload, cfg, { fileName: "test.bin", contentType: "application/octet-stream" });
  const oligos = enc.encoded.oligos;
  const metadata = enc.encoded.metadata;
  console.log(`${TAG} encoded: ${oligos.length} oligos`);
  console.log(`${TAG} oligo[0] sequence: ${oligos[0].sequence}`);
  console.log(`${TAG} oligo[0] length: ${oligos[0].sequence.length} nt`);

  // Save the encoded oligo for inspection
  fs.writeFileSync("/tmp/v53-oligo.txt", oligos[0].sequence);
  console.log(`${TAG} saved oligo to /tmp/v53-oligo.txt`);

  // Simulate with ZERO noise, 1x coverage (just 1 read of the 1 oligo)
  // Wait — we need at least 1 read per oligo. With 1 oligo and 1x coverage, we get 1 read.
  const simResult = simulate(oligos, { ...PRESET_CLEAN, coverage: 1 });
  const reads = simResult.reads;
  console.log(`${TAG} simulated: ${reads.length} reads`);
  console.log(`${TAG} read[0] sequence: ${reads[0].sequence}`);
  console.log(`${TAG} read[0] length: ${reads[0].sequence.length} nt`);

  // Decode
  const decoded = await decodeReadsUltra(reads, metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`${TAG} decoded: hashMatch=${decoded.hashMatches}, bytes=${decoded.data.length}`);
  console.log(`${TAG} decoded stats:`, JSON.stringify(decoded.stats, null, 2));

  // Compare byte-by-byte
  console.log(`${TAG} payload (first 50): ${Array.from(payload.slice(0, 50)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);
  console.log(`${TAG} decoded  (first 50): ${Array.from(decoded.data.slice(0, 50)).map(b => b.toString(16).padStart(2, '0')).join(' ')}`);

  // Print oligo inner DNA region (between primers)
  const primerLen = cfg.primerLength;
  const innerDna = oligos[0].sequence.slice(primerLen, oligos[0].sequence.length - primerLen);
  console.log(`${TAG} inner DNA (between primers): ${innerDna}`);
  console.log(`${TAG} inner DNA length: ${innerDna.length} nt`);
}

main().catch(e => { console.error(e); process.exit(1); });
