/**
 * Main DNA storage codec: encode files to synthetic DNA oligos and decode them back.
 *
 * ENCODING PIPELINE
 *   1. Read input file as bytes.
 *   2. (Optional) DEFLATE-compress via pako.
 *   3. Compute SHA-256 hash for integrity verification.
 *   4. Split into chunks of `payloadBytes` per oligo.
 *   5. Apply outer Reed-Solomon across oligos (parity oligos appended).
 *   6. For each oligo:
 *      a. Build inner block: address(4B) + payload + inner RS parity + CRC-16.
 *      b. Encode inner block to DNA via 2-bit mapping.
 *      c. Screen for GC content (40-60%) and homopolymer (max 3).
 *      d. If screen fails, XOR source bytes with seed and retry (up to maxRetries).
 *      e. Prepend forward primer, append reverse primer.
 *   7. Return EncodedFile with metadata + oligos.
 *
 * DECODING PIPELINE (see decode.ts for the recovery engine)
 *   1. Strip primers from each received read.
 *   2. Decode DNA -> bytes.
 *   3. Extract address (index + seed) from each oligo.
 *   4. Verify CRC-16; flag failures as erasures.
 *   5. Apply inner RS per oligo (corrects small substitution errors).
 *   6. Cluster oligos by index; take majority consensus if multiple copies.
 *   7. Apply outer RS across oligos (corrects missing strands via erasures).
 *   8. Concatenate payloads in index order; trim to fileSize.
 *   9. (Optional) DEFLATE-decompress.
 *   10. Verify SHA-256 hash matches metadata.
 */

import { deflate } from "pako";
import { deriveAddress, deriveHierarchicalAddress, deriveArchiveSalt, type AddressingConfig, type HierarchicalAddress } from './addressing';
import { ReedSolomon } from "./reedsolomon";
import { ReedSolomon216 } from "./reedsolomon216";
import { LDPCInnerCode } from "./ldpc-codec";
import { ConvolutionalInnerCode, DEFAULT_CONV_CONFIG } from "./convolutional";
import { NASA_K9_CONFIG } from "./convolutional-k9";
import { crc16, crc16Bytes, verifyCrc16 } from "./crc16";
import { bytesToDna, dnaToBytes, satisfiesConstraints, xorWithSeed, gcContent, maxHomopolymerRun, whitenAddress, Base } from "./mapping";
import { bytesToGoldmanDna, goldmanDnaToBytes } from "./goldman";
import { bytesToConstrainedDna, constrainedDnaToBytes, bytesToSplitConstrainedDna, splitConstrainedDnaToBytesWithErasure } from "./constrained-mapping";
import { bytesToSrtDna } from "./srt-constrained";
import { bytesToArithmeticDna, bytesToArithmeticDnaBlocked, bytesToArithmeticDnaCrc } from "./markov-arithmetic";
import { interleaveCodewords, deinterleaveCodewords } from "./interleaving";
import {
  CodecConfig,
  CodecMetadata,
  EncodedFile,
  Oligo,
  OligoLayout,
  computeLayout,
  computeLayoutAuto,
  DEFAULT_CONFIG,
} from "./types";

// --- Primer design ---

/**
 * Generate a balanced random primer of given length.
 * Uses a deterministic PRNG so the same config produces the same primers.
 */
function generatePrimer(length: number, seed: number): string {
  const bases = "ACGT";
  let state = seed >>> 0;
  let primer = "";
  let prev = "";
  // Aim for ~50% GC and no homopolymers.
  while (primer.length < length) {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    state = state >>> 0;
    const base = bases[state % 4];
    if (base !== prev) {
      primer += base;
      prev = base;
    }
  }
  return primer;
}

/** Get or generate primers from config. */
function getPrimers(cfg: CodecConfig): { fwd: string; rev: string } {
  const fwd = cfg.forwardPrimer ?? generatePrimer(cfg.primerLength, 0xC0FFEE);
  const rev = cfg.reversePrimer ?? generatePrimer(cfg.primerLength, 0xDEADBEEF);
  if (fwd.length !== cfg.primerLength || rev.length !== cfg.primerLength) {
    throw new Error(
      `Primer length ${fwd.length}/${rev.length} does not match config ${cfg.primerLength}`,
    );
  }
  return { fwd, rev };
}

// --- SHA-256 (using Node crypto with chunked hashing for large data) ---

async function sha256(data: Uint8Array): Promise<string> {
  // Always use Node's createHash which handles large data via streaming
  const { createHash } = await import("crypto");
  const h = createHash("sha256");
  // Hash in 64MB chunks to avoid buffer limits
  for (let i = 0; i < data.length; i += 64 * 1024 * 1024) {
    h.update(data.subarray(i, Math.min(i + 64 * 1024 * 1024, data.length)));
  }
  return h.digest("hex");
}

// --- Encode ---

export interface EncodeResult {
  encoded: EncodedFile;
  stats: {
    rawSize: number;
    compressedSize: number;
    oligoCount: number;
    payloadBytesPerOligo: number;
    netDensityBitsPerNt: number;
    overheadPercent: number;
    screeningRetries: number;
    encodeTimeMs: number;
  };
}

/**
 * Encode a file (as bytes) into synthetic DNA oligos.
 */
export async function encodeFile(
  data: Uint8Array,
  cfg: CodecConfig,
  meta: { fileName: string; contentType: string },
): Promise<EncodeResult> {
  const t0 = Date.now();
  const layout = computeLayoutAuto(cfg);
  const { fwd, rev } = getPrimers(cfg);
  const useConvInner = !!cfg.useConvolutionalInner;

  // 1) Compress (optional) — before encryption since encrypted data is incompressible
  let compressed: Uint8Array = data;
  let compression: "none" | "deflate" = "none";
  if (cfg.compress) {
    compressed = deflate(data, { level: 9 });
    compression = "deflate";
  }

  // 2) Encrypt (optional) — after compression since encrypted data is incompressible
  let processedData: Uint8Array = compressed;
  let encryptionSalt: Uint8Array | null = null;
  if (cfg.encryptPassword) {
    const { encrypt } = await import("./encryption");
    const salt = new Uint8Array(16);
    const { randomBytes } = await import("crypto");
    salt.set(randomBytes(16));
    const result = encrypt(compressed, { password: cfg.encryptPassword, salt });
    processedData = result.ciphertext;
    encryptionSalt = salt;
  }

  // 3) Hash the ORIGINAL data (for integrity verification at decode time)
  const fileHash = await sha256(data);

  // 3) Split into payload chunks.
  //    Last chunk may be padded with zeros; we record fileSize to trim later.
  const chunkSize = layout.payloadBytes;
  const dataOligoCount = Math.max(1, Math.ceil(compressed.length / chunkSize));
  // Pad compressed data to a multiple of chunkSize
  const paddedLen = dataOligoCount * chunkSize;
  const padded = new Uint8Array(paddedLen);
  padded.set(compressed, 0);
  // Record original compressed length so we can trim padding after decode
  // (we'll use fileSize from metadata + decompression to handle this)

  // 4) Outer RS: across oligos.
  //    Use GF(2^8) RS when oligo count <= 255, GF(2^16) RS when > 255.
  //    GF(2^16) treats pairs of bytes as single 16-bit symbols, supporting
  //    up to 65535 oligos per RS block.
  //
  //    For data requiring >65535 oligos, the CALLER must split into multiple
  //    shards and call encodeFile() for each shard independently.
  const parityCount = Math.max(2, Math.ceil(dataOligoCount * cfg.outerParityRatio));
  if (dataOligoCount + parityCount > 65535) {
    throw new Error(
      `Data too large for single RS block: ${dataOligoCount + parityCount} oligos > 65535 limit. ` +
      `Split the data into smaller shards (e.g., 64MB each) and call encodeFile() for each shard.`
    );
  }
  const useGF216 = dataOligoCount + parityCount > 255;

  let outerN: number;
  let outerK: number;
  let totalOligoCount: number;

  if (useGF216) {
    outerN = dataOligoCount + parityCount;
    outerK = dataOligoCount;
    totalOligoCount = outerN;
  } else {
    outerN = dataOligoCount + parityCount;
    outerK = dataOligoCount;
    totalOligoCount = outerN;
  }
  if (outerN <= outerK) {
    throw new Error(
      `Outer RS config invalid: n=${outerN} <= k=${outerK}. Increase outerParityRatio.`,
    );
  }

  // Create the appropriate RS encoder
  const outerRs8 = !useGF216 ? new ReedSolomon({ n: outerN, k: outerK }) : null;
  const outerRs216 = useGF216 ? new ReedSolomon216({ n: outerN, k: outerK }) : null;

  // 5) Outer RS parity computation
  const parityBytes = new Uint8Array(parityCount * chunkSize);
  let screeningRetries = 0;

  if (useGF216 && outerRs216) {
    // GF(2^16): process byte pairs
    // Pre-allocate reusable buffer to avoid GC pressure
    const dataSymbols16 = new Uint16Array(dataOligoCount);
    const numPairs = Math.floor(chunkSize / 2);
    for (let pairIdx = 0; pairIdx < numPairs; pairIdx++) {
      const j0 = pairIdx * 2;
      const j1 = pairIdx * 2 + 1;
      // Build 16-bit data symbols: symbol[i] = (byte[j0] << 8) | byte[j1]
      for (let i = 0; i < dataOligoCount; i++) {
        dataSymbols16[i] = (padded[i * chunkSize + j0] << 8) | padded[i * chunkSize + j1];
      }
      const parity = outerRs216.parity(dataSymbols16);
      // Unpack 16-bit parity back to bytes
      for (let i = 0; i < parity.length; i++) {
        parityBytes[i * chunkSize + j0] = (parity[i] >> 8) & 0xff;
        parityBytes[i * chunkSize + j1] = parity[i] & 0xff;
      }
    }
    // Handle odd byte if chunkSize is odd
    if (chunkSize % 2 === 1) {
      const j = chunkSize - 1;
      if (dataOligoCount + parityCount <= 255) {
        // Use GF(2^8) for the last odd byte (within 255 limit)
        // Pre-allocate reusable buffer to avoid GC pressure
        const dataSymbols8 = new Uint8Array(dataOligoCount);
        for (let i = 0; i < dataOligoCount; i++) {
          dataSymbols8[i] = padded[i * chunkSize + j];
        }
        const rs8 = new ReedSolomon({ n: outerN, k: outerK });
        const parity = rs8.parity(dataSymbols8);
        for (let i = 0; i < parity.length; i++) {
          parityBytes[i * chunkSize + j] = parity[i];
        }
      } else {
        // For n > 255, use GF(2^16) with the odd byte as the low byte (high=0)
        // Reuse the pre-allocated Uint16Array buffer
        for (let i = 0; i < dataOligoCount; i++) {
          dataSymbols16[i] = padded[i * chunkSize + j] & 0xff; // high byte = 0
        }
        const parity = outerRs216.parity(dataSymbols16);
        for (let i = 0; i < parity.length; i++) {
          parityBytes[i * chunkSize + j] = parity[i] & 0xff; // store low byte only
        }
      }
    }
  } else if (outerRs8) {
    // GF(2^8): process each byte position independently
    // Pre-allocate reusable buffer to avoid GC pressure
    const dataSymbols = new Uint8Array(dataOligoCount);
    for (let j = 0; j < chunkSize; j++) {
      // Reuse dataSymbols buffer
      for (let i = 0; i < dataOligoCount; i++) {
        dataSymbols[i] = padded[i * chunkSize + j];
      }
      const parity = outerRs8.parity(dataSymbols);
      for (let i = 0; i < parity.length; i++) {
        parityBytes[i * chunkSize + j] = parity[i];
      }
    }
  }

  // 6) Build each oligo's inner block and encode to DNA
  const useLDPC = (cfg.innerCode ?? "rs") === "ldpc";
  // v63: Use GF(2^16) inner RS when n > 255 (enables 1000nt+ oligos).
  const innerRsN = layout.addressBytes + layout.payloadBytes + layout.innerParityBytes;
  const innerRs = innerRsN > 255
    ? new ReedSolomon216({
        n: innerRsN,
        k: layout.addressBytes + layout.payloadBytes,
      })
    : new ReedSolomon({
        n: innerRsN,
        k: layout.addressBytes + layout.payloadBytes,
      });
  // Wait: inner RS should be over (address + payload) -> parity. CRC is over (address + payload + parity).
  // Let's re-derive:
  //   inner block = [address(4B), payload(chunkSize B), innerParity(8B), crc(2B)]
  //   inner RS encodes (address + payload) -> innerParity
  //   CRC is computed over (address + payload + innerParity) and appended.
  //   Total inner bytes = 4 + chunkSize + 8 + 2 = 14 + chunkSize.
  //
  //   With oligoLength=200, primerLen=20: innerNt = 160, innerBytes = 40.
  //   So chunkSize = 40 - 4 - 8 - 2 = 26 bytes per oligo.
  //
  //   inner RS: n = 4 + 26 + 8 = 38, k = 4 + 26 = 30. Good (n <= 255).

  // Inner RS k and n
  // v62: For arithmetic-v2, the LDPC codeword does NOT include the address.
  //   Normal mode:    innerK = addressBytes + payloadBytes, innerN = innerK + parity
  //   Arithmetic-v2:  innerK = payloadBytes (NO address),   innerN = innerK + parity
  const useArithmeticV2 = (cfg.mappingMode ?? "direct") === "arithmetic";
  const innerK = useArithmeticV2
    ? layout.payloadBytes
    : layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  // The CRC covers the entire inner block (address + payload + parity), but is
  // not part of the RS codeword. So we have:
  //   RS codeword = address(4B) + payload(chunkSize B) + parity(8B)  -> length innerN
  //   Inner block (DNA-encoded) = RS codeword + CRC(2B)              -> length innerN + 2 = totalInnerBytes

  // v63: Use GF(2^16) inner RS when innerN > 255 (enables 1000nt+ oligos).
  // LDPC is still the primary inner code; RS is only a fallback.
  const innerRsReal = innerN > 255
    ? new ReedSolomon216({ n: innerN, k: innerK })
    : new ReedSolomon({ n: innerN, k: innerK });

  // LDPC inner code (if configured). Same byte-oriented interface as ReedSolomon.
  // LDPC operates on bits internally: n_bits = innerN * 8, k_bits = innerK * 8.
  // Note: WASM LDPC encoder is slower than JS for encode due to boundary overhead
  // (3.5µs vs 2.7µs per op). WASM is only faster for batch decode.
  const innerLdpcReal = useLDPC ? new LDPCInnerCode({ n: innerN, k: innerK }) : null;

  // v52: Convolutional inner code (HEDGES-style). When useConvolutionalInner
  // is enabled, the LDPC codeword (innerN bytes) is conv-encoded at rate 1/2
  // before DNA mapping. The conv-encoded region occupies `convEncodedBytes`
  // between the address and CRC.
  //
  // v61: For nanopore channel, use NASA K=9 (memory=8, d_free=24) instead of
  // the default K=3 (memory=2, d_free=5). The encoder/decoder pair must use
  // matching configs. The K=9 code's 5× stronger correction capability is
  // what enables the indel-tolerant Viterbi to distinguish insertions from
  // substitutions at 9% IDS.
  const convCfg = (cfg.channel === "nanopore")
    ? NASA_K9_CONFIG
    : DEFAULT_CONV_CONFIG;
  const convInner = useConvInner ? new ConvolutionalInnerCode(innerN, convCfg) : null;

  const oligos: Oligo[] = [];

  for (let oligoIdx = 0; oligoIdx < totalOligoCount; oligoIdx++) {
    // Get payload bytes for this oligo
    const payload = new Uint8Array(chunkSize);
    if (oligoIdx < dataOligoCount) {
      payload.set(padded.slice(oligoIdx * chunkSize, (oligoIdx + 1) * chunkSize), 0);
    } else if (oligoIdx < dataOligoCount + parityCount) {
      // Parity oligo
      const parityIdx = oligoIdx - dataOligoCount;
      payload.set(
        parityBytes.slice(parityIdx * chunkSize, (parityIdx + 1) * chunkSize),
        0,
      );
    }

    // Raw address: 3 bytes index + 1 byte seed (seed set during retry)
    const rawAddress = new Uint8Array(4);
    rawAddress[0] = (oligoIdx >> 16) & 0xff;
    rawAddress[1] = (oligoIdx >> 8) & 0xff;
    rawAddress[2] = oligoIdx & 0xff;
    rawAddress[3] = 0; // seed placeholder

    // Try encoding with seed = 0, 1, 2, ... until constraints are satisfied.
    const constraints = {
      gcMin: cfg.constraints.gcMin,
      gcMax: cfg.constraints.gcMax,
      maxHomopolymer: cfg.constraints.maxHomopolymer,
    };
    const useGoldman = (cfg.mappingMode ?? "direct") === "goldman";
    const useConstrained = (cfg.mappingMode ?? "direct") === "constrained";
    const useSrt = (cfg.mappingMode ?? "direct") === "srt";
    const useArithmetic = (cfg.mappingMode ?? "direct") === "arithmetic";
    let seed = 0;
    let dna = "";
    let attempts = 0;
    let bestDna = "";
    let bestSeed = 0;
    let bestSatisfied = false;

    if (useConvInner) {
      // v52 HEDGES-style: LDPC codeword → conv-encoded → DNA.
      // Layout: [address(4) + conv_encoded(LDPC codeword) + CRC(2)]
      // The address is direct-mapped at the start (for clustering).
      // The LDPC codeword (which also contains the address + payload + parity)
      // is conv-encoded at rate 1/2 to provide indel tolerance.
      // CRC is over the original LDPC codeword (not the conv-encoded bytes).
      const address = rawAddress.slice();
      address[3] = 0;
      const whitenedAddress = whitenAddress(address);
      const rsData = new Uint8Array(innerK);
      rsData.set(whitenedAddress, 0);
      rsData.set(payload, layout.addressBytes);
      const rsCodeword = useLDPC && innerLdpcReal
        ? innerLdpcReal.encode(rsData)
        : innerRsReal.encode(rsData);
      // CRC over the original LDPC codeword (pre-conv)
      const crc = crc16Bytes(rsCodeword);
      // Conv-encode the LDPC codeword
      if (!convInner) throw new Error("convInner not initialized");
      let convEncoded = convInner.encode(rsCodeword);
      // For nanopore channel, insert periodic CRC-8 sync markers for resynchronization
      // This allows the decoder to resync after burst indels instead of losing the whole oligo
      let convEncodedWithMarkers: Uint8Array | null = null;
      if (cfg.channel === 'nanopore') {
        const { insertCRCMarkers, DEFAULT_CRC_MARKER_CONFIG } = await import('./crcmarker');
        // Insert markers every ~30nt in the conv-encoded region
        // 30nt at 4nt/byte = 7.5 bytes; use segmentSize=8 (~32nt spacing)
        const markerCfg = { ...DEFAULT_CRC_MARKER_CONFIG, segmentSize: 8 };
        convEncodedWithMarkers = insertCRCMarkers(convEncoded, markerCfg);
        convEncoded = convEncodedWithMarkers;
      }
      // Assemble inner block: address + conv-encoded + CRC
      // (pad to totalInnerBytes if there's any leftover space)
      const innerBlock = new Uint8Array(layout.totalInnerBytes);
      innerBlock.set(whitenedAddress, 0);
      innerBlock.set(convEncoded, layout.addressBytes);
      innerBlock.set(crc, layout.addressBytes + layout.convEncodedBytes);
      // Direct 2-bit mapping (conv inner only supports direct mapping in v52)
      dna = bytesToDna(innerBlock);
      bestDna = dna;
      bestSeed = 0;
      bestSatisfied = true;
    } else if (useGoldman) {
      // Goldman rotational mapping: GUARANTEES no homopolymers (max run = 1).
      // No screening needed — every oligo is valid. GC content is ~50% on average.
      // Seed is always 0 (no XOR re-encoding needed).
      const address = rawAddress.slice();
      address[3] = 0; // seed = 0
      const whitenedAddress = whitenAddress(address);
      const effectivePayload = payload; // no XOR (seed = 0)
      const rsData = new Uint8Array(innerK);
      rsData.set(whitenedAddress, 0);
      rsData.set(effectivePayload, layout.addressBytes);
      const rsCodeword = useLDPC && innerLdpcReal
        ? innerLdpcReal.encode(rsData)
        : innerRsReal.encode(rsData);
      const crc = crc16Bytes(rsCodeword);
      const innerBlock = new Uint8Array(totalNtBytes(layout));
      innerBlock.set(rsCodeword, 0);
      innerBlock.set(crc, rsCodeword.length);
      // Goldman mapping: bytes → trits → DNA (homopolymer-free)
      const goldmanMode = cfg.goldmanMode ?? "fast";
      dna = bytesToGoldmanDna(innerBlock, "A", goldmanMode);
      bestDna = dna;
      bestSeed = 0;
      bestSatisfied = true;
      // Goldman guarantees maxHomopolymer = 1, but GC may still be out of range.
      // For high-entropy data, GC is ~50% on average. If out of range, we could
      // retry with a different startBase, but for simplicity we accept it.
      // (GC violations are much less problematic than homopolymers for synthesis.)
    } else if (useConstrained) {
      // Split constrained mapping: direct for address (4B), constrained for rest.
      // Address uses direct mapping (no erasures → reliable clustering).
      // Payload+CRC+parity use constrained (homopolymer-free, 1.1% erasures).
      const address = rawAddress.slice();
      address[3] = 0; // seed = 0
      const whitenedAddress = whitenAddress(address);
      const effectivePayload = payload;
      const rsData = new Uint8Array(innerK);
      rsData.set(whitenedAddress, 0);
      rsData.set(effectivePayload, layout.addressBytes);
      const rsCodeword = useLDPC && innerLdpcReal
        ? innerLdpcReal.encode(rsData)
        : innerRsReal.encode(rsData);
      const crc = crc16Bytes(rsCodeword);
      const innerBlock = new Uint8Array(totalNtBytes(layout));
      innerBlock.set(rsCodeword, 0);
      innerBlock.set(crc, rsCodeword.length);
      // Split: address bytes use direct, rest uses constrained
      dna = bytesToSplitConstrainedDna(innerBlock, cfg.constraints.maxHomopolymer, layout.addressBytes);
      bestDna = dna;
      bestSeed = 0;
      bestSatisfied = true;
    } else if (useSrt) {
      // SRT constrained mapping: 2.0 bits/nt, homopolymer ≤ 3 guaranteed.
      // Breaks homopolymers by injecting 1-bit errors (LDPC corrects them).
      // Zero screening retries. Fully reversible via LDPC.
      const address = rawAddress.slice();
      address[3] = 0;
      const whitenedAddress = whitenAddress(address);
      const effectivePayload = payload;
      const rsData = new Uint8Array(innerK);
      rsData.set(whitenedAddress, 0);
      rsData.set(effectivePayload, layout.addressBytes);
      const rsCodeword = useLDPC && innerLdpcReal
        ? innerLdpcReal.encode(rsData)
        : innerRsReal.encode(rsData);
      const crc = crc16Bytes(rsCodeword);
      const innerBlock = new Uint8Array(totalNtBytes(layout));
      innerBlock.set(rsCodeword, 0);
      innerBlock.set(crc, rsCodeword.length);
      // SRT: 2-bit mapping with homopolymer breaking (injects errors, LDPC corrects)
      dna = bytesToSrtDna(innerBlock, cfg.constraints.maxHomopolymer, totalNtBytes(layout) * 4);
      bestDna = dna;
      bestSeed = 0;
      bestSatisfied = true;
    } else if (useArithmetic) {
      // v62: Arithmetic-v2 mode (address OUTSIDE the arithmetic stream).
      //
      // This fixes the v57-v61 "arithmetic mode still broken" issue where the
      // address was inside the arithmetic stream, causing:
      //   1. Address corruption from arithmetic termination (last bytes unreliable)
      //   2. Address corruption from IDS (indels shift arithmetic interval state)
      //
      // v62 layout (between primers):
      //   [Address (16 nt direct DNA)] [Arithmetic stream (payload+parity)]
      //
      // The LDPC codeword = payload + parity (NO address, NO CRC-16).
      // Per-block CRC-8 sync markers inside the arithmetic stream replace CRC-16.
      // The address is direct DNA → robust to indels via k-mer clustering.
      const address = rawAddress.slice();
      address[3] = 0; // seed = 0
      const whitenedAddress = whitenAddress(address);

      // Build LDPC codeword = payload + parity (NO address in LDPC)
      // For arithmetic-v2, the LDPC instance has k = payloadBytes (set up above)
      const ldpcCodeword = useLDPC && innerLdpcReal
        ? innerLdpcReal.encode(payload)
        : innerRsReal.encode(payload);

      // Split LDPC codeword into payload + parity for arithmetic encoding
      const payloadPart = ldpcCodeword.slice(0, layout.payloadBytes);
      const parityPart = ldpcCodeword.slice(layout.payloadBytes);

      // Encode address as direct DNA (16 nt)
      const addressDna = bytesToDna(whitenedAddress);

      // Arithmetic-encode the LDPC codeword (payload + parity) with per-block CRC-8
      const blockSize = cfg.arithmeticBlockSize ?? 80;
      const arithData = new Uint8Array(layout.payloadBytes + layout.innerParityBytes);
      arithData.set(payloadPart, 0);
      arithData.set(parityPart, layout.payloadBytes);
      // v62: targetLen is the DNA length (in nt), NOT byte count.
      // The arithmetic stream occupies: totalInnerNt - addressNt(16) nt.
      // Pass undefined to let the encoder auto-compute from data length.
      const arithmeticDna = bytesToArithmeticDnaCrc(
        arithData,
        cfg.constraints.maxHomopolymer,
        undefined, // auto-compute DNA length from data
        blockSize,
      );

      // Concatenate: address (direct) + arithmetic stream
      dna = addressDna + arithmeticDna;

      // Pad to expected DNA length (addressNt + arithmeticNt)
      const expectedArithDnaLen = cfg.oligoLength - 2 * cfg.primerLength;
      if (dna.length < expectedArithDnaLen) {
        const lastBase = dna[dna.length - 1] as Base;
        const padBase: Base = lastBase === "A" ? "C" : "A";
        dna += padBase.repeat(expectedArithDnaLen - dna.length);
      }

      bestDna = dna;
      bestSeed = 0;
      bestSatisfied = true;
    } else {
      // Direct 2-bit mapping with constraint screening + seed retries.
      // Optimized: precompute the non-changing parts, only re-encode the changed parts.
      const baseAddress = rawAddress.slice();
      const baseRsData = new Uint8Array(innerK);
      // Pre-set the address (will be overwritten with seed)
      baseRsData.set(whitenAddress(baseAddress), 0);
      baseRsData.set(payload, layout.addressBytes);

      while (attempts <= cfg.maxRetries) {
        // Only change the seed byte in the address
        baseAddress[3] = seed & 0xff;
        const whitenedAddress = whitenAddress(baseAddress);
        baseRsData.set(whitenedAddress, 0);
        // Only XOR payload if seed != 0
        const effectivePayload = seed === 0 ? payload : xorWithSeed(payload, seed);
        baseRsData.set(effectivePayload, layout.addressBytes);

        const rsCodeword = useLDPC && innerLdpcReal
          ? innerLdpcReal.encode(baseRsData)
          : innerRsReal.encode(baseRsData);

        const crc = crc16Bytes(rsCodeword);
        const innerBlock = new Uint8Array(totalNtBytes(layout));
        innerBlock.set(rsCodeword, 0);
        innerBlock.set(crc, rsCodeword.length);

        dna = bytesToDna(innerBlock);
        const ok = satisfiesConstraints(dna, constraints);

        if (ok) {
          bestDna = dna;
          bestSeed = seed;
          bestSatisfied = true;
          break;
        }
        if (!bestDna) {
          bestDna = dna;
          bestSeed = seed;
        }
        attempts++;
        seed = attempts;
      }
    }

    screeningRetries += attempts;
    dna = bestDna;
    seed = bestSeed;

    // Compute stats on the payload region (excluding primers)
    const gc = gcContent(dna);
    const maxHp = maxHomopolymerRun(dna);

    // Pad DNA to exactly the expected inner nt length.
    // v62: For arithmetic-v2, the DNA length = oligoLength - 2*primerLength
    // (address 16nt + arithmetic stream). For other modes, it's totalInnerBytes * 4.
    const expectedDnaLen = useArithmeticV2
      ? (cfg.oligoLength - 2 * cfg.primerLength)
      : layout.totalInnerBytes * 4;
    if (dna.length < expectedDnaLen) {
      // v62: For arithmetic-v2, use alternating bases to avoid homopolymers
      if (useArithmeticV2) {
        const lastBase = dna[dna.length - 1] || "A";
        const padBase: string = lastBase === "A" ? "C" : "A";
        dna += padBase.repeat(expectedDnaLen - dna.length);
      } else {
        dna = dna + "A".repeat(expectedDnaLen - dna.length);
      }
    } else if (dna.length > expectedDnaLen) {
      dna = dna.slice(0, expectedDnaLen);
    }

    const sequence = fwd + dna + rev;
    oligos.push({
      index: oligoIdx,
      sequence,
      gc,
      maxHomopolymer: maxHp,
      seed,
      payloadBytes: chunkSize,
      length: sequence.length,
    });
  }

  // 6b) Interleaving post-processing (Kim 2024)
  // If interleaveDepth > 1, interleave the PAYLOAD+PARITY bytes across groups
  // of oligos. The ADDRESS bytes (first 4) are NOT interleaved — they must
  // remain at fixed positions for clustering to work.
  //
  // This spreads burst errors: a burst in one oligo becomes 1 error per codeword
  // across `depth` oligos, each easily corrected by LDPC.
  const interleaveDepth = cfg.interleaveDepth ?? 0;
  if (interleaveDepth > 1 && oligos.length >= interleaveDepth) {
    const primerLen = fwd.length;
    const innerDnaLen = layout.totalInnerBytes * 4;
    const addressBytes = layout.addressBytes; // 4 bytes — NOT interleaved
    const interleaveRegionLen = layout.totalInnerBytes - addressBytes;

    for (let g = 0; g < oligos.length; g += interleaveDepth) {
      const groupSize = Math.min(interleaveDepth, oligos.length - g);
      if (groupSize < 2) continue;

      // Extract the interleavable region (bytes after address) from each oligo
      const regions: Uint8Array[] = [];
      for (let i = 0; i < groupSize; i++) {
        const oligo = oligos[g + i];
        // Pad the oligo sequence to exactly oligoLength (some oligos may be shorter
        // due to encoding quirks — padding with 'A' is safe)
        let seq = oligo.sequence;
        if (seq.length < cfg.oligoLength) {
          seq = seq + "A".repeat(cfg.oligoLength - seq.length);
        }
        const innerDna = seq.slice(primerLen, primerLen + innerDnaLen);
        const fullBlock = dnaToBytes(innerDna);
        // Skip the first `addressBytes` bytes, keep the rest
        regions.push(fullBlock.slice(addressBytes));
      }

      // Interleave the regions
      const interleaved = interleaveCodewords(regions);
      const blockLen = regions[0].length;

      // Assign interleaved bytes back to oligos
      for (let i = 0; i < groupSize; i++) {
        const oligo = oligos[g + i];
        // Pad the sequence to exactly oligoLength
        let seq = oligo.sequence;
        if (seq.length < cfg.oligoLength) {
          seq = seq + "A".repeat(cfg.oligoLength - seq.length);
        }
        const innerDna = seq.slice(primerLen, primerLen + innerDnaLen);
        const fullBlock = dnaToBytes(innerDna);
        // Replace the interleavable region with interleaved bytes
        const interleavedChunk = interleaved.slice(i * blockLen, (i + 1) * blockLen);
        fullBlock.set(interleavedChunk, addressBytes);
        const newInnerDna = bytesToDna(fullBlock);
        // Reconstruct the sequence with primers (preserve original length)
        const oldSuffix = oligo.sequence.slice(primerLen + innerDnaLen);
        const prefix = oligo.sequence.slice(0, primerLen);
        oligo.sequence = prefix + newInnerDna + oldSuffix;
      }
    }
  }

  // 7) Build metadata
  const metadata: CodecMetadata = {
    fileName: meta.fileName,
    fileSize: data.length,
    fileHash,
    contentType: meta.contentType,
    compression,
    rawSize: data.length,
    oligoCount: totalOligoCount,
    payloadBytesPerOligo: chunkSize,
    innerRS: { n: innerN, k: innerK },
    innerCode: useLDPC ? "ldpc" : "rs",
    ldpcDecoder: cfg.ldpcDecoder as "hard" | "osd" | "bp" | "auto" | undefined,
    mappingMode: (cfg.mappingMode ?? "direct") as "direct" | "goldman" | "constrained" | "srt" | "arithmetic",
    goldmanMode: (cfg.goldmanMode ?? "fast") as "fast" | "dense",
    outerRS: { n: outerN, k: outerK },
    parityOligos: parityCount,
    interleaveDepth: interleaveDepth,
    encryptionSalt: encryptionSalt ? Buffer.from(encryptionSalt).toString("base64") : undefined,
    channel: (cfg.channel ?? "illumina") as "illumina" | "nanopore",
    lowCoverageTrigger: cfg.lowCoverageTrigger ?? 5,
    useConvolutionalInner: useConvInner,
    version: 1,
    encodedAt: new Date().toISOString(),
  };

  const encoded: EncodedFile = {
    metadata,
    oligos,
    forwardPrimer: fwd,
    reversePrimer: rev,
  };

  const encodeTimeMs = Date.now() - t0;
  // Net density: (compressed_size * 8 bits) / (total_oligos * oligo_length nt)
  const totalNt = totalOligoCount * cfg.oligoLength;
  const netDensityBitsPerNt = (compressed.length * 8) / totalNt;
  const overheadPercent =
    ((totalNt - compressed.length * 4) / (totalNt)) * 100;

  return {
    encoded,
    stats: {
      rawSize: data.length,
      compressedSize: compressed.length,
      oligoCount: totalOligoCount,
      payloadBytesPerOligo: chunkSize,
      netDensityBitsPerNt,
      overheadPercent,
      screeningRetries,
      encodeTimeMs,
    },
  };
}

// ─── v3.0 Canonical Archive (portable, deterministic intermediate format) ─────

/**
 * A single block in a canonical archive — packed 2-bit bytes before DNA mapping.
 *
 * This is the deterministic, portable intermediate representation. No DNA
 * constraint screening has been applied yet, so the bytes may or may not
 * satisfy GC/homopolymer constraints when mapped to DNA.
 */
export interface CanonicalBlock {
  /** Oligo index (sequential or content-derived). */
  index: number;
  /** Inner block bytes (address + payload + parity + CRC). */
  innerBytes: Uint8Array;
  /** Seed used for XOR retry (0 if no retry needed). */
  seed: number;
  /** Whether this is a parity oligo. */
  isParity: boolean;
}

/**
 * Canonical archive — the portable, deterministic intermediate format.
 *
 * Produced by `encodeToCanonical()` (steps 1–6 of the pipeline). Consumed by
 * `canonicalToSynthesis()` (step 7: DNA mapping + constraint screening).
 *
 * Separating these two phases enables:
 *   1. Testing the codec without DNA mapping (faster, deterministic).
 *   2. Different synthesis strategies for the same canonical data.
 *   3. Content-derived addressing (address derived from canonical bytes,
 *      not sequential position).
 */
export interface CanonicalArchive {
  metadata: CodecMetadata;
  /** Packed 2-bit blocks (before DNA mapping). */
  blocks: CanonicalBlock[];
  /** Forward primer. */
  forwardPrimer: string;
  /** Reverse primer. */
  reversePrimer: string;
  /** Archive salt for content-derived addressing (if addressMode !== 'sequential'). */
  archiveSalt?: Uint8Array;
}

// ─── v3.0 encodeToCanonical ───────────────────────────────────────────────────

/**
 * Encode data to a canonical archive (steps 1–6 only).
 *
 * This is the deterministic, portable intermediate format. It performs:
 *   1. (Optional) DEFLATE-compress.
 *   2. Compute SHA-256 hash.
 *   3. Split into chunks.
 *   4. Outer RS parity.
 *   5. For each oligo: inner RS/LDPC + CRC + pack to 2-bit bytes.
 *   6. Return CanonicalArchive WITHOUT DNA mapping or constraint screening.
 *
 * Use `canonicalToSynthesis()` to map canonical bytes → DNA with constraint
 * screening (where seed retries happen).
 *
 * @param data - Raw file bytes.
 * @param cfg  - Codec configuration.
 * @param meta - File metadata (name, content type).
 * @returns A CanonicalArchive — the portable, deterministic intermediate format.
 */
export async function encodeToCanonical(
  data: Uint8Array,
  cfg: CodecConfig,
  meta: { fileName: string; contentType: string },
): Promise<CanonicalArchive> {
  const layout = computeLayoutAuto(cfg);
  const { fwd, rev } = getPrimers(cfg);
  const useConvInner = !!cfg.useConvolutionalInner;

  // 1) Compress (optional)
  let compressed: Uint8Array = data;
  let compression: "none" | "deflate" = "none";
  if (cfg.compress) {
    compressed = deflate(data, { level: 9 });
    compression = "deflate";
  }

  // 2) Encrypt (optional)
  let processedData: Uint8Array = compressed;
  let encryptionSalt: Uint8Array | null = null;
  if (cfg.encryptPassword) {
    const { encrypt } = await import("./encryption");
    const salt = new Uint8Array(16);
    const { randomBytes } = await import("crypto");
    salt.set(randomBytes(16));
    const result = encrypt(compressed, { password: cfg.encryptPassword, salt });
    processedData = result.ciphertext;
    encryptionSalt = salt;
  }

  // 3) Hash the ORIGINAL data
  const fileHash = await sha256(data);

  // 3) Split into payload chunks
  const chunkSize = layout.payloadBytes;
  const dataOligoCount = Math.max(1, Math.ceil(compressed.length / chunkSize));
  const paddedLen = dataOligoCount * chunkSize;
  const padded = new Uint8Array(paddedLen);
  padded.set(compressed, 0);

  // 4) Outer RS
  const parityCount = Math.max(2, Math.ceil(dataOligoCount * cfg.outerParityRatio));
  if (dataOligoCount + parityCount > 65535) {
    throw new Error(
      `Data too large for single RS block: ${dataOligoCount + parityCount} oligos > 65535 limit. ` +
      `Split the data into smaller shards and call encodeToCanonical() for each shard.`
    );
  }
  const useGF216 = dataOligoCount + parityCount > 255;
  const outerN = dataOligoCount + parityCount;
  const outerK = dataOligoCount;
  const totalOligoCount = outerN;

  if (outerN <= outerK) {
    throw new Error(
      `Outer RS config invalid: n=${outerN} <= k=${outerK}. Increase outerParityRatio.`,
    );
  }

  const outerRs8 = !useGF216 ? new ReedSolomon({ n: outerN, k: outerK }) : null;
  const outerRs216 = useGF216 ? new ReedSolomon216({ n: outerN, k: outerK }) : null;

  // 5) Outer RS parity computation
  const parityBytes = new Uint8Array(parityCount * chunkSize);

  if (useGF216 && outerRs216) {
    const dataSymbols16 = new Uint16Array(dataOligoCount);
    const numPairs = Math.floor(chunkSize / 2);
    for (let pairIdx = 0; pairIdx < numPairs; pairIdx++) {
      const j0 = pairIdx * 2;
      const j1 = pairIdx * 2 + 1;
      for (let i = 0; i < dataOligoCount; i++) {
        dataSymbols16[i] = (padded[i * chunkSize + j0] << 8) | padded[i * chunkSize + j1];
      }
      const parity = outerRs216.parity(dataSymbols16);
      for (let i = 0; i < parity.length; i++) {
        parityBytes[i * chunkSize + j0] = (parity[i] >> 8) & 0xff;
        parityBytes[i * chunkSize + j1] = parity[i] & 0xff;
      }
    }
    if (chunkSize % 2 === 1) {
      const j = chunkSize - 1;
      if (dataOligoCount + parityCount <= 255) {
        const dataSymbols8 = new Uint8Array(dataOligoCount);
        for (let i = 0; i < dataOligoCount; i++) {
          dataSymbols8[i] = padded[i * chunkSize + j];
        }
        const rs8 = new ReedSolomon({ n: outerN, k: outerK });
        const parity = rs8.parity(dataSymbols8);
        for (let i = 0; i < parity.length; i++) {
          parityBytes[i * chunkSize + j] = parity[i];
        }
      } else {
        for (let i = 0; i < dataOligoCount; i++) {
          dataSymbols16[i] = padded[i * chunkSize + j] & 0xff;
        }
        const parity = outerRs216.parity(dataSymbols16);
        for (let i = 0; i < parity.length; i++) {
          parityBytes[i * chunkSize + j] = parity[i] & 0xff;
        }
      }
    }
  } else if (outerRs8) {
    const dataSymbols = new Uint8Array(dataOligoCount);
    for (let j = 0; j < chunkSize; j++) {
      for (let i = 0; i < dataOligoCount; i++) {
        dataSymbols[i] = padded[i * chunkSize + j];
      }
      const parity = outerRs8.parity(dataSymbols);
      for (let i = 0; i < parity.length; i++) {
        parityBytes[i * chunkSize + j] = parity[i];
      }
    }
  }

  // 6) Build each oligo's inner block (address + payload + parity + CRC)
  const useLDPC = (cfg.innerCode ?? "rs") === "ldpc";
  const innerRsN = layout.addressBytes + layout.payloadBytes + layout.innerParityBytes;
  const innerRs = innerRsN > 255
    ? new ReedSolomon216({ n: innerRsN, k: layout.addressBytes + layout.payloadBytes })
    : new ReedSolomon({ n: innerRsN, k: layout.addressBytes + layout.payloadBytes });

  const useArithmeticV2 = (cfg.mappingMode ?? "direct") === "arithmetic";
  const innerK = useArithmeticV2
    ? layout.payloadBytes
    : layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;

  const innerRsReal = innerN > 255
    ? new ReedSolomon216({ n: innerN, k: innerK })
    : new ReedSolomon({ n: innerN, k: innerK });

  const innerLdpcReal = useLDPC ? new LDPCInnerCode({ n: innerN, k: innerK }) : null;

  // Generate archive salt for content-derived addressing
  const addressMode = cfg.addressMode ?? 'sequential';
  const archiveSalt = addressMode !== 'sequential'
    ? (cfg.archiveSalt ?? deriveArchiveSalt())
    : undefined;

  const blocks: CanonicalBlock[] = [];

  for (let oligoIdx = 0; oligoIdx < totalOligoCount; oligoIdx++) {
    const payload = new Uint8Array(chunkSize);
    if (oligoIdx < dataOligoCount) {
      payload.set(padded.slice(oligoIdx * chunkSize, (oligoIdx + 1) * chunkSize), 0);
    } else if (oligoIdx < dataOligoCount + parityCount) {
      const parityIdx = oligoIdx - dataOligoCount;
      payload.set(
        parityBytes.slice(parityIdx * chunkSize, (parityIdx + 1) * chunkSize),
        0,
      );
    }

    // Build address
    const rawAddress = new Uint8Array(4);
    if (addressMode === 'content-derived' || addressMode === 'hierarchical') {
      // Content-derived addressing: derive address from payload hash
      const addrConfig: AddressingConfig = {
        mode: addressMode,
        archiveSalt: archiveSalt!,
        addressBytes: layout.addressBytes,
        hierarchicalDepth: addressMode === 'hierarchical' ? 3 : undefined,
      };
      const derivedAddr = deriveAddress(payload, addrConfig);
      rawAddress.set(derivedAddr.subarray(0, Math.min(derivedAddr.length, 4)), 0);
      rawAddress[3] = 0; // seed placeholder
    } else {
      // Sequential addressing (legacy)
      rawAddress[0] = (oligoIdx >> 16) & 0xff;
      rawAddress[1] = (oligoIdx >> 8) & 0xff;
      rawAddress[2] = oligoIdx & 0xff;
      rawAddress[3] = 0; // seed placeholder
    }

    const whitenedAddress = whitenAddress(rawAddress);

    // Build RS data: whitened address + payload
    const rsData = new Uint8Array(innerK);
    rsData.set(whitenedAddress, 0);
    rsData.set(payload, layout.addressBytes);

    // Inner RS/LDPC encode
    const rsCodeword = useLDPC && innerLdpcReal
      ? innerLdpcReal.encode(rsData)
      : innerRsReal.encode(rsData);

    // CRC-16 over RS codeword
    const crc = crc16Bytes(rsCodeword);

    // Pack inner block: RS codeword + CRC
    const innerBlock = new Uint8Array(totalNtBytes(layout));
    innerBlock.set(rsCodeword, 0);
    innerBlock.set(crc, rsCodeword.length);

    const isParity = oligoIdx >= dataOligoCount;

    blocks.push({
      index: oligoIdx,
      innerBytes: innerBlock,
      seed: 0, // No screening yet — seed will be set in canonicalToSynthesis
      isParity,
    });
  }

  // Build metadata
  const interleaveDepth = cfg.interleaveDepth ?? 0;
  const metadata: CodecMetadata = {
    fileName: meta.fileName,
    fileSize: data.length,
    fileHash,
    contentType: meta.contentType,
    compression,
    rawSize: data.length,
    oligoCount: totalOligoCount,
    payloadBytesPerOligo: chunkSize,
    innerRS: { n: innerN, k: innerK },
    innerCode: useLDPC ? "ldpc" : "rs",
    ldpcDecoder: cfg.ldpcDecoder as "hard" | "osd" | "bp" | "auto" | undefined,
    mappingMode: (cfg.mappingMode ?? "direct") as "direct" | "goldman" | "constrained" | "srt" | "arithmetic",
    goldmanMode: (cfg.goldmanMode ?? "fast") as "fast" | "dense",
    outerRS: { n: outerN, k: outerK },
    parityOligos: parityCount,
    interleaveDepth,
    encryptionSalt: encryptionSalt ? Buffer.from(encryptionSalt).toString("base64") : undefined,
    channel: (cfg.channel ?? "illumina") as "illumina" | "nanopore",
    lowCoverageTrigger: cfg.lowCoverageTrigger ?? 5,
    useConvolutionalInner: useConvInner,
    version: 1,
    encodedAt: new Date().toISOString(),
  };

  return {
    metadata,
    blocks,
    forwardPrimer: fwd,
    reversePrimer: rev,
    archiveSalt,
  };
}

// ─── v3.0 canonicalToSynthesis ────────────────────────────────────────────────

/**
 * Map canonical archive bytes → DNA with constraint screening.
 *
 * This is where GC/homopolymer screening and seed retries happen. Separating
 * this step from the canonical encoding enables:
 *   1. Testing the codec without DNA mapping (faster, deterministic).
 *   2. Different synthesis strategies for the same canonical data.
 *   3. Content-derived addressing (address derived from canonical bytes).
 *
 * @param archive - Canonical archive from `encodeToCanonical()`.
 * @param cfg    - Codec configuration (constraints, mapping mode, etc.).
 * @returns EncodedFile with DNA oligos ready for synthesis.
 */
export function canonicalToSynthesis(
  archive: CanonicalArchive,
  cfg: CodecConfig,
): EncodedFile {
  const layout = computeLayoutAuto(cfg);
  const { fwd, rev } = getPrimers(cfg);
  const useArithmeticV2 = (cfg.mappingMode ?? "direct") === "arithmetic";
  const useGoldman = (cfg.mappingMode ?? "direct") === "goldman";
  const useConstrained = (cfg.mappingMode ?? "direct") === "constrained";
  const useSrt = (cfg.mappingMode ?? "direct") === "srt";
  const useArithmetic = (cfg.mappingMode ?? "direct") === "arithmetic";

  const constraints = {
    gcMin: cfg.constraints.gcMin,
    gcMax: cfg.constraints.gcMax,
    maxHomopolymer: cfg.constraints.maxHomopolymer,
  };

  const oligos: Oligo[] = [];

  for (const block of archive.blocks) {
    let dna = "";
    let seed = 0;

    if (useGoldman) {
      // Goldman rotational mapping: no screening needed
      dna = bytesToGoldmanDna(block.innerBytes, "A", cfg.goldmanMode ?? "fast");
      seed = 0;
    } else if (useConstrained) {
      // Split constrained mapping
      dna = bytesToSplitConstrainedDna(block.innerBytes, cfg.constraints.maxHomopolymer, layout.addressBytes);
      seed = 0;
    } else if (useSrt) {
      // SRT constrained mapping
      dna = bytesToSrtDna(block.innerBytes, cfg.constraints.maxHomopolymer, totalNtBytes(layout) * 4);
      seed = 0;
    } else {
      // Direct 2-bit mapping with constraint screening + seed retries
      let bestDna = "";
      let bestSeed = 0;
      let bestSatisfied = false;
      let attempts = 0;
      let currentSeed = 0;

      while (attempts <= cfg.maxRetries) {
        // Apply XOR with seed if retrying (seed > 0)
        const effectiveBytes = currentSeed === 0
          ? block.innerBytes
          : xorWithSeed(block.innerBytes, currentSeed);

        const candidateDna = bytesToDna(effectiveBytes);
        const ok = satisfiesConstraints(candidateDna, constraints);

        if (ok) {
          bestDna = candidateDna;
          bestSeed = currentSeed;
          bestSatisfied = true;
          break;
        }
        if (!bestDna) {
          bestDna = candidateDna;
          bestSeed = currentSeed;
        }
        attempts++;
        currentSeed = attempts;
      }

      dna = bestDna;
      seed = bestSeed;
    }

    // Compute stats on the payload region (excluding primers)
    const gc = gcContent(dna);
    const maxHp = maxHomopolymerRun(dna);

    // Pad DNA to exactly the expected inner nt length
    const expectedDnaLen = useArithmeticV2
      ? (cfg.oligoLength - 2 * cfg.primerLength)
      : layout.totalInnerBytes * 4;
    if (dna.length < expectedDnaLen) {
      if (useArithmeticV2) {
        const lastBase = dna[dna.length - 1] || "A";
        const padBase: string = lastBase === "A" ? "C" : "A";
        dna += padBase.repeat(expectedDnaLen - dna.length);
      } else {
        dna = dna + "A".repeat(expectedDnaLen - dna.length);
      }
    } else if (dna.length > expectedDnaLen) {
      dna = dna.slice(0, expectedDnaLen);
    }

    const sequence = fwd + dna + rev;
    oligos.push({
      index: block.index,
      sequence,
      gc,
      maxHomopolymer: maxHp,
      seed,
      payloadBytes: layout.payloadBytes,
      length: sequence.length,
    });
  }

  // Interleaving post-processing (same as encodeFile)
  const interleaveDepth = cfg.interleaveDepth ?? 0;
  if (interleaveDepth > 1 && oligos.length >= interleaveDepth) {
    const primerLen = fwd.length;
    const innerDnaLen = layout.totalInnerBytes * 4;
    const addressBytes = layout.addressBytes;

    for (let g = 0; g < oligos.length; g += interleaveDepth) {
      const groupSize = Math.min(interleaveDepth, oligos.length - g);
      if (groupSize < 2) continue;

      const regions: Uint8Array[] = [];
      for (let i = 0; i < groupSize; i++) {
        const oligo = oligos[g + i];
        let seq = oligo.sequence;
        if (seq.length < cfg.oligoLength) {
          seq = seq + "A".repeat(cfg.oligoLength - seq.length);
        }
        const innerDna = seq.slice(primerLen, primerLen + innerDnaLen);
        const fullBlock = dnaToBytes(innerDna);
        regions.push(fullBlock.slice(addressBytes));
      }

      const interleaved = interleaveCodewords(regions);
      const blockLen = regions[0].length;

      for (let i = 0; i < groupSize; i++) {
        const oligo = oligos[g + i];
        let seq = oligo.sequence;
        if (seq.length < cfg.oligoLength) {
          seq = seq + "A".repeat(cfg.oligoLength - seq.length);
        }
        const innerDna = seq.slice(primerLen, primerLen + innerDnaLen);
        const fullBlock = dnaToBytes(innerDna);
        const interleavedChunk = interleaved.slice(i * blockLen, (i + 1) * blockLen);
        fullBlock.set(interleavedChunk, addressBytes);
        const newInnerDna = bytesToDna(fullBlock);
        const oldSuffix = oligo.sequence.slice(primerLen + innerDnaLen);
        const prefix = oligo.sequence.slice(0, primerLen);
        oligo.sequence = prefix + newInnerDna + oldSuffix;
      }
    }
  }

  return {
    metadata: archive.metadata,
    oligos,
    forwardPrimer: archive.forwardPrimer,
    reversePrimer: archive.reversePrimer,
  };
}

// ─── v3.0 deriveAndSetAddresses ───────────────────────────────────────────────

/**
 * Replace sequential addresses with content-derived ones.
 *
 * For `addressMode === 'content-derived'`: compute BLAKE3-derived addresses
 * from each block's payload bytes and update the block's `innerBytes` in place.
 *
 * For `addressMode === 'hierarchical'`: also compute pool/well/oligo components
 * and store the hierarchical address in the block's `index` field (encoded as
 * a single integer from the 4-byte address for backward compatibility).
 *
 * This function mutates the archive's blocks in place.
 *
 * @param archive - Canonical archive (blocks will be mutated).
 * @param cfg    - Codec configuration (addressMode, archiveSalt, etc.).
 *
 * @example
 * ```ts
 * const archive = await encodeToCanonical(data, cfg, meta);
 * if (cfg.addressMode === 'content-derived') {
 *   deriveAndSetAddresses(archive, cfg);
 * }
 * const encoded = canonicalToSynthesis(archive, cfg);
 * ```
 */
export function deriveAndSetAddresses(
  archive: CanonicalArchive,
  cfg: CodecConfig,
): void {
  const addressMode = cfg.addressMode ?? 'sequential';
  if (addressMode === 'sequential') {
    // Nothing to do — sequential addresses are already set during encoding.
    return;
  }

  const layout = computeLayoutAuto(cfg);
  const archiveSalt = archive.archiveSalt ?? cfg.archiveSalt ?? deriveArchiveSalt();
  // Persist the salt if we just generated it
  if (!archive.archiveSalt) {
    archive.archiveSalt = archiveSalt;
  }

  const addrConfig: AddressingConfig = {
    mode: addressMode,
    archiveSalt,
    addressBytes: layout.addressBytes,
    hierarchicalDepth: addressMode === 'hierarchical' ? 3 : undefined,
  };

  for (const block of archive.blocks) {
    // Extract payload from innerBytes: skip addressBytes, take payloadBytes
    const payloadBytes = block.innerBytes.slice(
      layout.addressBytes,
      layout.addressBytes + layout.payloadBytes,
    );

    // Derive content-based address
    const derivedAddr = deriveAddress(payloadBytes, addrConfig);

    // Replace the address portion of innerBytes
    // Address is the first `layout.addressBytes` of the inner block
    for (let i = 0; i < layout.addressBytes && i < derivedAddr.length; i++) {
      block.innerBytes[i] = derivedAddr[i];
    }

    // For hierarchical mode, compute and store the hierarchical components
    if (addressMode === 'hierarchical') {
      const hAddr = deriveHierarchicalAddress(payloadBytes, addrConfig);
      // Re-encode the hierarchical address as a 4-byte address
      // pool (2 bytes) + well (1 byte) + oligo (1 byte)
      const poolBytes = base32Decode(hAddr.pool);
      const wellBytes = base32Decode(hAddr.well);
      const oligoBytes = base32Decode(hAddr.oligo);

      // Write into innerBytes (first 4 bytes = address)
      block.innerBytes[0] = poolBytes.length > 0 ? poolBytes[0] : 0;
      block.innerBytes[1] = poolBytes.length > 1 ? poolBytes[1] : 0;
      block.innerBytes[2] = wellBytes.length > 0 ? wellBytes[0] : 0;
      block.innerBytes[3] = oligoBytes.length > 0 ? oligoBytes[0] : 0;

      // Update index to a deterministic value derived from the address
      block.index = (block.innerBytes[0] << 24) |
                    (block.innerBytes[1] << 16) |
                    (block.innerBytes[2] << 8) |
                    block.innerBytes[3];
    } else {
      // Content-derived mode: update index from the 4-byte address
      block.index = (derivedAddr[0] << 24) |
                    (derivedAddr[1] << 16) |
                    (derivedAddr[2] << 8) |
                    (derivedAddr.length > 3 ? derivedAddr[3] : 0);
    }
  }
}

/**
 * Decode a Crockford Base32 string to bytes.
 *
 * Inverse of `base32Encode` from addressing.ts. Handles the same alphabet:
 * `0-9, A-V` (excludes I, L, O, U).
 *
 * @param str - Base32-encoded string.
 * @returns Decoded byte array.
 */
function base32Decode(str: string): Uint8Array {
  if (!str) return new Uint8Array(0);

  const BASE32_ALPHABET = '0123456789ABCDEFGHJKMNPQRSTVWXYZ';
  const lookup = new Map<string, number>();
  for (let i = 0; i < BASE32_ALPHABET.length; i++) {
    lookup.set(BASE32_ALPHABET[i], i);
    lookup.set(BASE32_ALPHABET[i].toLowerCase(), i); // case-insensitive
  }

  let buffer = 0;
  let bits = 0;
  const output: number[] = [];

  for (let i = 0; i < str.length; i++) {
    const val = lookup.get(str[i]);
    if (val === undefined) continue; // skip invalid chars

    buffer = (buffer << 5) | val;
    bits += 5;

    if (bits >= 8) {
      bits -= 8;
      output.push((buffer >> bits) & 0xff);
    }
  }

  return new Uint8Array(output);
}

function totalNtBytes(layout: OligoLayout): number {
  return layout.addressBytes + layout.payloadBytes + layout.innerParityBytes + layout.crcBytes;
}

// --- Convenience: serialize EncodedFile to JSON ---

export function serializeEncodedFile(encoded: EncodedFile): string {
  return JSON.stringify(encoded, null, 2);
}

export function deserializeEncodedFile(json: string): EncodedFile {
  return JSON.parse(json) as EncodedFile;
}

export { DEFAULT_CONFIG, computeLayout, computeLayoutAuto };
