/**
 * v61: Arithmetic Mode with Address OUTSIDE the Arithmetic Stream
 *
 * Problem (v57-v60): The address (4 bytes = 16 nt) was encoded INSIDE the
 * LDPC codeword, which was then arithmetic-coded. This caused two issues:
 *   1. Address corruption from arithmetic termination: the last few bytes of
 *      each arithmetic block are unreliable (termination corruption), and if
 *      the address happens to land in this region, clustering fails.
 *   2. Address corruption from IDS: indels in the address region shift the
 *      arithmetic interval state, causing cascading decode errors.
 *
 * Solution (v61): Split the oligo layout:
 *   [Primer][Address (16 nt direct DNA)][Arithmetic-encoded payload+parity]
 *                                          ↑
 *                          Address sits OUTSIDE the arithmetic stream
 *
 * Layout (between primers):
 *   - Address: 4 bytes = 16 nt of direct 2-bit mapping (uncoded, raw)
 *     - 3 bytes oligo index (24-bit, supports up to 16M oligos)
 *     - 1 byte seed (for constraint screening, usually 0)
 *   - Payload+Parity: arithmetic-coded at 1.9 b/nt
 *     - LDPC codeword (payload + parity, no address!)
 *     - Per-block CRC-8 sync markers (DNA-Aeon style)
 *
 * Density impact:
 *   - Address overhead: 16 nt per oligo (was inside arithmetic stream, ~8 nt
 *     effective). The arithmetic stream loses ~16 nt of capacity.
 *   - At 700nt oligo / 660 nt inner, address is 16 nt = 2.4% overhead.
 *   - Arithmetic rate 1.9 b/nt × (660-16) = 1224 bits = 153 bytes per oligo.
 *   - With 8B parity: 145B payload per oligo.
 *   - Net density: 145 × 8 / 700 = 1.657 b/nt (close to v60's 1.664).
 *   - But arithmetic mode now WORKS, unlocking the 1.85+ b/nt path.
 *
 * Reversibility:
 *   - Address is read directly (no arithmetic decode) → robust to indels via
 *     k-mer clustering.
 *   - Arithmetic stream is independent → termination corruption no longer
 *     affects address.
 *   - LDPC codeword no longer contains address → simpler erasure decoding.
 *
 * References:
 *   - Ding et al. (2024). arXiv:2410.04886 — separates index from payload.
 *   - Welzel et al. (2023). DNA-Aeon, Nature Comms 14:628.
 */

import { Base } from "./mapping";
import {
  bytesToArithmeticDnaCrc,
  arithmeticDnaToBytesCrc,
} from "./markov-arithmetic";
import { bytesToDna, dnaToBytes, whitenAddress } from "./mapping";

/**
 * Layout for arithmetic-v2 oligo (address outside arithmetic stream).
 *
 * Total inner nt = addressNt + arithmeticNt
 *   - addressNt = 16 (4 bytes × 4 nt/byte, direct 2-bit mapping)
 *   - arithmeticNt = totalInnerNt - 16
 *
 * Arithmetic stream encodes:
 *   [payload (payloadBytes)] + [LDPC parity (innerParityBytes)]
 *   NO address, NO CRC-16 (per-block CRC-8 is inside the arithmetic stream)
 *
 * Total bytes encoded by arithmetic stream = payloadBytes + innerParityBytes
 * Required arithmetic nt = ceil(totalBytes * 8 / ARITH_CAPACITY_RATE)
 */
export interface ArithmeticV2Layout {
  /** Total inner nt (between primers) */
  totalInnerNt: number;
  /** Address length in nt (always 16 for 4-byte address) */
  addressNt: number;
  /** Arithmetic stream length in nt */
  arithmeticNt: number;
  /** Address length in bytes (always 4) */
  addressBytes: number;
  /** Payload bytes per oligo */
  payloadBytes: number;
  /** LDPC parity bytes */
  innerParityBytes: number;
  /** Total bytes encoded by arithmetic stream = payloadBytes + innerParityBytes */
  arithmeticDataBytes: number;
  /** LDPC codeword length (k + parity) — same as arithmeticDataBytes */
  innerN: number;
  /** LDPC info length (k = payloadBytes, NO address) */
  innerK: number;
}

/**
 * Compute the arithmetic-v2 layout for a given config.
 *
 * Capacity formula:
 *   arithmeticDataBytes = floor(arithmeticNt * ARITH_CAPACITY_RATE / 8)
 *   innerParityBytes = cfg.innerParityBytes
 *   payloadBytes = arithmeticDataBytes - innerParityBytes
 *
 * @param oligoLength Total oligo length in nt
 * @param primerLength Primer length in nt (per side)
 * @param innerParityBytes LDPC parity bytes
 */
export function computeArithmeticV2Layout(
  oligoLength: number,
  primerLength: number,
  innerParityBytes: number,
  blockSize: number = 80,
): ArithmeticV2Layout {
  const totalInnerNt = oligoLength - 2 * primerLength;
  const addressNt = 16; // 4 bytes × 4 nt/byte

  // v61: Round arithmeticNt DOWN to a multiple of blockSize.
  // The arithmetic encoder pads targetLen up to the next multiple of blockSize,
  // so we need arithmeticNt to already be a multiple to avoid length mismatch.
  //
  // v61 uses blockSize=80 (vs markov-arithmetic.ts default of 20):
  //   - blockSize=20: 4 bytes/block, 1 CRC = 25% overhead → density 1.0 b/nt
  //   - blockSize=80: 19 bytes/block, 1 CRC = 5.3% overhead → density 1.55 b/nt
  //   - blockSize=160: 38 bytes/block, 1 CRC = 2.6% overhead → density 1.65 b/nt
  // 80 is a good balance of density and error confinement (one block error
  // affects ≤19 bytes, easily corrected by LDPC erasure decoder).
  const rawArithmeticNt = totalInnerNt - addressNt;
  const arithmeticNt = Math.floor(rawArithmeticNt / blockSize) * blockSize;

  if (arithmeticNt < blockSize * 2) {
    throw new Error(
      `arithmetic-v2: arithmeticNt=${arithmeticNt} too small (need ≥${blockSize * 2}). ` +
      `Use longer oligo or shorter primer.`,
    );
  }

  // Capacity: matches markov-arithmetic.ts bytesToArithmeticDnaCrc
  const ARITH_CAPACITY_RATE = 1.95; // matches markov-arithmetic.ts
  const bytesPerBlockTotal = Math.max(2, Math.floor((blockSize * ARITH_CAPACITY_RATE) / 8));
  const bytesPerBlockData = bytesPerBlockTotal - 1; // 1 byte for CRC-8
  const numBlocks = arithmeticNt / blockSize;
  const arithmeticDataBytes = numBlocks * bytesPerBlockData;

  const innerN = arithmeticDataBytes;
  const innerK = innerN - innerParityBytes;
  let payloadBytes = innerK;

  if (payloadBytes <= 0) {
    throw new Error(
      `arithmetic-v2: payloadBytes=${payloadBytes} ≤ 0. Reduce innerParityBytes ` +
      `or increase oligoLength.`,
    );
  }

  // Force even for GF(2^16) outer RS compatibility
  if (payloadBytes % 2 === 1) {
    payloadBytes -= 1;
  }

  return {
    totalInnerNt,
    addressNt,
    arithmeticNt,
    addressBytes: 4,
    payloadBytes,
    innerParityBytes,
    arithmeticDataBytes,
    innerN,
    innerK: payloadBytes, // k = payload (no address in LDPC codeword)
  };
}

/**
 * Encode an oligo using arithmetic-v2 layout.
 *
 * Input:
 *   - oligoIdx: Oligo index (0..N-1)
 *   - payload: Payload bytes (length = layout.payloadBytes)
 *   - ldpcParity: LDPC parity bytes (length = layout.innerParityBytes)
 *
 * Output:
 *   - inner DNA string (length = layout.totalInnerNt)
 *     - First 16 nt: address (direct 2-bit mapping of whitened address)
 *     - Next layout.arithmeticNt: arithmetic-coded payload+parity
 *
 * @returns Inner DNA string (between primers)
 */
export function encodeArithmeticV2Oligo(
  oligoIdx: number,
  payload: Uint8Array,
  ldpcParity: Uint8Array,
  layout: ArithmeticV2Layout,
  maxHomopolymer: number = 3,
  blockSize: number = 80,
): string {
  // 1. Build address: 3 bytes index + 1 byte seed (0)
  const address = new Uint8Array(4);
  address[0] = (oligoIdx >> 16) & 0xff;
  address[1] = (oligoIdx >> 8) & 0xff;
  address[2] = oligoIdx & 0xff;
  address[3] = 0; // seed placeholder (no screening in v61 — see note below)
  const whitenedAddress = whitenAddress(address);
  const addressDna = bytesToDna(whitenedAddress); // 16 nt

  // 2. Build arithmetic data: payload + LDPC parity (NO address, NO CRC-16)
  //    The per-block CRC-8 sync markers are added INSIDE the arithmetic
  //    encoder (bytesToArithmeticDnaCrc), not as separate bytes.
  const arithData = new Uint8Array(layout.arithmeticDataBytes);
  arithData.set(payload, 0);
  arithData.set(ldpcParity, layout.payloadBytes);

  // 3. Arithmetic-encode the data (with per-block CRC-8 sync markers)
  //    Pass blockSize explicitly — must match the layout's blockSize.
  const arithmeticDna = bytesToArithmeticDnaCrc(
    arithData,
    maxHomopolymer,
    layout.arithmeticNt,
    blockSize,
  );

  // 4. Concatenate: address (direct) + arithmetic stream
  let dna = addressDna + arithmeticDna;

  // 5. Pad to totalInnerNt (the arithmeticNt may be < totalInnerNt - addressNt
  //    due to blockSize rounding). Use a base that doesn't extend the last
  //    homopolymer run beyond maxHomopolymer.
  if (dna.length < layout.totalInnerNt) {
    const lastBase = dna[dna.length - 1] as Base;
    const padBase: Base = lastBase === "A" ? "C" : "A";
    dna += padBase.repeat(layout.totalInnerNt - dna.length);
  }

  return dna;
}

/**
 * Decode an oligo using arithmetic-v2 layout.
 *
 * Input:
 *   - innerDna: Inner DNA string (between primers, length = layout.totalInnerNt)
 *
 * Output:
 *   - oligoIdx: Decoded oligo index (or -1 if address undecodable)
 *   - payload: Decoded payload bytes
 *   - ldpcParity: Decoded LDPC parity bytes
 *
 * Notes:
 *   - Address is decoded by direct 2-bit mapping (no arithmetic decode).
 *     This is robust to indels via k-mer clustering.
 *   - Arithmetic stream is decoded independently. Termination corruption
 *     only affects the last byte of the stream, not the address.
 */
export function decodeArithmeticV2Oligo(
  innerDna: string,
  layout: ArithmeticV2Layout,
  maxHomopolymer: number = 3,
  blockSize: number = 80,
): {
  oligoIdx: number;
  seed: number;
  payload: Uint8Array;
  ldpcParity: Uint8Array;
  arithmeticSuccess: boolean;
} {
  // 1. Extract address (first 16 nt) and decode via direct mapping
  const addressDna = innerDna.slice(0, layout.addressNt);
  let oligoIdx = -1;
  let seed = 0;
  try {
    const addressBytes = dnaToBytes(addressDna);
    // Unwhiten: whitenAddress is its own inverse (XOR-based)
    const unwhitened = whitenAddress(addressBytes);
    oligoIdx = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
    seed = unwhitened[3];
  } catch {
    // Address contained 'N' or was corrupted — caller should fall back to k-mer clustering
    oligoIdx = -1;
  }

  // 2. Extract arithmetic stream and decode
  const arithmeticDna = innerDna.slice(layout.addressNt);
  let arithData: Uint8Array;
  let arithmeticSuccess = true;
  try {
    const result = arithmeticDnaToBytesCrc(
      arithmeticDna,
      maxHomopolymer,
      layout.arithmeticDataBytes,
      blockSize,
    );
    arithData = result.data;
  } catch {
    // Arithmetic decode failed — return zeros, LDPC may recover via erasure
    arithData = new Uint8Array(layout.arithmeticDataBytes);
    arithmeticSuccess = false;
  }

  // 3. Split arithmetic data into payload + parity
  const payload = arithData.slice(0, layout.payloadBytes);
  const ldpcParity = arithData.slice(layout.payloadBytes);

  return { oligoIdx, seed, payload, ldpcParity, arithmeticSuccess };
}

/**
 * Compute the theoretical density of arithmetic-v2 mode.
 *
 * Density = (payloadBits × outerEfficiency) / oligoLength
 *
 * @param layout Arithmetic-v2 layout
 * @param outerParityRatio Outer RS parity ratio (0.0 - 1.0)
 * @param oligoLength Total oligo length in nt
 */
export function computeArithmeticV2Density(
  layout: ArithmeticV2Layout,
  outerParityRatio: number,
  oligoLength: number,
): number {
  const payloadBits = layout.payloadBytes * 8;
  const outerEfficiency = 1 / (1 + outerParityRatio);
  const infoBits = payloadBits * outerEfficiency;
  return infoBits / oligoLength;
}
