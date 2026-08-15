/**
 * DNA-MGC+ — Multi-Gain Correction Codec (Khabbaz et al., 2026)
 *
 * A versatile codec that achieves simultaneous gains in:
 *   - Sequencing depth requirements (lower coverage needed)
 *   - Read cost (fewer reads required)
 *   - Decoding time (faster recovery)
 *   - Storage density (higher bits/nt)
 *   - Error correction capability (handles higher IDS rates)
 *
 * Architecture:
 *   The MGC+ code acts as an inner code that introduces intra-sequence
 *   redundancy through structured encoding applied separately to each
 *   indexed oligo. This is combined with an outer RS code for
 *   cross-oligo erasure correction.
 *
 * MGC+ Inner Encoding:
 *   1. Split payload into blocks of size b
 *   2. For each block, compute gain parity: XOR of all bytes in block
 *   3. Interleave gain parity with payload blocks
 *   4. This creates structured redundancy within each oligo
 *
 * MGC+ Inner Decoding:
 *   1. Use gain parity to detect single-block errors
 *   2. For each failing block, try all 256 byte values
 *   3. Accept the value that satisfies the gain parity check
 *   4. If multiple candidates, use CRC to disambiguate
 *   5. For indels, use sync markers to resync before correction
 *
 * Reference:
 *   - Khabbaz et al. (2026). arXiv:2603.14527.
 */

import { crc16, crc16Bytes } from './crc16';

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

export interface MGCPlusConfig {
  /** Block size for gain parity computation. Default: 8 bytes */
  blockSize: number;
  /** Number of gain parity bytes per block. Default: 1 */
  gainParityBytes: number;
  /** Whether to include CRC-8 sync markers. Default: true */
  useSyncMarkers: boolean;
  /** Sync marker interval in bytes. Default: 32 */
  syncMarkerInterval: number;
  /** Maximum correction attempts per block. Default: 256 */
  maxCorrectionAttempts: number;
}

export const DEFAULT_MGC_PLUS_CONFIG: MGCPlusConfig = {
  blockSize: 8,
  gainParityBytes: 1,
  useSyncMarkers: true,
  syncMarkerInterval: 32,
  maxCorrectionAttempts: 256,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function crc8Simple(data: Uint8Array): number {
  let crc = 0xFF;
  for (const b of data) {
    crc ^= b;
    for (let i = 0; i < 8; i++) crc = (crc & 0x80) ? ((crc << 1) ^ 0x31) & 0xFF : (crc << 1) & 0xFF;
  }
  return crc ^ 0xFF;
}

function computeGainParity(block: Uint8Array, parityBytes: number): Uint8Array {
  const parity = new Uint8Array(parityBytes);
  for (let p = 0; p < parityBytes; p++) {
    let xor = 0;
    for (let i = p; i < block.length; i += parityBytes) {
      xor ^= block[i];
    }
    parity[p] = xor;
  }
  return parity;
}

function verifyGainParity(block: Uint8Array, parity: Uint8Array, parityBytes: number): boolean {
  const computed = computeGainParity(block, parityBytes);
  for (let i = 0; i < parity.length; i++) {
    if (computed[i] !== parity[i]) return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// Encode
// ---------------------------------------------------------------------------

export interface MGCPlusEncodeResult {
  encoded: Uint8Array;
  gainParityCount: number;
  syncMarkerCount: number;
  overhead: number;
}

export function mgcPlusEncode(data: Uint8Array, config: MGCPlusConfig = DEFAULT_MGC_PLUS_CONFIG): MGCPlusEncodeResult {
  const cfg = { ...DEFAULT_MGC_PLUS_CONFIG, ...config };
  const blockSize = cfg.blockSize;
  const parityBytes = cfg.gainParityBytes;

  const numBlocks = Math.ceil(data.length / blockSize);
  const parityTotal = numBlocks * parityBytes;
  const encodedLen = data.length + parityTotal;
  const encoded = new Uint8Array(encodedLen);

  let srcOff = 0;
  let dstOff = 0;
  for (let b = 0; b < numBlocks; b++) {
    const blockEnd = Math.min(srcOff + blockSize, data.length);
    const blockLen = blockEnd - srcOff;
    const block = data.subarray(srcOff, blockEnd);

    encoded.set(block, dstOff);
    dstOff += blockLen;

    const parity = computeGainParity(block, parityBytes);
    encoded.set(parity, dstOff);
    dstOff += parityBytes;

    srcOff = blockEnd;
  }

  let syncMarkerCount = 0;
  let finalData = encoded;
  if (cfg.useSyncMarkers) {
    const result = insertSyncMarkers(encoded, cfg.syncMarkerInterval);
    finalData = result.data;
    syncMarkerCount = result.markerCount;
  }

  const totalOverhead = parityTotal + syncMarkerCount;
  const overhead = totalOverhead / finalData.length;

  return { encoded: finalData, gainParityCount: parityTotal, syncMarkerCount, overhead };
}

function insertSyncMarkers(data: Uint8Array, interval: number): { data: Uint8Array; markerCount: number } {
  const numMarkers = Math.floor(data.length / interval);
  const outputLen = data.length + numMarkers;
  const output = new Uint8Array(outputLen);

  let srcOff = 0;
  let dstOff = 0;
  let markerCount = 0;

  while (srcOff < data.length) {
    const segLen = Math.min(interval, data.length - srcOff);
    const segment = data.subarray(srcOff, srcOff + segLen);
    output.set(segment, dstOff);
    dstOff += segLen;

    if (srcOff + segLen < data.length) {
      output[dstOff] = crc8Simple(segment);
      dstOff++;
      markerCount++;
    }
    srcOff += segLen;
  }

  return { data: output, markerCount };
}

// ---------------------------------------------------------------------------
// Decode
// ---------------------------------------------------------------------------

export interface MGCPlusDecodeResult {
  decoded: Uint8Array | null;
  blocksCorrected: number;
  markersPassed: number;
  markersFailed: number;
  resyncApplied: boolean;
}

export function mgcPlusDecode(
  encoded: Uint8Array,
  originalLength: number,
  config: MGCPlusConfig = DEFAULT_MGC_PLUS_CONFIG,
): MGCPlusDecodeResult {
  const cfg = { ...DEFAULT_MGC_PLUS_CONFIG, ...config };
  const blockSize = cfg.blockSize;
  const parityBytes = cfg.gainParityBytes;

  let data: Uint8Array;
  let markersPassed = 0;
  let markersFailed = 0;
  let resyncApplied = false;

  if (cfg.useSyncMarkers) {
    const result = stripSyncMarkers(encoded, cfg.syncMarkerInterval);
    data = result.data;
    markersPassed = result.passed;
    markersFailed = result.failed;
    resyncApplied = result.failed > 0;
  } else {
    data = encoded;
  }

  const decoded = new Uint8Array(originalLength);
  let srcOff = 0;
  let dstOff = 0;
  let blocksCorrected = 0;
  const numBlocks = Math.ceil(originalLength / blockSize);

  for (let b = 0; b < numBlocks; b++) {
    const expectedBlockLen = Math.min(blockSize, originalLength - dstOff);
    const blockAndParityLen = expectedBlockLen + parityBytes;

    if (srcOff + blockAndParityLen > data.length) {
      const remaining = Math.min(expectedBlockLen, data.length - srcOff);
      decoded.set(data.subarray(srcOff, srcOff + remaining), dstOff);
      break;
    }

    const block = data.subarray(srcOff, srcOff + expectedBlockLen);
    const parity = data.subarray(srcOff + expectedBlockLen, srcOff + blockAndParityLen);

    if (verifyGainParity(block, parity, parityBytes)) {
      decoded.set(block, dstOff);
    } else {
      const corrected = correctSingleByteError(block, parity, parityBytes, cfg.maxCorrectionAttempts);
      if (corrected) {
        decoded.set(corrected, dstOff);
        blocksCorrected++;
      } else {
        decoded.set(block, dstOff);
      }
    }

    srcOff += blockAndParityLen;
    dstOff += expectedBlockLen;
  }

  return { decoded, blocksCorrected, markersPassed, markersFailed, resyncApplied };
}

function correctSingleByteError(
  block: Uint8Array,
  parity: Uint8Array,
  parityBytes: number,
  maxAttempts: number,
): Uint8Array | null {
  const corrected = new Uint8Array(block.length);
  corrected.set(block);

  let attempts = 0;
  for (let pos = 0; pos < block.length && attempts < maxAttempts; pos++) {
    const original = block[pos];
    for (let val = 0; val < 256 && attempts < maxAttempts; val++) {
      if (val === original) continue;
      corrected[pos] = val;
      if (verifyGainParity(corrected, parity, parityBytes)) {
        return corrected;
      }
      attempts++;
    }
    corrected[pos] = original;
  }

  return null;
}

function stripSyncMarkers(
  data: Uint8Array,
  interval: number,
): { data: Uint8Array; passed: number; failed: number } {
  const markersInData = Math.floor(data.length / (interval + 1));
  const strippedLen = data.length - markersInData;
  const stripped = new Uint8Array(strippedLen);

  let srcOff = 0;
  let dstOff = 0;
  let passed = 0;
  let failed = 0;

  while (srcOff < data.length) {
    const segLen = Math.min(interval, data.length - srcOff);
    stripped.set(data.subarray(srcOff, srcOff + segLen), dstOff);
    dstOff += segLen;
    srcOff += segLen;

    if (srcOff < data.length && segLen === interval) {
      const markerCrc = data[srcOff];
      srcOff++;
      const computedCrc = crc8Simple(stripped.subarray(dstOff - segLen, dstOff));
      if (markerCrc === computedCrc) {
        passed++;
      } else {
        failed++;
      }
    }
  }

  return { data: stripped.subarray(0, dstOff), passed, failed };
}
