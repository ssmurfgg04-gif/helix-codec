/**
 * CRC-Marker Inner Code (DNA-Aeon / HEDGES pattern)
 *
 * Inserts periodic CRC markers into the payload so the decoder can re-anchor
 * after insertions/deletions. When a marker CRC fails, the decoder knows an
 * indel occurred within the preceding marker window and can hypothesize
 * shifts to re-sync.
 *
 * Format (DNA-Aeon-inspired, sync=4):
 *   [payload_seg_0 (4B)][CRC-8 marker (1B)] [payload_seg_1 (4B)][CRC-8 marker] ...
 *   [payload_seg_N (4B)][CRC-8 marker] [inner RS parity]
 *
 * The CRC-8 marker is computed over the preceding 4-byte payload segment.
 * Marker every 4 bytes = 20% overhead, but enables indel re-synchronization.
 *
 * Decode strategy:
 *   1. Walk the read in 5-byte windows (4B payload + 1B marker).
 *   2. For each window, verify CRC-8 of the 4 payload bytes against the marker.
 *   3. If CRC passes: payload is intact (no indel in this window).
 *   4. If CRC fails: an indel occurred. Hypothesize ±1, ±2, ±3 base shifts
 *      and re-check subsequent markers until one validates.
 *   5. After 3 consecutive marker failures (checkpoint=3), declare the window
 *      as an erasure and let the outer RS code recover it.
 *
 * Reference:
 *   - Welzel et al. (2023). "DNA-Aeon." Nature Comms 14:433. (sync=4, checkpoint=3)
 *   - Press et al. (2020). "HEDGES." PNAS 117:31. (convolutional + re-sync)
 */

import { crc16 } from "./crc16";

export interface CRCMarkerConfig {
  /** Number of payload bytes between markers. Default: 4 (DNA-Aeon). */
  segmentSize: number;
  /** Marker size in bytes (CRC-8 = 1 byte, CRC-16 = 2 bytes). Default: 1. */
  markerSize: number;
  /** Number of consecutive failures before declaring erasure. Default: 3. */
  checkpoint: number;
}

export const DEFAULT_CRC_MARKER_CONFIG: CRCMarkerConfig = {
  segmentSize: 4,
  markerSize: 1,
  checkpoint: 3,
};

/**
 * Compute CRC-8 (polynomial 0x07, init 0x00, CCITT equivalent).
 * Used for marker bytes (1 byte each).
 */
export function crc8(data: Uint8Array): number {
  let crc = 0x00;
  for (let i = 0; i < data.length; i++) {
    crc ^= data[i];
    for (let j = 0; j < 8; j++) {
      if (crc & 0x80) {
        crc = ((crc << 1) ^ 0x07) & 0xff;
      } else {
        crc = (crc << 1) & 0xff;
      }
    }
  }
  return crc & 0xff;
}

/**
 * Insert CRC markers into a payload byte array.
 *
 * Input:  [p0, p1, p2, p3, p4, p5, p6, p7, ...] (raw payload)
 * Output: [p0, p1, p2, p3, CRC8(p0..p3), p4, p5, p6, p7, CRC8(p4..p7), ...]
 *
 * The output length = payload.length + ceil(payload.length / segmentSize) * markerSize.
 */
export function insertCRCMarkers(
  payload: Uint8Array,
  cfg: CRCMarkerConfig = DEFAULT_CRC_MARKER_CONFIG,
): Uint8Array {
  const numSegments = Math.ceil(payload.length / cfg.segmentSize);
  const outputLen = payload.length + numSegments * cfg.markerSize;
  const output = new Uint8Array(outputLen);

  let outIdx = 0;
  for (let seg = 0; seg < numSegments; seg++) {
    const segStart = seg * cfg.segmentSize;
    const segEnd = Math.min(segStart + cfg.segmentSize, payload.length);
    const segLen = segEnd - segStart;

    // Copy segment bytes
    output.set(payload.slice(segStart, segEnd), outIdx);
    outIdx += segLen;

    // Compute and insert marker
    const segBytes = payload.slice(segStart, segEnd);
    if (cfg.markerSize === 1) {
      output[outIdx] = crc8(segBytes);
      outIdx += 1;
    } else {
      // CRC-16 (2 bytes, big-endian)
      const crc = crc16(segBytes);
      output[outIdx] = (crc >> 8) & 0xff;
      output[outIdx + 1] = crc & 0xff;
      outIdx += 2;
    }
  }

  return output;
}

/**
 * Strip CRC markers from a payload byte array, returning the raw payload.
 * Also returns the list of segments where the CRC check FAILED (erasures).
 *
 * Used during decode to:
 *   1. Remove markers
 *   2. Detect which segments had errors (CRC mismatch)
 *   3. Pass failed segment indices as erasure hints to the RS decoder
 */
export function stripCRCMarkers(
  data: Uint8Array,
  cfg: CRCMarkerConfig = DEFAULT_CRC_MARKER_CONFIG,
): { payload: Uint8Array; failedSegments: number[]; segmentBoundaries: { start: number; end: number }[] } {
  const windowSize = cfg.segmentSize + cfg.markerSize;
  const numSegments = Math.ceil(data.length / windowSize);
  const payload = new Uint8Array(numSegments * cfg.segmentSize);
  const failedSegments: number[] = [];
  const segmentBoundaries: { start: number; end: number }[] = [];

  let payloadIdx = 0;
  for (let seg = 0; seg < numSegments; seg++) {
    const windowStart = seg * windowSize;
    const segBytes = data.slice(windowStart, windowStart + cfg.segmentSize);
    const markerBytes = data.slice(
      windowStart + cfg.segmentSize,
      windowStart + cfg.segmentSize + cfg.markerSize,
    );

    // Copy segment to payload
    payload.set(segBytes, payloadIdx);
    const segStart = payloadIdx;
    payloadIdx += cfg.segmentSize;
    segmentBoundaries.push({ start: segStart, end: payloadIdx });

    // Verify marker
    let expectedMarker: number;
    if (cfg.markerSize === 1) {
      expectedMarker = crc8(segBytes);
      if (markerBytes.length >= 1 && markerBytes[0] !== expectedMarker) {
        failedSegments.push(seg);
      }
    } else {
      const expectedCrc = crc16(segBytes);
      const expectedHi = (expectedCrc >> 8) & 0xff;
      const expectedLo = expectedCrc & 0xff;
      if (
        markerBytes.length >= 2 &&
        (markerBytes[0] !== expectedHi || markerBytes[1] !== expectedLo)
      ) {
        failedSegments.push(seg);
      }
    }
  }

  return { payload: payload.slice(0, payloadIdx), failedSegments, segmentBoundaries };
}

/**
 * Attempt to re-sync after an indel by hypothesizing base shifts.
 *
 * Given a read where a marker CRC failed, try shifting the read by ±1, ±2, ±3
 * positions and re-checking the NEXT marker. If a shift validates, the indel
 * is localized to that offset.
 *
 * This is a simplified version of DNA-Aeon's stack algorithm.
 *
 * @param read The full read DNA (after primer trimming)
 * @param failPos Position (in read coordinates) where the marker failed
 * @param cfg Marker config
 * @returns The hypothesized indel offset (positive = insertion, negative = deletion),
 *          or 0 if no valid shift found (declare erasure).
 */
export function resyncAfterIndel(
  read: Uint8Array,
  failPos: number,
  cfg: CRCMarkerConfig = DEFAULT_CRC_MARKER_CONFIG,
): number {
  const windowSize = cfg.segmentSize + cfg.markerSize;

  // Try shifts from -3 to +3 (skip 0, that's the original which failed)
  for (let shift = -3; shift <= 3; shift++) {
    if (shift === 0) continue;

    // Look at the NEXT window after the failed one, with the hypothesized shift
    const nextWindowStart = failPos + windowSize + shift;
    if (nextWindowStart < 0 || nextWindowStart + windowSize > read.length) continue;

    const nextSeg = read.slice(
      nextWindowStart,
      nextWindowStart + cfg.segmentSize,
    );
    const nextMarker = read.slice(
      nextWindowStart + cfg.segmentSize,
      nextWindowStart + cfg.segmentSize + cfg.markerSize,
    );

    let expected: number;
    if (cfg.markerSize === 1) {
      expected = crc8(nextSeg);
      if (nextMarker.length >= 1 && nextMarker[0] === expected) {
        return shift;
      }
    } else {
      const crc = crc16(nextSeg);
      if (
        nextMarker.length >= 2 &&
        nextMarker[0] === ((crc >> 8) & 0xff) &&
        nextMarker[1] === (crc & 0xff)
      ) {
        return shift;
      }
    }
  }

  return 0; // no valid shift found
}

/**
 * Compute the overhead of CRC markers.
 * Overhead = numMarkers * markerSize / (payload + numMarkers * markerSize)
 */
export function markerOverhead(
  payloadLength: number,
  cfg: CRCMarkerConfig = DEFAULT_CRC_MARKER_CONFIG,
): number {
  const numMarkers = Math.ceil(payloadLength / cfg.segmentSize);
  const markerBytes = numMarkers * cfg.markerSize;
  return markerBytes / (payloadLength + markerBytes);
}
