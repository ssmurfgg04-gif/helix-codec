/**
 * CRC-32 (ISO-HDLC) checksum.
 *
 * Polynomial: 0x04C11DB7 (reflected: 0xEDB88320)
 * Init: 0xFFFFFFFF, final XOR: 0xFFFFFFFF, reflected input/output.
 *
 * This is the same CRC-32 used in zlib, PNG, and Ethernet.
 * Used by the Mahoraga codec (Banal 2026) for per-sequence integrity.
 *
 * Reference:
 *   - ISO/IEC 13239:2002 (HDLC procedures)
 *   - Castagnoli (1993). "Optimization of cyclic redundancy-check codes
 *     with 24 and 32 parity bits."
 */

const POLY = 0xedb88320; // reflected polynomial

let table: Uint32Array | null = null;

function buildTable(): Uint32Array {
  const tbl = new Uint32Array(256);
  for (let i = 0; i < 256; i++) {
    let crc = i;
    for (let j = 0; j < 8; j++) {
      if (crc & 1) {
        crc = (crc >>> 1) ^ POLY;
      } else {
        crc = crc >>> 1;
      }
    }
    tbl[i] = crc >>> 0;
  }
  return tbl;
}

function getTable(): Uint32Array {
  if (!table) table = buildTable();
  return table;
}

/** Compute CRC-32 of a byte array. Returns a 32-bit unsigned integer. */
export function crc32(data: Uint8Array): number {
  const tbl = getTable();
  let crc = 0xffffffff;
  for (let i = 0; i < data.length; i++) {
    crc = (crc >>> 8) ^ tbl[(crc ^ data[i]) & 0xff];
  }
  return (crc ^ 0xffffffff) >>> 0;
}

/** Encode CRC-32 as 4 bytes (little-endian). */
export function crc32Bytes(data: Uint8Array): Uint8Array {
  const crc = crc32(data);
  return new Uint8Array([
    crc & 0xff,
    (crc >>> 8) & 0xff,
    (crc >>> 16) & 0xff,
    (crc >>> 24) & 0xff,
  ]);
}

/** Verify CRC-32: data must have 4-byte CRC appended (little-endian). */
export function verifyCrc32(data: Uint8Array): boolean {
  if (data.length < 4) return false;
  const payload = data.slice(0, data.length - 4);
  const expected =
    (data[data.length - 4]) |
    (data[data.length - 3] << 8) |
    (data[data.length - 2] << 16) |
    ((data[data.length - 1] << 24) >>> 0);
  const actual = crc32(payload) >>> 0;
  return actual === expected;
}
