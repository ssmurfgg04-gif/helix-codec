/**
 * CRC-16/CCITT-FALSE checksum.
 *
 * Polynomial: 0x1021 (x^16 + x^12 + x^5 + 1)
 * Init: 0xFFFF, no reflection, no xor-out.
 *
 * Used per-oligo to detect residual errors after RS decoding.
 * Reference: CRC-16/CCITT-FALSE from CRC RevEng catalogue.
 */

export function crc16(data: Uint8Array): number {
  let crc = 0xffff;
  for (let i = 0; i < data.length; i++) {
    crc ^= (data[i] << 8) & 0xffff;
    for (let j = 0; j < 8; j++) {
      if (crc & 0x8000) {
        crc = ((crc << 1) ^ 0x1021) & 0xffff;
      } else {
        crc = (crc << 1) & 0xffff;
      }
    }
  }
  return crc & 0xffff;
}

/** Encode CRC-16 as 2 bytes (big-endian). */
export function crc16Bytes(data: Uint8Array): Uint8Array {
  const crc = crc16(data);
  return new Uint8Array([(crc >> 8) & 0xff, crc & 0xff]);
}

/** Verify CRC-16: returns true if data (with appended 2-byte CRC) is valid. */
export function verifyCrc16(data: Uint8Array): boolean {
  if (data.length < 2) return false;
  const payload = data.slice(0, data.length - 2);
  const expected = (data[data.length - 2] << 8) | data[data.length - 1];
  return crc16(payload) === expected;
}
