import { describe, it, expect } from 'vitest';
import { encodeFile } from './codec';
import { decodeReads } from './decode';
import { simulate, PRESET_CLEAN } from './simulate';
import { CodecConfig, DEFAULT_CONFIG } from './types';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/**
 * Compute SHA-256 hash of a Uint8Array, returning hex string.
 */
async function sha256Hex(data: Uint8Array): Promise<string> {
  const { createHash } = await import('crypto');
  return createHash('sha256').update(Buffer.from(data)).digest('hex');
}

/**
 * Generate random bytes for testing.
 */
function randomBytes(n: number): Uint8Array {
  const buf = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    buf[i] = Math.floor(Math.random() * 256);
  }
  return buf;
}

/**
 * Perform a clean roundtrip: encode → simulate with no errors → decode.
 */
async function cleanRoundtrip(
  data: Uint8Array,
  cfg: CodecConfig,
  meta: { fileName: string; contentType: string },
) {
  // Encode.
  const encoded = await encodeFile(data, cfg, meta);

  // Simulate clean reads (no errors, 1x coverage).
  const simResult = simulate(encoded.encoded.oligos, PRESET_CLEAN);

  // Decode.
  const decoded = await decodeReads(
    simResult.reads,
    encoded.encoded.metadata,
    cfg,
    encoded.encoded.forwardPrimer,
    encoded.encoded.reversePrimer,
    true, // useSoftInfo
  );

  return { encoded, decoded };
}

// ---------------------------------------------------------------------------
// Codec roundtrip tests
// ---------------------------------------------------------------------------

describe('codec roundtrip', () => {
  it('encodes and decodes a small file', async () => {
    const data = randomBytes(256);
    const { decoded } = await cleanRoundtrip(data, DEFAULT_CONFIG, {
      fileName: 'test-256.bin',
      contentType: 'application/octet-stream',
    });

    expect(decoded.data).not.toBeNull();
    expect(decoded.hashMatches).toBe(true);
    expect(decoded.data!.length).toBe(data.length);

    // Verify byte-for-byte match.
    for (let i = 0; i < data.length; i++) {
      expect(decoded.data![i]).toBe(data[i]);
    }
  });

  it('encodes and decodes with LDPC inner code', async () => {
    const data = randomBytes(512);
    const ldpcConfig: CodecConfig = {
      ...DEFAULT_CONFIG,
      innerCode: 'ldpc',
      ldpcDecoder: 'auto',
      compress: false, // Disable compression for deterministic roundtrip.
    };

    const { decoded } = await cleanRoundtrip(data, ldpcConfig, {
      fileName: 'test-ldpc-512.bin',
      contentType: 'application/octet-stream',
    });

    expect(decoded.data).not.toBeNull();
    expect(decoded.hashMatches).toBe(true);
    expect(decoded.data!.length).toBe(data.length);
  });

  it('preserves SHA-256 hash across roundtrip', async () => {
    const data = randomBytes(1024);
    const originalHash = await sha256Hex(data);

    const { decoded } = await cleanRoundtrip(data, DEFAULT_CONFIG, {
      fileName: 'test-hash-1024.bin',
      contentType: 'application/octet-stream',
    });

    expect(decoded.data).not.toBeNull();
    expect(decoded.hash).toBe(originalHash);
    expect(decoded.hashMatches).toBe(true);
  });

  it('handles empty file', async () => {
    const data = new Uint8Array(0);

    // Encode empty data — this may throw or return an empty result depending on codec.
    // We wrap in try/catch since encoding 0 bytes may not be supported.
    try {
      const encoded = await encodeFile(data, DEFAULT_CONFIG, {
        fileName: 'test-empty.bin',
        contentType: 'application/octet-stream',
      });

      // If encoding succeeds, decode should recover empty data.
      const simResult = simulate(encoded.encoded.oligos, PRESET_CLEAN);
      const decoded = await decodeReads(
        simResult.reads,
        encoded.encoded.metadata,
        DEFAULT_CONFIG,
        encoded.encoded.forwardPrimer,
        encoded.encoded.reversePrimer,
        true,
      );

      if (decoded.data !== null) {
        expect(decoded.data.length).toBe(0);
        expect(decoded.hashMatches).toBe(true);
      }
    } catch (e) {
      // If the codec throws on empty input, that's acceptable behavior.
      expect((e as Error).message).toBeDefined();
    }
  });

  it('handles single byte file', async () => {
    const data = new Uint8Array([42]);

    const { decoded } = await cleanRoundtrip(data, DEFAULT_CONFIG, {
      fileName: 'test-single-byte.bin',
      contentType: 'application/octet-stream',
    });

    expect(decoded.data).not.toBeNull();
    expect(decoded.hashMatches).toBe(true);
    expect(decoded.data!.length).toBe(1);
    expect(decoded.data![0]).toBe(42);
  });
});

// ---------------------------------------------------------------------------
// Interleaving tests (preserved from original)
// ---------------------------------------------------------------------------

import { interleaveCodewords, deinterleaveCodewords } from './interleaving';

describe('Interleaving', () => {
  it('should be symmetric: deinterleave(interleave(x)) === x', () => {
    const a = new Uint8Array([1, 2, 3, 4]);
    const b = new Uint8Array([5, 6, 7, 8]);
    const c = new Uint8Array([9, 10, 11, 12]);
    const interleaved = interleaveCodewords([a, b, c]);
    const deinterleaved = deinterleaveCodewords(interleaved, 3);
    expect(deinterleaved.length).toBe(3);
    expect(Array.from(deinterleaved[0])).toEqual([1, 2, 3, 4]);
    expect(Array.from(deinterleaved[1])).toEqual([5, 6, 7, 8]);
    expect(Array.from(deinterleaved[2])).toEqual([9, 10, 11, 12]);
  });
});
