/**
 * Reed-Solomon encoder/decoder for DNA storage codec.
 *
 * Uses the proven `reedsolomon` npm package (ZXing port by cho45) for the core
 * error-correction logic. This library is battle-tested in production barcode
 * applications and supports up to RS(255, k) over GF(256).
 *
 * We add a thin wrapper that:
 *   1. Provides a clean TypeScript API matching our codec design.
 *   2. Tracks the number of corrections made.
 *   3. Adds a pure-erasure decoder for use in the outer across-strand code
 *      (where we know which strands are missing from the address index).
 *
 * Conventions:
 *   - Codeword is data-first, parity-last: out[0..k-1] = data, out[k..n-1] = parity.
 *   - Erasure positions are indices into the codeword array (0 = first byte).
 *
 * References:
 *   - ZXing ReedSolomonEncoder/Decoder (Apache 2.0)
 *   - Forney (1965) for error magnitudes in erasure decoding
 */

// reedsolomon is a CommonJS module (ZXing port). Use createRequire to avoid
// ESM/CJS interop issues when running with tsx or other ESM-first tools.
import { createRequire } from 'node:module';
const require2 = createRequire(import.meta.url ?? __filename);
const rsLib = require2("reedsolomon");
import { gfMul, gfDiv, gfInverse, gfPow } from "./gf256";
// v69: napi-rs native Reed-Solomon (FIRST PRIORITY for encode/parity)
import { getNativeAddon } from "./native/helix-napi";

export interface RSConfig {
  n: number; // total codeword length (data + parity), <= 255
  k: number; // data length, < n
}

export interface RSDecodeResult {
  data: Uint8Array;
  corrected: number; // number of unknown errors corrected
  erased: number; // number of erasures corrected (always 0 for library path)
}

export class ReedSolomon {
  readonly n: number;
  readonly k: number;
  readonly nsym: number; // number of parity symbols = n - k
  private encoder: any;
  private decoder: any;
  // Use our own GF(256) tables (matching reedsolomon's AZTEC_DATA_8: primitive 0x12D)
  // for the pure-erasure path. The reedsolomon lib uses primitive 0x12D too.
  readonly alpha = 2;
  readonly fcr = 1; // AZTEC_DATA_8 uses generatorBase=1

  constructor(cfg: RSConfig) {
    if (cfg.n > 255 || cfg.n <= 0) throw new Error(`RS n must be in 1..255, got ${cfg.n}`);
    if (cfg.k >= cfg.n || cfg.k <= 0) throw new Error(`RS k must be in 1..n-1, got ${cfg.k}`);
    this.n = cfg.n;
    this.k = cfg.k;
    this.nsym = cfg.n - cfg.k;
    const GF = rsLib.GenericGF.AZTEC_DATA_8();
    this.encoder = new rsLib.ReedSolomonEncoder(GF);
    this.decoder = new rsLib.ReedSolomonDecoder(GF);
  }

  /**
   * Encode k bytes -> n bytes (k data + nsym parity).
   * Output: out[0..k-1] = data, out[k..n-1] = parity.
   */
  encode(data: Uint8Array): Uint8Array {
    if (data.length !== this.k) {
      throw new Error(`RS encode expects ${this.k} bytes, got ${data.length}`);
    }
    // v69: napi-rs native FIRST PRIORITY — true FFI RS encode
    const addon = getNativeAddon();
    if (addon) {
      try {
        const result = addon.rsEncode(data, this.nsym);
        if (result && result.length === this.n) return result;
      } catch { /* fall through */ }
    }
    const msg = new Int32Array(this.n);
    for (let i = 0; i < this.k; i++) msg[i] = data[i];
    this.encoder.encode(msg, this.nsym);
    const out = new Uint8Array(this.n);
    for (let i = 0; i < this.n; i++) out[i] = msg[i] & 0xff;
    return out;
  }

  /** Encode only the parity bytes. */
  parity(data: Uint8Array): Uint8Array {
    const full = this.encode(data);
    return full.slice(this.k);
  }

  /**
   * Decode n bytes -> k bytes, correcting unknown errors.
   * Throws on uncorrectable error.
   *
   * NOTE: This path does NOT support erasures. Use decodeWithErasures() if you
   * know which positions are corrupted.
   */
  decode(recv: Uint8Array): RSDecodeResult {
    if (recv.length !== this.n) {
      throw new Error(`RS decode expects ${this.n} bytes, got ${recv.length}`);
    }
    const msg = new Int32Array(this.n);
    for (let i = 0; i < this.n; i++) msg[i] = recv[i];
    try {
      this.decoder.decode(msg, this.nsym);
    } catch (e) {
      throw new Error(`RS decode failed: ${(e as Error).message}`);
    }
    const out = new Uint8Array(this.k);
    for (let i = 0; i < this.k; i++) out[i] = msg[i] & 0xff;
    let corrected = 0;
    for (let i = 0; i < this.n; i++) if ((msg[i] & 0xff) !== recv[i]) corrected++;
    return { data: out, corrected, erased: 0 };
  }

  /**
   * Decode n bytes -> k bytes, using known erasure positions.
   *
   * This is a pure-erasure decoder: it can correct up to `nsym` erasures
   * (twice the error-correction capacity of the unknown-error decoder).
   * Used for the outer across-strand code where the address index tells us
   * which strands are missing.
   *
   * Algorithm:
   *   1. Compute syndromes S_i = r(alpha^(fcr+i)).
   *   2. Build erasure locator Lambda(x) = prod (1 - alpha^p * x).
   *   3. Compute Omega(x) = S(x) * Lambda(x) mod x^nsym.
   *   4. For each erasure position p, magnitude E_p = Omega(alpha^-p) / Lambda'(alpha^-p).
   *   5. Apply corrections.
   *
   * Returns: corrected data + count of erasures fixed.
   * Throws if too many erasures or any inconsistency.
   */
  decodeWithErasures(recv: Uint8Array, erasePos: number[]): RSDecodeResult {
    if (recv.length !== this.n) {
      throw new Error(`RS decode expects ${this.n} bytes, got ${recv.length}`);
    }
    const cleanPos = Array.from(new Set(erasePos)).filter((p) => p >= 0 && p < this.n);
    if (cleanPos.length > this.nsym) {
      throw new Error(`Too many erasures: ${cleanPos.length} > ${this.nsym}`);
    }
    if (cleanPos.length === 0) {
      // Fall back to standard decode (handles unknown errors)
      return this.decode(recv);
    }

    // Treat erased positions as 0 (so the "error" at position p is exactly the original value).
    // Then we can solve for the original values via Forney.
    const mutated = recv.slice();
    for (const p of cleanPos) mutated[p] = 0;

    // 1) Compute syndromes of the mutated codeword.
    // Codeword convention: msg[0] = highest power (BE).
    const synd = new Uint8Array(this.nsym);
    for (let i = 0; i < this.nsym; i++) {
      const x = gfPow(this.alpha, this.fcr + i);
      let s = 0;
      for (let j = 0; j < this.n; j++) s = gfMul(s, x) ^ mutated[j];
      synd[i] = s;
    }

    let hasError = false;
    for (let i = 0; i < this.nsym; i++) if (synd[i] !== 0) { hasError = true; break; }
    if (!hasError) {
      // No errors at all — recv was already a valid codeword (after zeroing erased positions).
      return { data: mutated.slice(0, this.k), corrected: 0, erased: cleanPos.length };
    }

    // 2) Build erasure locator Lambda(x) = prod (1 - alpha^p * x).
    // Note: convention. The error at BE position p contributes E_p * x^(n-1-p) to e(x).
    // So S_i = sum_p E_p * (alpha^(n-1-p))^i = sum_p E_p * X_p^i where X_p = alpha^(n-1-p).
    // Locator Lambda(x) = prod (1 - X_p * x), roots at x = X_p^-1 = alpha^-(n-1-p) = alpha^(p-n+1).
    // For BE position p, the "Chien index" is i = n-1-p, root alpha^-i, position return n-1-i = p.
    //
    // We'll use LE convention internally: poly[0] = constant.
    let LambdaLE = new Uint8Array([1]);
    for (const p of cleanPos) {
      const i = this.n - 1 - p; // Chien index
      const X = gfPow(this.alpha, i); // locator X = alpha^i
      // (1 - X*x) = (1 + X*x) in GF(2). In LE: [1, X]
      const factor = new Uint8Array([1, X]);
      LambdaLE = polyMulLE(LambdaLE, factor);
    }

    // 3) Compute Omega(x) = S(x) * Lambda(x) mod x^nsym. (LE convention)
    // S(x) = S_0 + S_1*x + S_2*x^2 + ...; LE: synd[i] = S_i.
    let OmegaLE = polyMulLE(synd, LambdaLE);
    if (OmegaLE.length > this.nsym) OmegaLE = OmegaLE.slice(0, this.nsym);

    // 4) Forney: E_p = Omega(X_p^-1) / Lambda'(X_p^-1)
    // Lambda'(x) is the formal derivative: drop even powers, shift odd powers down.
    const LambdaPrimeLE = formalDerivativeLE(LambdaLE);

    const corrected = mutated.slice();
    for (const p of cleanPos) {
      const i = this.n - 1 - p;
      const X = gfPow(this.alpha, i);
      const Xinv = gfInverse(X);
      const omegaVal = polyEvalLE(OmegaLE, Xinv);
      const lambdaPrimeVal = polyEvalLE(LambdaPrimeLE, Xinv);
      if (lambdaPrimeVal === 0) {
        throw new Error("Forney division by zero — erasure set is degenerate");
      }
      const magnitude = gfDiv(omegaVal, lambdaPrimeVal);
      corrected[p] = magnitude;
    }

    // 5) Verify by re-checking syndromes on the corrected codeword.
    for (let i = 0; i < this.nsym; i++) {
      const x = gfPow(this.alpha, this.fcr + i);
      let s = 0;
      for (let j = 0; j < this.n; j++) s = gfMul(s, x) ^ corrected[j];
      if (s !== 0) {
        throw new Error("Post-correction syndrome nonzero — erasure decoding failed");
      }
    }

    return { data: corrected.slice(0, this.k), corrected: 0, erased: cleanPos.length };
  }

  /** Check if a codeword has any errors (nonzero syndrome). */
  hasError(recv: Uint8Array): boolean {
    if (recv.length !== this.n) return true;
    for (let i = 0; i < this.nsym; i++) {
      const x = gfPow(this.alpha, this.fcr + i);
      let s = 0;
      for (let j = 0; j < this.n; j++) s = gfMul(s, x) ^ recv[j];
      if (s !== 0) return true;
    }
    return false;
  }
}

// --- LE polynomial helpers (for the pure-erasure path) ---

function polyMulLE(a: Uint8Array, b: Uint8Array): Uint8Array {
  const out = new Uint8Array(a.length + b.length - 1);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      out[i + j] ^= gfMul(a[i], b[j]);
    }
  }
  return out;
}

function polyEvalLE(poly: Uint8Array, x: number): number {
  let y = 0;
  let xPow = 1;
  for (let i = 0; i < poly.length; i++) {
    y ^= gfMul(poly[i], xPow);
    xPow = gfMul(xPow, x);
  }
  return y;
}

function formalDerivativeLE(poly: Uint8Array): Uint8Array {
  if (poly.length <= 1) return new Uint8Array(0);
  const out = new Uint8Array(poly.length - 1);
  // d/dx(sum c_i x^i) = sum c_i * i * x^(i-1). In GF(2^8), only odd i contribute (i mod 2 = 1).
  for (let i = 1; i < poly.length; i++) {
    if (i % 2 === 1) {
      out[i - 1] = poly[i];
    }
  }
  return out;
}

/** Pad/truncate data to fit an RS block of size k (shorter -> zero-padded). */
export function padToK(data: Uint8Array, k: number): Uint8Array {
  if (data.length === k) return data;
  if (data.length > k) throw new Error(`Data length ${data.length} exceeds k=${k}`);
  const padded = new Uint8Array(k);
  padded.set(data, 0);
  return padded;
}
