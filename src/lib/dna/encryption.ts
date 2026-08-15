/**
 * Encryption layer for BioArchive.
 *
 * Uses XChaCha20-Poly1305 (authenticated encryption with extended nonce)
 * via @noble/ciphers. Key derivation uses Argon2id → HKDF-SHA256:
 *   1. Argon2id (memory-hard KDF) for password stretching — makes brute-force
 *      attacks on low-entropy passwords expensive.
 *   2. HKDF-SHA256 for domain separation — binds the stretched key to this
 *      application so it cannot be reused elsewhere.
 * Argon2id is provided by @noble/hashes (pure-JS, no native binding needed).
 * If Argon2id is unavailable at runtime, falls back to PBKDF2-SHA256 with
 * 100,000 iterations.
 *
 * Pipeline:
 *   plaintext bytes
 *     → (optional) XChaCha20-Poly1305 encrypt with derived key
 *     → ciphertext bytes (24-byte nonce prepended + 16-byte auth tag)
 *
 * The encryption layer sits between compression and chunking, so:
 *   file → compress → encrypt → chunk → ECC → DNA encode
 *
 * Reference:
 *   - RFC 8439: ChaCha20-Poly1305
 *   - draft-irtf-cfrg-xchacha: XChaCha20-Poly1305 (extended nonce)
 *   - RFC 5869: HKDF
 *   - RFC 9106: Argon2
 */

import { xchacha20poly1305 } from "@noble/ciphers/chacha.js";
import { argon2id } from "@noble/hashes/argon2.js";
import { hkdf } from "@noble/hashes/hkdf.js";
import { pbkdf2 } from "@noble/hashes/pbkdf2.js";
import { sha256 } from "@noble/hashes/sha2.js";
import { randomBytes } from "@noble/hashes/utils.js";

export interface EncryptionResult {
  ciphertext: Uint8Array; // nonce(24) + ciphertext + tag(16)
  nonce: Uint8Array; // 24-byte nonce (also prepended to ciphertext)
  salt: Uint8Array; // 16-byte salt used for KDF
  keyId: string; // SHA-256 of derived key (hex, first 16 bytes)
}

export interface EncryptionConfig {
  password: string;
  // Optional fixed salt (for reproducibility). If not set, random.
  salt?: Uint8Array;
  // Optional fixed nonce (for reproducibility). If not set, random.
  nonce?: Uint8Array;
}

/**
 * Derive a 32-byte XChaCha20 key from a password using Argon2id → HKDF-SHA256.
 *
 * Step 1: Argon2id for password hardening (memory-hard KDF). This makes
 *         brute-force attacks on low-entropy passwords expensive.
 *         Parameters: m=65536 (64 MB), t=3 (3 passes), p=1 (single lane).
 *         Falls back to PBKDF2-SHA256 with 100,000 iterations if Argon2id
 *         throws at runtime (e.g. unsuitable environment).
 * Step 2: HKDF for domain separation (binds key to this application).
 */
export function deriveKey(password: string, salt: Uint8Array): Uint8Array {
  const ikm = new TextEncoder().encode(password);

  // Step 1: Password stretching
  let stretched: Uint8Array;
  try {
    stretched = argon2id(ikm, salt, {
      t: 3, // 3 iterations
      m: 65536, // 64 MB memory
      p: 1, // single lane
      dkLen: 32, // 32-byte output
    });
  } catch {
    // Fallback: PBKDF2-SHA256 with 100,000 iterations
    stretched = pbkdf2(sha256, ikm, salt, {
      c: 100_000,
      dkLen: 32,
    });
  }

  // Step 2: HKDF for domain separation (binds key to this application)
  const info = new TextEncoder().encode("bioarchive/v1/xchacha20poly1305");
  return hkdf(sha256, stretched, salt, info, 32);
}

/** Compute key ID (first 16 bytes of SHA-256 of the derived key, hex). */
export function computeKeyId(key: Uint8Array): string {
  const hash = sha256(key);
  return Array.from(hash.slice(0, 16))
    .map((b) => (b as number).toString(16).padStart(2, "0"))
    .join("");
}

/**
 * Encrypt plaintext using XChaCha20-Poly1305.
 * Returns ciphertext with nonce prepended (so the receiver needs only the key).
 */
export function encrypt(plaintext: Uint8Array, cfg: EncryptionConfig): EncryptionResult {
  const salt = cfg.salt ?? randomBytes(16);
  const key = deriveKey(cfg.password, salt);
  const nonce = cfg.nonce ?? randomBytes(24);
  const cipher = xchacha20poly1305(key, nonce);
  const sealed = cipher.encrypt(plaintext); // ciphertext + 16-byte tag
  // Prepend nonce to ciphertext
  const ciphertext = new Uint8Array(24 + sealed.length);
  ciphertext.set(nonce, 0);
  ciphertext.set(sealed, 24);
  return {
    ciphertext,
    nonce,
    salt,
    keyId: computeKeyId(key),
  };
}

/**
 * Decrypt ciphertext (with nonce prepended) using XChaCha20-Poly1305.
 * Throws if authentication fails (tag mismatch).
 */
export function decrypt(ciphertext: Uint8Array, password: string, salt: Uint8Array): Uint8Array {
  if (ciphertext.length < 24 + 16) {
    throw new Error("Ciphertext too short (need at least nonce + tag)");
  }
  const nonce = ciphertext.slice(0, 24);
  const sealed = ciphertext.slice(24);
  const key = deriveKey(password, salt);
  const cipher = xchacha20poly1305(key, nonce);
  const plaintext = cipher.decrypt(sealed);
  if (plaintext === null) {
    throw new Error("Decryption failed: authentication tag mismatch");
  }
  return plaintext;
}

/**
 * Encrypt a JSON-serializable object (used for encrypted metadata).
 * Returns a base64 string of (salt + ciphertext).
 */
export function encryptJson(obj: unknown, cfg: EncryptionConfig): string {
  const json = new TextEncoder().encode(JSON.stringify(obj));
  const result = encrypt(json, cfg);
  // Pack salt + ciphertext into a single buffer
  const packed = new Uint8Array(16 + result.ciphertext.length);
  packed.set(result.salt, 0);
  packed.set(result.ciphertext, 16);
  return bytesToBase64(packed);
}

/**
 * Decrypt a base64-encoded (salt + ciphertext) blob back to a JSON object.
 */
export function decryptJson<T>(packedB64: string, password: string): T {
  const packed = base64ToBytes(packedB64);
  const salt = packed.slice(0, 16);
  const ciphertext = packed.slice(16);
  const plaintext = decrypt(ciphertext, password, salt);
  const json = new TextDecoder().decode(plaintext);
  return JSON.parse(json) as T;
}

// --- Base64 helpers (for browser compatibility) ---

function bytesToBase64(bytes: Uint8Array): string {
  if (typeof Buffer !== "undefined") {
    return Buffer.from(bytes).toString("base64");
  }
  let binary = "";
  for (let i = 0; i < bytes.length; i++) binary += String.fromCharCode(bytes[i]);
  return btoa(binary);
}

function base64ToBytes(b64: string): Uint8Array {
  if (typeof Buffer !== "undefined") {
    return new Uint8Array(Buffer.from(b64, "base64"));
  }
  const binary = atob(b64);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return bytes;
}
