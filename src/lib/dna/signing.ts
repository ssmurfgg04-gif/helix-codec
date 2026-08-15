/**
 * Archive Signing (Ed25519)
 *
 * Signs the BioArchive manifest with an Ed25519 keypair so that recipients
 * can verify the archive was not tampered with and originated from the
 * claimed issuer.
 *
 * Uses @noble/curves for Ed25519 (RFC 8032). Keys are 32 bytes each;
 * signatures are 64 bytes.
 *
 * Pipeline:
 *   manifest.json → canonicalize (sorted keys) → SHA-256 → Ed25519 sign → signature
 *
 * Verification:
 *   signature + manifest + publicKey → valid/invalid
 *
 * Reference:
 *   - RFC 8032: Edwards-Curve Digital Signature Algorithm (EdDSA)
 *   - Bernstein et al. (2012). "High-speed high-security signatures."
 */

import { ed25519 } from "@noble/curves/ed25519.js";
import { sha256 } from "@noble/hashes/sha2.js";
import { randomBytes } from "@noble/hashes/utils.js";

export interface KeyPair {
  publicKey: Uint8Array; // 32 bytes
  privateKey: Uint8Array; // 32 bytes
}

/**
 * Generate a new Ed25519 keypair.
 */
export function generateKeyPair(): KeyPair {
  const privateKey = randomBytes(32);
  const publicKey = ed25519.getPublicKey(privateKey);
  return { publicKey, privateKey };
}

/**
 * Canonicalize a JSON object for signing (sorted keys, no whitespace).
 * This ensures the signature is reproducible regardless of key order.
 */
export function canonicalizeJson(obj: unknown): string {
  if (obj === null || typeof obj !== "object") {
    return JSON.stringify(obj);
  }
  if (Array.isArray(obj)) {
    return "[" + obj.map(canonicalizeJson).join(",") + "]";
  }
  const sortedKeys = Object.keys(obj as Record<string, unknown>).sort();
  const pairs = sortedKeys.map((k) => JSON.stringify(k) + ":" + canonicalizeJson((obj as Record<string, unknown>)[k]));
  return "{" + pairs.join(",") + "}";
}

/**
 * Sign a manifest (or any JSON-serializable object) with Ed25519.
 * Returns the 64-byte signature as hex.
 */
export function signManifest(obj: unknown, privateKey: Uint8Array): string {
  const canonical = canonicalizeJson(obj);
  const messageHash = sha256(new TextEncoder().encode(canonical));
  const signature = ed25519.sign(messageHash, privateKey);
  return Array.from(signature)
    .map((b) => (b as number).toString(16).padStart(2, "0"))
    .join("");
}

/**
 * Verify a manifest signature.
 * Returns true if the signature is valid for the given object + public key.
 */
export function verifyManifest(
  obj: unknown,
  signatureHex: string,
  publicKey: Uint8Array,
): boolean {
  try {
    const canonical = canonicalizeJson(obj);
    const messageHash = sha256(new TextEncoder().encode(canonical));
    const signature = new Uint8Array(
      signatureHex.match(/.{1,2}/g)?.map((b) => parseInt(b, 16)) ?? [],
    );
    if (signature.length !== 64) return false;
    return ed25519.verify(signature, messageHash, publicKey);
  } catch {
    return false;
  }
}

/**
 * Export a public key as hex string.
 */
export function publicKeyToHex(pk: Uint8Array): string {
  return Array.from(pk)
    .map((b) => (b as number).toString(16).padStart(2, "0"))
    .join("");
}

/**
 * Import a public key from hex string.
 */
export function hexToPublicKey(hex: string): Uint8Array {
  return new Uint8Array(
    hex.match(/.{1,2}/g)?.map((b) => parseInt(b, 16)) ?? [],
  );
}

/**
 * Sign a BioArchive manifest and attach the signature + public key.
 * Returns a new manifest object with `signature` and `signedBy` fields.
 */
export function signArchive<T extends Record<string, unknown>>(
  manifest: T,
  privateKey: Uint8Array,
): T & { signature: string; signedBy: string } {
  const publicKey = ed25519.getPublicKey(privateKey);
  const signature = signManifest(manifest, privateKey);
  return {
    ...manifest,
    signature,
    signedBy: publicKeyToHex(publicKey),
  };
}

/**
 * Verify a signed BioArchive manifest.
 * Strips `signature` and `signedBy` fields before verifying.
 */
export function verifyArchive<T extends Record<string, unknown>>(
  manifest: T & { signature?: string; signedBy?: string },
): boolean {
  if (!manifest.signature || !manifest.signedBy) return false;
  const { signature, signedBy, ...rest } = manifest;
  const publicKey = hexToPublicKey(signedBy!);
  return verifyManifest(rest, signature!, publicKey);
}
