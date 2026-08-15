/**
 * Post-Quantum Signatures — ML-DSA (Dilithium)
 *
 * ML-DSA (Module-Lattice-Based Digital Signature Algorithm) is a post-quantum
 * signature scheme standardized by NIST in FIPS 204. It provides security
 * against quantum computer attacks, making it suitable for DNA archives that
 * must remain secure for centuries.
 *
 * DNA archives stored for 100+ years will outlive the expected arrival of
 * quantum computers (2030-2050). Ed25519 (used in v3.0) is broken by Shor's
 * algorithm, but ML-DSA remains secure.
 *
 * This module provides an adapter interface for ML-DSA. The actual
 * implementation uses @noble/post-quantum (pure JS) or falls back to Ed25519
 * if post-quantum libs are unavailable.
 *
 * Reference:
 *   - NIST FIPS 204 (2024). "Module-Lattice-Based Digital Signature Standard."
 *   - Lyubashevsky et al. (2018). "CRYSTALS-Dilithium."
 *     github.com/pq-crystals/dilithium
 *   - @noble/post-quantum: github.com/paulmillr/noble-post-quantum
 */

import { ed25519 } from "@noble/curves/ed25519.js";
import { sha256 } from "@noble/hashes/sha2.js";
import { randomBytes } from "@noble/hashes/utils.js";

export type SignatureAlgorithm = "ed25519" | "ml-dsa-65";

export interface PostQuantumKeyPair {
  publicKey: Uint8Array;
  privateKey: Uint8Array;
  algorithm: SignatureAlgorithm;
}

export interface SignedArchive {
  data: unknown;
  signature: string; // hex
  signedBy: string; // public key hex
  algorithm: SignatureAlgorithm;
}

/**
 * Generate a keypair for the specified algorithm.
 * Falls back to Ed25519 if ML-DSA is unavailable.
 */
export function generatePQKeyPair(
  algorithm: SignatureAlgorithm = "ml-dsa-65",
): PostQuantumKeyPair {
  if (algorithm === "ml-dsa-65") {
    try {
      // Try to load @noble/post-quantum
       
      const mlDsa = require("@noble/post-quantum/ml-dsa");
      const keys = mlDsa.mlDsa65.keypair(randomBytes(32));
      return {
        publicKey: keys.publicKey,
        privateKey: keys.secretKey,
        algorithm: "ml-dsa-65",
      };
    } catch {
      // Fall back to Ed25519
      console.warn("ML-DSA unavailable, falling back to Ed25519");
      const privateKey = randomBytes(32);
      const publicKey = ed25519.getPublicKey(privateKey);
      return { publicKey, privateKey, algorithm: "ed25519" };
    }
  }

  // Ed25519
  const privateKey = randomBytes(32);
  const publicKey = ed25519.getPublicKey(privateKey);
  return { publicKey, privateKey, algorithm: "ed25519" };
}

/**
 * Canonicalize JSON for signing (sorted keys).
 */
function canonicalizeJson(obj: unknown): string {
  if (obj === null || typeof obj !== "object") return JSON.stringify(obj);
  if (Array.isArray(obj)) return "[" + obj.map(canonicalizeJson).join(",") + "]";
  const sortedKeys = Object.keys(obj as Record<string, unknown>).sort();
  const pairs = sortedKeys.map(
    (k) => JSON.stringify(k) + ":" + canonicalizeJson((obj as Record<string, unknown>)[k]),
  );
  return "{" + pairs.join(",") + "}";
}

/**
 * Sign data with post-quantum signature.
 */
export function pqSign(
  data: unknown,
  privateKey: Uint8Array,
  algorithm: SignatureAlgorithm = "ml-dsa-65",
): string {
  const canonical = canonicalizeJson(data);
  const messageHash = sha256(new TextEncoder().encode(canonical));

  if (algorithm === "ml-dsa-65") {
    try {
       
      const mlDsa = require("@noble/post-quantum/ml-dsa");
      const signature = mlDsa.mlDsa65.sign(privateKey, messageHash);
      return bytesToHex(signature);
    } catch {
      // Fall back to Ed25519
      const signature = ed25519.sign(messageHash, privateKey);
      return bytesToHex(signature);
    }
  }

  // Ed25519
  const signature = ed25519.sign(messageHash, privateKey);
  return bytesToHex(signature);
}

/**
 * Verify a post-quantum signature.
 */
export function pqVerify(
  data: unknown,
  signatureHex: string,
  publicKey: Uint8Array,
  algorithm: SignatureAlgorithm = "ml-dsa-65",
): boolean {
  try {
    const canonical = canonicalizeJson(data);
    const messageHash = sha256(new TextEncoder().encode(canonical));
    const signature = hexToBytes(signatureHex);

    if (algorithm === "ml-dsa-65") {
      try {
         
        const mlDsa = require("@noble/post-quantum/ml-dsa");
        return mlDsa.mlDsa65.verify(publicKey, messageHash, signature);
      } catch {
        // Fall back to Ed25519 verification
        return ed25519.verify(signature, messageHash, publicKey);
      }
    }

    return ed25519.verify(signature, messageHash, publicKey);
  } catch {
    return false;
  }
}

/**
 * Sign an archive with post-quantum signature.
 */
export function signArchivePQ<T extends Record<string, unknown>>(
  data: T,
  keyPair: PostQuantumKeyPair,
): SignedArchive {
  const signature = pqSign(data, keyPair.privateKey, keyPair.algorithm);
  return {
    data,
    signature,
    signedBy: bytesToHex(keyPair.publicKey),
    algorithm: keyPair.algorithm,
  };
}

/**
 * Verify a signed archive.
 */
export function verifyArchivePQ(signed: SignedArchive): boolean {
  if (!signed.signature || !signed.signedBy) return false;
  const { signature, signedBy, algorithm, data } = signed;
  const publicKey = hexToBytes(signedBy);
  return pqVerify(data, signature, publicKey, algorithm);
}

function bytesToHex(bytes: Uint8Array): string {
  return Array.from(bytes)
    .map((b) => b.toString(16).padStart(2, "0"))
    .join("");
}

function hexToBytes(hex: string): Uint8Array {
  return new Uint8Array(
    hex.match(/.{1,2}/g)?.map((b) => parseInt(b, 16)) ?? [],
  );
}

/**
 * Get the security level of an algorithm.
 */
export function getSecurityLevel(algorithm: SignatureAlgorithm): {
  classicalBits: number;
  quantumBits: number;
  publicKeyBytes: number;
  signatureBytes: number;
  quantumSafe: boolean;
} {
  if (algorithm === "ml-dsa-65") {
    return {
      classicalBits: 128,
      quantumBits: 128,
      publicKeyBytes: 1952,
      signatureBytes: 3293,
      quantumSafe: true,
    };
  }
  return {
    classicalBits: 128,
    quantumBits: 0, // Ed25519 is NOT quantum-safe
    publicKeyBytes: 32,
    signatureBytes: 64,
    quantumSafe: false,
  };
}
