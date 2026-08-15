/**
 * Self-verifying address-payload binding for the decode pipeline.
 *
 * This module provides the decode-side verification that integrates with `decode.ts`.
 * When an oligo's address doesn't match its payload (after inner RS decode),
 * it is marked as an **erasure** for the outer RS decoder.
 *
 * ## Why this matters
 *
 * CRC-16 has a ~1/65536 false-positive rate — meaning roughly 1 in 65,536
 * corrupted payloads will pass the CRC check and appear valid. For large archives
 * with thousands of oligos, this is a real risk.
 *
 * Content-derived addressing eliminates this risk: the address is a BLAKE3 hash
 * of the payload. A corrupted payload would need to produce the **same BLAKE3
 * hash**, which has a collision probability of ~1/2^256 — negligible even for
 * petabyte-scale archives.
 *
 * ## Integration with decode pipeline
 *
 * ```
 * Sequencing reads → Inner RS/LDPC decode → [THIS MODULE: verify addresses]
 *   → Augment erasure list → Outer RS decode → Recover original data
 * ```
 *
 * Call `verifyAllAddressBindings()` or `verifyAndAugmentErasures()` **after**
 * inner RS/LDPC decode but **before** outer RS decode.
 *
 * @module addressing-verify
 */

import {
  type AddressingConfig,
  verifyAddressBinding,
  deriveAddress,
} from './addressing';

// ─── Types ────────────────────────────────────────────────────────────────────

/** Result of address-payload verification for a batch of oligos. */
export interface AddressVerificationResult {
  /** Whether ALL address-payload bindings are valid. */
  valid: boolean;
  /**
   * Indices of oligos where the address does NOT match the payload.
   * These should be treated as erasures in the outer RS decoder.
   */
  erasureIndices: number[];
  /**
   * Expected addresses (re-derived from payloads), indexed by oligo position
   * in the input array.
   */
  expectedAddresses: Uint8Array[];
  /**
   * Claimed addresses (from the oligos' address fields), indexed by oligo
   * position in the input array.
   */
  claimedAddresses: Uint8Array[];
}

// ─── Core verification ────────────────────────────────────────────────────────

/**
 * Verify address-payload bindings for all recovered oligos.
 *
 * For each oligo, this function:
 *   1. Re-derives the expected address from the recovered payload using BLAKE3.
 *   2. Compares it with the claimed address (from the oligo's address field)
 *      using constant-time comparison.
 *   3. If mismatch → marks the oligo as an erasure (its payload is corrupted
 *      despite passing CRC-16).
 *
 * **When to call**: AFTER inner RS/LDPC decode, BEFORE outer RS decode.
 *
 * **Why this catches what CRC-16 misses**:
 *   - CRC-16 false-positive rate: ~1/65536 ≈ 1.5×10⁻⁵
 *   - BLAKE3 collision rate: ~1/2^256 ≈ 8.6×10⁻⁷⁸
 *   - Improvement factor: ~2²⁴⁰ ≈ 10⁷²
 *
 * @param oligos - Array of recovered oligos, each with an index, payload, and address.
 * @param config - The addressing configuration used during encoding.
 * @returns Verification result with validity flag, erasure indices, and addresses.
 *
 * @example
 * ```ts
 * const result = verifyAllAddressBindings(recoveredOligos, config);
 * if (!result.valid) {
 *   console.warn(`${result.erasureIndices.length} oligos failed address verification`);
 *   // Pass result.erasureIndices to the outer RS decoder as erasures
 * }
 * ```
 */
export function verifyAllAddressBindings(
  oligos: { index: number; payload: Uint8Array; address: Uint8Array }[],
  config: AddressingConfig,
): AddressVerificationResult {
  if (oligos.length === 0) {
    return {
      valid: true,
      erasureIndices: [],
      expectedAddresses: [],
      claimedAddresses: [],
    };
  }

  const erasureIndices: number[] = [];
  const expectedAddresses: Uint8Array[] = [];
  const claimedAddresses: Uint8Array[] = [];

  for (const oligo of oligos) {
    const expectedAddr = deriveAddress(oligo.payload, config);
    const claimedAddr = oligo.address;

    expectedAddresses.push(expectedAddr);
    claimedAddresses.push(claimedAddr);

    const bindingValid = verifyAddressBinding(oligo.payload, claimedAddr, config);

    if (!bindingValid) {
      erasureIndices.push(oligo.index);
    }
  }

  return {
    valid: erasureIndices.length === 0,
    erasureIndices,
    expectedAddresses,
    claimedAddresses,
  };
}

/**
 * Integrate with `decode.ts`: given a map of recovered payloads and the original
 * oligo addresses, verify bindings and return an updated erasure list.
 *
 * This is the primary integration point for the decode pipeline. It:
 *   1. Iterates over all recovered payloads.
 *   2. Looks up the corresponding claimed address.
 *   3. Re-derives the expected address from the payload.
 *   4. If mismatch, adds the oligo index to the erasure list.
 *   5. Merges with any existing erasures (e.g., from CRC failures).
 *
 * The returned erasure list is deduplicated and sorted in ascending order,
 * which is the format expected by the outer RS decoder.
 *
 * @param payloads        - Map of oligo index → recovered payload bytes.
 * @param addresses       - Map of oligo index → claimed address bytes (from header).
 * @param existingErasures - Erasure indices already identified (e.g., CRC failures).
 * @param config          - The addressing configuration used during encoding.
 * @returns Augmented erasure list (sorted, deduplicated).
 *
 * @example
 * ```ts
 * // After inner RS decode:
 * const erasures = verifyAndAugmentErasures(
 *   recoveredPayloads,   // Map<number, Uint8Array>
 *   originalAddresses,   // Map<number, Uint8Array>
 *   crcErasures,        // number[]
 *   config,             // AddressingConfig
 * );
 * // Pass erasures to outer RS decoder
 * outerRsDecode(dataShards, parityShards, erasures);
 * ```
 */
export function verifyAndAugmentErasures(
  payloads: Map<number, Uint8Array>,
  addresses: Map<number, Uint8Array>,
  existingErasures: number[],
  config: AddressingConfig,
): number[] {
  // Use a Set for O(1) deduplication.
  const erasureSet = new Set<number>(existingErasures);

  // Verify each recovered payload against its claimed address.
  payloads.forEach((payload, index) => {
    const claimedAddress = addresses.get(index);

    if (claimedAddress === undefined) {
      // No address on record for this oligo — can't verify.
      // This can happen if the address field itself was lost to a harder erasure.
      // Trust the payload (it passed inner RS), but log a warning in debug mode.
      return;
    }

    const bindingValid = verifyAddressBinding(payload, claimedAddress, config);

    if (!bindingValid) {
      // Address mismatch → payload is corrupted despite passing inner RS + CRC.
      // Mark as erasure for the outer RS decoder to correct.
      erasureSet.add(index);
    }
  });

  // Also check for oligos that have addresses but no recovered payload.
  // These are already erasures (missing data), but ensure they're in the set.
  addresses.forEach((_addr, index) => {
    if (!payloads.has(index)) {
      erasureSet.add(index);
    }
  });

  // Sort in ascending order (expected by RS decoder).
  return Array.from(erasureSet).sort((a, b) => a - b);
}
