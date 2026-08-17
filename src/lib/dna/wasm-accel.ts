// WASM-accelerated DNA codec — Rust + wasm-pack
// Drop-in replacement for pure-JS GF(256) + RS + DNA mapping.

let wasm: any = null;

export async function initWasm(): Promise<void> {
  if (wasm) return;
  const wasmModule = await import("./pkg/helix_dna_wasm.js");
  await wasmModule.default();
  wasm = wasmModule;
  wasm.init_gf();
}

export function isWasmAvailable(): boolean {
  return wasm !== null;
}

// --- GF(256) ---

export function rsEncode(data: Uint8Array, n: number, k: number): Uint8Array {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.rs_encode(data, n, k);
}

export function rsParity(data: Uint8Array, n: number, k: number): Uint8Array {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.rs_parity(data, n, k);
}

export function rsDecodeErasure(
  recv: Uint8Array,
  erasePos: number[] | Uint32Array,
  n: number,
  k: number,
): Uint8Array {
  if (!wasm) throw new Error("WASM not initialized");
  const eraseArr = erasePos instanceof Uint32Array ? Array.from(erasePos) : erasePos;
  return wasm.rs_decode_erasure(recv, eraseArr, n, k);
}

// --- DNA Mapping ---

export function bytesToDnaWasm(data: Uint8Array): string {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.bytes_to_dna(data);
}

export function dnaToBytesWasm(dna: string): Uint8Array {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.dna_to_bytes(dna);
}

export function gcContentWasm(dna: string): number {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.gc_content(dna);
}

export function maxHomopolymerRunWasm(dna: string): number {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.max_homopolymer_run(dna);
}

export function satisfiesConstraintsWasm(
  dna: string,
  gcMin: number,
  gcMax: number,
  maxHp: number,
): boolean {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.satisfies_constraints(dna, gcMin, gcMax, maxHp);
}

// --- Viterbi Decoder (Rust hot path) ---

/**
 * K=9 Indel-Tolerant Viterbi decode via Rust WASM (~5ms vs ~800ms in JS).
 *
 * @param receivedBytes  Received bytes (hard decisions)
 * @param llrF64  Packed LLR values as IEEE 754 f64 LE bytes (8 bytes per LLR, one per bit).
 *               Pass empty Uint8Array for hard-decision decoding.
 * @param numInfoBits  Number of information bits to decode
 * @param maxDrift  Maximum net indel drift to track (default 15)
 * @param insertionPenaltyX10  Insertion penalty × 10 (e.g. 15 for 1.5)
 * @param deletionPenaltyX10  Deletion penalty × 10 (e.g. 10 for 1.0)
 * @returns Decoded bytes
 */
export function viterbiK9Decode(
  receivedBytes: Uint8Array,
  llrF64: Uint8Array,
  numInfoBits: number,
  maxDrift: number,
  insertionPenaltyX10: number,
  deletionPenaltyX10: number,
): Uint8Array {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.viterbi_k9_decode(receivedBytes, llrF64, numInfoBits, maxDrift, insertionPenaltyX10, deletionPenaltyX10);
}

/**
 * K=7 (Voyager) Indel-Tolerant Viterbi decode via Rust WASM.
 */
export function viterbiK7Decode(
  receivedBytes: Uint8Array,
  llrF64: Uint8Array,
  numInfoBits: number,
  maxDrift: number,
  insertionPenaltyX10: number,
  deletionPenaltyX10: number,
): Uint8Array {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.viterbi_k7_decode(receivedBytes, llrF64, numInfoBits, maxDrift, insertionPenaltyX10, deletionPenaltyX10);
}

/**
 * K=9 NASA convolutional encode via Rust WASM.
 */
export function convK9Encode(infoBytes: Uint8Array): Uint8Array {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.conv_k9_encode(infoBytes);
}

/**
 * K=7 Voyager convolutional encode via Rust WASM.
 */
export function convK7Encode(infoBytes: Uint8Array): Uint8Array {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.conv_k7_encode(infoBytes);
}

// --- Benchmarks ---

export function benchRsEncode(iterations: number, n: number, k: number): number {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.bench_rs_encode(iterations, n, k);
}

export function benchDnaMapping(size: number, iterations: number): number {
  if (!wasm) throw new Error("WASM not initialized");
  return wasm.bench_dna_mapping(size, iterations);
}
