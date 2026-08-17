/**
 * Node.js WASM loader for helix-dna-wasm — for use in tsx scripts and Node.js CLI tools.
 *
 * Loads the wasm-pack --target nodejs build which uses CommonJS require()
 * and doesn't need a bundler or fetch(). This avoids the ES module compatibility
 * issues that the --target web build has in Node.js.
 *
 * Usage:
 *   import { initWasmNode, viterbiK9Decode } from './wasm-node-loader';
 *   await initWasmNode();
 *   const decoded = viterbiK9Decode(recv, llr, numBits, 15, 15, 10);
 */

let wasm: any = null;

export async function initWasmNode(): Promise<void> {
  if (wasm) return;
  try {
    // wasm-pack --target nodejs produces a CJS module with a synchronous __wbg_init
    const wasmModule = await import("./wasm-node/helix_dna_wasm.js");
    // The nodejs target auto-initializes on require; no need to call default()
    wasm = wasmModule;
  } catch (e) {
    console.error("Failed to load Node.js WASM module:", e);
    throw e;
  }
}

export function isWasmNodeAvailable(): boolean {
  return wasm !== null;
}

// --- Viterbi Decoder ---

export function viterbiK9Decode(
  receivedBytes: Uint8Array,
  llrF64: Uint8Array,
  numInfoBits: number,
  maxDrift: number,
  insertionPenaltyX10: number,
  deletionPenaltyX10: number,
): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized — call initWasmNode() first");
  return wasm.viterbi_k9_decode(receivedBytes, llrF64, numInfoBits, maxDrift, insertionPenaltyX10, deletionPenaltyX10);
}

export function viterbiK7Decode(
  receivedBytes: Uint8Array,
  llrF64: Uint8Array,
  numInfoBits: number,
  maxDrift: number,
  insertionPenaltyX10: number,
  deletionPenaltyX10: number,
): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized — call initWasmNode() first");
  return wasm.viterbi_k7_decode(receivedBytes, llrF64, numInfoBits, maxDrift, insertionPenaltyX10, deletionPenaltyX10);
}

export function convK9Encode(infoBytes: Uint8Array): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.conv_k9_encode(infoBytes);
}

export function convK7Encode(infoBytes: Uint8Array): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.conv_k7_encode(infoBytes);
}

// --- GF(256) / RS ---

export function rsEncode(handle: number, data: Uint8Array): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.rs_encode(handle, data);
}

export function rsDecode(handle: number, recv: Uint8Array): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.rs_decode(handle, recv);
}

export function rsCreate(n: number, k: number): number {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.rs_create(n, k);
}

export function rsFree(handle: number): void {
  if (!wasm) throw new Error("Node WASM not initialized");
  wasm.rs_free(handle);
}

// --- Pack/Unpack ---

export function packDnaToBits(dna: string): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.pack_dna_to_bits(dna);
}

export function unpackBitsToDna(packed: Uint8Array, numBases: number): string {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.unpack_bits_to_dna(packed, numBases);
}

// --- Compress ---

export function arithCompress(data: Uint8Array): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.arith_compress(data);
}

export function arithDecompress(data: Uint8Array): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.arith_decompress(data);
}

// --- Simulate ---

export function simulateSingle(
  oligo: Uint8Array, synthSub: number, synthIns: number, synthDel: number,
  seqSub: number, seqIns: number, seqDel: number, seed: number,
): Uint8Array {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.simulate_single(oligo, synthSub, synthIns, synthDel, seqSub, seqIns, seqDel, seed);
}

// --- Version ---

export function version(): string {
  if (!wasm) throw new Error("Node WASM not initialized");
  return wasm.version();
}
