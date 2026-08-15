/**
 * Decode Worker — runs batch_decode_all on a chunk of reads
 *
 * This script runs in a Node.js worker_thread. It:
 *   1. Loads the WASM module
 *   2. Receives read data (as transferable ArrayBuffers)
 *   3. Runs batch_decode_all on its chunk
 *   4. Returns the per-oligo result array
 */

import { parentPort, workerData } from "worker_threads";

// @ts-ignore — require the WASM module (works in worker_threads)
const wasm = require("./wasm-pkg/helix_dna_wasm.js");

interface WorkerData {
  allReads: ArrayBuffer;
  readOffsets: ArrayBuffer;
  readLengths: ArrayBuffer;
  fwdPrimer: ArrayBuffer;
  revPrimer: ArrayBuffer;
  oligoCount: number;
  innerN: number;
  innerK: number;
  totalInnerBytes: number;
}

async function main() {
  const data = workerData as WorkerData;

  const allReads = new Uint8Array(data.allReads);
  const readOffsets = new Uint8Array(data.readOffsets);
  const readLengths = new Uint8Array(data.readLengths);
  const fwdPrimer = new Uint8Array(data.fwdPrimer);
  const revPrimer = new Uint8Array(data.revPrimer);

  try {
    const result = wasm.batch_decode_all(
      allReads,
      readOffsets,
      readLengths,
      fwdPrimer,
      revPrimer,
      data.oligoCount,
      data.innerN,
      data.innerK,
      data.totalInnerBytes,
    );

    // Send result back — transfer the buffer for zero-copy
    const resultBuf = result.buffer.slice(0);
    parentPort?.postMessage(
      { result: new Uint8Array(resultBuf) },
      [resultBuf],
    );
  } catch (e: any) {
    parentPort?.postMessage({ error: e.message || String(e) });
  }
}

main().catch((e) => {
  parentPort?.postMessage({ error: e.message || String(e) });
});
