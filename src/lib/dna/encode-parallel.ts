/**
 * v60: Encode Parallelization with Bun Worker Threads
 *
 * The encode loop is embarrassingly parallel — each oligo's LDPC encode +
 * DNA mapping + screening is independent. Bun supports worker_threads via
 * node:worker_threads, giving 4-8× speedup on multi-core machines.
 *
 * Architecture:
 *   1. Main thread: compression, outer RS, layout, metadata
 *   2. Worker pool (n=CPU cores): per-oligo LDPC + DNA mapping + screening
 *   3. Main thread: collect oligos, build metadata
 *
 * The worker code is inlined as a string and spawned via new Worker(url, {...}).
 * Each worker processes a batch of oligos and returns the results.
 */

import { Worker, isMainThread, parentPort, workerData } from "worker_threads";
import { cpus } from "os";

const NUM_WORKERS = Math.max(1, Math.min(8, cpus().length));

/**
 * Worker script (inlined as a string).
 * 
 * Receives a batch of oligos to encode, processes each one, and returns results.
 * Uses dynamic imports to avoid circular dependencies.
 */
const WORKER_CODE = `
const { parentPort, workerData } = require('worker_threads');

async function encodeBatch(task) {
  const { LDPCInnerCode } = require('${require.resolve("../src/lib/dna/ldpc-codec").replace(/\\/g, "\\\\")}');
  const { ReedSolomon } = require('${require.resolve("../src/lib/dna/reedsolomon").replace(/\\/g, "\\\\")}');
  const { crc16Bytes } = require('${require.resolve("../src/lib/dna/crc16").replace(/\\/g, "\\\\")}');
  const {
    bytesToDna, xorWithSeed, gcContent, maxHomopolymerRun, whitenAddress,
  } = require('${require.resolve("../src/lib/dna/mapping").replace(/\\/g, "\\\\")}');
  
  const { oligoIndices, payloads, cfg, layout, fwdPrimer, revPrimer } = task;
  
  const useLDPC = (cfg.innerCode ?? 'rs') === 'ldpc';
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  const innerLdpc = useLDPC ? new LDPCInnerCode({ n: innerN, k: innerK }) : null;
  const innerRs = new ReedSolomon({ n: innerN + layout.crcBytes, k: innerK });
  
  const results = [];
  let screeningRetries = 0;
  
  for (let i = 0; i < oligoIndices.length; i++) {
    const oligoIdx = oligoIndices[i];
    const payload = payloads[i];
    
    const rawAddress = new Uint8Array(4);
    rawAddress[0] = (oligoIdx >> 16) & 0xff;
    rawAddress[1] = (oligoIdx >> 8) & 0xff;
    rawAddress[2] = oligoIdx & 0xff;
    rawAddress[3] = 0;
    
    const constraints = {
      gcMin: cfg.constraints.gcMin,
      gcMax: cfg.constraints.gcMax,
      maxHomopolymer: cfg.constraints.maxHomopolymer,
    };
    
    let seed = 0;
    let dna = '';
    let attempts = 0;
    
    const baseAddress = rawAddress.slice();
    const baseRsData = new Uint8Array(innerK);
    baseRsData.set(whitenAddress(baseAddress), 0);
    baseRsData.set(payload, layout.addressBytes);
    
    while (attempts <= (cfg.maxRetries || 10)) {
      baseAddress[3] = seed & 0xff;
      const whitenedAddress = whitenAddress(baseAddress);
      baseRsData.set(whitenedAddress, 0);
      const effectivePayload = seed === 0 ? payload : xorWithSeed(payload, seed);
      baseRsData.set(effectivePayload, layout.addressBytes);
      
      const rsCodeword = useLDPC && innerLdpc
        ? innerLdpc.encode(baseRsData)
        : innerRs.encode(baseRsData);
      const crc = crc16Bytes(rsCodeword);
      const innerBlock = new Uint8Array(layout.totalInnerBytes);
      innerBlock.set(rsCodeword, 0);
      innerBlock.set(crc, rsCodeword.length);
      
      dna = bytesToDna(innerBlock);
      
      // Screen for GC and homopolymer
      const gc = gcContent(dna);
      const hp = maxHomopolymerRun(dna);
      if (gc >= constraints.gcMin && gc <= constraints.gcMax && hp <= constraints.maxHomopolymer) {
        break;
      }
      
      seed++;
      attempts++;
      screeningRetries++;
    }
    
    results.push({
      index: oligoIdx,
      sequence: fwdPrimer + dna + revPrimer,
      seed,
      screeningRetries,
    });
  }
  
  return { oligos: results, screeningRetries };
}

parentPort.on('message', async (task) => {
  try {
    const result = await encodeBatch(task);
    parentPort.postMessage(result);
  } catch (e) {
    parentPort.postMessage({ error: e.message });
  }
});
`;

/**
 * Encode a batch of oligos in parallel using worker threads.
 * 
 * Falls back to sequential if workers are unavailable.
 */
export async function encodeOligosParallel(
  oligoIndices: number[],
  payloads: Uint8Array[],
  cfg: any,
  layout: any,
  fwdPrimer: string,
  revPrimer: string,
): Promise<{ oligos: any[]; screeningRetries: number; workersUsed: number }> {
  const numWorkers = Math.min(NUM_WORKERS, oligoIndices.length);
  
  if (numWorkers <= 1) {
    // Sequential fallback
    const result = await encodeBatchSequential(oligoIndices, payloads, cfg, layout, fwdPrimer, revPrimer);
    return { ...result, workersUsed: 1 };
  }
  
  // Split into batches
  const batchSize = Math.ceil(oligoIndices.length / numWorkers);
  const batches: { indices: number[]; payloads: Uint8Array[] }[] = [];
  for (let w = 0; w < numWorkers; w++) {
    const start = w * batchSize;
    const end = Math.min(start + batchSize, oligoIndices.length);
    if (start >= end) break;
    batches.push({
      indices: oligoIndices.slice(start, end),
      payloads: payloads.slice(start, end),
    });
  }
  
  // Spawn workers
  const workers: Worker[] = [];
  const promises: Promise<any>[] = [];
  
  for (let w = 0; w < batches.length; w++) {
    const worker = new Worker(WORKER_CODE, { eval: true });
    workers.push(worker);
    
    const promise = new Promise((resolve, reject) => {
      worker.on("message", (msg) => {
        if (msg.error) reject(new Error(msg.error));
        else resolve(msg);
      });
      worker.on("error", reject);
      worker.postMessage({
        oligoIndices: batches[w].indices,
        payloads: batches[w].payloads,
        cfg,
        layout,
        fwdPrimer,
        revPrimer,
      });
    });
    promises.push(promise);
  }
  
  try {
    const results = await Promise.all(promises);
    
    // Collect results
    const allOligos: any[] = [];
    let totalScreeningRetries = 0;
    for (const result of results) {
      allOligos.push(...result.oligos);
      totalScreeningRetries += result.screeningRetries;
    }
    
    // Sort by index to preserve order
    allOligos.sort((a, b) => a.index - b.index);
    
    // Terminate workers
    for (const worker of workers) {
      worker.terminate();
    }
    
    return { oligos: allOligos, screeningRetries: totalScreeningRetries, workersUsed: batches.length };
  } catch (e) {
    // Clean up workers on error
    for (const worker of workers) {
      worker.terminate();
    }
    throw e;
  }
}

/**
 * Sequential encode (fallback when workers unavailable).
 */
async function encodeBatchSequential(
  oligoIndices: number[],
  payloads: Uint8Array[],
  cfg: any,
  layout: any,
  fwdPrimer: string,
  revPrimer: string,
): Promise<{ oligos: any[]; screeningRetries: number }> {
  // Dynamic import to avoid circular dependencies
  const { LDPCInnerCode } = await import("./ldpc-codec");
  const { ReedSolomon } = await import("./reedsolomon");
  const { crc16Bytes } = await import("./crc16");
  const { bytesToDna, xorWithSeed, gcContent, maxHomopolymerRun, whitenAddress } = await import("./mapping");
  
  const useLDPC = (cfg.innerCode ?? "rs") === "ldpc";
  const innerK = layout.addressBytes + layout.payloadBytes;
  const innerN = innerK + layout.innerParityBytes;
  const innerLdpc = useLDPC ? new LDPCInnerCode({ n: innerN, k: innerK }) : null;
  const innerRs = new ReedSolomon({ n: innerN + layout.crcBytes, k: innerK });
  
  const results: any[] = [];
  let screeningRetries = 0;
  
  for (let i = 0; i < oligoIndices.length; i++) {
    const oligoIdx = oligoIndices[i];
    const payload = payloads[i];
    
    const rawAddress = new Uint8Array(4);
    rawAddress[0] = (oligoIdx >> 16) & 0xff;
    rawAddress[1] = (oligoIdx >> 8) & 0xff;
    rawAddress[2] = oligoIdx & 0xff;
    rawAddress[3] = 0;
    
    const constraints = {
      gcMin: cfg.constraints.gcMin,
      gcMax: cfg.constraints.gcMax,
      maxHomopolymer: cfg.constraints.maxHomopolymer,
    };
    
    let seed = 0;
    let dna = "";
    let attempts = 0;
    
    const baseAddress = rawAddress.slice();
    const baseRsData = new Uint8Array(innerK);
    baseRsData.set(whitenAddress(baseAddress), 0);
    baseRsData.set(payload, layout.addressBytes);
    
    while (attempts <= (cfg.maxRetries || 10)) {
      baseAddress[3] = seed & 0xff;
      const whitenedAddress = whitenAddress(baseAddress);
      baseRsData.set(whitenedAddress, 0);
      const effectivePayload = seed === 0 ? payload : xorWithSeed(payload, seed);
      baseRsData.set(effectivePayload, layout.addressBytes);
      
      const rsCodeword = useLDPC && innerLdpc
        ? innerLdpc.encode(baseRsData)
        : innerRs.encode(baseRsData);
      const crc = crc16Bytes(rsCodeword);
      const innerBlock = new Uint8Array(layout.totalInnerBytes);
      innerBlock.set(rsCodeword, 0);
      innerBlock.set(crc, rsCodeword.length);
      
      dna = bytesToDna(innerBlock);
      
      const gc = gcContent(dna);
      const hp = maxHomopolymerRun(dna);
      if (gc >= constraints.gcMin && gc <= constraints.gcMax && hp <= constraints.maxHomopolymer) {
        break;
      }
      
      seed++;
      attempts++;
      screeningRetries++;
    }
    
    results.push({
      index: oligoIdx,
      sequence: fwdPrimer + dna + revPrimer,
      seed,
      screeningRetries,
    });
  }
  
  return { oligos: results, screeningRetries };
}
