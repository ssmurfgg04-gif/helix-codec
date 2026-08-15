// Benchmark WASM vs JS LDPC decode speed.
import * as wasm from "../src/lib/dna/wasm-decode";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";

async function main() {
  console.log("=== WASM vs JS LDPC Benchmark ===\n");

  wasm.ensureWasm();
  console.log("WASM initialized\n");

  const n = 62, k = 58; // Same as production config

  // Create JS and WASM LDPC codes
  const jsLdpc = new LDPCInnerCode({ n, k });
  const wasmLdpc = new wasm.WasmLdpcCode(n, k);

  // Create test data
  const info = new Uint8Array(k);
  for (let i = 0; i < k; i++) info[i] = (i * 31 + 17) & 0xff;
  const codeword = jsLdpc.encode(info);

  // Verify WASM gives same result
  const wasmCw = wasmLdpc.encode(info);
  let match = true;
  for (let i = 0; i < n; i++) {
    if (codeword[i] !== wasmCw[i]) { match = false; break; }
  }
  console.log("Encode match:", match ? "✅" : "❌");

  // Benchmark JS encode
  const N = 100000;
  const t0 = Date.now();
  for (let i = 0; i < N; i++) jsLdpc.encode(info);
  const jsEncMs = Date.now() - t0;
  console.log(`JS encode:   ${N} ops in ${jsEncMs}ms (${(N / (jsEncMs / 1000)).toFixed(0)} ops/s)`);

  // Benchmark WASM encode
  const t1 = Date.now();
  for (let i = 0; i < N; i++) wasmLdpc.encode(info);
  const wasmEncMs = Date.now() - t1;
  console.log(`WASM encode: ${N} ops in ${wasmEncMs}ms (${(N / (wasmEncMs / 1000)).toFixed(0)} ops/s)`);
  console.log(`Encode speedup: ${(jsEncMs / wasmEncMs).toFixed(1)}x\n`);

  // Benchmark JS decode (zero-syndrome = fast path)
  const t2 = Date.now();
  for (let i = 0; i < N; i++) { try { jsLdpc.decode(codeword); } catch {} }
  const jsDecMs = Date.now() - t2;
  console.log(`JS decode:   ${N} ops in ${jsDecMs}ms (${(N / (jsDecMs / 1000)).toFixed(0)} ops/s)`);

  // Benchmark WASM decode
  const t3 = Date.now();
  for (let i = 0; i < N; i++) { try { wasmLdpc.decode(codeword); } catch {} }
  const wasmDecMs = Date.now() - t3;
  console.log(`WASM decode: ${N} ops in ${wasmDecMs}ms (${(N / (wasmDecMs / 1000)).toFixed(0)} ops/s)`);
  console.log(`Decode speedup: ${(jsDecMs / wasmDecMs).toFixed(1)}x\n`);

  // Benchmark DNA→bytes
  const { bytesToDna } = await import("../src/lib/dna/mapping");
  const dna = bytesToDna(codeword);
  const t4 = Date.now();
  for (let i = 0; i < N; i++) wasm.wasmDnaToBytes(dna);
  const wasmDnaMs = Date.now() - t4;
  console.log(`WASM dnaToBytes: ${N} ops in ${wasmDnaMs}ms`);

  const { dnaToBytes } = await import("../src/lib/dna/mapping");
  const t5 = Date.now();
  for (let i = 0; i < N; i++) dnaToBytes(dna);
  const jsDnaMs = Date.now() - t5;
  console.log(`JS dnaToBytes:   ${N} ops in ${jsDnaMs}ms`);
  console.log(`DNA speedup: ${(jsDnaMs / wasmDnaMs).toFixed(1)}x\n`);

  // Benchmark CRC-16
  const { crc16Bytes } = await import("../src/lib/dna/crc16");
  const t6 = Date.now();
  for (let i = 0; i < N; i++) wasm.wasmCrc16Bytes(codeword);
  const wasmCrcMs = Date.now() - t6;
  console.log(`WASM CRC-16: ${N} ops in ${wasmCrcMs}ms`);

  const t7 = Date.now();
  for (let i = 0; i < N; i++) crc16Bytes(codeword);
  const jsCrcMs = Date.now() - t7;
  console.log(`JS CRC-16:   ${N} ops in ${jsCrcMs}ms`);
  console.log(`CRC speedup: ${(jsCrcMs / wasmCrcMs).toFixed(1)}x\n`);

  // Summary
  console.log("=== SUMMARY ===");
  console.log(`Encode: ${jsEncMs}ms → ${wasmEncMs}ms = ${(jsEncMs / wasmEncMs).toFixed(1)}x faster`);
  console.log(`Decode: ${jsDecMs}ms → ${wasmDecMs}ms = ${(jsDecMs / wasmDecMs).toFixed(1)}x faster`);
  console.log(`DNA:    ${jsDnaMs}ms → ${wasmDnaMs}ms = ${(jsDnaMs / wasmDnaMs).toFixed(1)}x faster`);
  console.log(`CRC:    ${jsCrcMs}ms → ${wasmCrcMs}ms = ${(jsCrcMs / wasmCrcMs).toFixed(1)}x faster`);
}

main().catch((e) => { console.error(e); process.exit(1); });
