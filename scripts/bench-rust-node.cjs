/**
 * Direct Rust WASM benchmark — Node.js target.
 */

const fs = require("fs");
const path = require("path");

async function main() {
  console.log("Helix Rust WASM Direct Benchmark (Node.js)");
  console.log("═".repeat(50));

  // Load the Node.js WASM module
  const wasmPkgPath = path.join(__dirname, "..", "src", "lib", "dna", "wasm-pkg-rust-node", "helix_dna_wasm.js");
  if (!fs.existsSync(wasmPkgPath)) {
    console.log("WASM Node.js module not found. Building...");
    console.log("Run: cd rust/helix-dna-wasm && wasm-pack build --target nodejs --release --out-dir ../../src/lib/dna/wasm-pkg-rust-node");
    return;
  }

  const wasmPkg = await import(wasmPkgPath);
  console.log(`WASM version: ${wasmPkg.version()}`);

  const wasmPath = path.join(__dirname, "..", "src", "lib", "dna", "wasm-pkg-rust-node", "helix_dna_wasm_bg.wasm");
  if (fs.existsSync(wasmPath)) {
    console.log(`WASM size: ${(fs.statSync(wasmPath).size / 1024).toFixed(1)} KB`);
  }

  // -----------------------------------------------------------------------
  // Pack/Unpack benchmark
  // -----------------------------------------------------------------------
  console.log("\n=== Pack/Unpack ===");

  const sizes = [1000, 10000, 100000, 1000000];
  for (const numBases of sizes) {
    const bases = "ACGT";
    let dna = "";
    for (let i = 0; i < numBases; i++) dna += bases[i & 3];

    // Pack
    const packStart = performance.now();
    const packed = wasmPkg.pack_dna_to_bits(dna);
    const packMs = performance.now() - packStart;

    // Unpack to string
    const unpackStart = performance.now();
    const unpacked = wasmPkg.unpack_bits_to_dna(packed, numBases);
    const unpackMs = performance.now() - unpackStart;

    // Verify
    const verified = unpacked === dna;

    const packThroughput = numBases / packMs / 1000;
    const unpackThroughput = numBases / unpackMs / 1000;
    console.log(`  ${numBases.toLocaleString()} bases: pack=${packMs.toFixed(2)}ms (${packThroughput.toFixed(1)} MB/s) unpack=${unpackMs.toFixed(2)}ms (${unpackThroughput.toFixed(1)} MB/s) verified=${verified}`);
  }

  // -----------------------------------------------------------------------
  // BHE benchmark
  // -----------------------------------------------------------------------
  console.log("\n=== BHE Encoding (u128 FSM) ===");

  const bheSizes = [15, 30, 60, 100];
  for (const size of bheSizes) {
    const data = new Uint8Array(size);
    for (let i = 0; i < size; i++) data[i] = (i * 7 + 13) & 0xFF;

    // k=1 encode
    const k1Start = performance.now();
    const k1Dna = wasmPkg.bhe_encode_k1(data);
    const k1Ms = performance.now() - k1Start;

    // k=1 decode
    const k1DecStart = performance.now();
    const k1Decoded = wasmPkg.bhe_decode_k1(k1Dna, size);
    const k1DecMs = performance.now() - k1DecStart;

    // k=1 roundtrip
    let k1Roundtrip = true;
    if (k1Decoded.length !== size) k1Roundtrip = false;
    else for (let i = 0; i < size; i++) if (k1Decoded[i] !== data[i]) { k1Roundtrip = false; break; }

    // k=1 validate
    const k1Valid = wasmPkg.bhe_validate(k1Dna, 1);

    // k=3 encode
    const k3Start = performance.now();
    const k3Dna = wasmPkg.bhe_encode_fsm(data, 3);
    const k3Ms = performance.now() - k3Start;

    // k=3 decode
    const k3DecStart = performance.now();
    const k3Decoded = wasmPkg.bhe_decode_fsm(k3Dna, size, 3);
    const k3DecMs = performance.now() - k3DecStart;

    // k=3 roundtrip
    let k3Roundtrip = true;
    if (k3Decoded.length !== size) k3Roundtrip = false;
    else for (let i = 0; i < size; i++) if (k3Decoded[i] !== data[i]) { k3Roundtrip = false; break; }

    const k3Valid = wasmPkg.bhe_validate(k3Dna, 3);

    console.log(`  ${size}B: k=1 enc=${k1Ms.toFixed(2)}ms dec=${k1DecMs.toFixed(2)}ms valid=${k1Valid} rt=${k1Roundtrip}  k=3 enc=${k3Ms.toFixed(2)}ms dec=${k3DecMs.toFixed(2)}ms valid=${k3Valid} rt=${k3Roundtrip}`);
  }

  // -----------------------------------------------------------------------
  // Arithmetic compression benchmark
  // -----------------------------------------------------------------------
  console.log("\n=== Arithmetic Compression ===");

  const compSizes = [1000, 10000, 50000];
  for (const size of compSizes) {
    const data = new Uint8Array(size);
    for (let i = 0; i < size; i++) data[i] = (i * 37 + 17) & 0xFF;

    const compStart = performance.now();
    const compressed = wasmPkg.arith_compress(data);
    const compMs = performance.now() - compStart;

    const decStart = performance.now();
    const decompressed = wasmPkg.arith_decompress(compressed);
    const decMs = performance.now() - decStart;

    let roundtrip = true;
    if (decompressed.length !== size) roundtrip = false;
    else for (let i = 0; i < Math.min(size, decompressed.length); i++) if (decompressed[i] !== data[i]) { roundtrip = false; break; }

    const ratio = (compressed.length / size).toFixed(3);
    console.log(`  ${size}B: compress=${compMs.toFixed(2)}ms decompress=${decMs.toFixed(2)}ms ratio=${ratio} roundtrip=${roundtrip}`);
  }

  // -----------------------------------------------------------------------
  // Simulation benchmark
  // -----------------------------------------------------------------------
  console.log("\n=== Simulation ===");

  const oligoLen = 200;
  const oligo = new Uint8Array(oligoLen);
  for (let i = 0; i < oligoLen; i++) oligo[i] = "ACGT".charCodeAt(i & 3);

  const simStart = performance.now();
  const numSims = 1000;
  for (let i = 0; i < numSims; i++) {
    wasmPkg.simulate_single(oligo, 0.002, 0.0005, 0.001, 0.001, 0.0001, 0.0001, i);
  }
  const simMs = performance.now() - simStart;
  console.log(`  ${numSims} oligos (${oligoLen}nt) Illumina: ${simMs.toFixed(1)}ms total, ${(simMs / numSims).toFixed(2)}ms/oligo`);

  const nanoStart = performance.now();
  for (let i = 0; i < numSims; i++) {
    wasmPkg.simulate_single(oligo, 0.05, 0.04, 0.04, 0.05, 0.03, 0.03, i + 10000);
  }
  const nanoMs = performance.now() - nanoStart;
  console.log(`  ${numSims} oligos (${oligoLen}nt) Nanopore: ${nanoMs.toFixed(1)}ms total, ${(nanoMs / numSims).toFixed(2)}ms/oligo`);

  // -----------------------------------------------------------------------
  // Hamming distance benchmark
  // -----------------------------------------------------------------------
  console.log("\n=== Bit-parallel Hamming ===");

  const hamSize = 100000;
  const a = new Uint8Array(hamSize);
  const b = new Uint8Array(hamSize);
  for (let i = 0; i < hamSize; i++) { a[i] = i & 0xFF; b[i] = (i + 1) & 0xFF; }

  const hamStart = performance.now();
  const hamDist = wasmPkg.bit_parallel_hamming(a, b);
  const hamMs = performance.now() - hamStart;
  console.log(`  ${hamSize} bytes: ${hamMs.toFixed(2)}ms distance=${hamDist}`);

  // -----------------------------------------------------------------------
  // Rolling hash benchmark
  // -----------------------------------------------------------------------
  console.log("\n=== Rolling Hash (Rabin-Karp) ===");

  const rhPacked = wasmPkg.pack_dna_to_bits("ACGT".repeat(25000)); // 100K bases
  const rhStart = performance.now();
  const rhResult = wasmPkg.rolling_hash(rhPacked, 21);
  const rhMs = performance.now() - rhStart;
  console.log(`  100K bases, k=21: ${rhMs.toFixed(2)}ms → ${rhResult.length.toLocaleString()} hashes`);

  // -----------------------------------------------------------------------
  // RS benchmark via WASM
  // -----------------------------------------------------------------------
  console.log("\n=== Reed-Solomon (Rust GF(256)) ===");

  const rsHandle = wasmPkg.rs_create(255, 223);
  const rsData = new Uint8Array(223);
  for (let i = 0; i < 223; i++) rsData[i] = i & 0xFF;

  const rsEncStart = performance.now();
  const rsEncoded = wasmPkg.rs_encode(rsHandle, rsData);
  const rsEncMs = performance.now() - rsEncStart;

  const rsDecStart = performance.now();
  const rsDecoded = wasmPkg.rs_decode(rsHandle, rsEncoded);
  const rsDecMs = performance.now() - rsDecStart;

  let rsRoundtrip = true;
  for (let i = 0; i < 223; i++) if (rsDecoded[i] !== rsData[i]) { rsRoundtrip = false; break; }

  console.log(`  RS(255,223): encode=${rsEncMs.toFixed(2)}ms decode=${rsDecMs.toFixed(2)}ms roundtrip=${rsRoundtrip}`);

  // RS with errors
  const rsCorrupted = rsEncoded.slice();
  rsCorrupted[10] ^= 0x55;
  rsCorrupted[100] ^= 0xAA;
  const rsErrDecStart = performance.now();
  const rsErrDecoded = wasmPkg.rs_decode(rsHandle, rsCorrupted);
  const rsErrDecMs = performance.now() - rsErrDecStart;
  let rsErrRoundtrip = true;
  for (let i = 0; i < 223; i++) if (rsErrDecoded[i] !== rsData[i]) { rsErrRoundtrip = false; break; }
  console.log(`  RS(255,223) with 2 errors: decode=${rsErrDecMs.toFixed(2)}ms roundtrip=${rsErrRoundtrip}`);

  wasmPkg.rs_free(rsHandle);

  // -----------------------------------------------------------------------
  // Summary
  // -----------------------------------------------------------------------
  console.log("\n" + "═".repeat(50));
  console.log("All Rust WASM benchmarks completed!");
  console.log("═".repeat(50));
}

main().catch(e => {
  console.error("Benchmark failed:", e);
  process.exit(1);
});
