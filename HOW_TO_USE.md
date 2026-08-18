# How to Use Helix Codec

This guide explains how to install, build, run, and use the Helix Codec DNA storage system. It is written for both human developers and AI coding agents.

---

## Quick Start

### Prerequisites

- **Node.js** >= 20 (or **Bun** >= 1.0)
- **Rust** toolchain (for native Viterbi addon — `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y`)
- **Python 3** + **make** + **g++** (for `@ronomon/reed-solomon` native addon, optional)

### Install

```bash
# Clone the repo
git clone https://github.com/ssmurfgg04-gif/helix-codec.git
cd helix-codec

# Install dependencies (skip native build scripts if build tools unavailable)
npm install --ignore-scripts

# Or with Bun (faster)
bun install --ignore-scripts
```

### Run the Dev Server

```bash
# Start Next.js dev server on port 3000
npm run dev
# Or
bun run dev
```

Open `http://localhost:3000` in your browser to use the web UI.

---

## Core API: Encode and Decode Files

### Encode a File to DNA Oligos

```typescript
import { HelixCodec } from "./src/lib/dna/codec";

// Create codec instance
const codec = new HelixCodec({
  // Optional: encrypt with password
  // encryptPassword: "my-secret-key",
});

// Encode a file (Uint8Array) to DNA oligos
const fileData = new TextEncoder().encode("Hello, DNA storage!");
const encodeResult = await codec.encode(fileData);

// encodeResult contains:
//   .oligos: string[]  — DNA oligo sequences (A/C/G/T only)
//   .metadata: object  — decoding metadata (must be saved!)
//   .hash: string      — SHA-256 of original file
console.log(`Encoded to ${encodeResult.oligos.length} oligos`);
```

### Decode DNA Oligos Back to a File

```typescript
// Simulate sequencing errors (substitutions, insertions, deletions)
const noisyOligos = simulateChannel(encodeResult.oligos, {
  substitutionRate: 0.01,
  insertionRate: 0.005,
  deletionRate: 0.005,
});

// Decode noisy reads back to the original file
const decodeResult = await codec.decode(noisyOligos, encodeResult.metadata);

// decodeResult.data is a Uint8Array matching the original file
console.log(new TextDecoder().decode(decodeResult.data));
// → "Hello, DNA storage!"
```

---

## Running Validation Scripts

The `scripts/` directory contains many validation and benchmark scripts. All are run with `npx tsx`:

### Cascade Validation (Full Pipeline Test)

```bash
# Full Viterbi + OSD + LDPC + RS cascade validation
npx tsx scripts/cascade-validation.ts

# Quick cascade (fewer oligos, faster)
npx tsx scripts/cascade-quick.ts
```

### Individual Component Tests

```bash
# Reed-Solomon roundtrip
npx tsx scripts/test-rs.ts

# LDPC inner code
npx tsx scripts/test-ldpc.ts

# LDPC erasure decoding
npx tsx scripts/test-ldpc-erasure.ts

# Holographic sharding
npx tsx scripts/test-holographic.ts

# OSD (Ordered Statistics Decoding)
npx tsx scripts/test-osd2.ts

# Full codec roundtrip
npx tsx scripts/test-codec.ts

# Zero-noise channel
npx tsx scripts/test-zero-noise.ts
```

### Benchmarking

```bash
# Quick benchmark
node scripts/quick-bench.js

# Full benchmark suite
npx tsx scripts/bench-final.ts

# Speed benchmark
npx tsx scripts/bench-speed.ts

# Viterbi K=9 timing
npx tsx scripts/v64-k9-timing.ts
```

### Nanopore Channel Tests

```bash
# 9% IDS Nanopore simulation
npx tsx scripts/nanopore-9pct-validation.ts

# Full Nanopore validation
npx tsx scripts/nanopore-validation.ts

# Indel-tolerant Viterbi test
npx tsx scripts/v60-indel-viterbi-test.ts
```

---

## Rust Viterbi Native Addon

The K=9 Viterbi decoder is the #1 hot path (~800ms/oligo in JS, ~5ms in Rust). It's built as a napi-rs native addon.

### Build the Rust Addon

```bash
cd rust/helix-dna-napi

# Build release
cargo build --release

# Run Viterbi unit tests (12 tests)
cargo test

# Run diagnostic example
cargo run --example diag
```

### Use the Native Addon from TypeScript

```typescript
import { viterbi_k9_decode, conv_k9_encode, napi_version } from "helix-dna-napi";

// Check version
console.log(napi_version());
// → "helix-dna-napi v0.4.2 — Viterbi v4.2 ..."

// Encode with K=9 convolutional code
const data = new Uint8Array([0, 1, 2, 3, 4]);
const encoded = conv_k9_encode(data);

// Decode with indel-tolerant Viterbi
const decoded = viterbi_k9_decode(encoded, {
  maxDrift: 15,
  insertionPenalty: 1.5,
  deletionPenalty: 1.5,
});
```

### Rust WASM Viterbi (Alternative)

The WASM build is in `rust/helix-dna-wasm/`:

```bash
cd rust/helix-dna-wasm
wasm-pack build --target nodejs --release
```

The TypeScript code automatically falls back from native → WASM → pure JS.

---

## Decode Pipeline Stages

The full decode pipeline is a cascade of error-correction stages:

```
Noisy Reads → MSA Consensus → Viterbi (K=9) → OSD → LDPC → Reed-Solomon → Original Data
```

| Stage | Corrects | Code |
|-------|----------|------|
| MSA Consensus | Insertions & deletions (position-level) | `viterbi-preprocess.ts` |
| Viterbi K=9 | Convolutional decoding + indel tolerance | `convolutional-indel.ts` / Rust napi |
| OSD-0/1/2 | Residual errors after Viterbi | `osd.ts` / `osd-full.ts` |
| LDPC | Substitution errors (BP + OSD decoder) | `ldpc-codec.ts` |
| Reed-Solomon | Cross-oligo erasure correction | `reedsolomon.ts` / `fast-rs.ts` |

### Key Configuration Parameters

```typescript
// Viterbi config
const viterbiConfig = {
  maxDrift: 15,            // Max insertion-deletion drift (15 covers 99.99% at 9% IDS)
  insertionPenalty: 1.5,   // Cost of insertion transition (MUST equal deletionPenalty)
  deletionPenalty: 1.5,    // Cost of deletion transition
};

// LDPC config: getCachedLDPCInner(n, k)
// n = total codeword bytes (payload + parity)
// k = data bytes
const ldpcCode = getCachedLDPCInner(payloadBytes + ldpcParityBytes, payloadBytes);

// Reed-Solomon config
const rs = new ReedSolomon({ n: 255, k: 223 }); // Can correct up to 32 erasures
```

---

## Project Structure

```
helix-codec/
├── src/
│   ├── lib/dna/           # Core codec library
│   │   ├── codec.ts       # Main HelixCodec class (encode/decode)
│   │   ├── encode-parallel.ts  # Parallel encoding
│   │   ├── decode.ts      # Decode orchestrator
│   │   ├── convolutional-indel.ts  # Indel-tolerant Viterbi (TS)
│   │   ├── convolutional-k9.ts     # NASA K=9 config
│   │   ├── ldpc-codec.ts  # LDPC inner code
│   │   ├── reedsolomon.ts # Reed-Solomon outer code
│   │   ├── osd.ts         # Ordered Statistics Decoding
│   │   ├── viterbi-preprocess.ts  # MSA + HMM indel correction
│   │   ├── simulate.ts    # Channel simulation
│   │   ├── msa-consensus.ts  # Multiple sequence alignment
│   │   ├── native/viterbi-napi.ts  # Rust napi bridge
│   │   └── ...            # Many more modules
│   ├── components/dna/    # Web UI panels
│   └── app/               # Next.js app
├── rust/
│   ├── helix-dna-napi/    # Rust napi-rs addon (Viterbi hot path)
│   │   ├── src/lib.rs     # Viterbi v4.2 decoder
│   │   └── tests/         # 12 unit tests
│   └── helix-dna-wasm/    # Rust WASM build (alternative)
│       ├── src/viterbi.rs # WASM Viterbi decoder
│       └── src/pack.rs    # 2-bit pack/unpack
├── scripts/               # Validation & benchmark scripts
├── test-data/             # Test datasets
├── .github/workflows/     # CI pipeline
└── package.json           # Node.js project
```

---

## Running Tests

### Unit Tests (vitest)

```bash
# Run all unit tests
npm test

# Watch mode
npm run test:watch

# With coverage
npm run test:coverage
```

### Lint

```bash
# ESLint (flat config v9)
npm run lint
```

### Rust Tests

```bash
cd rust/helix-dna-napi
cargo test           # 12 Viterbi tests
cargo test --release # Optimized build
```

---

## CI Pipeline

The GitHub Actions CI (`.github/workflows/ci.yml`) runs three jobs:

1. **lint**: `bun install --ignore-scripts` + `bun run lint`
2. **test**: `npm ci --ignore-scripts` + `npm test` + codec scripts + Rust Viterbi tests
3. **build**: `bun install --ignore-scripts` + `bun run build` (Next.js production build)

The `--ignore-scripts` flag is used because `@ronomon/reed-solomon` has a native C++ addon that requires node-gyp build tools. The addon is rebuilt separately when possible (non-critical — the JS fallback works).

---

## Common Tasks for AI Agents

### Adding a new error-correction stage

1. Implement the stage in `src/lib/dna/your-stage.ts`
2. Wire it into the decode pipeline in `src/lib/dna/decode.ts`
3. Add a validation script in `scripts/test-your-stage.ts`
4. Add unit tests in a `.test.ts` file
5. Run `npm test` and `npx tsx scripts/test-your-stage.ts` to verify

### Tuning Viterbi penalties

- `insertionPenalty` and `deletionPenalty` MUST be equal (1.5 is default)
- `maxDrift` should be 15 for ≤9% IDS, 20 for higher rates
- Run `cd rust/helix-dna-napi && cargo test test_k9_sweep -- --nocapture` to see error rates

### Debugging decode failures

1. Run `npx tsx scripts/cascade-validation.ts` for per-stage breakdown
2. Run `npx tsx scripts/v60-indel-viterbi-test.ts` for Viterbi-specific debugging
3. Run `cd rust/helix-dna-napi && cargo run --example diag` for Rust-level diagnostics

---

## Important Notes

- **Standalone Viterbi cannot fully correct at >3% IDS** — the full MSA + Viterbi + LDPC + RS cascade is required for Nanopore-level noise (5-9% IDS)
- **All density/reliability figures are simulation-only** — this codec has never been synthesized in a wet lab
- **Encryption is optional but recommended** — without it, DNA data is vulnerable to unauthorized read/write
- **The `package-lock.json` must be regenerated if dependencies change** — empty version strings in the lockfile cause `npm ci` to fail
