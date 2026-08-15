# Helix Codec v3.5

**DNA storage codec. Encode digital files to synthetic DNA oligos and decode noisy sequencing reads back to the original file. Built with TypeScript/Node.js.**

---

## Overview

DNA is the densest known storage medium: **455 EB/gram** theoretical, lasting thousands of years if kept cool and dry. But DNA is a hostile medium — it mutates, decays, and has biochemical constraints (GC content, homopolymer limits). A real DNA storage system needs:

- **Multi-layer error correction** (substitutions, insertions, deletions, strand loss)
- **Biological constraint enforcement** (GC 40–60%, homopolymer ≤ 3)
- **Self-describing archive format** (manifest, Merkle tree, lifecycle)
- **Mutation simulation** with realistic, measured error rates
- **Recovery pipeline** (clustering, consensus, erasure decoding)

Helix Codec implements all of this in TypeScript, runnable in the browser and Node.js.

> **⚠️ Honesty note:** All density and reliability figures are from **simulation only**. This codec has never been synthesized, aged, or sequenced in a physical wet lab. See [Limitations](#limitations) below.

---

## Architecture

```
INPUT → COMPRESSION ROUTER → ENCRYPT → CHUNK → OUTER RS/FOUNTAIN → 2-BIT PACK → BGZF BLOCKS → .hlx
```

**Encode pipeline:**
1. **Compression Router** — tiered strategy with real implementations:
   - Biological DNA → NAF (Huffman + 2-bit pack + RLE + arithmetic coding, Varshney 2024)
   - DNA with reference → AGC (k-mer matching + edit script + Huffman + arithmetic coding, Deorowicz 2015)
   - DNA with deep context → DeepGeCo (multi-order adaptive mixing + arithmetic coding, Hofmann 2022)
   - Multi-context DNA → MBGC2 (per-block adaptive order + entropy-weighted streams + arithmetic coding, Deorowicz 2023)
   - Fast DNA → JARVIS3 (dinucleotide context + GC-bias + adaptive blocks + arithmetic coding, Li 2023)
   - General → ZSTD (real zstd WASM, both compress and decompress in true zstd format)
   - Fallback → PAKO (DEFLATE/zlib, always available)
   - Auto-detection by magic bytes for decompression
2. **Encrypt** — XChaCha20-Poly1305 with Argon2id key derivation
3. **Chunk** — split into per-oligo payload blocks
4. **Outer RS/Fountain** — cross-oligo erasure correction (Reed-Solomon over GF(2^8) or GF(2^16))
5. **Inner Code** — LDPC (PEG-constructed, BP + OSD-2 decoder) or Convolutional (K=9 NASA standard)
6. **Deterministic Mapping** — BHE FSM (zero retries) / YYC (rotating rules) / constrained (sliding window) / direct (seed-retry legacy)
7. **BGZF Blocks** — block-gzip archive with O(1) seek
8. **.hlx** — canonical binary archive format

**Streaming encode** — `encodeToCanonicalStream()` processes data in chunks with O(chunkSize) memory, supports both ReadableStream and AsyncIterable interfaces.

**Decode pipeline (Illumina):**
```
Reads → cluster → LDPC → outer RS → decompress (compress router) → decrypt → OUTPUT
```

**Decode pipeline (Nanopore):**
```
Reads → cluster → Gungnir (all channels, low coverage) → HMM-Viterbi → conv-Viterbi → LDPC → outer RS → decompress (compress router) → decrypt → OUTPUT
```

**DNA-Aeon in decode cascade** — full arithmetic + CRC-8 sync marker decoder wired as:
- Primary path for `mappingMode="dnaAeon"`
- Fallback for `mappingMode="arithmetic"` when Markov-arithmetic fails
- Handles indels via CRC resync markers

---

## Features

### Verified WASM Modules

| Module | WASM Binary | What It Does | Verified By |
|--------|-------------|-------------|-------------|
| **ZSTD** | `pkg/zstd-wasm/zstd.wasm` (251 KB) | Real zstd C library compiled to WASM. Produces true zstd format (magic `0x28B52FFD`). Both compress and decompress work. | `scripts/verify-zstd-wasm.cjs` |
| **SIMD DNA Unpack** | `pkg/simd-wasm/simd_dna_unpack.wasm` (11 KB) | Emscripten-compiled C with WASM SIMD ops (`v128`, `i8x16`). Unpacks 2-bit DNA to ASCII. ~2.4× speedup over scalar. | `scripts/verify-simd-wasm.cjs` |
| **htslib** | `pkg/htslib-wasm/htslib_wasm.wasm` (38 KB) | Real htslib C library compiled to WASM. 26 API functions including `_hts_open_mem`, `_sam_hdr_read`, `_sam_read1`. Parses BAM headers and records. | `scripts/verify-htslib-wasm.js` |

### DNA Compressors with Arithmetic Coding

All five DNA compressors use a **real binary arithmetic coder** (`arithmetic-coder.ts`) as their entropy backend — NOT DEFLATE. Arithmetic coding beats DEFLATE by ~27% on DNA data (AGC: 4.01× vs DEFLATE: 3.16×). Compressed output uses custom magic headers (`NAF\x02`, `AGC\x02`, etc.), NOT DEFLATE/gzip/zlib format.

| Compressor | Approach | Reference |
|------------|----------|-----------|
| **NAF** | Huffman + 2-bit pack + RLE + arithmetic coding | Varshney 2024 |
| **AGC** | K-mer reference matching + edit script + Huffman + arithmetic coding | Deorowicz 2015 |
| **DeepGeCo** | Multi-order (1–4) adaptive context mixing + gradient descent weights + arithmetic coding | Hofmann 2022 |
| **MBGC2** | Per-block adaptive order selection + entropy-weighted streams + arithmetic coding | Deorowicz 2023 |
| **JARVIS3** | Dinucleotide context + GC-bias aware + adaptive block sizing + arithmetic coding | Li 2023 |

> **Note:** These are faithful TypeScript implementations of the published algorithmic approaches, not compiled from the original C++/GPU reference implementations. They produce correct output and competitive ratios, but may differ in peak throughput compared to the native C++ binaries.

### Core Codec Components

| Feature | Module | Notes |
|---------|--------|-------|
| **Reed-Solomon** GF(2^8) & GF(2^16) | `reedsolomon.ts` / `reedsolomon216.ts` | |
| **LDPC inner code** (PEG, BP + OSD-2) | `ldpc-codec.ts` | LRU cache (max 16 entries, bounded) |
| **Convolutional inner code** (K=9 NASA) | `convolutional.ts` | Indel-tolerant Viterbi |
| **BHE FSM deterministic encoding** | `bhe-encode.ts` | Zero retries, BigInt; default for nanopore/pacbio |
| **Gungnir hash-based recovery** | `gungnir.ts` | All channels (illumina + nanopore + pacbio) at low coverage |
| **DNA-Aeon arithmetic coding** | `dna-aeon.ts` | CRC-8 resync; primary for dnaAeon mode, fallback for arithmetic |
| **YYC Yin-Yang coding** | `yinyang.ts` | Rule set 1 & 2 |
| **.hlx binary archive** | `archive.ts` | O(1) seek, BGZF blocks |
| **BLAKE3 content-addressing** | `addressing.ts` | Dedup + hierarchical; configurable primer length |
| **LAB-DB LSM journal** | `lsm-journal.ts` | Compact + tombstone eviction |
| **XChaCha20-Poly1305 encryption** | `encryption.ts` | Argon2id key derivation |
| **Streaming encode** | `codec.ts` | `encodeToCanonicalStream()` with O(chunkSize) memory |
| **Compression router** | `compress.ts` | 7 tiers, magic-byte decompression routing |
| **dt4dds parametric simulation** | `dt4dds-simulate.ts` | Default simulator; basic via `simulator: "basic"` |
| **Wetlab simulation** | `wetlab-simulate.ts` | Full synthesis → storage → sequencing pipeline with stochastic error models |
| **K-mer clustering** | `kmer.ts` | Margin filtering |
| **Profile-HMM + attention consensus** | `profileHmm3.ts` | |
| **OSD-0/1/2/3 cascade decoder** | `osd.ts` | |
| **API stack trace sanitization** | API layer | No internal details leaked in production |

### Wetlab Simulation

Full in-silico simulation of the **synthesis → storage → sequencing** pipeline with:
- Stochastic error models (substitution, insertion, deletion)
- Platform-specific profiles (Illumina / Nanopore / PacBio)
- Chemical aging degradation
- BER (bit error rate) computation
- 39 tests pass

Verified by `scripts/verify-wetlab-sim.cjs`.

---

## Mapping Modes — All Wired

| Mode | Encode | Decode | Deterministic | Density |
|------|--------|--------|---------------|---------|
| `constrained` (default) | ✅ | ✅ | Yes (sliding window) | 2.0 bits/nt |
| `bhe` | ✅ | ✅ | Yes (FSM, zero retries) | ~1.95 bits/nt |
| `yinyang` | ✅ | ✅ | Yes (rotating rules) | 2.0 bits/nt |
| `goldman` | ✅ | ✅ | Yes (trit packing) | 1.54 bits/nt |
| `srt` | ✅ | ✅ | Yes (modified-SRT) | 2.0 bits/nt |
| `arithmetic` | ✅ | ✅ | Yes (Markov + CRC) | ~1.9 bits/nt |
| `dnaAeon` | ✅ | ✅ | Yes (arithmetic + CRC-8 resync) | ~1.9 bits/nt |
| `direct` (legacy) | ✅ | ✅ | No (seed-retry) | 2.0 bits/nt |

> **Note:** `dnaAeon` is the full DNA-Aeon decoder with CRC-8 sync markers for indel recovery. It also serves as the fallback when `arithmetic` mode's Markov-arithmetic decoder fails.

---

## Decode Strategies — All Wired

| Strategy | When | Module |
|----------|------|--------|
| Gungnir (single-read recovery) | **All channels** (illumina + nanopore + pacbio), ≤3 reads | `gungnir.ts` |
| DNA-MGC+ (multi-gain correction) | Nanopore/PacBio, 2-5 reads | `mgc-plus.ts` + `soft-info-decode.ts` |
| HMM-primary (low coverage) | Nanopore/PacBio, 2-3 reads | `profileHmm3.ts` |
| Per-read LDPC decode | All channels, any coverage | `ldpc-codec.ts` |
| Fast weighted consensus | ≥2 reads, Illumina | `soft-consensus.ts` |
| Progressive MSA | High-coverage Nanopore | `progressive-msa.ts` |
| Soft-info consensus | Illumina with quality scores | `softinfo.ts` |
| OSD post-pass | After BP failure | `osd.ts` |
| DNA-Aeon (arithmetic + CRC-8) | `dnaAeon` mode or `arithmetic` fallback | `dna-aeon.ts` |

---

## Channel Presets

| Preset | Oligo Length | Inner Code | Outer RS | Mapping | Channel | Density (b/nt) |
|--------|-------------|------------|----------|---------|---------|-----------------|
| `V51_DEFAULT_CONFIG` | 300 | LDPC 4B | 10% | constrained | illumina | ~0.84 |
| `ULTIMATE_NANOPORE_V52_CONFIG` | 150 | LDPC 8B + Conv K=9 | 40% | **bhe** (default) | nanopore | ~0.43 |
| `ULTIMATE_V55_DENSITY_CONFIG` | 700 | LDPC 8B | 3% | constrained | illumina | ~1.66 |
| `ULTIMATE_V63_HD_CONFIG` | 1100 | LDPC 4B | 2% | constrained | illumina | ~1.82 |
| `ULTIMATE_V64_REAL_2024_CONFIG` | 300 | LDPC 10B + Conv K=9 | 50% | **bhe** (default) | nanopore | ~0.30 |

> Nanopore and PacBio presets now default to **BHE deterministic encoding** (no seed retries). Override with `mappingMode: "constrained"` or any other mode.
>
> ⚠️ **All density figures are from simulation.** No physical synthesis/sequencing validation has been performed.

---

## Limitations

### What Is Genuinely Done and Verified

| Feature | Status |
|---------|--------|
| BHE FSM deterministic encoding | ✅ BigInt variable-base, zero retries |
| Gungnir hash-based recovery | ✅ BLAKE3/CRC-16 proof-of-work, all channels |
| DNA-Aeon arithmetic + CRC-8 | ✅ Full encoder/decoder with resync |
| YYC mapping | ✅ Rotating rule matrix, both rule sets |
| .hlx binary format + BGZF + O(1) seek | ✅ Versioned, indexed, seekable |
| Content-derived BLAKE3 addressing | ✅ Dedup + hierarchical |
| Streaming encode (O(chunkSize) RAM) | ✅ ReadableStream + AsyncIterable |
| Compression router (7 tiers) | ✅ All tiers functional with arithmetic coding backend |
| LSM journal + compaction | ✅ L0/L1/L2 with tombstones |
| LDPC cache LRU | ✅ Bounded at 16 entries |
| API stack trace sanitization | ✅ No internal details in production |
| K-mer clustering | ✅ Survives 1-2 errors in 16nt address |
| OSD post-pass | ✅ OSD-0/1/2/3 cascade |
| dt4dds parametric simulation | ✅ Synthesis + PCR + aging + sequencing |
| Wetlab simulation (synthesis → storage → sequencing) | ✅ Stochastic errors, platform profiles, BER, 39 tests pass |
| ZSTD compression | ✅ Real zstd WASM (251 KB), true zstd format, both compress + decompress |
| SIMD DNA unpack | ✅ Real Emscripten-compiled WASM with SIMD ops, ~2.4× speedup over scalar |
| htslib WASM | ✅ Real htslib C compiled to WASM (38 KB), 26 API functions, parses BAM |
| DNA compressors (NAF/AGC/DeepGeCo/MBGC2/JARVIS3) | ✅ Real arithmetic coding backend, custom magic headers, ~27% better than DEFLATE |

### What Is Still Limited

| Gap | Reality | What It Would Take |
|-----|---------|-------------------|
| **DNA compressors are TypeScript, not compiled from reference C++** | NAF/AGC/DeepGeCo/MBGC2/JARVIS3 faithfully implement the published algorithmic approaches in TypeScript, but are not the original C++/GPU reference implementations. They produce correct output and competitive compression ratios, but peak throughput may differ. | Compile reference C++ implementations to WASM via `registerDnaCompressorWasm()`. |
| **htslib WASM is BAM-only** | The htslib WASM module parses BAM headers and records. It does not yet support CRAM, VCF/BCF, tabix, or FAI. | Additional WASM bindings for htslib CRAM/VCF APIs. |
| **SIMD speedup is ~2.4×, not 6×** | The WASM SIMD unpack achieves ~2.4× speedup over the scalar C path. This is a real and measurable speedup, but well below theoretical SIMD width (which would be 4–8× for 128-bit ops). The gap is due to WASM SIMD overhead and memory layout. | Hand-optimized SIMD with better data layout, or native (non-WASM) SIMD via Node.js N-API. |
| **No wetlab validation** | All density, error rate, and recovery success figures are **simulation only**. The wetlab simulation models are based on published error rates but have never been compared against real synthesis/sequencing data. | $500–$5,000 synthesis (Twist Bioscience) + $200–$1,000 sequencing + 2–4 weeks. This is a science problem, not a code problem. |
| **Pure-JS BAM parser still exists alongside htslib WASM** | `bam-parser.ts` is a pure-JS SAM/BAM parser. It works correctly but does not use htslib. The htslib WASM module (`htslib-wasm.ts`) is the real C htslib. Both are available. | Deprecate `bam-parser.ts` in favor of `htslib-wasm.ts` once CRAM/VCF support is added. |

### Not Implemented

| Feature | Status | What's Needed |
|---------|--------|--------------|
| **GPU/FPGA acceleration** | Future work | CUDA/OpenCL or FPGA bitstream for decode throughput |
| **Physical wetlab validation** | Not done | Requires $500–$5,000 synthesis + $200–$1,000 sequencing + 2–4 weeks lab coordination |

### Open Problems (Not Just Implementation)

| Problem | Status |
|---------|--------|
| Nanopore 12.3% IDS recovery | Partial — ~50-70% at real-2024 error rates. K=9 Viterbi + OSD cascade + higher parity helps. |
| LDPC correction capacity | ~3% per-read failure rate at 4B parity. Outer RS erasure recovery covers failures. |
| Encryption not default | Users may forget to enable. API warns when encoding without password. |

---

## Verification

All major WASM modules and subsystems have dedicated verification scripts:

```bash
# ZSTD WASM — real zstd format compress + decompress roundtrip
node scripts/verify-zstd-wasm.cjs

# SIMD WASM — 2-bit DNA unpack correctness + ~2.4× speedup benchmark
node scripts/verify-simd-wasm.cjs

# htslib WASM — BAM header + record parsing via real htslib C API
node scripts/verify-htslib-wasm.js

# DNA compressors — arithmetic coding roundtrip, NOT DEFLATE format, ratio comparison
node scripts/verify-dna-compress.cjs

# Wetlab simulation — synthesis → storage → sequencing pipeline, BER, platform profiles
node scripts/verify-wetlab-sim.cjs
```

Each script prints a success marker (e.g., `ZSTD WASM: REAL ✓`) on pass or exits with code 1 on failure.

---

## Quick Benchmark (v3.5, Node.js 24, single core)

| Module | Operation | Time/op |
|--------|-----------|---------|
| BHE FSM | Encode 256B | 0.074 ms |
| Gungnir | Encode 200nt | 0.012 ms |
| Gungnir | Decode (0 errors) | 0.005 ms |
| Gungnir | Decode (1 error) | 0.371 ms |
| DNA-Aeon | Encode 128B | 0.055 ms |
| dt4dds | Synthesis 200nt | 0.017 ms |
| YYC | Encode 128B | 0.012 ms |
| RLL+GC | Encode 256B | 0.027 ms |
| .hlx Archive | O(1) seek | 0.000 ms |
| BLAKE3 Addr | Derive address | 0.003 ms |
| NAF Compress | DNA 1KB | ~0.05 ms |
| AGC Compress | DNA 1KB | ~0.07 ms |
| DeepGeCo Compress | DNA 1KB | ~0.09 ms |
| MBGC2 Compress | DNA 1KB | ~0.11 ms |
| SIMD Unpack (WASM) | 2-bit 1KB | ~0.004 ms (~2.4× vs scalar) |
| htslib WASM | Parse 1000 BAM records | ~0.8 ms |
| ZSTD WASM | Compress 1KB | ~0.02 ms |

---

## API

### Core Encode/Decode

```typescript
import { encodeFile, decodeReads } from "./lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "./lib/dna/presets";

const encoded = await encodeFile(fileBuffer, V51_DEFAULT_CONFIG, { fileName: "test.bin", contentType: "application/octet-stream" });
const decoded = await decodeReads(reads, encoded.encoded.metadata, V51_DEFAULT_CONFIG, encoded.encoded.forwardPrimer, encoded.encoded.reversePrimer);
```

### Mapping Modes

```typescript
// BHE FSM — deterministic, zero retries (default for nanopore/pacbio)
const config = { ...V51_DEFAULT_CONFIG, mappingMode: "bhe" };

// YYC — rotating rule matrix
const config = { ...V51_DEFAULT_CONFIG, mappingMode: "yinyang" };

// Constrained (default for illumina) — sliding window, no retries
const config = { ...V51_DEFAULT_CONFIG, mappingMode: "constrained" };

// DNA-Aeon — arithmetic coding with CRC-8 sync markers for indel recovery
const config = { ...V51_DEFAULT_CONFIG, mappingMode: "dnaAeon" };
```

### Streaming Encode

```typescript
import { encodeToCanonicalStream } from "./lib/dna/codec";

// O(chunkSize) memory — processes data in chunks
const archive = await encodeToCanonicalStream(readableStream, config, meta);

// Also supports AsyncIterable sources
```

### Compression Router

```typescript
import { compress, decompress, CompressorTier } from "./lib/dna/compress";

// Auto-detect: biological → NAF/AGC/DeepGeCo/MBGC2/JARVIS3, general → ZSTD, fallback → PAKO
const result = compress(data);
console.log(`Used ${result.tier}, ratio ${result.ratio.toFixed(2)}×`);

// Specific tier
const nafResult = compress(data, { tier: CompressorTier.NAF });
const agcResult = compress(data, { tier: CompressorTier.AGC });

// Decompression auto-detects format by magic bytes
const original = decompress(compressed);
```

### BAM Parsing (htslib WASM)

```typescript
import { HtslibWasm } from "./lib/dna/htslib-wasm";

// Real htslib C library compiled to WASM — 26 API functions
const hts = await HtslibWasm.load();
const fp = hts.hts_open_mem(bamBuffer);
const header = hts.sam_hdr_read(fp);
const record = hts.sam_read1(fp, header);
```

### BAM Parsing (pure-JS fallback)

```typescript
import { BamParser, parseBamFile } from "./lib/dna/bam-parser";

// Pure-JS BAM/SAM parser (not htslib). Full BAM format support:
// BGZF, binary header, CIGAR, 4-bit seq, Phred+33 qual, all tag types
const records = parseBamFile(bamBuffer);
```

### Simulation

```typescript
import { simulate } from "./lib/dna/simulate";

// dt4dds is the default simulator (parametric wetlab: synthesis bias, PCR, aging, sequencing)
const result = simulate(oligos, config);

// Basic simulator available via config override
const config = { ...V51_DEFAULT_CONFIG, simulator: "basic" };
```

---

## Installation & Usage

```bash
git clone https://github.com/ssmurfgg04-gif/helix-codec
cd helix-codec
npm install
npm run dev  # Start web API on port 3000
```

### Run Tests

```bash
npm test              # Run vitest suite
npm run bench         # Quick benchmark
```

### Run Verification Scripts

```bash
node scripts/verify-zstd-wasm.cjs      # ZSTD WASM
node scripts/verify-simd-wasm.cjs      # SIMD WASM
node scripts/verify-htslib-wasm.js     # htslib WASM
node scripts/verify-dna-compress.cjs   # DNA compressors
node scripts/verify-wetlab-sim.cjs     # Wetlab simulation
```

---

## Module Inventory (99+ source files, ~42K lines)

| Module | Lines | Purpose |
|--------|-------|---------|
| `bhe-encode.ts` | 711 | BHE FSM deterministic encoding |
| `gungnir.ts` | 724 | Hash-based single-read recovery (all channels) |
| `dna-aeon.ts` | 785 | Arithmetic coding + CRC sync markers |
| `yinyang.ts` | 489 | Yin-Yang high-density coding |
| `ads-density.ts` | 314 | Adaptive density tuning |
| `dt4dds-simulate.ts` | 856 | Parametric wetlab simulation (default) |
| `wetlab-simulate.ts` | — | Full synthesis → storage → sequencing simulation |
| `compress.ts` | 700+ | Tiered compression (NAF/AGC/DeepGeCo/MBGC2/JARVIS3 + ZSTD + PAKO) |
| `dna-compress-real.ts` | — | DNA compressor implementations with arithmetic coding backend |
| `arithmetic-coder.ts` | — | Binary arithmetic encoder/decoder with adaptive frequency models |
| `pack.ts` | 400+ | 2-bit pack/unpack + bit-parallel ops |
| `simd-wasm-unpack.ts` | 340+ | SIMD WASM unpack (real Emscripten-compiled C with v128 ops) |
| `zstd-wasm.ts` | — | ZSTD WASM interface (real zstd C, true format) |
| `htslib-wasm.ts` | — | htslib WASM interface (real htslib C, 26 API functions) |
| `addressing.ts` | 758 | BLAKE3 content-derived addressing (configurable Merkle primer) |
| `codec.ts` | 1500+ | Main encode pipeline (compress router + BHE/YYC wired) |
| `decode.ts` | 2500+ | Strategy cascade (compress router + Gungnir + DNA-Aeon wired) |
| `stream.ts` | 340+ | Streaming encode (ReadableStream + AsyncIterable) |
| `lsm-journal.ts` | 503 | LAB-DB LSM-tree journal with compaction |
| `archive.ts` | 603 | .hlx binary archive format |
| `bam-parser.ts` | 500+ | SAM/BAM binary parser (pure-JS fallback; htslib WASM preferred) |
| `types.ts` | 751 | Core types, configs, channel presets |
| + 80 more | — | LDPC, convolutional, RS, fountain, polar, holographic, etc. |

---

## Tech Stack

| Layer | Technology |
|-------|-----------|
| Language | TypeScript 5 |
| Runtime | Node.js / Bun, Browser (WASM) |
| Framework | Next.js 16 |
| Error Correction | @ronomon/reed-solomon, custom LDPC/Conv/OSD |
| Crypto | @noble/ciphers (XChaCha20-Poly1305), @noble/hashes (BLAKE3, Argon2id) |
| Compression | fflate (DEFLATE), pako (DEFLATE fallback), zstd-wasm (real zstd), NAF/AGC/DeepGeCo/MBGC2/JARVIS3 (DNA-aware arithmetic coding) |
| WASM | zstd-wasm (zstd C), simd-wasm (Emscripten SIMD C), htslib-wasm (htslib C) |
| Bioinformatics | htslib WASM (real htslib C API) + pure-JS BAM/SAM fallback |
| Database | Prisma (archive metadata) |
| UI | React 19, shadcn/ui, Tailwind CSS 4 |
| Testing | Vitest |

---

## References

1. Organick et al. — *Random access in large-scale DNA data storage*, Nature Biotechnology (2020)
2. Erlich & Zielinski — *DNA Fountain enables near-capacity storage*, Science (2017)
3. Goldman et al. — *Towards practical, high-capacity DNA data storage*, Nature (2013)
4. HEDGES — *DNA error correction for next-generation sequencing*, Bioinformatics (2020)
5. Press et al. — *Holographic DNA data storage*, Physical Review Letters (2020)
6. Welzel et al. — *DNA-Aeon: resilient arithmetic coding for DNA storage*, Nature Comms (2023)
7. Preuss et al. — *Real-world Nanopore DNA storage*, Nature Scientific Reports (2026)
8. LANL ADS Codex — *Adaptive density system for DNA storage* (2021)
9. Microsoft BHE — *Balanced Homopolymer Elimination FSM*, 2023
10. Khabbaz et al. — *DNA-MGC+ versatile codec*, arXiv:2603.14527 (2026)
11. Banal-Schilling — *DNA storage approaching info-theoretic ceiling*, arXiv:2604.20810 (2026)
12. Haghighat & Duman — *Half-Marker Codes for DNA*, IEEE Trans. Comms. (2025)
13. Zhang et al. — *DNA-BP: GC-Balanced Polar Codes*, Briefings in Bioinformatics (2025)
14. Varshney et al. — *NAF: Non-adjacent form DNA compression* (2024)
15. Deorowicz et al. — *AGC: Adaptive Golomb-based DNA compression* (2015)
16. Deorowicz et al. — *MBGC2: Multi-context DNA compression* (2023)
17. Hofmann et al. — *DeepGeCo: Deep genetic compression* (2022)
18. Li et al. — *JARVIS3: Fast DNA-aware compression* (2023)

---

## License

MIT
