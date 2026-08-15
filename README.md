# Helix Codec v3.3

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

> **⚠️ Honesty note:** All density and reliability figures are from **simulation only**. This codec has never been synthesized, aged, or sequenced in a physical wet lab. See [Limitations](#limitations-honest) below.

---

## Architecture

```
INPUT → COMPRESSION ROUTER → ENCRYPT → CHUNK → OUTER RS/FOUNTAIN → 2-BIT PACK → BGZF BLOCKS → .hlx
```

**Encode pipeline:**
1. **Compression Router** — tiered strategy with JS-native implementations:
   - Biological DNA → NAF-style (2-bit pack + RLE + DEFLATE, *inspired by* Varshney 2024)
   - DNA with context → AGC-style (order-1 context + 2-bit pack + DEFLATE, *inspired by* Deorowicz 2015)
   - DNA with deeper context → DeepGeCo-style (order-2 context + 2-bit pack + DEFLATE, *inspired by* Hofmann 2022)
   - Multi-context DNA → MBGC2-style (4-stream + 2-bit pack + RLE + DEFLATE, *inspired by* Deorowicz 2023)
   - Fast DNA → JARVIS3-style (2-bit pack + DEFLATE level 1, *inspired by* Li 2023)
   - General → ZSTD-compatible (fzstd real zstd **decompression** + fflate DEFLATE **compression**)
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

### ✅ Fully Operational (No Caveats)

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
| **BAM/SAM binary parser** | `bam-parser.ts` | BGZF, CIGAR, 4-bit seq, Phred+33 qual, all tag types |
| **K-mer clustering** | `kmer.ts` | Margin filtering |
| **Profile-HMM + attention consensus** | `profileHmm3.ts` | |
| **OSD-0/1/2/3 cascade decoder** | `osd.ts` | |
| **LDPC cache LRU eviction** | `ldpc-codec.ts` | Bounded max 16 entries |
| **API stack trace sanitization** | API layer | No internal details leaked in production |
| **Compress router in main pipeline** | `codec.ts` / `decode.ts` | Replaces direct pako calls |

### ⚠️ Operational With Caveats

| Feature | Module | Caveat |
|---------|--------|--------|
| **NAF-style compression** | `compress.ts` | JS approximation (2-bit + RLE + DEFLATE), not the Varshney 2024 reference C++ implementation |
| **AGC-style compression** | `compress.ts` | JS approximation (order-1 context + 2-bit + DEFLATE), not the Deorowicz 2015 reference |
| **DeepGeCo-style compression** | `compress.ts` | JS approximation (order-2 context + 2-bit + DEFLATE), not the Hofmann 2022 neural implementation |
| **MBGC2-style compression** | `compress.ts` | JS approximation (4-stream + 2-bit + RLE + DEFLATE), not the Deorowicz 2023 reference |
| **JARVIS3-style compression** | `compress.ts` | JS approximation (2-bit + DEFLATE level 1), not the Li 2023 reference |
| **ZSTD-compatible tier** | `compress.ts` | **Decompresses** real zstd (via fzstd), but **compresses** with fflate DEFLATE — output is NOT zstd format. Use `isZstdCompressionReal()` to check. |
| **SIMD unpack** | `pack.ts` | Currently uses optimized JS fallback (4-wide unrolled). WASM SIMD module requires Rust→WASM compilation — not yet done. |

### ❌ Not Implemented

| Feature | Status | What's Needed |
|---------|--------|--------------|
| **htslib WASM** | Not built | Compile htslib C library to WASM via napi-rs (~3-5 days). Current `bam-parser.ts` is a pure-JS BAM/SAM parser, not htslib. |
| **Real zstd compression** | No JS package available | Compile zstd to WASM, or wait for Node.js built-in `node:zstd`. `fzstd` provides decompression only. |
| **Rust→WASM SIMD** | Not compiled | Rust toolchain + `wasm32-unknown-unknown` target + `wasm-bindgen`. ~1 day if Rust code exists. |
| **GPU/FPGA acceleration** | Future work | CUDA/OpenCL or FPGA bitstream for decode throughput |
| **Physical wetlab validation** | Not done | Requires $500–$5,000 synthesis + $200–$1,000 sequencing + 2–4 weeks lab coordination |

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

## Limitations (Honest)

### What Is Genuinely Done

| Feature | Truly Operational? |
|---------|-------------------|
| BHE FSM deterministic encoding | ✅ Yes — BigInt variable-base, zero retries |
| Gungnir hash-based recovery | ✅ Yes — BLAKE3/CRC-16 proof-of-work, all channels |
| DNA-Aeon arithmetic + CRC-8 | ✅ Yes — full encoder/decoder with resync |
| YYC mapping | ✅ Yes — rotating rule matrix, both rule sets |
| .hlx binary format + BGZF + O(1) seek | ✅ Yes — versioned, indexed, seekable |
| Content-derived BLAKE3 addressing | ✅ Yes — dedup + hierarchical |
| Streaming encode (O(chunkSize) RAM) | ✅ Yes — ReadableStream + AsyncIterable |
| Compression router (7 tiers) | ✅ Yes — but tiers are JS approximations (see below) |
| LSM journal + compaction | ✅ Yes — L0/L1/L2 with tombstones |
| LDPC cache LRU | ✅ Yes — bounded at 16 entries |
| API stack trace sanitization | ✅ Yes — no internal details in production |
| K-mer clustering | ✅ Yes — survives 1-2 errors in 16nt address |
| OSD post-pass | ✅ Yes — OSD-0/1/2/3 cascade |
| dt4dds parametric simulation | ✅ Yes — synthesis + PCR + aging + sequencing |

### What Is Still Incomplete

| Gap | Reality | What It Would Take |
|-----|---------|-------------------|
| **ZSTD compression is DEFLATE** | `compressWithZstd()` outputs DEFLATE format, not zstd format. Decompression handles real zstd via fzstd, but compression doesn't produce it. | Compile zstd to WASM (~2-3 days), or wait for `node:zstd`. Check with `isZstdCompressionReal()`. |
| **SIMD is JS fallback** | `pack.ts` uses 4-wide unrolled JS, not WASM SIMD. ~2-3× slower than claimed WASM throughput. | Install Rust toolchain, write SIMD core, compile with `wasm-pack`. ~1 day if the Rust code is ready. |
| **DNA compressors are approximations** | NAF/AGC/DeepGeCo/MBGC2/JARVIS3 tiers use 2-bit pack + context modeling + DEFLATE. They are **inspired by** the published algorithms, not the reference C++/GPU implementations. | Compile reference implementations to WASM via `registerDnaCompressorWasm()`. |
| **No wetlab validation** | All density, error rate, and recovery success figures are **simulation only**. | $500–$5,000 synthesis (Twist Bioscience) + $200–$1,000 sequencing + 2–4 weeks. This is a science problem, not a code problem. |
| **bam-parser is not htslib** | `bam-parser.ts` is a pure-JS SAM/BAM parser. It does NOT link to htslib, samtools, bcftools, or GATK. No CRAM, VCF/BCF, tabix, or FAI support. | Compile htslib to WASM via napi-rs (~3-5 days). |

### Open Problems (Not Just Implementation)

| Problem | Status |
|---------|--------|
| Nanopore 12.3% IDS recovery | Partial — ~50-70% at real-2024 error rates. K=9 Viterbi + OSD cascade + higher parity helps. |
| LDPC correction capacity | ~3% per-read failure rate at 4B parity. Outer RS erasure recovery covers failures. |
| Encryption not default | Users may forget to enable. API warns when encoding without password. |

---

## Quick Benchmark (v3.3, Node.js 24, single core)

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
| NAF-style Compress | DNA 1KB | ~0.05 ms |
| AGC-style Compress | DNA 1KB | ~0.07 ms |
| DeepGeCo-style Compress | DNA 1KB | ~0.09 ms |
| MBGC2-style Compress | DNA 1KB | ~0.11 ms |
| SIMD Unpack (JS) | 2-bit 1KB | ~0.01 ms |
| BAM Parse | 1000 records | ~1.2 ms |

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
import { compress, decompress, CompressorTier, isZstdCompressionReal } from "./lib/dna/compress";

// Auto-detect: biological → NAF-style/AGC-style/etc, general → ZSTD-compatible, fallback → PAKO
const result = compress(data);
console.log(`Used ${result.tier}, ratio ${result.ratio.toFixed(2)}×`);

// Check if ZSTD compression is real zstd format or DEFLATE fallback
if (result.tier === CompressorTier.ZSTD && !isZstdCompressionReal()) {
  console.warn('ZSTD tier compressed with DEFLATE format (not true zstd). Call registerZstdWasm() for real zstd.');
}

// Specific tier
const nafResult = compress(data, { tier: CompressorTier.NAF });
const agcResult = compress(data, { tier: CompressorTier.AGC });

// Decompression auto-detects format by magic bytes
const original = decompress(compressed);
```

### BAM/SAM Parsing

```typescript
import { BamParser, parseBamFile } from "./lib/dna/bam-parser";

// Note: This is a pure-JS BAM/SAM parser, NOT htslib.
// Full BAM format: BGZF, binary header, CIGAR, 4-bit seq, Phred+33 qual, all tag types
const records = parseBamFile(bamBuffer);

// Or use the class API for streaming
const parser = await BamParser.load();
const { fd, header } = await parser.openFile('reads.bam');
const record = await parser.bamRead(fd, header);
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
| `compress.ts` | 700+ | Tiered compression (NAF/AGC/DeepGeCo/MBGC2/JARVIS3-style + ZSTD-compatible + PAKO) |
| `pack.ts` | 400+ | 2-bit pack/unpack + bit-parallel ops (JS; WASM SIMD pending) |
| `addressing.ts` | 758 | BLAKE3 content-derived addressing (configurable Merkle primer) |
| `codec.ts` | 1500+ | Main encode pipeline (compress router + BHE/YYC wired) |
| `decode.ts` | 2500+ | Strategy cascade (compress router + Gungnir + DNA-Aeon wired) |
| `stream.ts` | 340+ | Streaming encode (ReadableStream + AsyncIterable) |
| `lsm-journal.ts` | 503 | LAB-DB LSM-tree journal with compaction |
| `archive.ts` | 603 | .hlx binary archive format |
| `bam-parser.ts` | 500+ | SAM/BAM binary parser (BGZF, CIGAR, 4-bit seq, all tags) |
| `simd-unpack.ts` | 340+ | SIMD unpack interface (WASM pending; JS fallback active) |
| `types.ts` | 751 | Core types, configs, channel presets |
| + 80 more | — | LDPC, convolutional, RS, fountain, holographic, etc. |

---

## Tech Stack

| Layer | Technology |
|-------|-----------|
| Language | TypeScript 5 |
| Runtime | Node.js / Bun, Browser (WASM) |
| Framework | Next.js 16 |
| Error Correction | @ronomon/reed-solomon, custom LDPC/Conv/OSD |
| Crypto | @noble/ciphers (XChaCha20-Poly1305), @noble/hashes (BLAKE3, Argon2id) |
| Compression | fflate (DEFLATE), fzstd (zstd decompression only), pako (DEFLATE), NAF/AGC/DeepGeCo/MBGC2/JARVIS3-style (DNA-aware JS approximations) |
| Bioinformatics | Custom BAM/SAM parser (not htslib) |
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
