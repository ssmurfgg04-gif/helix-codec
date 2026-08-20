---
Task ID: v69-rust-rewrites-pr84-fix
Agent: main
Task: Complete 5 Rust rewrites (pack/bhe/compress/ecc/simulate), fix PR #84 CI, attempt Yeast full 12MB

Work Log:
- PR #84 CI fix:
  - Investigated root cause via subagent — Goldman branch (codec.ts:653-678) set bestSatisfied=true unconditionally without calling satisfiesConstraints. With 4096B payload + goldman + default constraints, oligo 19 specifically drifted to GC=0.39.
  - Applied Option A (durable encoder fix): added seed-based retry loop to Goldman branch mirroring useSrt/useConstrained (codec.ts:656-697).
  - Also applied verification to deterministic branches that previously set bestSatisfied=true blindly:
    - useConvInner (codec.ts:652-655): now calls satisfiesConstraints
    - useBHE (codec.ts:866-870): now calls satisfiesConstraints
    - useArithmetic (codec.ts:845-848): now calls satisfiesConstraints
  - Added console.warn for silent screening failures (codec.ts:958-972) so future regressions surface in logs rather than only at downstream tests.
  - Verified scripts/test-codec.ts: "PASS: All oligos satisfy constraints (GC + homopolymer)" — CI test no longer fails.
  - The "Clean decode hash matches" failure is pre-existing (verified via git stash).

- Rust rewrites (5 modules) at rust/helix-dna-napi/src/:
  - pack.rs (~210 lines): packDnaToBits, unpackBitsToDna, complementPacked, reverseComplementPacked, bitParallelHamming, rollingHash, gcContent, maxHomopolymerRun. All tested PASS.
  - bhe.rs (~270 lines): bheEncode (k=1 fast path + k>1 FSM fallback), bheDecode. Tested PASS roundtrip.
  - compress.rs (~225 lines): compressZstd, decompressZstd, isAlreadyCompressed. Custom LZ77+RLE format with magic 0x48 0x4E 0x41 0x50 0x01 (HNAPv1). Tested PASS roundtrip.
  - ecc.rs (~345 lines): rsEncode, rsDecode, rsDecodeErasures, rsParity, rsVersion. RS encode/parity PASS. RS decode has Berlekamp-Massey bug (err_loc update edge case); JS fallback used in production for decode.
  - simulate.rs (~270 lines): simulateOligoReads, simulateBasic, readStats, simulateVersion. Uses Mulberry32 PRNG. Returns flat [coverage_u32, r0_len_u32, r0_bytes, ...] format. Tested PASS.
  - lib.rs: declared pub mod pack/bhe/compress/ecc/simulate.
  - Built successfully with cargo build --release. Output: target/release/libhelix_dna_napi.so (488KB).
  - All 28 napi exports verified loaded via direct dlopen test.

- TS wrappers and wiring:
  - Created src/lib/dna/native/helix-napi.ts — single loader for all 5 new modules + the existing Viterbi module. Reuses the same .so file.
  - Wired pack.ts: packDnaToBits, unpackBitsToDna, complement, reverseComplement, bitParallelHamming all use native FIRST PRIORITY, fall through to WASM then JS.
  - Wired compress.ts: decompress() detects native HNAP magic 0x48 0x4E 0x41 0x50 0x01 and dispatches to native decompressZstd. (compress() not yet wired to write native format — kept current behavior for backward compatibility with existing data; native decompress is a new path.)
  - Wired simulate.ts: simulate() (basic mode only — dt4dds still JS) uses native simulateOligoReads for bulk read simulation. Parses flat output into SequencingRead objects.
  - Wired reedsolomon.ts: encode() uses native FIRST PRIORITY; parity() inherits from encode(); decode() still uses JS (Rust BM has bug).

- Yeast full 12MB test:
  - Host has only 3.9GB total RAM (verified via `free -h`).
  - Tried with NODE_OPTIONS=--max-old-space-size=3072 — process killed (OOM).
  - Tried with --max-old-space-size=3960 --max-semi-space-size=128 — also OOMs.
  - E. coli K-12 (4.4MB) encodes successfully in ~73s in this env, but the 12MB Yeast triggers Node heap growth beyond physical RAM.
  - Conclusion: cannot run full 12MB Yeast in this 4GB environment. Validated via 4×1MB chunks in prior session (all 4 PASSED with roundtrip=true, hash=true).
  - To run full 12MB Yeast, recommend a host with ≥8GB RAM.

Stage Summary:
- PR #84 CI fix: COMPLETE — encoder retry loop + verification + warning; test passes.
- Rust rewrites: 5/5 modules built and exported. pack/bhe/compress/simulate fully working; RS encode working; RS decode has BM bug (JS fallback).
- TS wiring: 4/5 source files (pack.ts, compress.ts, simulate.ts, reedsolomon.ts) wired to use native FIRST PRIORITY. BHE encoder is already JS-only and wasn't rewired (kept separate).
- Yeast full: cannot run in 4GB env; chunked validation already PASS.
- Outstanding: cargo binary wipes on session restart (had to reinstall rustup this session). Consider adding a Makefile/script that builds the addon automatically.
