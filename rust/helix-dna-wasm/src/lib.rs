//! Helix DNA Codec — Rust Hot Paths
//!
//! This crate implements the CPU-bound modules of helix-codec in Rust,
//! compiled to WASM with SIMD for browser and Node.js via napi-rs.
//!
//! Modules:
//!   - viterbi:  K=9/K=7 indel-tolerant Viterbi decoder (THE #1 hot path)
//!   - pack:    2-bit DNA pack/unpack with WASM SIMD (v128/i8x16)
//!   - ecc:     Reed-Solomon GF(256) + LDPC belief propagation with SIMD
//!   - compress: Arithmetic coding with context models
//!   - simulate: Parallel per-oligo stochastic simulation
//!   - bhe:     Bit-parallel FSM encoding (u128 instead of BigInt)

pub mod viterbi;
pub mod pack;
pub mod ecc;
pub mod compress;
pub mod bhe;
pub mod simulate;

use wasm_bindgen::prelude::*;

/// Package version
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}
