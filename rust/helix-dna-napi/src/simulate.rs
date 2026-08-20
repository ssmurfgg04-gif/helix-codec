//! simulate.rs — Wetlab sequencing simulation
//!
//! Native Rust port of the hot paths from `src/lib/dna/simulate.ts` and
//! `src/lib/dna/dt4dds-simulate.ts`.
//!
//! Models:
//!   1. Substitution / insertion / deletion errors per position
//!   2. Coverage: each oligo is read N times with independent errors
//!   3. Dropout: fraction of oligos completely lost
//!   4. Position-dependent synthesis errors (5'/3' end effects)

use napi::bindgen_prelude::*;
use napi_derive::napi;

// Mulberry32 PRNG for reproducibility
struct Rng {
    state: u32,
}
impl Rng {
    fn new(seed: u32) -> Self {
        Self { state: seed }
    }
    fn next(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x6D2B79F5);
        let mut t = self.state;
        t = (t ^ (t >> 15)).wrapping_mul(t | 1);
        t ^= t.wrapping_add(t.wrapping_mul(t ^ (t >> 7)));
        ((t ^ (t >> 14)) >> 0) as f64 / 4294967296.0
    }
    fn next_u32(&mut self) -> u32 {
        self.state = self.state.wrapping_add(0x6D2B79F5);
        let mut t = self.state;
        t = (t ^ (t >> 15)).wrapping_mul(t | 1);
        t ^= t.wrapping_add(t.wrapping_mul(t ^ (t >> 7)));
        t ^ (t >> 14)
    }
    fn range(&mut self, lo: u32, hi: u32) -> u32 {
        if hi <= lo {
            return lo;
        }
        lo + (self.next_u32() % (hi - lo))
    }
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

#[napi(object)]
pub struct SimulationConfig {
    pub substitution_rate: f64,
    pub insertion_rate: f64,
    pub deletion_rate: f64,
    pub coverage: u32,
    pub dropout_rate: f64,
    pub seed: u32,
    /// Optional: position-dependent multiplier at 5' end (default 1.5)
    pub five_prime_mult: Option<f64>,
    /// Optional: position-dependent multiplier at 3' end (default 2.0)
    pub three_prime_mult: Option<f64>,
    /// Optional: enable position-dependent error scaling
    pub position_dependent: Option<bool>,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            substitution_rate: 0.001,
            insertion_rate: 0.0005,
            deletion_rate: 0.001,
            coverage: 20,
            dropout_rate: 0.0,
            seed: 0,
            five_prime_mult: Some(1.5),
            three_prime_mult: Some(2.0),
            position_dependent: Some(false),
        }
    }
}

/// Apply per-position errors (substitution/insertion/deletion) to a single
/// read. Returns the corrupted read.
///
/// `position_mult` is a function that takes (pos, len) and returns the error
/// multiplier for that position. Used for 5'/3' end effects.
fn apply_errors(seq: &[u8], cfg: &SimulationConfig, rng: &mut Rng, position_mult: impl Fn(usize, usize) -> f64) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(seq.len() + 8);
    let len = seq.len();
    for (i, &base) in seq.iter().enumerate() {
        let pm = position_mult(i, len);
        let sub_rate = cfg.substitution_rate * pm;
        let ins_rate = cfg.insertion_rate * pm;
        let del_rate = cfg.deletion_rate * pm;
        let r = rng.next();
        if r < del_rate {
            // Deletion: skip this base
            continue;
        }
        // Insertion (before this base)
        let r2 = rng.next();
        if r2 < ins_rate {
            // Insert a random base (not the current one, to be more realistic)
            let ins_base = BASES[rng.range(0, 4) as usize];
            out.push(ins_base);
        }
        // Substitution
        let r3 = rng.next();
        if r3 < sub_rate {
            // Substitute with a different base
            let mut new_base;
            loop {
                new_base = BASES[rng.range(0, 4) as usize];
                if new_base != base {
                    break;
                }
            }
            out.push(new_base);
        } else {
            out.push(base);
        }
    }
    out
}

/// Simulate sequencing reads for a single oligo.
///
/// Returns a flat array of (read_idx, read_seq) pairs:
///   [r0_len, r0_seq_bytes, r1_len, r1_seq_bytes, ...]
///
/// This flat format avoids the overhead of nested arrays across the FFI.
#[napi]
pub fn simulate_oligo_reads(
    oligo_seq: String,
    config: Option<SimulationConfig>,
) -> Uint8Array {
    let cfg = config.unwrap_or_default();
    let seq = oligo_seq.into_bytes();
    if seq.is_empty() {
        return Uint8Array::new(vec![]);
    }

    // Check dropout
    let mut rng = Rng::new(if cfg.seed == 0 { 1 } else { cfg.seed });
    if cfg.dropout_rate > 0.0 && rng.next() < cfg.dropout_rate {
        // Oligo completely lost — return empty
        return Uint8Array::new(vec![]);
    }

    let coverage = cfg.coverage.max(1);
    let pos_dep = cfg.position_dependent.unwrap_or(false);
    let five_m = cfg.five_prime_mult.unwrap_or(1.5);
    let three_m = cfg.three_prime_mult.unwrap_or(2.0);

    // Build flat output: [coverage_count_u32, r0_len_u32, r0_bytes, r1_len_u32, r1_bytes, ...]
    let mut out: Vec<u8> = Vec::with_capacity(4 + coverage as usize * (4 + seq.len() + 8));
    // First 4 bytes: coverage count (little-endian u32)
    out.extend_from_slice(&coverage.to_le_bytes());

    let oligo_idx_hash = rng.next_u32(); // for per-oligo seed variation
    let seq_len = seq.len();
    for c in 0..coverage {
        let seed = if cfg.seed == 0 {
            // Non-deterministic: use time + counter
            (oligo_idx_hash.wrapping_add(c)).wrapping_mul(0x9E3779B9)
        } else {
            cfg.seed.wrapping_add(c * 1000).wrapping_add(oligo_idx_hash)
        };
        let mut read_rng = Rng::new(seed);
        let read = if pos_dep {
            apply_errors(&seq, &cfg, &mut read_rng, |pos, len| {
                // 5' end: first 10% of read has higher error
                // 3' end: last 20% of read has higher error
                let pos_frac = pos as f64 / len.max(1) as f64;
                if pos_frac < 0.1 {
                    five_m
                } else if pos_frac > 0.8 {
                    three_m
                } else {
                    1.0
                }
            })
        } else {
            apply_errors(&seq, &cfg, &mut read_rng, |_, _| 1.0)
        };
        // Write read length (u32 LE) then bytes
        let rlen = read.len() as u32;
        out.extend_from_slice(&rlen.to_le_bytes());
        out.extend_from_slice(&read);
    }

    Uint8Array::new(out)
}

/// Simulate a basic substitution-only channel (for quick testing).
/// Returns a single corrupted sequence.
#[napi]
pub fn simulate_basic(
    oligo_seq: String,
    sub_rate: f64,
    ins_rate: f64,
    del_rate: f64,
    seed: u32,
) -> String {
    let cfg = SimulationConfig {
        substitution_rate: sub_rate,
        insertion_rate: ins_rate,
        deletion_rate: del_rate,
        coverage: 1,
        dropout_rate: 0.0,
        seed,
        five_prime_mult: Some(1.0),
        three_prime_mult: Some(1.0),
        position_dependent: Some(false),
    };
    let mut rng = Rng::new(if seed == 0 { 1 } else { seed });
    let seq = oligo_seq.into_bytes();
    let out = apply_errors(&seq, &cfg, &mut rng, |_, _| 1.0);
    String::from_utf8(out).unwrap_or_default()
}

/// Compute summary statistics for a set of reads.
/// Returns [num_reads, total_bases, mean_len, gc_fraction].
#[napi]
pub fn read_stats(reads: Uint8Array) -> Float64Array {
    let data: &[u8] = &reads;
    if data.len() < 4 {
        return Float64Array::new(vec![0.0, 0.0, 0.0, 0.0]);
    }
    let mut num_reads: u32 = 0;
    let mut total_bases: u64 = 0;
    let mut gc: u64 = 0;
    let mut i = 4; // skip coverage count
    while i + 4 <= data.len() {
        let rlen = u32::from_le_bytes([data[i], data[i + 1], data[i + 2], data[i + 3]]) as usize;
        i += 4;
        if i + rlen > data.len() {
            break;
        }
        num_reads += 1;
        total_bases += rlen as u64;
        for j in 0..rlen {
            let b = data[i + j];
            if b == b'G' || b == b'C' {
                gc += 1;
            }
        }
        i += rlen;
    }
    let mean_len = if num_reads > 0 {
        total_bases as f64 / num_reads as f64
    } else {
        0.0
    };
    let gc_frac = if total_bases > 0 {
        gc as f64 / total_bases as f64
    } else {
        0.0
    };
    Float64Array::new(vec![num_reads as f64, total_bases as f64, mean_len, gc_frac])
}

/// Return the napi-rs version string for simulate.
#[napi]
pub fn simulate_version() -> String {
    String::from("helix-dna-napi simulate.rs v1.0 — wetlab simulation (sub/ins/del + coverage + dropout + position-dep)")
}
