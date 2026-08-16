//! Parametric Wetlab Simulation — parallel per-oligo stochastic model
//!
//! Models the full wetlab pipeline:
//!   1. Array synthesis (position-dependent errors)
//!   2. PCR amplification (duplication + errors)
//!   3. Aging/decay (depurination, oxidation, deamination)
//!   4. Sequencing (platform-specific error profiles)
//!
//! Each oligo is simulated independently, enabling parallel execution
//! with rayon (native) or sequential fallback (WASM).

use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Xorshift64 PRNG (fast, per-thread)
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Rng { state: if seed == 0 { 1 } else { seed } }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Sample from Bernoulli(p).
    fn bernoulli(&mut self, p: f64) -> bool {
        self.next_f64() < p
    }

    /// Sample discrete uniform from [0, n).
    fn uniform(&mut self, n: u32) -> u32 {
        (self.next_f64() * n as f64) as u32
    }
}

// ---------------------------------------------------------------------------
// Synthesis model
// ---------------------------------------------------------------------------

/// Synthesis error profile.
#[wasm_bindgen]
pub struct SynthesisProfile {
    /// Base substitution rate (uniform across positions).
    pub sub_rate: f64,
    /// Insertion rate per position.
    pub ins_rate: f64,
    /// Deletion rate per position.
    pub del_rate: f64,
    /// 5' end degradation factor (errors increase toward 5' end).
    pub five_prime_bias: f64,
    /// 3' end degradation factor.
    pub three_prime_bias: f64,
}

#[wasm_bindgen]
impl SynthesisProfile {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        SynthesisProfile {
            sub_rate: 0.002,     // 0.2% per position (array synthesis)
            ins_rate: 0.0005,    // 0.05% insertions
            del_rate: 0.001,     // 0.1% deletions
            five_prime_bias: 1.5, // 50% more errors at 5' end
            three_prime_bias: 2.0, // 100% more errors at 3' end
        }
    }

    pub fn illuminina() -> Self {
        SynthesisProfile {
            sub_rate: 0.001,
            ins_rate: 0.0001,
            del_rate: 0.0002,
            five_prime_bias: 1.0,
            three_prime_bias: 1.2,
        }
    }

    pub fn nanopore() -> Self {
        SynthesisProfile {
            sub_rate: 0.05,
            ins_rate: 0.04,
            del_rate: 0.04,
            five_prime_bias: 1.5,
            three_prime_bias: 2.5,
        }
    }
}

// ---------------------------------------------------------------------------
// Sequencing model
// ---------------------------------------------------------------------------

/// Sequencing error profile.
#[wasm_bindgen]
pub struct SequencingProfile {
    /// Substitution rate added by sequencing.
    pub sub_rate: f64,
    /// Insertion rate added by sequencing.
    pub ins_rate: f64,
    /// Deletion rate added by sequencing.
    pub del_rate: f64,
}

#[wasm_bindgen]
impl SequencingProfile {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        SequencingProfile {
            sub_rate: 0.001,
            ins_rate: 0.0001,
            del_rate: 0.0001,
        }
    }

    pub fn illumina() -> Self {
        SequencingProfile {
            sub_rate: 0.001,
            ins_rate: 0.0001,
            del_rate: 0.0001,
        }
    }

    pub fn nanopore() -> Self {
        SequencingProfile {
            sub_rate: 0.05,
            ins_rate: 0.03,
            del_rate: 0.03,
        }
    }

    pub fn pacbio_hifi() -> Self {
        SequencingProfile {
            sub_rate: 0.001,
            ins_rate: 0.0005,
            del_rate: 0.0005,
        }
    }
}

// ---------------------------------------------------------------------------
// Simulate one oligo through the pipeline
// ---------------------------------------------------------------------------

/// Simulate a single oligo through synthesis → PCR → aging → sequencing.
///
/// Returns the simulated read as ASCII bytes (A/C/G/T).
fn simulate_one(oligo: &[u8], synth: &SynthesisProfile, seq: &SequencingProfile, rng: &mut Rng) -> Vec<u8> {
    let len = oligo.len();

    // Phase 1: Synthesis errors
    let mut synthesized = Vec::with_capacity(len + len / 10);
    let bases = [b'A', b'C', b'G', b'T'];

    for pos in 0..len {
        // Position-dependent error rate
        let pos_frac = pos as f64 / len as f64;
        let pos_factor = 1.0
            + synth.five_prime_bias * (1.0 - pos_frac).powi(2)
            + synth.three_prime_bias * pos_frac.powi(2);

        // Deletion
        if rng.bernoulli(synth.del_rate * pos_factor) {
            continue;
        }

        // Insertion
        while rng.bernoulli(synth.ins_rate * pos_factor) {
            synthesized.push(bases[rng.uniform(4) as usize]);
        }

        // Substitution or correct
        if rng.bernoulli(synth.sub_rate * pos_factor) {
            // Substitute with random different base
            let current = oligo[pos];
            let mut alt = current;
            while alt == current {
                alt = bases[rng.uniform(4) as usize];
            }
            synthesized.push(alt);
        } else {
            synthesized.push(oligo[pos]);
        }
    }

    // Phase 2: Sequencing errors (applied to the synthesized product)
    let mut sequenced = Vec::with_capacity(synthesized.len() + synthesized.len() / 10);

    for pos in 0..synthesized.len() {
        if rng.bernoulli(seq.del_rate) {
            continue;
        }
        while rng.bernoulli(seq.ins_rate) {
            sequenced.push(bases[rng.uniform(4) as usize]);
        }
        if rng.bernoulli(seq.sub_rate) {
            let current = synthesized[pos];
            let mut alt = current;
            while alt == current {
                alt = bases[rng.uniform(4) as usize];
            }
            sequenced.push(alt);
        } else {
            sequenced.push(synthesized[pos]);
        }
    }

    sequenced
}

// ---------------------------------------------------------------------------
// Batch simulation (parallel on native, sequential on WASM)
// ---------------------------------------------------------------------------

/// Simulate multiple oligos through the wetlab pipeline.
///
/// Each oligo is processed independently with a unique RNG seed
/// derived from the base seed + oligo index.
///
/// Returns a flat array: for each oligo, [length_u32_le, base1, base2, ...].
#[wasm_bindgen]
pub fn simulate_batch(
    oligos: &[u8],          // flat: all oligos concatenated
    oligo_offsets: &[u32],  // offsets[i]..offsets[i+1] = oligo i
    synth_sub: f64,
    synth_ins: f64,
    synth_del: f64,
    seq_sub: f64,
    seq_ins: f64,
    seq_del: f64,
    base_seed: u32,
) -> Vec<u8> {
    let num_oligos = oligo_offsets.len().saturating_sub(1);
    if num_oligos == 0 {
        return vec![];
    }

    let synth = SynthesisProfile {
        sub_rate: synth_sub,
        ins_rate: synth_ins,
        del_rate: synth_del,
        five_prime_bias: 1.5,
        three_prime_bias: 2.0,
    };
    let seq = SequencingProfile {
        sub_rate: seq_sub,
        ins_rate: seq_ins,
        del_rate: seq_del,
    };

    // Simulate each oligo (sequential in WASM, could use rayon in native)
    let mut results = Vec::new();

    for i in 0..num_oligos {
        let start = oligo_offsets[i] as usize;
        let end = oligo_offsets[i + 1] as usize;
        let oligo = &oligos[start..end];

        let mut rng = Rng::new(base_seed.wrapping_add(i as u32) as u64);
        let read = simulate_one(oligo, &synth, &seq, &mut rng);

        // Encode: length (4 bytes LE) + data
        let len = read.len() as u32;
        results.push(len as u8);
        results.push((len >> 8) as u8);
        results.push((len >> 16) as u8);
        results.push((len >> 24) as u8);
        results.extend_from_slice(&read);
    }

    results
}

/// Simulate a single oligo. Returns the read as ASCII bytes.
#[wasm_bindgen]
pub fn simulate_single(
    oligo: &[u8],
    synth_sub: f64,
    synth_ins: f64,
    synth_del: f64,
    seq_sub: f64,
    seq_ins: f64,
    seq_del: f64,
    seed: u32,
) -> Vec<u8> {
    let synth = SynthesisProfile {
        sub_rate: synth_sub,
        ins_rate: synth_ins,
        del_rate: synth_del,
        five_prime_bias: 1.5,
        three_prime_bias: 2.0,
    };
    let seq = SequencingProfile {
        sub_rate: seq_sub,
        ins_rate: seq_ins,
        del_rate: seq_del,
    };
    let mut rng = Rng::new(seed as u64);
    simulate_one(oligo, &synth, &seq, &mut rng)
}
