//! helix-dna-napi: Native Node.js addon for K=9 Viterbi decoder
//!
//! v3: Production-hardened indel-tolerant Viterbi with trellis-step iteration.
//!
//! Key design: matches the WASM/TypeScript reference implementation exactly.
//!   - Receives num_info_bits explicitly (like WASM) OR estimates from received length
//!   - Iterates for exactly totalChannelUses = (numInfoBits + memory) * rate steps
//!   - I-chain propagation after each step
//!   - Correct traceback: I stays at same step, M/D advance
//!   - Zero-tail + drift penalty in final state selection
//!   - LLR soft-decision support

use napi::bindgen_prelude::*;
use napi_derive::napi;
use std::sync::OnceLock;

// ===========================================================================
// Convolutional Code Configuration
// ===========================================================================

const K9_MEMORY: usize = 8;
const K9_GENERATORS: [u16; 2] = [0o561, 0o753]; // NASA standard
const K7_MEMORY: usize = 6;
const K7_GENERATORS: [u16; 2] = [0o171, 0o133]; // Voyager

static K9_TABLE: OnceLock<TransitionTable> = OnceLock::new();
static K7_TABLE: OnceLock<TransitionTable> = OnceLock::new();
fn k9_table() -> &'static TransitionTable {
    K9_TABLE.get_or_init(|| TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2))
}
fn k7_table() -> &'static TransitionTable {
    K7_TABLE.get_or_init(|| TransitionTable::build(K7_MEMORY, &K7_GENERATORS, 2))
}

struct TransitionTable {
    outputs: Vec<u16>,
    next_states: Vec<u16>,
    num_states: usize,
    memory: usize,
    rate: usize,
}

impl TransitionTable {
    fn build(memory: usize, generators: &[u16], rate: usize) -> Self {
        let ns = 1usize << memory;
        let mut outputs = vec![0u16; ns * 2];
        let mut next_states = vec![0u16; ns * 2];
        for s in 0..ns {
            for ib in 0..2 {
                next_states[s * 2 + ib] = ((s >> 1) | ((ib as usize) << (memory - 1))) as u16;
                let mut ob = 0u16;
                for g in 0..rate {
                    let mut b = 0u16;
                    if (generators[g] & 1) != 0 { b ^= ib as u16; }
                    for i in 0..memory {
                        if (generators[g] & (1 << (i + 1))) != 0 { b ^= ((s >> i) & 1) as u16; }
                    }
                    ob |= b << g;
                }
                outputs[s * 2 + ib] = ob;
            }
        }
        TransitionTable { outputs, next_states, num_states: ns, memory, rate }
    }
    #[inline] fn output(&self, s: usize, i: usize) -> u16 { self.outputs[s * 2 + i] }
    #[inline] fn next_state(&self, s: usize, i: usize) -> usize { self.next_states[s * 2 + i] as usize }
}

// ===========================================================================
// Encoder
// ===========================================================================

fn conv_encode(data: &[u8], table: &TransitionTable) -> Vec<u8> {
    let ob = (data.len() * 8 + table.memory) * table.rate;
    let mut r = vec![0u8; (ob + 7) / 8];
    let mut st = 0usize; let mut p = 0usize;
    for &b in data {
        for bi in 0..8 {
            let ib = ((b >> (7 - bi)) & 1) as usize;
            let o = table.output(st, ib); st = table.next_state(st, ib);
            for g in 0..table.rate { if (o >> (table.rate - 1 - g)) & 1 != 0 { r[p / 8] |= 1 << (7 - (p % 8)); } p += 1; }
        }
    }
    for _ in 0..table.memory {
        let o = table.output(st, 0); st = table.next_state(st, 0);
        for g in 0..table.rate { if (o >> (table.rate - 1 - g)) & 1 != 0 && p / 8 < r.len() { r[p / 8] |= 1 << (7 - (p % 8)); } p += 1; }
    }
    r.truncate((p + 7) / 8); r
}

#[inline]
fn get_bit(d: &[u8], p: usize) -> u8 {
    if p / 8 >= d.len() { 0 } else { (d[p / 8] >> (7 - (p % 8))) & 1 }
}

// ===========================================================================
// napi config
// ===========================================================================

#[napi(object)]
#[derive(Clone, Copy)]
pub struct ViterbiConfig {
    pub max_drift: Option<u32>,
    pub insertion_penalty: Option<f64>,
    pub deletion_penalty: Option<f64>,
    pub use_llr: Option<bool>,
    pub expected_length: Option<u32>,
    /// Number of information bits (including tail). If not provided,
    /// estimated from received length as (received_bits / rate) - memory.
    pub num_info_bits: Option<u32>,
}

impl Default for ViterbiConfig {
    fn default() -> Self {
        ViterbiConfig {
            max_drift: Some(10),
            insertion_penalty: Some(1.5),
            deletion_penalty: Some(1.0),
            use_llr: Some(false),
            expected_length: None,
            num_info_bits: None,
        }
    }
}

fn resolve_config(cfg: Option<&ViterbiConfig>) -> (u32, f64, f64, Option<usize>, Option<usize>) {
    let c = cfg.copied().unwrap_or_default();
    (
        c.max_drift.unwrap_or(10),
        c.insertion_penalty.unwrap_or(1.5),
        c.deletion_penalty.unwrap_or(1.0),
        c.expected_length.map(|v| v as usize),
        c.num_info_bits.map(|v| v as usize),
    )
}

// ===========================================================================
// Indel-Tolerant Viterbi Decoder — Production v3
// ===========================================================================

/// Indel-tolerant Viterbi decoder matching WASM/TypeScript reference exactly.
///
/// State layout (same as TS/WASM):
///   phase 0:  index = cs * nd + (drift + maxDrift)           [0 .. ns*nd)
///   phase 1:  index = ns*nd + (cs*2+pending)*nd + (drift+maxDrift)
///
/// Trellis iteration:
///   - totalInfoSteps = numInfoBits + memory (zero tail)
///   - totalChannelUses = totalInfoSteps * rate
///   - Iterate step 0..totalChannelUses
///   - At each step: M/D transitions → next step, I propagation at next step
///
/// Traceback:
///   - Walk from totalChannelUses back to 0
///   - I transitions: stay at same step
///   - M/D transitions: advance step
fn viterbi_decode_indel(
    received: &[u8],
    table: &TransitionTable,
    max_drift: u32,
    ins_pen: f64,
    del_pen: f64,
    llr_data: Option<&[f32]>,
    num_info_bits_override: Option<usize>,
) -> Vec<u8> {
    let ns = table.num_states;
    let memory = table.memory;
    let rate = table.rate;
    let md = max_drift as i32;
    let nd = 2 * max_drift as usize + 1;
    let rb = received.len() * 8;
    if rb == 0 { return Vec::new(); }

    // Compute numInfoBits: either from override or estimated from received length
    let num_info_bits = if let Some(nib) = num_info_bits_override {
        nib
    } else {
        // Estimate: for rate-1/2, received_bits ≈ (numInfoBits + memory) * rate
        // So numInfoBits ≈ (received_bits / rate) - memory
        // This is correct for clean channels; for noisy channels with indels,
        // the drift is bounded by maxDrift so the estimate is close.
        if rb >= memory * rate { (rb / rate) - memory } else { rb / rate }
    };

    let total_info_steps = num_info_bits + memory;
    let total_channel_uses = total_info_steps * rate;

    let p0 = ns * nd;           // phase-1 region starts here
    let ta = ns * 3 * nd;       // total augmented states per step
    let inf = f64::INFINITY;

    // Unpack received bits to flat array for faster access
    let mut recv_bits: Vec<u8> = vec![0u8; rb];
    for i in 0..rb { recv_bits[i] = (received[i / 8] >> (7 - (i % 8))) & 1; }

    // Unpack LLR if provided
    let has_llr = llr_data.is_some() && llr_data.unwrap().len() >= rb;
    let llr_f64: Vec<f64> = if has_llr {
        llr_data.unwrap().iter().map(|&v| v as f64).collect()
    } else {
        Vec::new()
    };

    // Allocate step-indexed trellis storage
    // max_steps: total_channel_uses + extra for insertions at end
    let max_steps = total_channel_uses + max_drift as usize + 10;
    let total_cells = (max_steps + 1) * ta;
    let mut path_metric = vec![inf; total_cells];
    let mut back_ptr = vec![0u32; total_cells];     // aug index of predecessor
    let mut trans_type = vec![0u8; total_cells];     // 0=M, 1=I, 2=D
    let mut input_bit = vec![-1i8; total_cells];     // -1=none, 0/1=input committed

    // Initialize: (cs=0, phase=0, pending=0, drift=0) at step 0
    let start_aug = max_drift as usize; // cs=0, phase=0, drift=0
    path_metric[start_aug] = 0.0;

    // --- I-chain propagation helper ---
    // Processes insertion transitions within a single step's offset.
    // Iterates drifts in INCREASING order so I→I→I… chains propagate in one pass.
    let propagate_insertions = |step_offset: usize, step_for_recv: i32,
        path_metric: &mut Vec<f64>, back_ptr: &mut Vec<u32>,
        trans_type: &mut Vec<u8>, input_bit: &mut Vec<i8>| {
        for cs in 0..ns {
            // Phase 0 insertions
            for di in 0..nd - 1 {
                let drift = di as i32 - md;
                let recv_pos = step_for_recv + drift;
                if recv_pos < 0 || recv_pos as usize >= rb { continue; }
                let aug = cs * nd + di;
                let m = path_metric[step_offset + aug];
                if m >= 1e29 { continue; }
                let next_aug = cs * nd + di + 1;
                let new_metric = m + ins_pen;
                if new_metric < path_metric[step_offset + next_aug] {
                    path_metric[step_offset + next_aug] = new_metric;
                    back_ptr[step_offset + next_aug] = aug as u32;
                    trans_type[step_offset + next_aug] = 1; // I
                    input_bit[step_offset + next_aug] = -1;
                }
            }
            // Phase 1 insertions (both pending values)
            for pend in 0..2usize {
                for di in 0..nd - 1 {
                    let drift = di as i32 - md;
                    let recv_pos = step_for_recv + drift;
                    if recv_pos < 0 || recv_pos as usize >= rb { continue; }
                    let aug = p0 + (cs * 2 + pend) * nd + di;
                    let m = path_metric[step_offset + aug];
                    if m >= 1e29 { continue; }
                    let next_aug = p0 + (cs * 2 + pend) * nd + di + 1;
                    let new_metric = m + ins_pen;
                    if new_metric < path_metric[step_offset + next_aug] {
                        path_metric[step_offset + next_aug] = new_metric;
                        back_ptr[step_offset + next_aug] = aug as u32;
                        trans_type[step_offset + next_aug] = 1; // I
                        input_bit[step_offset + next_aug] = -1;
                    }
                }
            }
        }
    };

    // Run insertions at step 0 (handles insertions BEFORE the first M/D transition)
    propagate_insertions(0, 0, &mut path_metric, &mut back_ptr, &mut trans_type, &mut input_bit);

    // --- Main trellis recursion ---
    // Iterate for exactly totalChannelUses steps (matching WASM/TS reference).
    for step in 0..total_channel_uses {
        let step_offset = step * ta;
        let next_offset = (step + 1) * ta;

        // M and D transitions
        for cs in 0..ns {
            // Phase 0: encoder is about to emit G1. Try both input bits.
            for di in 0..nd {
                let aug = cs * nd + di;
                let m = path_metric[step_offset + aug];
                if m >= 1e29 { continue; }
                let drift = di as i32 - md;
                let recv_pos = step as i32 + drift;

                for ib in 0..2usize {
                    let out = table.output(cs, ib) as usize;
                    // Phase 0 emits the FIRST output bit = (out >> (rate-1)) & 1
                    // Encoder emits: for g in 0..rate { emit (out >> (rate-1-g)) & 1 }
                    // So g=0 (first) = MSB = G1 = (out >> (rate-1)) & 1
                    let emitted = (out >> (rate - 1)) & 1;

                    // M transition: emit first output bit, consume 1 received bit, move to phase 1
                    if recv_pos >= 0 && (recv_pos as usize) < rb {
                        let rp = recv_pos as usize;
                        let rb_bit = recv_bits[rp] as usize;
                        let dist = if has_llr {
                            llr_f64[rp].abs() * if emitted == rb_bit { -1.0 } else { 1.0 }
                        } else if emitted != rb_bit { 1.0 } else { 0.0 };
                        let next_aug = p0 + (cs * 2 + ib) * nd + di;
                        let new_metric = m + dist;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = aug as u32;
                            trans_type[next_offset + next_aug] = 0; // M
                            input_bit[next_offset + next_aug] = -1;
                        }
                    }

                    // D transition: emit G1, consume 0 bits, move to phase 1, drift-1
                    if di > 0 {
                        let next_aug = p0 + (cs * 2 + ib) * nd + di - 1;
                        let new_metric = m + del_pen;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = aug as u32;
                            trans_type[next_offset + next_aug] = 2; // D
                            input_bit[next_offset + next_aug] = -1;
                        }
                    }
                }
            }

            // Phase 1: encoder is about to emit G2. Use pending_input.
            for pend in 0..2usize {
                for di in 0..nd {
                    let aug = p0 + (cs * 2 + pend) * nd + di;
                    let m = path_metric[step_offset + aug];
                    if m >= 1e29 { continue; }
                    let drift = di as i32 - md;
                    let recv_pos = step as i32 + drift;

                    let out = table.output(cs, pend) as usize;
                    let next_state = table.next_state(cs, pend);
                    // Phase 1 emits the SECOND output bit = (out >> (rate-2)) & 1
                    // Encoder emits: g=0 (first at phase 0) = G1, g=1 (second at phase 1) = G2
                    let emitted = (out >> (rate - 2)) & 1;

                    // M transition: emit second output bit, consume 1 received bit, move to phase 0, commit input
                    if recv_pos >= 0 && (recv_pos as usize) < rb {
                        let rp = recv_pos as usize;
                        let rb_bit = recv_bits[rp] as usize;
                        let dist = if has_llr {
                            llr_f64[rp].abs() * if emitted == rb_bit { -1.0 } else { 1.0 }
                        } else if emitted != rb_bit { 1.0 } else { 0.0 };
                        let next_aug = next_state * nd + di;
                        let new_metric = m + dist;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = aug as u32;
                            trans_type[next_offset + next_aug] = 0; // M
                            input_bit[next_offset + next_aug] = pend as i8;
                        }
                    }

                    // D transition: emit G2, consume 0 bits, move to phase 0, drift-1, commit input
                    if di > 0 {
                        let next_aug = next_state * nd + di - 1;
                        let new_metric = m + del_pen;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = aug as u32;
                            trans_type[next_offset + next_aug] = 2; // D
                            input_bit[next_offset + next_aug] = pend as i8;
                        }
                    }
                }
            }
        }

        // I-chain propagation at next step
        propagate_insertions(next_offset, step as i32 + 1,
            &mut path_metric, &mut back_ptr, &mut trans_type, &mut input_bit);
    }

    // === Find best final state ===
    // Only consider phase-0 states at totalChannelUses.
    // Apply penalties: conv_state≠0 → +50 (zero tail), |drift|→0 → +0.5*|drift|
    let final_offset = total_channel_uses * ta;
    let mut best_aug = 0usize;
    let mut best_metric = inf;

    for cs in 0..ns {
        for di in 0..nd {
            let aug = cs * nd + di; // phase 0
            let m = path_metric[final_offset + aug];
            if m >= 1e29 { continue; }
            let drift = (di as i32 - md).abs() as f64;
            let penalty = if cs != 0 { 50.0 } else { 0.0 } + drift * 0.5;
            let total = m + penalty;
            if total < best_metric {
                best_metric = total;
                best_aug = aug;
            }
        }
    }

    if best_metric >= 1e29 {
        // No valid path found — return zeros of expected length
        return vec![0u8; num_info_bits / 8];
    }

    // === Traceback ===
    // Walk from totalChannelUses back to step 0.
    // I transitions stay at same step; M/D transitions advance step.
    let mut decoded_bits: Vec<u8> = Vec::with_capacity(num_info_bits);
    let mut step = total_channel_uses;
    let mut aug = best_aug;
    let mut safety = total_channel_uses * 4;

    while step > 0 && safety > 0 {
        safety -= 1;
        let offset = step * ta + aug;
        let tt = trans_type[offset];
        let ib = input_bit[offset];
        let prev_aug = back_ptr[offset] as usize;

        match tt {
            0 | 2 => {
                // M or D: step advances by 1
                if ib >= 0 {
                    decoded_bits.push(ib as u8);
                }
                step -= 1;
                aug = prev_aug;
            }
            1 => {
                // I: step stays the same
                aug = prev_aug;
            }
            _ => {
                step -= 1;
                aug = prev_aug;
            }
        }
    }

    decoded_bits.reverse();

    // Remove tail bits (last `memory` bits are zero-tail)
    if decoded_bits.len() > memory {
        decoded_bits.truncate(decoded_bits.len() - memory);
    }

    // Pack bits to bytes
    let nb = decoded_bits.len() / 8;
    let mut result = vec![0u8; nb];
    for i in 0..nb * 8 {
        if i < decoded_bits.len() && decoded_bits[i] != 0 {
            result[i / 8] |= 1 << (7 - (i % 8));
        }
    }
    result
}

// ===========================================================================
// napi Public API
// ===========================================================================

#[napi]
pub fn viterbi_k9_decode(received: Buffer, config: Option<ViterbiConfig>) -> Result<Buffer> {
    let (md, ip, dp, el, nib) = resolve_config(config.as_ref());
    let mut r = viterbi_decode_indel(received.as_ref(), k9_table(), md, ip, dp, None, nib);
    if let Some(l) = el { r.truncate(l); }
    Ok(r.into())
}

#[napi]
pub fn viterbi_k7_decode(received: Buffer, config: Option<ViterbiConfig>) -> Result<Buffer> {
    let (md, ip, dp, el, nib) = resolve_config(config.as_ref());
    let mut r = viterbi_decode_indel(received.as_ref(), k7_table(), md, ip, dp, None, nib);
    if let Some(l) = el { r.truncate(l); }
    Ok(r.into())
}

#[napi]
pub fn viterbi_k9_decode_with_llr(
    received: Buffer,
    llr: Float32Array,
    config: Option<ViterbiConfig>,
) -> Result<Buffer> {
    let (md, ip, dp, el, nib) = resolve_config(config.as_ref());
    let ls: Vec<f32> = llr.to_vec();
    let mut r = viterbi_decode_indel(received.as_ref(), k9_table(), md, ip, dp, Some(&ls), nib);
    if let Some(l) = el { r.truncate(l); }
    Ok(r.into())
}

#[napi]
pub fn conv_k9_encode(data: Buffer) -> Buffer {
    conv_encode(data.as_ref(), k9_table()).into()
}

#[napi]
pub fn conv_k7_encode(data: Buffer) -> Buffer {
    conv_encode(data.as_ref(), k7_table()).into()
}

/// Standard (non-indel) Viterbi decoder for K=9. ~0.5ms per oligo.
#[napi]
pub fn viterbi_k9_decode_standard(received: Buffer) -> Buffer {
    let t = k9_table();
    let ns = t.num_states;
    let rb = received.len() * 8;
    if rb == 0 { return Vec::<u8>::new().into(); }
    let ni = f64::INFINITY; let nst = rb / 2;
    let mut c = vec![ni; ns]; let mut p = vec![ni; ns];
    let mut tr: Vec<(u16, u8)> = vec![(0, 0); nst * ns]; c[0] = 0.0;
    let mut rp = 0usize;
    for s in 0..nst {
        std::mem::swap(&mut c, &mut p); c.fill(ni);
        let b0 = get_bit(received.as_ref(), rp); let b1 = get_bit(received.as_ref(), rp + 1); rp += 2;
        for st in 0..ns { if p[st] >= 1e29 { continue; } for ib in 0..2 {
            let o = t.output(st, ib);
            let co = p[st] + if ((o >> (t.rate - 1)) & 1) as u8 != b0 { 1.0 } else { 0.0 } + if ((o >> (t.rate - 2)) & 1) as u8 != b1 { 1.0 } else { 0.0 };
            let n = t.next_state(st, ib); if co < c[n] { c[n] = co; tr[s * ns + n] = (st as u16, ib as u8); }
        }}
    }
    let mut bs = 0usize; let mut bc = c[0]; for s in 1..ns { if c[s] < bc { bc = c[s]; bs = s; } }
    let mut db = Vec::with_capacity(nst); let mut st = bs;
    for s in (0..nst).rev() { let (ps, ib) = tr[s * ns + st]; db.push(ib); st = ps as usize; }
    db.reverse(); if db.len() > K9_MEMORY { db.truncate(db.len() - K9_MEMORY); }
    let nb = db.len() / 8; let mut r = vec![0u8; nb];
    for i in 0..nb * 8 { if i < db.len() && db[i] != 0 { r[i / 8] |= 1 << (7 - (i % 8)); } }
    r.into()
}

#[napi]
pub fn napi_version() -> String {
    format!("helix-dna-napi v0.3.0 — Viterbi v3 (trellis-step iteration, I-chain propagation, zero-tail penalty, numInfoBits, LLR)")
}
