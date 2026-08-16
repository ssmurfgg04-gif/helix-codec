//! K=9 Indel-Tolerant Viterbi Decoder — Rust hot path
//!
//! NASA standard convolutional code (memory=8, K=9, d_free=24):
//!   G1 = 561 (octal) = 0x171
//!   G2 = 753 (octal) = 0x1EB
//!
//! The augmented trellis tracks (conv_state, bit_phase, pending_input, drift).
//! For K=9: 256 conv states × 3 phase layouts × (2*maxDrift+1) drifts.
//!
//! This is the #1 hot path of the codec — ~800ms/oligo in JS, ~5ms in Rust.

use wasm_bindgen::prelude::*;

// ===========================================================================
// Convolutional Code Configuration
// ===========================================================================

/// NASA K=9 rate-1/2 convolutional code (d_free=24)
const K9_MEMORY: usize = 8;
const K9_NUM_STATES: usize = 1 << K9_MEMORY; // 256
const K9_GENERATORS: [u16; 2] = [0o561, 0o753]; // NASA standard
const K9_RATE: usize = 2;

/// Voyager K=7 rate-1/2 convolutional code (d_free=10)
const K7_MEMORY: usize = 6;
const K7_NUM_STATES: usize = 1 << K7_MEMORY; // 64
const K7_GENERATORS: [u16; 2] = [0o171, 0o133];

/// Transition table: precomputed outputs and next states for all (state, input).
struct TransitionTable {
    /// outputs[state * 2 + input] = output bits (rate bits packed in low bits)
    outputs: Vec<u16>,
    /// next_states[state * 2 + input] = next state
    next_states: Vec<u16>,
    num_states: usize,
    memory: usize,
    rate: usize,
}

impl TransitionTable {
    fn build(memory: usize, generators: &[u16], rate: usize) -> Self {
        let num_states = 1usize << memory;
        let mut outputs = vec![0u16; num_states * 2];
        let mut next_states = vec![0u16; num_states * 2];

        for state in 0..num_states {
            for input in 0..2usize {
                let reg = (input << memory) | state;
                let mut output = 0u16;
                for g in 0..rate {
                    let gen = generators[g];
                    let mut bit = 0u16;
                    for b in 0..=memory {
                        if (gen >> b) & 1 != 0 {
                            bit ^= (reg >> b) as u16 & 1;
                        }
                    }
                    output = (output << 1) | bit;
                }
                let next_state = (reg >> 1) & (num_states - 1);
                let idx = state * 2 + input;
                outputs[idx] = output;
                next_states[idx] = next_state as u16;
            }
        }

        TransitionTable { outputs, next_states, num_states, memory, rate }
    }
}

// Lazy-initialized transition tables
use std::sync::LazyLock;
static K9_TABLE: LazyLock<TransitionTable> = LazyLock::new(|| {
    TransitionTable::build(K9_MEMORY, &K9_GENERATORS, K9_RATE)
});
static K7_TABLE: LazyLock<TransitionTable> = LazyLock::new(|| {
    TransitionTable::build(K7_MEMORY, &K7_GENERATORS, K9_RATE)
});

// ===========================================================================
// Indel-Tolerant Viterbi Decoder
// ===========================================================================

/// Augmented state index layout:
///   phase 0: (cs * 3 + 0) * W + (drift + max_drift)
///   phase 1: (cs * 3 + 1 + pending) * W + (drift + max_drift)
///
/// States per conv state: 3 * W where W = 2*max_drift+1
/// Total augmented states: num_conv_states * 3 * W

#[inline(always)]
fn aug_index(cs: usize, phase: usize, pending: usize, drift: i32, max_drift: i32, w: usize) -> usize {
    let drift_idx = (drift + max_drift) as usize;
    if phase == 0 {
        (cs * 3 + 0) * w + drift_idx
    } else {
        (cs * 3 + 1 + pending) * w + drift_idx
    }
}

/// Decode soft-decision LLR-weighted received bits through the indel-tolerant Viterbi trellis.
///
/// # Arguments
/// * `received_bits` — hard-decision received bits (0 or 1)
/// * `received_llr` — per-bit log-likelihood ratios (|LLR| = confidence; sign indicates hard decision)
/// * `num_info_bits` — number of information bits to decode
/// * `memory` — convolutional code memory (6 for K=7, 8 for K=9)
/// * `max_drift` — maximum net indel drift to track
/// * `insertion_penalty` — cost of an insertion transition
/// * `deletion_penalty` — cost of a deletion transition
///
/// # Returns
/// Decoded information bits as Vec<u8>
fn viterbi_decode(
    received_bits: &[u8],
    received_llr: &[f64],
    num_info_bits: usize,
    memory: usize,
    max_drift: i32,
    insertion_penalty: f64,
    deletion_penalty: f64,
) -> Vec<u8> {
    // Select transition table based on memory
    let tbl = if memory == 8 { &*K9_TABLE } else { &*K7_TABLE };
    let num_conv_states = tbl.num_states;
    let rate = tbl.rate;

    let total_info_steps = num_info_bits + memory;
    let total_channel_uses = total_info_steps * rate;

    let w = (2 * max_drift + 1) as usize;
    let num_aug = num_conv_states * 3 * w;
    let max_steps = total_channel_uses + max_drift as usize + 10;
    let total_cells = max_steps * num_aug;

    // Allocate trellis storage
    let mut path_metric = vec![f64::INFINITY; total_cells];
    let mut back_ptr = vec![0u32; total_cells]; // packed: step * num_aug + aug
    let mut trans_type = vec![0u8; total_cells]; // 0=M, 1=I, 2=D
    let mut input_bit = vec![0i8; total_cells]; // -1=none, 0/1=input committed

    // Initialize start state
    let start_aug = aug_index(0, 0, 0, 0, max_drift, w);
    path_metric[start_aug] = 0.0;

    let received_len = received_bits.len();
    let has_llr = received_llr.len() >= received_len;

    // Propagate insertion transitions within a step
    let propagate_insertions = |step_offset: usize, step_for_recv: usize,
        path_metric: &mut Vec<f64>, back_ptr: &mut Vec<u32>,
        trans_type: &mut Vec<u8>, input_bit: &mut Vec<i8>| {
        for cs in 0..num_conv_states {
            for drift in -max_drift..max_drift { // drift < max_drift (can increase by 1)
                let recv_pos = step_for_recv as i32 + drift;
                if recv_pos < 0 || recv_pos as usize >= received_len { continue; }

                // Phase 0
                let aug = aug_index(cs, 0, 0, drift, max_drift, w);
                let m = path_metric[step_offset + aug];
                if m.is_finite() {
                    let next_aug = aug_index(cs, 0, 0, drift + 1, max_drift, w);
                    let new_metric = m + insertion_penalty;
                    if new_metric < path_metric[step_offset + next_aug] {
                        path_metric[step_offset + next_aug] = new_metric;
                        back_ptr[step_offset + next_aug] = (aug) as u32;
                        trans_type[step_offset + next_aug] = 1; // I
                        input_bit[step_offset + next_aug] = -1;
                    }
                }

                // Phase 1 (both pending values)
                for pending in 0..2usize {
                    let aug = aug_index(cs, 1, pending, drift, max_drift, w);
                    let m = path_metric[step_offset + aug];
                    if !m.is_finite() { continue; }
                    let next_aug = aug_index(cs, 1, pending, drift + 1, max_drift, w);
                    let new_metric = m + insertion_penalty;
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

    // Run insertions at step 0
    propagate_insertions(0, 0, &mut path_metric, &mut back_ptr, &mut trans_type, &mut input_bit);

    // Main trellis recursion
    for step in 0..total_channel_uses {
        let step_offset = step * num_aug;
        let next_offset = (step + 1) * num_aug;

        // M and D transitions
        for cs in 0..num_conv_states {
            // Phase 0: try both inputs
            for drift in -max_drift..=max_drift {
                let aug = aug_index(cs, 0, 0, drift, max_drift, w);
                let m = path_metric[step_offset + aug];
                if !m.is_finite() { continue; }
                let recv_pos = (step as i32 + drift) as usize;

                for input in 0..2usize {
                    let idx = cs * 2 + input;
                    let output = tbl.outputs[idx] as usize;
                    let g1 = (output >> (rate - 1)) & 1;

                    // M transition: emit G1, consume 1 received bit
                    if recv_pos < received_len {
                        let rb = received_bits[recv_pos] as usize;
                        let dist = if has_llr {
                            received_llr[recv_pos].abs() * if g1 == rb { -1.0 } else { 1.0 }
                        } else {
                            if g1 != rb { 1.0 } else { 0.0 }
                        };
                        let next_aug = aug_index(cs, 1, input, drift, max_drift, w);
                        let new_metric = m + dist;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = (step_offset + aug) as u32;
                            trans_type[next_offset + next_aug] = 0; // M
                            input_bit[next_offset + next_aug] = -1;
                        }
                    }

                    // D transition: emit G1, consume 0 bits
                    if drift - 1 >= -max_drift {
                        let next_aug = aug_index(cs, 1, input, drift - 1, max_drift, w);
                        let new_metric = m + deletion_penalty;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = (step_offset + aug) as u32;
                            trans_type[next_offset + next_aug] = 2; // D
                            input_bit[next_offset + next_aug] = -1;
                        }
                    }
                }
            }

            // Phase 1: use pending_input
            for pending in 0..2usize {
                for drift in -max_drift..=max_drift {
                    let aug = aug_index(cs, 1, pending, drift, max_drift, w);
                    let m = path_metric[step_offset + aug];
                    if !m.is_finite() { continue; }
                    let recv_pos = (step as i32 + drift) as usize;

                    let idx = cs * 2 + pending;
                    let output = tbl.outputs[idx] as usize;
                    let next_state = tbl.next_states[idx] as usize;
                    let g2 = (output >> (rate - 2)) & 1;

                    // M transition: emit G2
                    if recv_pos < received_len {
                        let rb = received_bits[recv_pos] as usize;
                        let dist = if has_llr {
                            received_llr[recv_pos].abs() * if g2 == rb { -1.0 } else { 1.0 }
                        } else {
                            if g2 != rb { 1.0 } else { 0.0 }
                        };
                        let next_aug = aug_index(next_state, 0, 0, drift, max_drift, w);
                        let new_metric = m + dist;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = (step_offset + aug) as u32;
                            trans_type[next_offset + next_aug] = 0; // M
                            input_bit[next_offset + next_aug] = pending as i8;
                        }
                    }

                    // D transition: emit G2
                    if drift - 1 >= -max_drift {
                        let next_aug = aug_index(next_state, 0, 0, drift - 1, max_drift, w);
                        let new_metric = m + deletion_penalty;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = (step_offset + aug) as u32;
                            trans_type[next_offset + next_aug] = 2; // D
                            input_bit[next_offset + next_aug] = pending as i8;
                        }
                    }
                }
            }
        }

        // Propagate insertions at next step
        propagate_insertions(next_offset, step + 1, &mut path_metric, &mut back_ptr, &mut trans_type, &mut input_bit);
    }

    // === Traceback ===
    let final_offset = total_channel_uses * num_aug;
    let mut best_aug = 0usize;
    let mut best_metric = f64::INFINITY;

    for cs in 0..num_conv_states {
        for drift in -max_drift..=max_drift {
            let aug = aug_index(cs, 0, 0, drift, max_drift, w);
            let m = path_metric[final_offset + aug];
            if !m.is_finite() { continue; }
            let penalty = if cs != 0 { 50.0 } else { 0.0 } + drift.abs() as f64 * 0.5;
            let total = m + penalty;
            if total < best_metric {
                best_metric = total;
                best_aug = aug;
            }
        }
    }

    if best_metric.is_infinite() {
        return vec![0u8; num_info_bits];
    }

    // Traceback: walk from final step back to step 0.
    // back_ptr stores absolute cell index (step_offset + aug) for M/D,
    // and just aug (within the same step) for I transitions.
    // We reconstruct step from the stored value.
    let mut decoded = vec![0u8; num_info_bits];
    let _tail_channel_uses = memory * rate;
    let total_input_commit_points = total_info_steps;
    let mut step = total_channel_uses;
    let mut aug = best_aug;
    let mut commit_pos = total_input_commit_points - 1;

    let mut safety = total_channel_uses * 4;
    while step > 0 && safety > 0 {
        safety -= 1;
        let offset = step * num_aug + aug;
        let tt = trans_type[offset];
        let ib = input_bit[offset];
        let prev_packed = back_ptr[offset] as usize;

        // Recover predecessor step and augmented state
        let (prev_step, prev_aug) = if tt == 1 {
            // I transition: back_ptr stores just aug (same step)
            (step, prev_packed)
        } else {
            // M or D: back_ptr stores step_offset + aug
            let ps = prev_packed / num_aug;
            let pa = prev_packed % num_aug;
            (ps, pa)
        };

        if tt == 0 || tt == 2 {
            // M or D: step advances by 1
            if ib >= 0 {
                if commit_pos < num_info_bits {
                    decoded[commit_pos] = ib as u8;
                }
                if commit_pos > 0 { commit_pos -= 1; }
            }
            step -= 1;
        }
        // For I transitions, step stays the same
        aug = prev_aug;
    }

    decoded
}

// ===========================================================================
// WASM Public API
// ===========================================================================

/// K=9 Indel-Tolerant Viterbi decode (byte-oriented).
///
/// # Arguments
/// * `received_bytes` — received bytes (hard decisions)
/// * `llr_bytes` — packed LLR values as IEEE 754 f64 little-endian bytes (8 bytes per LLR, one per bit)
/// * `num_info_bits` — number of information bits
/// * `max_drift` — maximum drift to track (default 15)
/// * `insertion_penalty_x10` — insertion penalty × 10 (e.g. 15 for 1.5)
/// * `deletion_penalty_x10` — deletion penalty × 10 (e.g. 10 for 1.0)
///
/// # Returns
/// Decoded bytes
#[wasm_bindgen]
pub fn viterbi_k9_decode(
    received_bytes: &[u8],
    llr_f64: &[u8],
    num_info_bits: usize,
    max_drift: usize,
    insertion_penalty_x10: usize,
    deletion_penalty_x10: usize,
) -> Vec<u8> {
    let max_drift = max_drift as i32;
    let ins_pen = insertion_penalty_x10 as f64 / 10.0;
    let del_pen = deletion_penalty_x10 as f64 / 10.0;

    // Unpack received bytes to bits
    let total_bits = received_bytes.len() * 8;
    let mut received_bits = vec![0u8; total_bits];
    for i in 0..total_bits {
        received_bits[i] = (received_bytes[i / 8] >> (7 - (i % 8))) & 1;
    }

    // Unpack LLR f64 values
    let num_llr = llr_f64.len() / 8;
    let mut llr = vec![0.0f64; num_llr];
    for i in 0..num_llr {
        let bytes = &llr_f64[i * 8..(i + 1) * 8];
        llr[i] = f64::from_le_bytes(bytes.try_into().unwrap_or([0u8; 8]));
    }

    let decoded_bits = viterbi_decode(
        &received_bits, &llr, num_info_bits,
        K9_MEMORY, max_drift, ins_pen, del_pen,
    );

    // Pack bits to bytes
    let num_bytes = (decoded_bits.len() + 7) / 8;
    let mut decoded_bytes = vec![0u8; num_bytes];
    for i in 0..decoded_bits.len() {
        decoded_bytes[i / 8] |= decoded_bits[i] << (7 - (i % 8));
    }
    decoded_bytes
}

/// K=7 (Voyager) Indel-Tolerant Viterbi decode (byte-oriented).
#[wasm_bindgen]
pub fn viterbi_k7_decode(
    received_bytes: &[u8],
    llr_f64: &[u8],
    num_info_bits: usize,
    max_drift: usize,
    insertion_penalty_x10: usize,
    deletion_penalty_x10: usize,
) -> Vec<u8> {
    let max_drift = max_drift as i32;
    let ins_pen = insertion_penalty_x10 as f64 / 10.0;
    let del_pen = deletion_penalty_x10 as f64 / 10.0;

    let total_bits = received_bytes.len() * 8;
    let mut received_bits = vec![0u8; total_bits];
    for i in 0..total_bits {
        received_bits[i] = (received_bytes[i / 8] >> (7 - (i % 8))) & 1;
    }

    let num_llr = llr_f64.len() / 8;
    let mut llr = vec![0.0f64; num_llr];
    for i in 0..num_llr {
        let bytes = &llr_f64[i * 8..(i + 1) * 8];
        llr[i] = f64::from_le_bytes(bytes.try_into().unwrap_or([0u8; 8]));
    }

    let decoded_bits = viterbi_decode(
        &received_bits, &llr, num_info_bits,
        K7_MEMORY, max_drift, ins_pen, del_pen,
    );

    let num_bytes = (decoded_bits.len() + 7) / 8;
    let mut decoded_bytes = vec![0u8; num_bytes];
    for i in 0..decoded_bits.len() {
        decoded_bytes[i / 8] |= decoded_bits[i] << (7 - (i % 8));
    }
    decoded_bytes
}

/// Convolutional encode with K=9 NASA code.
#[wasm_bindgen]
pub fn conv_k9_encode(info_bytes: &[u8]) -> Vec<u8> {
    let tbl = &*K9_TABLE;
    let num_info_bits = info_bytes.len() * 8;

    // Unpack info bytes to bits
    let mut info_bits = vec![0u8; num_info_bits + tbl.memory]; // + zero tail
    for i in 0..num_info_bits {
        info_bits[i] = (info_bytes[i / 8] >> (7 - (i % 8))) & 1;
    }

    // Encode: for each info bit, compute output and advance state
    let total_bits = info_bits.len() * tbl.rate;
    let mut output_bits = vec![0u8; total_bits];
    let mut state = 0usize;

    for (i, &input) in info_bits.iter().enumerate() {
        let idx = state * 2 + input as usize;
        let output = tbl.outputs[idx] as usize;
        state = tbl.next_states[idx] as usize;

        // Emit rate output bits: G1 (MSB) first, then G2 (LSB)
        // output = (G1 << 1) | G2, so bit 1 = G1, bit 0 = G2
        // For rate=2: position 0 = G1 = (output >> 1) & 1, position 1 = G2 = output & 1
        for r in 0..tbl.rate {
            output_bits[i * tbl.rate + r] = ((output >> (tbl.rate - 1 - r)) & 1) as u8;
        }
    }

    // Pack to bytes
    let num_bytes = (output_bits.len() + 7) / 8;
    let mut out = vec![0u8; num_bytes];
    for i in 0..output_bits.len() {
        out[i / 8] |= output_bits[i] << (7 - (i % 8));
    }
    out
}

/// Convolutional encode with K=7 Voyager code.
#[wasm_bindgen]
pub fn conv_k7_encode(info_bytes: &[u8]) -> Vec<u8> {
    let tbl = &*K7_TABLE;
    let num_info_bits = info_bytes.len() * 8;

    let mut info_bits = vec![0u8; num_info_bits + tbl.memory];
    for i in 0..num_info_bits {
        info_bits[i] = (info_bytes[i / 8] >> (7 - (i % 8))) & 1;
    }

    let total_bits = info_bits.len() * tbl.rate;
    let mut output_bits = vec![0u8; total_bits];
    let mut state = 0usize;

    for (i, &input) in info_bits.iter().enumerate() {
        let idx = state * 2 + input as usize;
        let output = tbl.outputs[idx] as usize;
        state = tbl.next_states[idx] as usize;
        for r in 0..tbl.rate {
            output_bits[i * tbl.rate + r] = ((output >> (tbl.rate - 1 - r)) & 1) as u8;
        }
    }

    let num_bytes = (output_bits.len() + 7) / 8;
    let mut out = vec![0u8; num_bytes];
    for i in 0..output_bits.len() {
        out[i / 8] |= output_bits[i] << (7 - (i % 8));
    }
    out
}

// ===========================================================================
// Standalone test (cargo test)
// ===========================================================================

// ===========================================================================
// Standard (non-indel) Viterbi — simpler, used as reference and for no-noise case
// ===========================================================================

/// Standard Viterbi decode without indel tolerance.
/// Each info step produces `rate` channel bits. Traceback from state 0.
fn viterbi_decode_standard(
    received_bits: &[u8],
    num_info_bits: usize,
    memory: usize,
) -> Vec<u8> {
    let tbl = if memory == 8 { &*K9_TABLE } else { &*K7_TABLE };
    let num_states = tbl.num_states;
    let rate = tbl.rate;
    let total_info_steps = num_info_bits + memory;
    let total_channel_uses = total_info_steps * rate;

    // PM[step * num_states + state], step 0..total_info_steps
    let mut pm = vec![f64::INFINITY; (total_info_steps + 1) * num_states];
    let mut prev_state = vec![0u16; (total_info_steps + 1) * num_states];
    let mut inp_bit = vec![0u8; (total_info_steps + 1) * num_states];

    pm[0] = 0.0; // start at state 0

    for info_step in 0..total_info_steps {
        let ch_base = info_step * rate;

        for s in 0..num_states {
            let m = pm[info_step * num_states + s];
            if !m.is_finite() { continue; }

            for input in 0..2usize {
                let idx = s * 2 + input;
                let output = tbl.outputs[idx] as usize;
                let ns = tbl.next_states[idx] as usize;

                // Compute distance: sum over rate bits
                let mut dist = 0.0f64;
                for r in 0..rate {
                    let emitted = (output >> (rate - 1 - r)) & 1;
                    let ch_pos = ch_base + r;
                    if ch_pos < received_bits.len() && emitted != received_bits[ch_pos] as usize {
                        dist += 1.0;
                    }
                }

                let new_m = m + dist;
                let target = (info_step + 1) * num_states + ns;
                if new_m < pm[target] {
                    pm[target] = new_m;
                    prev_state[target] = s as u16;
                    inp_bit[target] = input as u8;
                }
            }
        }
    }

    // Traceback from state 0 (zero-tail forcing)
    let mut decoded = vec![0u8; num_info_bits];
    let mut s = 0usize;
    for step in (1..=total_info_steps).rev() {
        let idx = step * num_states + s;
        if step - 1 < num_info_bits {
            decoded[step - 1] = inp_bit[idx];
        }
        s = prev_state[idx] as usize;
    }

    decoded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_k9_encode_decode_standard() {
        // First test with 2 bytes (known good from standalone)
        let info_small = vec![0xABu8, 0xCD];
        let encoded_small = conv_k9_encode(&info_small);
        let mut bits_small = vec![0u8; encoded_small.len() * 8];
        for i in 0..bits_small.len() { bits_small[i] = (encoded_small[i / 8] >> (7 - (i % 8))) & 1; }
        let decoded_small = viterbi_decode_standard(&bits_small, info_small.len() * 8, K9_MEMORY);
        let mut dec_bytes_small = vec![0u8; info_small.len()];
        for i in 0..decoded_small.len().min(info_small.len() * 8) {
            dec_bytes_small[i / 8] |= decoded_small[i] << (7 - (i % 8));
        }
        eprintln!("2-byte test: decoded={:?}, expected={:?}", dec_bytes_small, info_small);
        assert_eq!(dec_bytes_small, info_small, "K=9 2-byte roundtrip");

        // Now test with 8 bytes
        let info = vec![0xABu8, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89];
        let encoded = conv_k9_encode(&info);
        eprintln!("8-byte encoded {} bytes", encoded.len());

        // Unpack encoded bytes to bits
        let mut bits = vec![0u8; encoded.len() * 8];
        for i in 0..bits.len() { bits[i] = (encoded[i / 8] >> (7 - (i % 8))) & 1; }
        eprintln!("bits.len()={}, expected={}", bits.len(), (info.len() * 8 + K9_MEMORY) * 2);

        // Standard Viterbi decode (no indel)
        let num_info_bits = info.len() * 8;
        let decoded_bits = viterbi_decode_standard(&bits, num_info_bits, K9_MEMORY);
        eprintln!("decoded {} bits", decoded_bits.len());

        let mut decoded_bytes = vec![0u8; info.len()];
        for i in 0..decoded_bits.len().min(info.len() * 8) {
            decoded_bytes[i / 8] |= decoded_bits[i] << (7 - (i % 8));
        }

        assert_eq!(decoded_bytes, info, "K=9 standard Viterbi roundtrip should be exact");
    }

    #[test]
    fn test_k7_encode_decode_standard() {
        let info = vec![0xABu8, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89];
        let encoded = conv_k7_encode(&info);

        let mut bits = vec![0u8; encoded.len() * 8];
        for i in 0..bits.len() { bits[i] = (encoded[i / 8] >> (7 - (i % 8))) & 1; }

        let num_info_bits = info.len() * 8;
        let decoded_bits = viterbi_decode_standard(&bits, num_info_bits, K7_MEMORY);

        let mut decoded_bytes = vec![0u8; info.len()];
        for i in 0..decoded_bits.len().min(info.len() * 8) {
            decoded_bytes[i / 8] |= decoded_bits[i] << (7 - (i % 8));
        }

        assert_eq!(decoded_bytes, info, "K=7 standard Viterbi roundtrip should be exact");
    }

    #[test]
    fn test_k9_indel_no_noise() {
        let info = vec![0xABu8, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89];
        let encoded = conv_k9_encode(&info);

        let mut bits = vec![0u8; encoded.len() * 8];
        for i in 0..bits.len() { bits[i] = (encoded[i / 8] >> (7 - (i % 8))) & 1; }

        let num_info_bits = info.len() * 8;
        let decoded_bits = viterbi_decode(&bits, &[], num_info_bits, K9_MEMORY, 15, 1.5, 1.0);

        let mut decoded_bytes = vec![0u8; info.len()];
        for i in 0..decoded_bits.len().min(info.len() * 8) {
            decoded_bytes[i / 8] |= decoded_bits[i] << (7 - (i % 8));
        }

        assert_eq!(decoded_bytes, info, "K=9 indel Viterbi roundtrip with no noise should be exact");
    }

    #[test]
    fn test_k9_indel_with_noise() {
        let info = vec![0x42u8, 0x37, 0xA5, 0xC3];
        let encoded = conv_k9_encode(&info);

        let mut bits = vec![0u8; encoded.len() * 8];
        for i in 0..bits.len() { bits[i] = (encoded[i / 8] >> (7 - (i % 8))) & 1; }

        // Add 3% substitution noise
        let mut noisy = bits.clone();
        let mut rng = 42u32;
        for i in 0..noisy.len() {
            rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
            if (rng % 100) < 3 { noisy[i] ^= 1; }
        }

        let num_info_bits = info.len() * 8;
        let decoded_bits = viterbi_decode(&noisy, &[], num_info_bits, K9_MEMORY, 5, 1.5, 1.0);

        let mut decoded_bytes = vec![0u8; info.len()];
        for i in 0..decoded_bits.len().min(info.len() * 8) {
            decoded_bytes[i / 8] |= decoded_bits[i] << (7 - (i % 8));
        }

        assert_eq!(decoded_bytes, info, "K=9 indel Viterbi should correct 3% sub noise on 4 bytes");
    }
}
