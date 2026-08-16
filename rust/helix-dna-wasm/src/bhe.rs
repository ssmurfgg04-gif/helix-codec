//! Bounded Homopolymer Encoding (BHE) — u128 bit-parallel FSM
//!
//! Replaces JS BigInt FSM with native u128 arithmetic:
//!   - k=1: base-2 → base-3 conversion using u128 (50× faster than BigInt)
//!   - k>1: FSM with precomputed transition table and u128 path counts
//!
//! The FSM transition table is built once and cached. Encoding/decoding
//! walks the FSM using native integer arithmetic instead of BigInt.
//!
//! For inputs > 16 bytes, chunked encoding is used:
//!   - k=1: each 15-byte chunk is encoded independently with a "previous base"
//!     constraint at chunk boundaries to prevent junction homopolymers.
//!   - k>1 (FSM): each 15-byte chunk is encoded independently; a small header
//!     stores each chunk's DNA length so the decoder can split correctly.

use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Codebook tables
// ---------------------------------------------------------------------------

/// For each previous base (0..3), the 3 non-previous bases in order.
const CODEBOOK: [[usize; 3]; 4] = [
    [1, 2, 3], // prev=A(0) → C(1), G(2), T(3)
    [0, 2, 3], // prev=C(1) → A(0), G(2), T(3)
    [0, 1, 3], // prev=G(2) → A(0), C(1), T(3)
    [0, 1, 2], // prev=T(3) → A(0), C(1), G(2)
];

/// Inverse codebook: inv_codebook[prev][base] = trit (0..2), or 255 if invalid.
const INV_CODEBOOK: [[u8; 4]; 4] = {
    let mut table = [[255u8; 4]; 4];
    let mut prev = 0;
    while prev < 4 {
        let mut t = 0;
        while t < 3 {
            table[prev][CODEBOOK[prev][t]] = t as u8;
            t += 1;
        }
        prev += 1;
    }
    table
};

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const BASE_TO_IDX: [u8; 256] = {
    let mut table = [255u8; 256];
    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table
};

// ---------------------------------------------------------------------------
// u128 helpers
// ---------------------------------------------------------------------------

/// Convert bytes to u128 (big-endian). Works for up to 16 bytes.
fn bytes_to_u128(bytes: &[u8]) -> u128 {
    let mut v = 0u128;
    for &b in bytes {
        v = (v << 8) | b as u128;
    }
    v
}

/// Convert u128 to exactly num_bytes bytes (big-endian).
fn u128_to_bytes(value: u128, num_bytes: usize) -> Vec<u8> {
    let mut out = vec![0u8; num_bytes];
    let mut v = value;
    for i in (0..num_bytes).rev() {
        out[i] = (v & 0xFF) as u8;
        v >>= 8;
    }
    out
}

// ---------------------------------------------------------------------------
// k=1 Fast Path: base-2 → base-3 using u128
// ---------------------------------------------------------------------------

/// Compute minimum DNA length for k=1 encoding without previous base constraint.
/// Value space: 4 * 3^(n-1). We need 4 * 3^(n-1) > 2^total_bits.
fn k1_dna_length(num_bytes: usize) -> usize {
    if num_bytes == 0 { return 0; }
    let total_bits = num_bytes * 8;
    let log2_3: f64 = 1.584_962_500_721_156;
    std::cmp::max(1, ((total_bits as f64 - 2.0) / log2_3).ceil() as usize + 1)
}

/// Compute minimum DNA length for k=1 encoding WITH previous base constraint.
/// Value space: 3^n. We need 3^n > 2^total_bits.
fn k1_dna_length_with_prev(num_bytes: usize) -> usize {
    if num_bytes == 0 { return 0; }
    let total_bits = num_bytes * 8;
    let log2_3: f64 = 1.584_962_500_721_156;
    std::cmp::max(1, (total_bits as f64 / log2_3).ceil() as usize + 1)
}

/// k=1 encode for a single chunk (≤16 bytes), no previous base constraint.
/// First base has 4 choices, subsequent have 3 each.
fn bhe_encode_k1_direct(bytes: &[u8]) -> Vec<u8> {
    debug_assert!(bytes.len() <= 16);
    let dna_length = k1_dna_length(bytes.len());
    let value = bytes_to_u128(bytes);

    let mut bases = vec![0u8; dna_length];
    let mut v = value;
    for i in (1..dna_length).rev() {
        bases[i] = (v % 3) as u8;
        v /= 3;
    }
    bases[0] = (v % 4) as u8;

    let mut out = vec![0u8; dna_length];
    let mut prev_idx = bases[0] as usize;
    out[0] = BASES[prev_idx];
    for i in 1..dna_length {
        let trit = bases[i] as usize;
        let next_idx = CODEBOOK[prev_idx][trit];
        out[i] = BASES[next_idx];
        prev_idx = next_idx;
    }
    out
}

/// k=1 encode for a single chunk (≤16 bytes), WITH previous base constraint.
/// First base has 3 choices (not prev_base), subsequent have 3 each.
fn bhe_encode_k1_with_prev(bytes: &[u8], prev_base: u8) -> Vec<u8> {
    debug_assert!(bytes.len() <= 16);
    let dna_length = k1_dna_length_with_prev(bytes.len());
    let value = bytes_to_u128(bytes);
    let prev_idx = BASE_TO_IDX[prev_base as usize] as usize;

    // Decompose value into all base-3 trits
    let mut trits = vec![0u8; dna_length];
    let mut v = value;
    for i in (0..dna_length).rev() {
        trits[i] = (v % 3) as u8;
        v /= 3;
    }

    // Map to DNA: first base uses codebook relative to prev_base
    let mut out = vec![0u8; dna_length];
    let first_idx = CODEBOOK[prev_idx][trits[0] as usize];
    out[0] = BASES[first_idx];
    let mut cur_idx = first_idx;
    for i in 1..dna_length {
        let next_idx = CODEBOOK[cur_idx][trits[i] as usize];
        out[i] = BASES[next_idx];
        cur_idx = next_idx;
    }
    out
}

/// k=1 decode for a single chunk, no previous base constraint.
fn bhe_decode_k1_direct(dna: &[u8], num_bytes: usize) -> Vec<u8> {
    let n = dna.len();
    let first_idx = BASE_TO_IDX[dna[0] as usize] as usize;

    let mut value = first_idx as u128;
    let mut prev_idx = first_idx;
    for i in 1..n {
        let base_idx = BASE_TO_IDX[dna[i] as usize] as usize;
        let trit = INV_CODEBOOK[prev_idx][base_idx];
        if trit == 255 { return vec![]; }
        value = value * 3 + trit as u128;
        prev_idx = base_idx;
    }
    u128_to_bytes(value, num_bytes)
}

/// k=1 decode for a single chunk, WITH previous base constraint.
fn bhe_decode_k1_with_prev(dna: &[u8], num_bytes: usize, prev_base: u8) -> Vec<u8> {
    let n = dna.len();
    let prev_idx = BASE_TO_IDX[prev_base as usize] as usize;
    let first_idx = BASE_TO_IDX[dna[0] as usize] as usize;
    let first_trit = INV_CODEBOOK[prev_idx][first_idx];
    if first_trit == 255 { return vec![]; }

    let mut value = first_trit as u128;
    let mut cur_idx = first_idx;
    for i in 1..n {
        let base_idx = BASE_TO_IDX[dna[i] as usize] as usize;
        let trit = INV_CODEBOOK[cur_idx][base_idx];
        if trit == 255 { return vec![]; }
        value = value * 3 + trit as u128;
        cur_idx = base_idx;
    }
    u128_to_bytes(value, num_bytes)
}

// ---------------------------------------------------------------------------
// k=1 public API with chunking
// ---------------------------------------------------------------------------

const K1_CHUNK_SIZE: usize = 15; // 15 bytes = 120 bits, safe for u128

/// BHE k=1 encode: convert bytes to DNA with no homopolymers.
///
/// Uses u128 arithmetic for ~50× speedup over JS BigInt.
/// For inputs ≤ 16 bytes: direct single-value encoding.
/// For inputs > 16 bytes: chunked encoding with junction homopolymer avoidance.
#[wasm_bindgen]
pub fn bhe_encode_k1(bytes: &[u8]) -> Vec<u8> {
    if bytes.is_empty() {
        return vec![];
    }

    if bytes.len() <= 16 {
        return bhe_encode_k1_direct(bytes);
    }

    // Chunked encoding: split into 15-byte chunks
    let mut out = Vec::new();
    let mut offset = 0;
    let mut prev_base: Option<u8> = None;

    while offset < bytes.len() {
        let end = std::cmp::min(offset + K1_CHUNK_SIZE, bytes.len());
        let chunk = &bytes[offset..end];

        let dna = match prev_base {
            Some(pb) => bhe_encode_k1_with_prev(chunk, pb),
            None => bhe_encode_k1_direct(chunk),
        };

        out.extend_from_slice(&dna);
        if !dna.is_empty() {
            prev_base = Some(dna[dna.len() - 1]);
        }
        offset = end;
    }
    out
}

/// BHE k=1 decode: convert DNA back to bytes.
#[wasm_bindgen]
pub fn bhe_decode_k1(dna: &[u8], num_bytes: usize) -> Vec<u8> {
    if dna.is_empty() || num_bytes == 0 {
        return vec![];
    }

    if num_bytes <= 16 {
        return bhe_decode_k1_direct(dna, num_bytes);
    }

    // Chunked decoding: split DNA at known chunk boundaries
    let mut out = Vec::with_capacity(num_bytes);
    let mut dna_offset = 0;
    let mut byte_offset = 0;
    let mut prev_base: Option<u8> = None;

    while byte_offset < num_bytes {
        let chunk_bytes = std::cmp::min(K1_CHUNK_SIZE, num_bytes - byte_offset);
        let chunk_dna_len = match prev_base {
            Some(_) => k1_dna_length_with_prev(chunk_bytes),
            None => k1_dna_length(chunk_bytes),
        };

        if dna_offset + chunk_dna_len > dna.len() {
            // Not enough DNA data
            return vec![];
        }
        let chunk_dna = &dna[dna_offset..dna_offset + chunk_dna_len];

        let decoded = match prev_base {
            Some(pb) => bhe_decode_k1_with_prev(chunk_dna, chunk_bytes, pb),
            None => bhe_decode_k1_direct(chunk_dna, chunk_bytes),
        };

        out.extend_from_slice(&decoded);
        if !chunk_dna.is_empty() {
            prev_base = Some(chunk_dna[chunk_dna.len() - 1]);
        }
        dna_offset += chunk_dna_len;
        byte_offset += chunk_bytes;
    }
    out
}

// ---------------------------------------------------------------------------
// k>1 FSM-based encode/decode using u128 path counts
// ---------------------------------------------------------------------------

/// FSM state count cache.
struct BheFsm {
    max_run: usize,
    num_states: usize,
    /// transitions[state * 4 + base] -> next_state (-1 = invalid)
    transitions: Vec<i16>,
}

impl BheFsm {
    fn new(max_run: usize) -> Self {
        let num_states = 1 + 4 * max_run;
        let mut transitions = vec![-1i16; num_states * 4];

        // From initial state (0): all bases valid
        for b in 0..4 {
            transitions[0 * 4 + b] = (1 + b * max_run) as i16;
        }

        // From state (b, r): encoded as 1 + b * maxRun + (r - 1)
        for b in 0..4 {
            for r in 1..=max_run {
                let state = 1 + b * max_run + (r - 1);
                for b2 in 0..4 {
                    let next_state = if b2 == b && r < max_run {
                        (1 + b * max_run + r) as i16
                    } else if b2 != b {
                        (1 + b2 * max_run) as i16
                    } else {
                        -1 // would exceed maxRun
                    };
                    transitions[state * 4 + b2] = next_state;
                }
            }
        }

        BheFsm { max_run, num_states, transitions }
    }

    /// Count valid strings of length `len` from each state.
    fn count_paths(&self, len: usize) -> Vec<u128> {
        let mut prev = vec![1u128; self.num_states];
        for _ in 1..=len {
            let mut curr = vec![0u128; self.num_states];
            for s in 0..self.num_states {
                let mut total = 0u128;
                for b in 0..4 {
                    let ns = self.transitions[s * 4 + b] as usize;
                    if ns < self.num_states {
                        total += prev[ns];
                    }
                }
                curr[s] = total;
            }
            prev = curr;
        }
        prev
    }

    /// Total valid strings of length `len` from initial state.
    fn total_count(&self, len: usize) -> u128 {
        if len == 0 { return 1; }
        let counts = self.count_paths(len - 1);
        let mut total = 0u128;
        for b in 0..4 {
            let ns = self.transitions[0 * 4 + b] as usize;
            if ns < self.num_states {
                total += counts[ns];
            }
        }
        total
    }

    /// Encode a value (fits in u128) as DNA of the given length.
    fn encode_value(&self, value: u128, dna_length: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(dna_length);
        let mut remaining = value;
        let mut state = 0usize;

        for pos in 0..dna_length {
            let remaining_len = dna_length - pos - 1;
            let counts = self.count_paths(remaining_len);

            let mut chosen = false;
            for b in 0..4 {
                let ns = self.transitions[state * 4 + b] as usize;
                if ns >= self.num_states { continue; }
                let count = counts[ns];
                if remaining < count {
                    out.push(BASES[b]);
                    state = ns;
                    chosen = true;
                    break;
                }
                remaining -= count;
            }
            if !chosen { break; }
        }
        out
    }

    /// Decode DNA of the given length back to a u128 value.
    fn decode_value(&self, dna: &[u8]) -> u128 {
        let dna_length = dna.len();
        let mut value = 0u128;
        let mut state = 0usize;

        for pos in 0..dna_length {
            let base_idx = BASE_TO_IDX[dna[pos] as usize] as usize;
            let remaining_len = dna_length - pos - 1;
            let counts = self.count_paths(remaining_len);

            for b2 in 0..base_idx {
                let ns = self.transitions[state * 4 + b2] as usize;
                if ns < self.num_states {
                    value += counts[ns];
                }
            }

            let ns = self.transitions[state * 4 + base_idx] as usize;
            if ns >= self.num_states { return 0; }
            state = ns;
        }
        value
    }
}

/// FSM encode for a single chunk (≤16 bytes).
/// Returns the encoded DNA.
fn bhe_encode_fsm_direct(bytes: &[u8], max_run: usize) -> Vec<u8> {
    debug_assert!(bytes.len() <= 16);
    let fsm = BheFsm::new(max_run);
    let value = bytes_to_u128(bytes);

    let log2_avg: f64 = 1.584_962_500_721_156;
    let estimated_dna_len = ((bytes.len() * 8) as f64 / log2_avg).ceil() as usize + 8;
    let max_search = estimated_dna_len * 2 + 64;

    let mut dna_length = 1;
    while dna_length < max_search {
        if fsm.total_count(dna_length) > value {
            break;
        }
        dna_length += 1;
    }

    fsm.encode_value(value, dna_length)
}

/// FSM decode for a single chunk.
fn bhe_decode_fsm_direct(dna: &[u8], num_bytes: usize, max_run: usize) -> Vec<u8> {
    let fsm = BheFsm::new(max_run);
    let value = fsm.decode_value(dna);
    u128_to_bytes(value, num_bytes)
}

const FSM_CHUNK_SIZE: usize = 15;

/// BHE k>1 encode using FSM with u128 arithmetic.
///
/// For inputs ≤ 16 bytes: direct single-value encoding.
/// For inputs > 16 bytes: chunked encoding with a header storing chunk DNA lengths.
#[wasm_bindgen]
pub fn bhe_encode_fsm(bytes: &[u8], max_run: usize) -> Vec<u8> {
    if bytes.is_empty() {
        return vec![];
    }
    if max_run <= 1 {
        return bhe_encode_k1(bytes);
    }

    if bytes.len() <= 16 {
        return bhe_encode_fsm_direct(bytes, max_run);
    }

    // Chunked encoding
    let mut chunks_dna: Vec<Vec<u8>> = Vec::new();
    let mut offset = 0;
    while offset < bytes.len() {
        let end = std::cmp::min(offset + FSM_CHUNK_SIZE, bytes.len());
        let chunk = &bytes[offset..end];
        let dna = bhe_encode_fsm_direct(chunk, max_run);
        chunks_dna.push(dna);
        offset = end;
    }

    // Build output with header:
    //   [1 byte: num_chunks] [2 bytes per chunk: DNA length as u16 BE] [concatenated DNA]
    let num_chunks = chunks_dna.len();
    let header_len = 1 + 2 * num_chunks;
    let total_dna_len: usize = chunks_dna.iter().map(|c| c.len()).sum();

    let mut out = Vec::with_capacity(header_len + total_dna_len);
    out.push(num_chunks as u8);
    for chunk in &chunks_dna {
        let len = chunk.len() as u16;
        out.push((len >> 8) as u8);
        out.push((len & 0xFF) as u8);
    }
    for chunk in &chunks_dna {
        out.extend_from_slice(chunk);
    }
    out
}

/// BHE k>1 decode using FSM with u128 arithmetic.
#[wasm_bindgen]
pub fn bhe_decode_fsm(dna: &[u8], num_bytes: usize, max_run: usize) -> Vec<u8> {
    if dna.is_empty() || num_bytes == 0 {
        return vec![];
    }
    if max_run <= 1 {
        return bhe_decode_k1(dna, num_bytes);
    }

    // Check if this is chunked output (has a header)
    // Chunked header: [1 byte num_chunks] [2 bytes per chunk DNA length]
    // A non-chunked output is pure DNA (all bytes in A,C,G,T range)
    // We detect chunked output by checking if the first byte is NOT a valid DNA base
    let first_byte = dna[0];
    let is_chunked = first_byte != b'A' && first_byte != b'C' && first_byte != b'G' && first_byte != b'T';

    if !is_chunked {
        // Direct (non-chunked) decode
        return bhe_decode_fsm_direct(dna, num_bytes, max_run);
    }

    // Chunked decode: read header
    if dna.len() < 1 {
        return vec![];
    }
    let num_chunks = dna[0] as usize;
    let header_len = 1 + 2 * num_chunks;
    if dna.len() < header_len {
        return vec![];
    }

    // Read chunk DNA lengths
    let mut chunk_dna_lengths: Vec<usize> = Vec::with_capacity(num_chunks);
    for i in 0..num_chunks {
        let hi = dna[1 + 2 * i] as usize;
        let lo = dna[1 + 2 * i + 1] as usize;
        chunk_dna_lengths.push((hi << 8) | lo);
    }

    // Compute chunk byte counts
    let mut chunk_byte_counts: Vec<usize> = Vec::with_capacity(num_chunks);
    let mut remaining = num_bytes;
    for i in 0..num_chunks {
        let cb = if i < num_chunks - 1 {
            std::cmp::min(FSM_CHUNK_SIZE, remaining)
        } else {
            remaining
        };
        chunk_byte_counts.push(cb);
        remaining -= cb;
    }

    // Decode each chunk
    let mut out = Vec::with_capacity(num_bytes);
    let mut dna_offset = header_len;
    for i in 0..num_chunks {
        let chunk_dna_len = chunk_dna_lengths[i];
        if dna_offset + chunk_dna_len > dna.len() {
            return vec![];
        }
        let chunk_dna = &dna[dna_offset..dna_offset + chunk_dna_len];
        let decoded = bhe_decode_fsm_direct(chunk_dna, chunk_byte_counts[i], max_run);
        out.extend_from_slice(&decoded);
        dna_offset += chunk_dna_len;
    }
    out
}

/// Validate BHE-encoded DNA satisfies homopolymer constraint.
#[wasm_bindgen]
pub fn bhe_validate(dna: &[u8], max_run: usize) -> bool {
    if dna.is_empty() {
        return true;
    }
    let mut run = 1usize;
    for i in 1..dna.len() {
        if dna[i] == dna[i - 1] {
            run += 1;
            if run > max_run {
                return false;
            }
        } else {
            run = 1;
        }
    }
    true
}
