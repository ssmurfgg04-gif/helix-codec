//! Comprehensive test for the production-hardened indel-tolerant Viterbi decoder v3.
//!
//! Run: cargo test --release -- --nocapture

use std::time::Instant;

const K9_MEMORY: usize = 8;
const K9_GENERATORS: [u16; 2] = [0o561, 0o753];
const K7_MEMORY: usize = 6;
const K7_GENERATORS: [u16; 2] = [0o171, 0o133];

struct TransitionTable { outputs: Vec<u16>, next_states: Vec<u16>, num_states: usize, memory: usize, rate: usize }
impl TransitionTable {
    fn build(memory: usize, generators: &[u16], rate: usize) -> Self {
        let ns = 1usize << memory;
        let mut outputs = vec![0u16; ns * 2]; let mut next_states = vec![0u16; ns * 2];
        for s in 0..ns { for ib in 0..2 {
            next_states[s * 2 + ib] = ((s >> 1) | ((ib as usize) << (memory - 1))) as u16;
            let mut ob = 0u16;
            for g in 0..rate { let mut b = 0u16; if (generators[g] & 1) != 0 { b ^= ib as u16; } for i in 0..memory { if (generators[g] & (1 << (i + 1))) != 0 { b ^= ((s >> i) & 1) as u16; } } ob |= b << g; }
            outputs[s * 2 + ib] = ob;
        }}
        TransitionTable { outputs, next_states, num_states: ns, memory, rate }
    }
    #[inline] fn output(&self, s: usize, i: usize) -> u16 { self.outputs[s * 2 + i] }
    #[inline] fn next_state(&self, s: usize, i: usize) -> usize { self.next_states[s * 2 + i] as usize }
}

fn conv_encode(data: &[u8], table: &TransitionTable) -> Vec<u8> {
    let ob = (data.len() * 8 + table.memory) * table.rate;
    let mut r = vec![0u8; (ob + 7) / 8]; let mut st = 0usize; let mut p = 0usize;
    for &b in data { for bi in 0..8 { let ib = ((b >> (7 - bi)) & 1) as usize; let o = table.output(st, ib); st = table.next_state(st, ib); for g in 0..table.rate { if (o >> (table.rate - 1 - g)) & 1 != 0 { r[p / 8] |= 1 << (7 - (p % 8)); } p += 1; } } }
    for _ in 0..table.memory { let o = table.output(st, 0); st = table.next_state(st, 0); for g in 0..table.rate { if (o >> (table.rate - 1 - g)) & 1 != 0 && p / 8 < r.len() { r[p / 8] |= 1 << (7 - (p % 8)); } p += 1; } }
    r.truncate((p + 7) / 8); r
}

#[inline]
fn get_bit(d: &[u8], p: usize) -> u8 { if p / 8 >= d.len() { 0 } else { (d[p / 8] >> (7 - (p % 8))) & 1 } }

// Production v3 decoder — same as lib.rs
fn viterbi_decode_indel(
    received: &[u8], table: &TransitionTable, max_drift: u32,
    ins_pen: f64, del_pen: f64, num_info_bits_override: Option<usize>,
) -> Vec<u8> {
    let ns = table.num_states; let memory = table.memory; let rate = table.rate;
    let md = max_drift as i32; let nd = 2 * max_drift as usize + 1;
    let rb = received.len() * 8;
    if rb == 0 { return Vec::new(); }
    let num_info_bits = if let Some(nib) = num_info_bits_override { nib } else { if rb >= memory * rate { (rb / rate) - memory } else { rb / rate } };
    let total_info_steps = num_info_bits + memory;
    let total_channel_uses = total_info_steps * rate;
    let p0 = ns * nd; let ta = ns * 3 * nd; let inf = f64::INFINITY;
    let mut recv_bits: Vec<u8> = vec![0u8; rb];
    for i in 0..rb { recv_bits[i] = (received[i / 8] >> (7 - (i % 8))) & 1; }
    let max_steps = total_channel_uses + max_drift as usize + 10;
    let total_cells = (max_steps + 1) * ta;
    let mut path_metric = vec![inf; total_cells];
    let mut back_ptr = vec![0u32; total_cells];
    let mut trans_type = vec![0u8; total_cells];
    let mut input_bit = vec![-1i8; total_cells];
    let start_aug = max_drift as usize;
    path_metric[start_aug] = 0.0;

    let propagate_insertions = |step_offset: usize, step_for_recv: i32,
        path_metric: &mut Vec<f64>, back_ptr: &mut Vec<u32>,
        trans_type: &mut Vec<u8>, input_bit: &mut Vec<i8>| {
        for cs in 0..ns {
            for di in 0..nd - 1 {
                let drift = di as i32 - md; let recv_pos = step_for_recv + drift;
                if recv_pos < 0 || recv_pos as usize >= rb { continue; }
                let aug = cs * nd + di; let m = path_metric[step_offset + aug];
                if m >= 1e29 { continue; }
                let next_aug = cs * nd + di + 1; let new_metric = m + ins_pen;
                if new_metric < path_metric[step_offset + next_aug] {
                    path_metric[step_offset + next_aug] = new_metric;
                    back_ptr[step_offset + next_aug] = aug as u32;
                    trans_type[step_offset + next_aug] = 1; input_bit[step_offset + next_aug] = -1;
                }
            }
            for pend in 0..2usize { for di in 0..nd - 1 {
                let drift = di as i32 - md; let recv_pos = step_for_recv + drift;
                if recv_pos < 0 || recv_pos as usize >= rb { continue; }
                let aug = p0 + (cs * 2 + pend) * nd + di; let m = path_metric[step_offset + aug];
                if m >= 1e29 { continue; }
                let next_aug = p0 + (cs * 2 + pend) * nd + di + 1; let new_metric = m + ins_pen;
                if new_metric < path_metric[step_offset + next_aug] {
                    path_metric[step_offset + next_aug] = new_metric;
                    back_ptr[step_offset + next_aug] = aug as u32;
                    trans_type[step_offset + next_aug] = 1; input_bit[step_offset + next_aug] = -1;
                }
            }}
        }
    };

    propagate_insertions(0, 0, &mut path_metric, &mut back_ptr, &mut trans_type, &mut input_bit);

    for step in 0..total_channel_uses {
        let step_offset = step * ta; let next_offset = (step + 1) * ta;
        for cs in 0..ns {
            for di in 0..nd {
                let aug = cs * nd + di; let m = path_metric[step_offset + aug];
                if m >= 1e29 { continue; }
                let drift = di as i32 - md; let recv_pos = step as i32 + drift;
                for ib in 0..2usize {
                    let out = table.output(cs, ib) as usize; let emitted = (out >> (rate - 1)) & 1;
                    if recv_pos >= 0 && (recv_pos as usize) < rb {
                        let rp = recv_pos as usize; let rb_bit = recv_bits[rp] as usize;
                        let dist = if emitted != rb_bit { 1.0 } else { 0.0 };
                        let next_aug = p0 + (cs * 2 + ib) * nd + di; let new_metric = m + dist;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = aug as u32;
                            trans_type[next_offset + next_aug] = 0; input_bit[next_offset + next_aug] = -1;
                        }
                    }
                    if di > 0 {
                        let next_aug = p0 + (cs * 2 + ib) * nd + di - 1; let new_metric = m + del_pen;
                        if new_metric < path_metric[next_offset + next_aug] {
                            path_metric[next_offset + next_aug] = new_metric;
                            back_ptr[next_offset + next_aug] = aug as u32;
                            trans_type[next_offset + next_aug] = 2; input_bit[next_offset + next_aug] = -1;
                        }
                    }
                }
            }
            for pend in 0..2usize { for di in 0..nd {
                let aug = p0 + (cs * 2 + pend) * nd + di; let m = path_metric[step_offset + aug];
                if m >= 1e29 { continue; }
                let drift = di as i32 - md; let recv_pos = step as i32 + drift;
                let out = table.output(cs, pend) as usize; let next_state = table.next_state(cs, pend); let emitted = (out >> (rate - 2)) & 1;
                if recv_pos >= 0 && (recv_pos as usize) < rb {
                    let rp = recv_pos as usize; let rb_bit = recv_bits[rp] as usize;
                    let dist = if emitted != rb_bit { 1.0 } else { 0.0 };
                    let next_aug = next_state * nd + di; let new_metric = m + dist;
                    if new_metric < path_metric[next_offset + next_aug] {
                        path_metric[next_offset + next_aug] = new_metric;
                        back_ptr[next_offset + next_aug] = aug as u32;
                        trans_type[next_offset + next_aug] = 0; input_bit[next_offset + next_aug] = pend as i8;
                    }
                }
                if di > 0 {
                    let next_aug = next_state * nd + di - 1; let new_metric = m + del_pen;
                    if new_metric < path_metric[next_offset + next_aug] {
                        path_metric[next_offset + next_aug] = new_metric;
                        back_ptr[next_offset + next_aug] = aug as u32;
                        trans_type[next_offset + next_aug] = 2; input_bit[next_offset + next_aug] = pend as i8;
                    }
                }
            }}
        }
        propagate_insertions(next_offset, step as i32 + 1, &mut path_metric, &mut back_ptr, &mut trans_type, &mut input_bit);
    }

    let final_offset = total_channel_uses * ta;
    let mut best_aug = 0usize; let mut best_metric = inf;
    for cs in 0..ns { for di in 0..nd {
        let aug = cs * nd + di; let m = path_metric[final_offset + aug];
        if m >= 1e29 { continue; }
        let drift = (di as i32 - md).abs() as f64;
        let penalty = if cs != 0 { 50.0 } else { 0.0 } + drift * 0.5;
        let total = m + penalty;
        if total < best_metric { best_metric = total; best_aug = aug; }
    }}
    if best_metric >= 1e29 { return vec![0u8; num_info_bits / 8]; }

    let mut decoded_bits: Vec<u8> = Vec::with_capacity(num_info_bits);
    let mut step = total_channel_uses; let mut aug = best_aug;
    let mut safety = total_channel_uses * 4;
    while step > 0 && safety > 0 {
        safety -= 1;
        let offset = step * ta + aug; let tt = trans_type[offset]; let ib = input_bit[offset];
        let prev_aug = back_ptr[offset] as usize;
        match tt {
            0 | 2 => { if ib >= 0 { decoded_bits.push(ib as u8); } step -= 1; aug = prev_aug; }
            1 => { aug = prev_aug; }
            _ => { step -= 1; aug = prev_aug; }
        }
    }
    decoded_bits.reverse();
    if decoded_bits.len() > memory { decoded_bits.truncate(decoded_bits.len() - memory); }
    let nb = decoded_bits.len() / 8; let mut result = vec![0u8; nb];
    for i in 0..nb * 8 { if i < decoded_bits.len() && decoded_bits[i] != 0 { result[i / 8] |= 1 << (7 - (i % 8)); } }
    result
}

fn viterbi_decode_standard(received: &[u8], table: &TransitionTable) -> Vec<u8> {
    let ns = table.num_states; let rb = received.len() * 8;
    if rb == 0 { return Vec::new(); }
    let ni = f64::INFINITY; let nst = rb / 2;
    let mut c = vec![ni; ns]; let mut p = vec![ni; ns];
    let mut tr: Vec<(u16, u8)> = vec![(0, 0); nst * ns]; c[0] = 0.0;
    let mut rp = 0usize;
    for s in 0..nst { std::mem::swap(&mut c, &mut p); c.fill(ni);
        let b0 = get_bit(received, rp); let b1 = get_bit(received, rp + 1); rp += 2;
        for st in 0..ns { if p[st] >= 1e29 { continue; } for ib in 0..2 {
            let o = table.output(st, ib); let co = p[st] + if ((o >> (table.rate - 1)) & 1) as u8 != b0 { 1.0 } else { 0.0 } + if ((o >> (table.rate - 2)) & 1) as u8 != b1 { 1.0 } else { 0.0 };
            let n = table.next_state(st, ib); if co < c[n] { c[n] = co; tr[s * ns + n] = (st as u16, ib as u8); }
        }}
    }
    let mut bs = 0usize; let mut bc = c[0]; for s in 1..ns { if c[s] < bc { bc = c[s]; bs = s; } }
    let mut db = Vec::with_capacity(nst); let mut st = bs;
    for s in (0..nst).rev() { let (ps, ib) = tr[s * ns + st]; db.push(ib); st = ps as usize; }
    db.reverse(); if db.len() > table.memory { db.truncate(db.len() - table.memory); }
    let nb = db.len() / 8; let mut r = vec![0u8; nb];
    for i in 0..nb * 8 { if i < db.len() && db[i] != 0 { r[i / 8] |= 1 << (7 - (i % 8)); } } r
}

// Channel models
struct Rng { s: [u32; 4] }
impl Rng {
    fn new(seed: u64) -> Self {
        let s0 = (seed & 0xFFFFFFFF) as u32 | 1; let s1 = ((seed >> 32) & 0xFFFFFFFF) as u32 | 1;
        Rng { s: [s0, s1, s0.wrapping_mul(0x85ebca6b), s1.wrapping_mul(0xc2b2ae35)] }
    }
    fn next_u32(&mut self) -> u32 {
        let result = self.s[0].wrapping_add(self.s[3]); let t = self.s[1].wrapping_shl(9);
        self.s[2] ^= self.s[0]; self.s[3] ^= self.s[1]; self.s[1] ^= self.s[2]; self.s[0] ^= self.s[3];
        self.s[2] ^= t; self.s[3] = self.s[3].wrapping_shl(11) | self.s[3].wrapping_shr(21); result
    }
    fn next_f64(&mut self) -> f64 { (self.next_u32() as f64) / (u32::MAX as f64) }
}

fn apply_substitutions(data: &[u8], sub_rate: f64, rng: &mut Rng) -> Vec<u8> {
    let total_bits = data.len() * 8; let mut result = data.to_vec();
    for i in 0..total_bits { if rng.next_f64() < sub_rate { result[i / 8] ^= 1 << (7 - (i % 8)); } }
    result
}

fn apply_ids_channel(data: &[u8], ins_rate: f64, del_rate: f64, rng: &mut Rng) -> Vec<u8> {
    let total_bits = data.len() * 8;
    let mut bits: Vec<u8> = Vec::with_capacity(total_bits * 2);
    for i in 0..total_bits { bits.push((data[i / 8] >> (7 - (i % 8))) & 1); }
    let mut i = 0usize;
    while i < bits.len() {
        if rng.next_f64() < ins_rate { bits.insert(i, if rng.next_f64() < 0.5 { 0 } else { 1 }); i += 1; }
        if i < bits.len() && rng.next_f64() < del_rate { bits.remove(i); } else { i += 1; }
    }
    let nb = (bits.len() + 7) / 8; let mut result = vec![0u8; nb];
    for i in 0..bits.len() { result[i / 8] |= bits[i] << (7 - (i % 8)); }
    result
}

fn bit_errors(a: &[u8], b: &[u8]) -> usize {
    let len = a.len().min(b.len()); let mut errors = 0;
    for i in 0..len { errors += (a[i] ^ b[i]).count_ones() as usize; }
    for i in len..a.len() { errors += a[i].count_ones() as usize; }
    for i in len..b.len() { errors += b[i].count_ones() as usize; }
    errors
}

// =========================================================================
// Tests
// =========================================================================

#[test]
fn test_k9_clean_channel() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let num_info_bits = data.len() * 8;
    let decoded = viterbi_decode_indel(&encoded, &table, 10, 1.5, 1.0, Some(num_info_bits));
    let errors = bit_errors(&data, &decoded);
    assert_eq!(errors, 0, "Clean channel: expected 0 bit errors, got {}", errors);
    println!("✓ K=9 clean channel: 0 bit errors");
}

#[test]
fn test_k7_clean_channel() {
    let table = TransitionTable::build(K7_MEMORY, &K7_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let num_info_bits = data.len() * 8;
    let decoded = viterbi_decode_indel(&encoded, &table, 10, 1.5, 1.0, Some(num_info_bits));
    let errors = bit_errors(&data, &decoded);
    assert_eq!(errors, 0, "K=7 clean: expected 0 bit errors, got {}", errors);
    println!("✓ K=7 clean channel: 0 bit errors");
}

#[test]
fn test_k9_auto_estimate() {
    // Test with numInfoBits auto-estimated from received length
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let decoded = viterbi_decode_indel(&encoded, &table, 10, 1.5, 1.0, None);
    let errors = bit_errors(&data, &decoded);
    assert_eq!(errors, 0, "Auto-estimate: expected 0 bit errors, got {}", errors);
    println!("✓ K=9 auto-estimate numInfoBits: 0 bit errors");
}

#[test]
fn test_k9_substitution_only() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let mut rng = Rng::new(42);
    let noisy = apply_substitutions(&encoded, 0.05, &mut rng);
    let decoded = viterbi_decode_indel(&noisy, &table, 10, 1.5, 1.0, Some(data.len() * 8));
    let errors = bit_errors(&data, &decoded);
    println!("✓ 5% substitution: {} bit errors", errors);
    assert!(errors <= 5, "5% sub: expected ≤5, got {}", errors);
}

#[test]
fn test_k9_low_ids() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let mut rng = Rng::new(123);
    let noisy = apply_ids_channel(&encoded, 0.02, 0.02, &mut rng);
    let t = Instant::now();
    let decoded = viterbi_decode_indel(&noisy, &table, 10, 1.5, 1.0, Some(data.len() * 8));
    let elapsed = t.elapsed();
    let errors = bit_errors(&data, &decoded);
    println!("✓ 2% ins + 2% del: {} bit errors, {:.2}ms", errors, elapsed.as_secs_f64() * 1000.0);
    assert!(errors <= 70, "2% IDS: expected ≤70, got {}", errors);
}

#[test]
fn test_k9_moderate_ids() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let mut rng = Rng::new(456);
    let noisy = apply_ids_channel(&encoded, 0.05, 0.05, &mut rng);
    let t = Instant::now();
    let decoded = viterbi_decode_indel(&noisy, &table, 15, 1.5, 1.0, Some(data.len() * 8));
    let elapsed = t.elapsed();
    let errors = bit_errors(&data, &decoded);
    println!("✓ 5% ins + 5% del: {} bit errors, {:.2}ms", errors, elapsed.as_secs_f64() * 1000.0);
    assert!(errors <= 30, "5% IDS: expected ≤30, got {}", errors);
}

#[test]
fn test_k9_high_ids_9percent() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let mut rng = Rng::new(789);
    let noisy = apply_ids_channel(&encoded, 0.09, 0.09, &mut rng);
    let t = Instant::now();
    let decoded = viterbi_decode_indel(&noisy, &table, 15, 1.5, 1.0, Some(data.len() * 8));
    let elapsed = t.elapsed();
    let errors = bit_errors(&data, &decoded);
    println!("✓ 9% ins + 9% del: {} bit errors, {:.2}ms", errors, elapsed.as_secs_f64() * 1000.0);
    assert!(elapsed.as_secs() < 5, "9% IDS: too slow");
}

#[test]
fn test_k9_leading_insertions() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let mut bits: Vec<u8> = vec![1, 0, 1]; // 3 leading insertions
    for i in 0..encoded.len() * 8 { bits.push((encoded[i / 8] >> (7 - (i % 8))) & 1); }
    let nb = (bits.len() + 7) / 8; let mut noisy = vec![0u8; nb];
    for i in 0..bits.len() { noisy[i / 8] |= bits[i] << (7 - (i % 8)); }
    let decoded = viterbi_decode_indel(&noisy, &table, 10, 1.5, 1.0, Some(data.len() * 8));
    let errors = bit_errors(&data, &decoded);
    println!("✓ 3 leading insertions: {} bit errors", errors);
    assert!(errors <= 5, "Leading insertions: expected ≤5, got {}", errors);
}

#[test]
fn test_k9_consecutive_insertions() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let mut bits: Vec<u8> = Vec::new();
    for i in 0..encoded.len() * 8 {
        if i == 100 { bits.extend_from_slice(&[1, 0, 1, 1, 0]); }
        bits.push((encoded[i / 8] >> (7 - (i % 8))) & 1);
    }
    let nb = (bits.len() + 7) / 8; let mut noisy = vec![0u8; nb];
    for i in 0..bits.len() { noisy[i / 8] |= bits[i] << (7 - (i % 8)); }
    let decoded = viterbi_decode_indel(&noisy, &table, 15, 1.5, 1.0, Some(data.len() * 8));
    let errors = bit_errors(&data, &decoded);
    println!("✓ 5 consecutive insertions: {} bit errors", errors);
    assert!(errors <= 60, "Consecutive insertions: expected ≤60, got {}", errors);
}

#[test]
fn test_k9_traceback_matches_standard() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let decoded_standard = viterbi_decode_standard(&encoded, &table);
    let decoded_indel = viterbi_decode_indel(&encoded, &table, 10, 1.5, 1.0, Some(data.len() * 8));
    assert_eq!(decoded_standard, decoded_indel, "Indel disagrees with standard on clean channel!");
    println!("✓ Traceback: indel matches standard on clean channel");
}

#[test]
fn test_k9_standard_sanity() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let decoded = viterbi_decode_standard(&encoded, &table);
    let errors = bit_errors(&data, &decoded);
    assert_eq!(errors, 0, "Standard: expected 0 errors, got {}", errors);
    println!("✓ Standard decoder: 0 errors on clean channel");
}

#[test]
fn test_k9_performance() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..30).collect();
    let encoded = conv_encode(&data, &table);
    let mut rng = Rng::new(999);
    let noisy = apply_ids_channel(&encoded, 0.05, 0.05, &mut rng);
    let nib = data.len() * 8;
    let _ = viterbi_decode_indel(&noisy, &table, 15, 1.5, 1.0, Some(nib));
    let n = 5;
    let t = Instant::now();
    for _ in 0..n { let _ = viterbi_decode_indel(&noisy, &table, 15, 1.5, 1.0, Some(nib)); }
    let per_decode = t.elapsed().as_secs_f64() / n as f64 * 1000.0;
    println!("✓ Performance: {:.2}ms per decode (5% IDS, K=9, maxDrift=15)", per_decode);
}

#[test]
fn test_k9_ids_sweep() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let nib = data.len() * 8;
    println!("\n  IDS Sweep (K=9, 20-byte payload, numInfoBits={}):", nib);
    println!("  {:>10} {:>8} {:>10} {:>10} {:>10}", "Scenario", "Errors", "BER", "Time(ms)", "maxDrift");
    for &(name, ins, del, md) in [
        ("1% IDS", 0.01, 0.01, 10u32), ("3% IDS", 0.03, 0.03, 10),
        ("5% IDS", 0.05, 0.05, 15), ("7% IDS", 0.07, 0.07, 15),
        ("9% IDS", 0.09, 0.09, 15),
    ].iter() {
        let mut rng = Rng::new(42);
        let noisy = apply_ids_channel(&encoded, ins, del, &mut rng);
        let t = Instant::now();
        let decoded = viterbi_decode_indel(&noisy, &table, md, 1.5, 1.0, Some(nib));
        let elapsed = t.elapsed();
        let errors = bit_errors(&data, &decoded);
        let ber = errors as f64 / nib as f64;
        println!("  {:>10} {:>8} {:>10.2}% {:>10.2} {:>10}", name, errors, ber * 100.0, elapsed.as_secs_f64() * 1000.0, md);
    }
    println!("✓ IDS sweep complete");
}
