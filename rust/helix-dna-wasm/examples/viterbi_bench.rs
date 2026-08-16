use std::time::Instant;

fn main() {
    let memory = 8usize;
    let generators: [u16; 2] = [0o561, 0o753];
    let rate = 2usize;
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
                    if (gen >> b) & 1 != 0 { bit ^= (reg >> b) as u16 & 1; }
                }
                output = (output << 1) | bit;
            }
            let next_state = (reg >> 1) & (num_states - 1);
            outputs[state * 2 + input] = output;
            next_states[state * 2 + input] = next_state as u16;
        }
    }

    let num_info_bytes = 40usize;
    let num_info_bits = num_info_bytes * 8;
    let total_info_steps = num_info_bits + memory;
    let total_channel_uses = total_info_steps * rate;

    let mut info_bits = vec![0u8; total_info_steps];
    for i in 0..num_info_bits { info_bits[i] = ((i * 7 + 13) % 2) as u8; }

    let mut state = 0usize;
    let mut encoded_bits = Vec::new();
    for &input in info_bits.iter() {
        let idx = state * 2 + input as usize;
        let output = outputs[idx] as usize;
        state = next_states[idx] as usize;
        for r in 0..rate {
            encoded_bits.push(((output >> (rate - 1 - r)) & 1) as u8);
        }
    }

    let mut noisy = encoded_bits.clone();
    for i in (0..noisy.len()).step_by(20) { if i < noisy.len() { noisy[i] ^= 1; } }

    // Benchmark standard Viterbi
    let iters = 100;
    let t0 = Instant::now();
    for _ in 0..iters {
        let mut pm = vec![f64::INFINITY; (total_info_steps + 1) * num_states];
        let mut ps = vec![0u16; (total_info_steps + 1) * num_states];
        let mut ib = vec![0u8; (total_info_steps + 1) * num_states];
        pm[0] = 0.0;
        for step in 0..total_info_steps {
            let ch = step * rate;
            for s in 0..num_states {
                let m = pm[step * num_states + s];
                if !m.is_finite() { continue; }
                for inp in 0..2usize {
                    let idx = s * 2 + inp;
                    let out = outputs[idx] as usize;
                    let ns = next_states[idx] as usize;
                    let mut d = 0.0f64;
                    for r in 0..rate {
                        let e = (out >> (rate - 1 - r)) & 1;
                        let p = ch + r;
                        if p < noisy.len() && e != noisy[p] as usize { d += 1.0; }
                    }
                    let nm = m + d;
                    let t = (step + 1) * num_states + ns;
                    if nm < pm[t] { pm[t] = nm; ps[t] = s as u16; ib[t] = inp as u8; }
                }
            }
        }
    }
    let std_ms = t0.elapsed().as_millis() as f64 / iters as f64;
    println!("K=9 standard Viterbi: {:.1} ms/decode ({} info bits)", std_ms, num_info_bits);

    // Indel Viterbi stats
    let max_drift: i32 = 15;
    let w = (2 * max_drift + 1) as usize;
    let num_aug = num_states * 3 * w;
    let total_cells = (total_channel_uses + max_drift as usize + 10) * num_aug;
    println!("Indel augmented states: {}, total cells: {}, ~{} MB", 
        num_aug, total_cells, total_cells * 14 / 1_000_000);
    println!("JS Viterbi: ~800 ms/decode | Rust estimate: ~5-15 ms/decode (50-160x speedup)");
}
