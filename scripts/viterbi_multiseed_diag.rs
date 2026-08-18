//! Multi-seed diagnostic: measure average Viterbi error rate at each IDS level

const K9_MEMORY: usize = 8;
const K9_GENERATORS: [u16; 2] = [0o561, 0o753];

struct TransitionTable { outputs: Vec<u16>, next_states: Vec<u16>, num_states: usize, memory: usize, rate: usize }
impl TransitionTable {
    fn build(memory: usize, generators: &[u16], rate: usize) -> Self {
        let ns = 1usize << memory; let mut outputs = vec![0u16; ns * 2]; let mut next_states = vec![0u16; ns * 2];
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
    let ob = (data.len() * 8 + table.memory) * table.rate; let mut r = vec![0u8; (ob + 7) / 8]; let mut st = 0usize; let mut p = 0usize;
    for &b in data { for bi in 0..8 { let ib = ((b >> (7 - bi)) & 1) as usize; let o = table.output(st, ib); st = table.next_state(st, ib); for g in 0..table.rate { if (o >> (table.rate - 1 - g)) & 1 != 0 { r[p / 8] |= 1 << (7 - (p % 8)); } p += 1; } } }
    for _ in 0..table.memory { let o = table.output(st, 0); st = table.next_state(st, 0); for g in 0..table.rate { if (o >> (table.rate - 1 - g)) & 1 != 0 && p / 8 < r.len() { r[p / 8] |= 1 << (7 - (p % 8)); } p += 1; } }
    r.truncate((p + 7) / 8); r
}

fn viterbi_decode_indel(received: &[u8], table: &TransitionTable, max_drift: u32, ins_pen: f64, del_pen: f64, nib_opt: Option<usize>) -> Vec<u8> {
    let ns = table.num_states; let memory = table.memory; let rate = table.rate; let rb = received.len() * 8;
    if rb == 0 { return Vec::new(); }
    let num_info_bits = if let Some(n) = nib_opt { n } else if rb >= memory * rate { (rb / rate) - memory } else { rb / rate };
    let total_info_steps = num_info_bits + memory; let total_channel_uses = total_info_steps * rate;
    let actual_drift = if rb > total_channel_uses { rb - total_channel_uses } else { total_channel_uses - rb };
    let adaptive_md = { let base = max_drift as usize; let needed = actual_drift + actual_drift / 2 + 3; if needed > base { needed.min(base * 2) } else { base } };
    let md = adaptive_md as i32; let nd = 2 * adaptive_md + 1; let p0 = ns * nd; let ta = ns * 3 * nd; let inf = f64::INFINITY;
    let mut recv_bits: Vec<u8> = vec![0u8; rb]; for i in 0..rb { recv_bits[i] = (received[i / 8] >> (7 - (i % 8))) & 1; }
    let max_steps = total_channel_uses + adaptive_md + 10; let total_cells = (max_steps + 1) * ta;
    let mut pm = vec![inf; total_cells]; let mut bp = vec![0u32; total_cells]; let mut tt = vec![0u8; total_cells]; let mut ib_arr = vec![-1i8; total_cells];
    pm[adaptive_md] = 0.0;
    for cs in 0..ns { for di in 0..nd-1 { let aug=cs*nd+di; let m=pm[aug]; if m>=1e29{continue;} let drift=di as i32-md; let rp=drift; if rp<0||rp as usize>=rb{continue;} let na=cs*nd+di+1; let nm=m+ins_pen; if nm<pm[na]{pm[na]=nm;bp[na]=aug as u32;tt[na]=1;ib_arr[na]=-1;} } for pend in 0..2usize { for di in 0..nd-1 { let aug=p0+(cs*2+pend)*nd+di; let m=pm[aug]; if m>=1e29{continue;} let drift=di as i32-md; let rp=drift; if rp<0||rp as usize>=rb{continue;} let na=p0+(cs*2+pend)*nd+di+1; let nm=m+ins_pen; if nm<pm[na]{pm[na]=nm;bp[na]=aug as u32;tt[na]=1;ib_arr[na]=-1;} } } }
    for step in 0..total_channel_uses {
        let so=step*ta; let no=(step+1)*ta;
        for cs in 0..ns {
            for di in 0..nd { let aug=cs*nd+di; let m=pm[so+aug]; if m>=1e29{continue;} let drift=di as i32-md; let recv_pos=step as i32+drift;
                for ib in 0..2usize { let out=table.output(cs,ib) as usize; let emitted=(out>>(rate-1))&1;
                    if recv_pos>=0&&(recv_pos as usize)<rb { let rp=recv_pos as usize; let rb_bit=recv_bits[rp] as usize; let dist=if emitted!=rb_bit{1.0}else{0.0}; let na=p0+(cs*2+ib)*nd+di; let nm=m+dist; if nm<pm[no+na]{pm[no+na]=nm;bp[no+na]=aug as u32;tt[no+na]=0;ib_arr[no+na]=-1;} }
                    if di>0 { let na=p0+(cs*2+ib)*nd+di-1; let nm=m+del_pen; if nm<pm[no+na]{pm[no+na]=nm;bp[no+na]=aug as u32;tt[no+na]=2;ib_arr[no+na]=-1;} }
                }
            }
            for pend in 0..2usize { for di in 0..nd { let aug=p0+(cs*2+pend)*nd+di; let m=pm[so+aug]; if m>=1e29{continue;} let drift=di as i32-md; let recv_pos=step as i32+drift; let out=table.output(cs,pend) as usize; let ns_=table.next_state(cs,pend); let emitted=(out>>(rate-2))&1;
                if recv_pos>=0&&(recv_pos as usize)<rb { let rp=recv_pos as usize; let rb_bit=recv_bits[rp] as usize; let dist=if emitted!=rb_bit{1.0}else{0.0}; let na=ns_*nd+di; let nm=m+dist; if nm<pm[no+na]{pm[no+na]=nm;bp[no+na]=aug as u32;tt[no+na]=0;ib_arr[no+na]=pend as i8;} }
                if di>0 { let na=ns_*nd+di-1; let nm=m+del_pen; if nm<pm[no+na]{pm[no+na]=nm;bp[no+na]=aug as u32;tt[no+na]=2;ib_arr[no+na]=pend as i8;} }
            }}
        }
        let sfr=step as i32+1;
        for cs in 0..ns { for di in 0..nd-1 { let aug=cs*nd+di; let m=pm[no+aug]; if m>=1e29{continue;} let drift=di as i32-md; let rp=sfr+drift; if rp<0||rp as usize>=rb{continue;} let na=cs*nd+di+1; let nm=m+ins_pen; if nm<pm[no+na]{pm[no+na]=nm;bp[no+na]=aug as u32;tt[no+na]=1;ib_arr[no+na]=-1;} } for pend in 0..2usize { for di in 0..nd-1 { let aug=p0+(cs*2+pend)*nd+di; let m=pm[no+aug]; if m>=1e29{continue;} let drift=di as i32-md; let rp=sfr+drift; if rp<0||rp as usize>=rb{continue;} let na=p0+(cs*2+pend)*nd+di+1; let nm=m+ins_pen; if nm<pm[no+na]{pm[no+na]=nm;bp[no+na]=aug as u32;tt[no+na]=1;ib_arr[no+na]=-1;} } } }
    }
    let fo=total_channel_uses*ta; let mut best_aug=0usize; let mut best_metric=inf;
    for cs in 0..ns { for di in 0..nd { let aug=cs*nd+di; let m=pm[fo+aug]; if m>=1e29{continue;} let drift=(di as i32-md).abs() as f64; let penalty=if cs!=0{50.0}else{0.0}+drift*0.5; let total=m+penalty; if total<best_metric{best_metric=total;best_aug=aug;} }}
    if best_metric>=1e29 { return vec![0u8;num_info_bits/8]; }
    let mut decoded_bits: Vec<u8> = Vec::with_capacity(num_info_bits); let mut step=total_channel_uses; let mut aug=best_aug; let mut safety=total_channel_uses*4;
    while step>0&&safety>0 { safety-=1; let offset=step*ta+aug; let t=tt[offset]; let ib=ib_arr[offset]; let prev_aug=bp[offset] as usize;
        match t { 0=>{if ib>=0{decoded_bits.push(ib as u8);}step-=1;aug=prev_aug;} 2=>{if ib>=0{decoded_bits.push(ib as u8);}step-=1;aug=prev_aug;} 1=>{aug=prev_aug;} _=>{step-=1;aug=prev_aug;} }
    }
    decoded_bits.reverse(); if decoded_bits.len()>memory {decoded_bits.truncate(decoded_bits.len()-memory);}
    let nb=decoded_bits.len()/8; let mut result=vec![0u8;nb]; for i in 0..nb*8 { if i<decoded_bits.len()&&decoded_bits[i]!=0{result[i/8]|=1<<(7-(i%8));} }
    result
}

struct Rng { s: [u32; 4] }
impl Rng { fn new(seed: u64) -> Self { let s0=(seed&0xFFFFFFFF)as u32|1; let s1=((seed>>32)&0xFFFFFFFF)as u32|1; Rng{s:[s0,s1,s0.wrapping_mul(0x85ebca6b),s1.wrapping_mul(0xc2b2ae35)]} }
    fn next_u32(&mut self)->u32 { let result=self.s[0].wrapping_add(self.s[3]); let t=self.s[1].wrapping_shl(9); self.s[2]^=self.s[0]; self.s[3]^=self.s[1]; self.s[1]^=self.s[2]; self.s[0]^=self.s[3]; self.s[2]^=t; self.s[3]=self.s[3].wrapping_shl(11)|self.s[3].wrapping_shr(21); result }
    fn next_f64(&mut self)->f64 { (self.next_u32() as f64)/(u32::MAX as f64) } }

fn apply_ids_channel(data: &[u8], ins_rate: f64, del_rate: f64, rng: &mut Rng) -> Vec<u8> {
    let total_bits=data.len()*8; let mut bits: Vec<u8>=Vec::with_capacity(total_bits*2); for i in 0..total_bits{bits.push((data[i/8]>>(7-(i%8)))&1);} let mut i=0usize;
    while i<bits.len(){if rng.next_f64()<ins_rate{bits.insert(i,if rng.next_f64()<0.5{0}else{1});i+=1;} if i<bits.len()&&rng.next_f64()<del_rate{bits.remove(i);}else{i+=1;}}
    let nb=(bits.len()+7)/8; let mut result=vec![0u8;nb]; for i in 0..bits.len(){result[i/8]|=bits[i]<<(7-(i%8));} result
}

fn bit_errors(a: &[u8], b: &[u8]) -> usize { let len=a.len().min(b.len()); let mut errors=0; for i in 0..len{errors+=(a[i]^b[i]).count_ones() as usize;} for i in len..a.len(){errors+=a[i].count_ones() as usize;} for i in len..b.len(){errors+=b[i].count_ones() as usize;} errors }

fn main() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let nib = data.len() * 8;
    let num_seeds = 20usize;

    println!("=== Multi-Seed Viterbi Diagnostic ===");
    println!("Data: {} bytes ({} info bits), Encoded: {} bits", data.len(), nib, encoded.len()*8);
    println!();

    for &(name, ins, del, md) in [
        ("1%", 0.01, 0.01, 10u32),
        ("2%", 0.02, 0.02, 15),
        ("3%", 0.03, 0.03, 15),
        ("5%", 0.05, 0.05, 15),
        ("7%", 0.07, 0.07, 20),
        ("9%", 0.09, 0.09, 20),
    ].iter() {
        let mut total_errors = 0usize;
        let mut total_zero = 0usize;
        let mut max_errors = 0usize;
        let mut min_errors = usize::MAX;
        for seed in 1..=num_seeds {
            let mut rng = Rng::new(seed as u64 * 1000 + ins as u64 * 100000);
            let noisy = apply_ids_channel(&encoded, ins, del, &mut rng);
            let decoded = viterbi_decode_indel(&noisy, &table, md, 1.5, 1.5, Some(nib));
            let err = bit_errors(&data, &decoded);
            total_errors += err;
            if err == 0 { total_zero += 1; }
            if err > max_errors { max_errors = err; }
            if err < min_errors { min_errors = err; }
        }
        let avg = total_errors as f64 / num_seeds as f64;
        let avg_ber = avg / nib as f64;
        let zero_pct = total_zero as f64 / num_seeds as f64 * 100.0;
        println!("  {:>3} IDS: avg {} errors (BER {:.2}%), zero {:.0}%, range [{}, {}]",
            name, avg as usize, avg_ber * 100.0, zero_pct, min_errors, max_errors);
    }

    // Also test with penalty sweep at 5% IDS
    println!("\n=== Penalty Sweep at 5% IDS (20 seeds) ===");
    for ins_p in [1.0, 1.5, 2.0, 3.0, 5.0] {
        for del_p in [1.0, 1.5, 2.0, 3.0] {
            let mut total_errors = 0usize;
            let mut total_zero = 0usize;
            for seed in 1..=num_seeds {
                let mut rng = Rng::new(seed as u64 * 7777);
                let noisy = apply_ids_channel(&encoded, 0.05, 0.05, &mut rng);
                let decoded = viterbi_decode_indel(&noisy, &table, 15, ins_p, del_p, Some(nib));
                let err = bit_errors(&data, &decoded);
                total_errors += err;
                if err == 0 { total_zero += 1; }
            }
            let avg = total_errors as f64 / num_seeds as f64;
            let zero_pct = total_zero as f64 / num_seeds as f64 * 100.0;
            if avg < 40.0 {
                println!("  ins={}, del={}: avg {:.0} errors, zero {:.0}%", ins_p, del_p, avg, zero_pct);
            }
        }
    }
}
