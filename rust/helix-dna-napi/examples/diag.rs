//! Diagnostic test: trace the Viterbi path for 5 consecutive insertions

use std::time::Instant;

const K9_MEMORY: usize = 8;
const K9_GENERATORS: [u16; 2] = [0o561, 0o753];

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

fn bit_errors(a: &[u8], b: &[u8]) -> usize {
    let len = a.len().min(b.len()); let mut errors = 0;
    for i in 0..len { errors += (a[i] ^ b[i]).count_ones() as usize; }
    for i in len..a.len() { errors += a[i].count_ones() as usize; }
    for i in len..b.len() { errors += b[i].count_ones() as usize; }
    errors
}

fn insert_bits_at(encoded: &[u8], pos: usize, ins_bits: &[u8]) -> Vec<u8> {
    let mut bits: Vec<u8> = Vec::new();
    for i in 0..encoded.len() * 8 { if i == pos { bits.extend_from_slice(ins_bits); } bits.push((encoded[i / 8] >> (7 - (i % 8))) & 1); }
    let nb = (bits.len() + 7) / 8; let mut result = vec![0u8; nb];
    for i in 0..bits.len() { result[i / 8] |= bits[i] << (7 - (i % 8)); }
    result
}

// Full trellis with diagnostics
fn viterbi_decode_indel_diag(
    received: &[u8], table: &TransitionTable, max_drift: u32,
    ins_pen: f64, del_pen: f64, num_info_bits: usize, verbose: bool,
) -> (Vec<u8>, usize, usize, usize, f64) {
    let ns = table.num_states; let memory = table.memory; let rate = table.rate;
    let rb = received.len() * 8;
    if rb == 0 { return (Vec::new(), 0, 0, 0, 0.0); }
    let total_info_steps = num_info_bits + memory;
    let total_channel_uses = total_info_steps * rate;
    let md = max_drift as i32; let nd = 2 * max_drift as usize + 1;
    let p0 = ns * nd; let ta = ns * 3 * nd; let inf = f64::INFINITY;
    let mut recv_bits: Vec<u8> = vec![0u8; rb];
    for i in 0..rb { recv_bits[i] = (received[i / 8] >> (7 - (i % 8))) & 1; }
    let max_steps = total_channel_uses + max_drift as usize + 10;
    let total_cells = (max_steps + 1) * ta;
    let mut pm = vec![inf; total_cells];
    let mut bp = vec![0u32; total_cells];
    let mut tt = vec![0u8; total_cells];
    let mut ib_arr = vec![-1i8; total_cells];
    let start_aug = max_drift as usize;
    pm[start_aug] = 0.0;

    // I-chain step 0
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
        // I-chain at step+1
        let sfr=step as i32+1;
        for cs in 0..ns { for di in 0..nd-1 { let aug=cs*nd+di; let m=pm[no+aug]; if m>=1e29{continue;} let drift=di as i32-md; let rp=sfr+drift; if rp<0||rp as usize>=rb{continue;} let na=cs*nd+di+1; let nm=m+ins_pen; if nm<pm[no+na]{pm[no+na]=nm;bp[no+na]=aug as u32;tt[no+na]=1;ib_arr[no+na]=-1;} } for pend in 0..2usize { for di in 0..nd-1 { let aug=p0+(cs*2+pend)*nd+di; let m=pm[no+aug]; if m>=1e29{continue;} let drift=di as i32-md; let rp=sfr+drift; if rp<0||rp as usize>=rb{continue;} let na=p0+(cs*2+pend)*nd+di+1; let nm=m+ins_pen; if nm<pm[no+na]{pm[no+na]=nm;bp[no+na]=aug as u32;tt[no+na]=1;ib_arr[no+na]=-1;} } } }
    }

    let fo=total_channel_uses*ta;
    let mut best_aug=0usize; let mut best_metric=inf; let mut best_cs=0usize; let mut best_drift=0i32;
    for cs in 0..ns { for di in 0..nd { let aug=cs*nd+di; let m=pm[fo+aug]; if m>=1e29{continue;} let drift=(di as i32-md).abs() as f64; let penalty=if cs!=0{50.0}else{0.0}+drift*0.5; let total=m+penalty; if total<best_metric{best_metric=total;best_aug=aug;best_cs=cs;best_drift=di as i32-md;} }}

    if verbose {
        println!("  Final: cs={}, drift={}, metric={:.2}, total={:.2}", best_cs, best_drift, pm[fo+best_aug], best_metric);
        let mut finals: Vec<(usize,i32,f64,f64)>=Vec::new();
        for cs in 0..ns { for di in 0..nd { let aug=cs*nd+di; let m=pm[fo+aug]; if m>=1e29{continue;} let drift=(di as i32-md).abs() as f64; let penalty=if cs!=0{50.0}else{0.0}+drift*0.5; finals.push((cs,di as i32-md,m,m+penalty)); }}
        finals.sort_by(|a,b|a.3.partial_cmp(&b.3).unwrap());
        println!("  Top 5 finals:");
        for (i,(cs,drift,m,total)) in finals.iter().take(5).enumerate() { println!("    #{}: cs={}, drift={}, m={:.2}, total={:.2}", i+1,cs,drift,m,total); }
    }

    if best_metric>=1e29 { return (vec![0u8;num_info_bits/8],0,0,0,best_metric); }

    let mut decoded_bits: Vec<u8> = Vec::with_capacity(num_info_bits);
    let mut step=total_channel_uses; let mut aug=best_aug; let mut safety=total_channel_uses*4;
    let mut nm=0usize; let mut ni=0usize; let mut nd=0usize;
    while step>0&&safety>0 { safety-=1; let offset=step*ta+aug; let t=tt[offset]; let ib=ib_arr[offset]; let prev_aug=bp[offset] as usize;
        match t { 0=>{if ib>=0{decoded_bits.push(ib as u8);}step-=1;aug=prev_aug;nm+=1;} 2=>{if ib>=0{decoded_bits.push(ib as u8);}step-=1;aug=prev_aug;nd+=1;} 1=>{aug=prev_aug;ni+=1;} _=>{step-=1;aug=prev_aug;} }
    }
    if verbose { println!("  Traceback: M={}, I={}, D={}, bits={}, step_end={}", nm,ni,nd,decoded_bits.len(),step); }
    decoded_bits.reverse();
    if decoded_bits.len()>memory {decoded_bits.truncate(decoded_bits.len()-memory);}
    let nb=decoded_bits.len()/8; let mut result=vec![0u8;nb];
    for i in 0..nb*8 { if i<decoded_bits.len()&&decoded_bits[i]!=0{result[i/8]|=1<<(7-(i%8));} }
    (result,nm,ni,nd,best_metric)
}

fn main() {
    let table = TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2);
    let data: Vec<u8> = (0u8..20).collect();
    let encoded = conv_encode(&data, &table);
    let nib = data.len() * 8;

    println!("=== Viterbi Diagnostic ===");
    println!("Data: {} bytes, Encoded: {} bits, totalChannelUses: {}", data.len(), encoded.len()*8, (nib+K9_MEMORY)*2);

    for n_ins in 1..=5 {
        let ins_bits: Vec<u8> = (0..n_ins).map(|i| (i % 2) as u8).collect();
        println!("\n--- {} insertion(s) at bit 100 ---", n_ins);
        let noisy = insert_bits_at(&encoded, 100, &ins_bits);
        let (decoded, nm, ni, nd, bm) = viterbi_decode_indel_diag(&noisy, &table, 15, 1.5, 1.5, nib, true);
        let errors = bit_errors(&data, &decoded);
        println!("  => {} errors, decoded={} bytes (expected {})", errors, decoded.len(), data.len());
    }

    println!("\n=== Penalty sweep (5 ins at 100) ===");
    let noisy5 = insert_bits_at(&encoded, 100, &[1,0,1,1,0]);
    for ins_p in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0] {
        for del_p in [0.5, 1.0, 1.5, 2.0, 3.0] {
            let (d,_,_,_,_) = viterbi_decode_indel_diag(&noisy5, &table, 15, ins_p, del_p, nib, false);
            let e = bit_errors(&data, &d);
            if e <= 5 { println!("  ins={}, del={}: {} errors **", ins_p, del_p, e); }
        }
    }

    println!("\n=== maxDrift sweep (5 ins at 100) ===");
    for md in [5, 10, 15, 20, 25, 30] {
        let (d,nm,ni,nd,bm) = viterbi_decode_indel_diag(&noisy5, &table, md, 1.5, 1.5, nib, true);
        let e = bit_errors(&data, &d);
        println!("  maxDrift={}: {} errors, M={}, I={}, D={}", md, e, nm, ni, nd);
    }
}
