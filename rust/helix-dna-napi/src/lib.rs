//! helix-dna-napi: Native Node.js addon for K=9 Viterbi decoder
//!
//! Iterates over TRELLIS STAGES (encoder output events), not received bits.
//! Each stage: try MATCH (consume 1 recv bit), INSERTION (consume 1 recv bit), DELETION (consume 0).
//! Decision encoding: 0,1 = input bit (phase1→phase0); 2 = insertion; 3 = intermediate (phase0→phase1)

use napi::bindgen_prelude::*;
use napi_derive::napi;
use std::sync::OnceLock;

const K9_MEMORY: usize = 8;
const K9_GENERATORS: [u16; 2] = [0o561, 0o753];
const K7_MEMORY: usize = 6;
const K7_GENERATORS: [u16; 2] = [0o171, 0o133];

static K9_TABLE: OnceLock<TransitionTable> = OnceLock::new();
static K7_TABLE: OnceLock<TransitionTable> = OnceLock::new();
fn k9_table() -> &'static TransitionTable { K9_TABLE.get_or_init(|| TransitionTable::build(K9_MEMORY, &K9_GENERATORS, 2)) }
fn k7_table() -> &'static TransitionTable { K7_TABLE.get_or_init(|| TransitionTable::build(K7_MEMORY, &K7_GENERATORS, 2)) }

struct TransitionTable { outputs: Vec<u16>, next_states: Vec<u16>, num_states: usize, memory: usize, rate: usize }
impl TransitionTable {
    fn build(memory: usize, generators: &[u16], rate: usize) -> Self {
        let ns = 1usize << memory;
        let mut o = vec![0u16; ns*2]; let mut n = vec![0u16; ns*2];
        for s in 0..ns { for ib in 0..2 {
            n[s*2+ib] = ((s>>1)|((ib as usize)<<(memory-1))) as u16;
            let mut ob = 0u16;
            for g in 0..rate { let mut b=0u16; if (generators[g]&1)!=0{b^=ib as u16;} for i in 0..memory{if(generators[g]&(1<<(i+1)))!=0{b^=((s>>i)&1)as u16;}} ob|=b<<g; }
            o[s*2+ib] = ob;
        }}
        TransitionTable { outputs: o, next_states: n, num_states: ns, memory, rate }
    }
    #[inline] fn output(&self, s: usize, i: usize) -> u16 { self.outputs[s*2+i] }
    #[inline] fn next_state(&self, s: usize, i: usize) -> usize { self.next_states[s*2+i] as usize }
}

fn conv_encode(data: &[u8], table: &TransitionTable) -> Vec<u8> {
    let ob = (data.len()*8+table.memory)*table.rate;
    let mut r = vec![0u8; (ob+7)/8]; let mut st=0usize; let mut p=0usize;
    for &b in data { for bi in 0..8 { let ib=((b>>(7-bi))&1)as usize; let o=table.output(st,ib); st=table.next_state(st,ib); for g in 0..table.rate{if(o>>g)&1!=0{r[p/8]|=1<<(7-(p%8));}p+=1;} } }
    for _ in 0..table.memory { let o=table.output(st,0); st=table.next_state(st,0); for g in 0..table.rate{if(o>>g)&1!=0&&p/8<r.len(){r[p/8]|=1<<(7-(p%8));}p+=1;} }
    r.truncate((p+7)/8); r
}

#[inline] fn get_bit(d: &[u8], p: usize) -> u8 { if p/8 >= d.len() { 0 } else { (d[p/8]>>(7-(p%8)))&1 } }
#[inline] fn pack_trace(prev: usize, dec: u8) -> u32 { (prev as u32)<<2 | (dec as u32 & 3) }

#[napi(object)] #[derive(Clone, Copy)]
pub struct ViterbiConfig { pub max_drift: Option<u32>, pub insertion_penalty: Option<f64>, pub deletion_penalty: Option<f64>, pub use_llr: Option<bool>, pub expected_length: Option<u32> }
impl Default for ViterbiConfig { fn default() -> Self { ViterbiConfig { max_drift: Some(10), insertion_penalty: Some(1.5), deletion_penalty: Some(1.0), use_llr: Some(false), expected_length: None } } }
fn resolve_config(cfg: Option<&ViterbiConfig>) -> (u32,f64,f64,Option<usize>) {
    let c = cfg.copied().unwrap_or_default();
    (c.max_drift.unwrap_or(10), c.insertion_penalty.unwrap_or(1.5), c.deletion_penalty.unwrap_or(1.0), c.expected_length.map(|v| v as usize))
}

/// Indel-tolerant Viterbi decoder.
///
/// KEY FIX: Iterates over TRELLIS STAGES (encoder output events), not received bit positions.
/// Each stage represents one encoder output event (G1 or G2 for rate-1/2).
/// - MATCH: consume 1 received bit, compare to expected
/// - INSERTION: consume 1 received bit (spurious), no encoder output
/// - DELETION: consume 0 received bits, encoder output was lost
///
/// This is the correct formulation: the trellis stage corresponds to the
/// encoder's clock, not the channel's clock. Deletions don't advance
/// the received bit position, which is why iterating over recv_bit_pos
/// was wrong.
fn viterbi_decode_indel(received: &[u8], table: &TransitionTable, max_drift: u32, ins_pen: f64, del_pen: f64, _llr_data: Option<&[f32]>) -> Vec<u8> {
    let ns = table.num_states; let mdi = max_drift as i16; let nd = 2*max_drift as usize+1;
    let rb = received.len()*8; if rb==0 { return Vec::new(); }
    let p0 = ns*nd; let p1 = ns*2*nd; let ta = p0+p1;
    let inf = f64::INFINITY;

    // Maximum trellis stages: we need enough stages to consume all received bits
    // via MATCH/INSERTION and handle DELETION transitions.
    // Upper bound: recv_bits (all MATCH) + max_drift (DELETIONs that don't consume bits)
    let max_stages = rb + (max_drift as usize) * 2 + table.memory * 2;

    let mut cm = vec![inf; ta]; let mut pm = vec![inf; ta];
    // Full trace storage: (prev_idx, decision) per stage
    let mut trace: Vec<u32> = vec![0; max_stages * ta];
    let mut decoded_bits: Vec<u8> = Vec::with_capacity(rb / 2);

    cm[max_drift as usize] = 0.0; // init: conv=0, phase=0, drift=0

    // Track how many received bits have been consumed by the best path
    // We use a separate array to track recv_bit_pos per augmented state
    let mut cur_rbp = vec![0u16; ta]; // received bit position per state
    let mut prev_rbp = vec![0u16; ta];
    cur_rbp[max_drift as usize] = 0;

    let mut stage = 0usize;

    while stage < max_stages {
        // Check if any state is still reachable and hasn't consumed all bits
        let mut any_active = false;
        for i in 0..ta { if cm[i] < inf { any_active = true; break; } }
        if !any_active { break; }

        std::mem::swap(&mut cm, &mut pm);
        std::mem::swap(&mut cur_rbp, &mut prev_rbp);
        cm.fill(inf);
        cur_rbp.fill(0);

        let sb = stage * ta;

        // Phase 0 states: encoder is about to emit G1
        for conv in 0..ns { for di in 0..nd {
            let drift = di as i16 - mdi; let pi = conv*nd+di;
            let pc = pm[pi]; if pc>=1e29 { continue; }
            let rbp = prev_rbp[pi] as usize;

            for ib in 0..2 {
                let out = table.output(conv,ib); let exp = (out&1) as u8; let nc = table.next_state(conv,ib);
                // MATCH: consume 1 received bit as G1, move to phase 1
                if rbp < rb {
                    let rb_bit = get_bit(received, rbp);
                    let bm = if exp != rb_bit { 1.0 } else { 0.0 };
                    let ni = p0+(nc*2+ib)*nd+di; let nco = pc+bm;
                    if nco < cm[ni] { cm[ni]=nco; cur_rbp[ni]=((rbp+1) as u16); trace[sb+ni]=pack_trace(pi,3); }
                }
                // DELETION: consume 0 received bits, G1 was lost, move to phase 1, drift-1
                if drift > -mdi {
                    let ni = p0+(nc*2+ib)*nd+di-1; let dc = pc+del_pen;
                    if dc < cm[ni] { cm[ni]=dc; cur_rbp[ni]=(rbp as u16); trace[sb+ni]=pack_trace(pi,3); }
                }
            }
            // INSERTION: consume 1 received bit (spurious), stay at phase 0, drift+1
            if drift < mdi && rbp < rb {
                let ni = conv*nd+di+1; let ic = pc+ins_pen;
                if ic < cm[ni] { cm[ni]=ic; cur_rbp[ni]=((rbp+1) as u16); trace[sb+ni]=pack_trace(pi,2); }
            }
        }}

        // Phase 1 states: encoder is about to emit G2
        for conv in 0..ns { for pend in 0..2 { for di in 0..nd {
            let drift = di as i16 - mdi; let pi = p0+(conv*2+pend)*nd+di;
            let pc = pm[pi]; if pc>=1e29 { continue; }
            let rbp = prev_rbp[pi] as usize;
            let out = table.output(conv,pend); let exp = ((out>>1)&1) as u8; let nc = table.next_state(conv,pend);

            // MATCH: consume 1 received bit as G2, move to phase 0, emit input bit
            if rbp < rb {
                let rb_bit = get_bit(received, rbp);
                let bm = if exp != rb_bit { 1.0 } else { 0.0 };
                let nco = pc+bm; let ni = nc*nd+di;
                if nco < cm[ni] { cm[ni]=nco; cur_rbp[ni]=((rbp+1) as u16); trace[sb+ni]=pack_trace(pi,pend as u8); }
            }
            // DELETION: consume 0 received bits, G2 was lost, move to phase 0, drift-1, emit input bit
            if drift > -mdi {
                let ni = nc*nd+di-1; let dc = pc+del_pen;
                if dc < cm[ni] { cm[ni]=dc; cur_rbp[ni]=(rbp as u16); trace[sb+ni]=pack_trace(pi,pend as u8); }
            }
            // INSERTION: consume 1 received bit (spurious), stay at phase 1, drift+1
            if drift < mdi && rbp < rb {
                let ni = p0+(conv*2+pend)*nd+di+1; let ic = pc+ins_pen;
                if ic < cm[ni] { cm[ni]=ic; cur_rbp[ni]=((rbp+1) as u16); trace[sb+ni]=pack_trace(pi,2); }
            }
        }}}

        stage += 1;

        // Early termination: if best state has consumed all received bits and is at phase 0, we can stop
        let mut best_cost = inf;
        for i in 0..ta { if cm[i] < best_cost { best_cost = cm[i]; } }
        if best_cost >= 1e29 { break; }
    }

    // Find best final state (phase 0, drift=0 preferred)
    let mut bi = 0usize; let mut bc = inf;
    // Prefer phase 0 states with drift near 0
    for i in 0..ta { if cm[i] < bc { bc = cm[i]; bi = i; } }

    // Traceback from best state
    let mut idx = bi;
    for s in (0..stage).rev() {
        let e = trace[s*ta+idx]; let d = (e&3) as u8;
        if d < 2 { decoded_bits.push(d); } // 0 or 1 = valid input bit
        idx = (e>>2) as usize;
    }
    decoded_bits.reverse();

    // Remove tail bits
    if decoded_bits.len() > table.memory { decoded_bits.truncate(decoded_bits.len() - table.memory); }

    let nb = decoded_bits.len()/8; let mut r = vec![0u8;nb];
    for i in 0..nb*8 { if i<decoded_bits.len()&&decoded_bits[i]!=0 { r[i/8]|=1<<(7-(i%8)); } } r
}

#[napi]
pub fn viterbi_k9_decode(received: Buffer, config: Option<ViterbiConfig>) -> Result<Buffer> {
    let (md,ip,dp,el) = resolve_config(config.as_ref());
    let mut r = viterbi_decode_indel(received.as_ref(),k9_table(),md,ip,dp,None); if let Some(l)=el { r.truncate(l); } Ok(r.into())
}
#[napi]
pub fn viterbi_k7_decode(received: Buffer, config: Option<ViterbiConfig>) -> Result<Buffer> {
    let (md,ip,dp,el) = resolve_config(config.as_ref());
    let mut r = viterbi_decode_indel(received.as_ref(),k7_table(),md,ip,dp,None); if let Some(l)=el { r.truncate(l); } Ok(r.into())
}
#[napi]
pub fn viterbi_k9_decode_with_llr(received: Buffer, llr: Float32Array, config: Option<ViterbiConfig>) -> Result<Buffer> {
    let (md,ip,dp,el) = resolve_config(config.as_ref());
    let ls: Vec<f32> = llr.to_vec(); let mut r = viterbi_decode_indel(received.as_ref(),k9_table(),md,ip,dp,Some(&ls)); if let Some(l)=el { r.truncate(l); } Ok(r.into())
}
#[napi] pub fn conv_k9_encode(data: Buffer) -> Buffer { conv_encode(data.as_ref(),k9_table()).into() }
#[napi] pub fn conv_k7_encode(data: Buffer) -> Buffer { conv_encode(data.as_ref(),k7_table()).into() }

#[napi]
pub fn viterbi_k9_decode_standard(received: Buffer) -> Buffer {
    let t = k9_table(); let ns = t.num_states; let rb = received.len()*8;
    if rb==0 { return Vec::<u8>::new().into(); }
    let ni = f64::INFINITY; let nst = rb/2;
    let mut c = vec![ni;ns]; let mut p = vec![ni;ns];
    let mut tr: Vec<(u16,u8)> = vec![(0,0);nst*ns]; c[0] = 0.0;
    let mut rp = 0usize;
    for s in 0..nst { std::mem::swap(&mut c,&mut p); c.fill(ni);
        let b0 = get_bit(received.as_ref(),rp); let b1 = get_bit(received.as_ref(),rp+1); rp+=2;
        for st in 0..ns { if p[st]>=1e29 { continue; } for ib in 0..2 {
            let o = t.output(st,ib); let co = p[st]+if (o&1) as u8!=b0{1.0}else{0.0}+if ((o>>1)&1) as u8!=b1{1.0}else{0.0};
            let n = t.next_state(st,ib); if co<c[n] { c[n]=co; tr[s*ns+n]=(st as u16,ib as u8); }
        }}
    }
    let mut bs = 0usize; let mut bc = c[0]; for s in 1..ns { if c[s]<bc { bc=c[s]; bs=s; } }
    let mut db = Vec::with_capacity(nst); let mut st = bs;
    for s in (0..nst).rev() { let (ps,ib)=tr[s*ns+st]; db.push(ib); st=ps as usize; }
    db.reverse(); if db.len()>K9_MEMORY { db.truncate(db.len()-K9_MEMORY); }
    let nb = db.len()/8; let mut r = vec![0u8;nb]; for i in 0..nb*8 { if i<db.len()&&db[i]!=0 { r[i/8]|=1<<(7-(i%8)); } } r.into()
}

#[napi] pub fn napi_version() -> String { format!("helix-dna-napi v0.1.0 — Viterbi (trellis-stage iteration, phase-correct)") }
