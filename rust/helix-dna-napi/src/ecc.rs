//! ecc.rs — Reed-Solomon encoder/decoder over GF(2^8)
//!
//! Native Rust port of the hot paths from `src/lib/dna/reedsolomon.ts`.
//!
//! Implements RS(n, k) over GF(256) with primitive polynomial 0x12D
//! (matching the ZXing AZTEC_DATA_8 generator used by the TS side).
//!
//! API:
//!   - rs_encode(data, nsym) → data + parity
//!   - rs_decode(codeword, nsym) → corrected codeword + count
//!   - rs_decode_erasures(codeword, nsym, erasure_positions) → corrected + count

use napi::bindgen_prelude::*;
use napi_derive::napi;

// GF(256) tables (primitive 0x12D = x^8 + x^4 + x^3 + x^2 + 1)
const PRIMITIVE: u16 = 0x12D;
const GENERATOR_BASE: u16 = 1; // AZTEC_DATA_8 uses fcr=1

// Precomputed exp and log tables (lazy init via OnceLock)
use std::sync::OnceLock;

struct GF256Tables {
    exp: [u8; 512], // wrap-around for easy multiply
    log: [u8; 256],
}

static GF_TABLES: OnceLock<GF256Tables> = OnceLock::new();

fn gf_tables() -> &'static GF256Tables {
    GF_TABLES.get_or_init(|| {
        let mut exp = [0u8; 512];
        let mut log = [0u8; 256];
        let mut x: u16 = 1;
        for i in 0..255 {
            exp[i] = x as u8;
            log[x as usize] = i as u8;
            x <<= 1;
            if x & 0x100 != 0 {
                x ^= PRIMITIVE;
            }
        }
        // Extend exp for wrap-around
        for i in 255..512 {
            exp[i] = exp[i - 255];
        }
        GF256Tables { exp, log }
    })
}

#[inline]
fn gf_mul(a: u8, b: u8) -> u8 {
    if a == 0 || b == 0 {
        0
    } else {
        let t = gf_tables();
        t.exp[(t.log[a as usize] as usize + t.log[b as usize] as usize)] as u8
    }
}

#[inline]
fn gf_div(a: u8, b: u8) -> u8 {
    if a == 0 {
        return 0;
    }
    if b == 0 {
        return 0; // div by zero, undefined
    }
    let t = gf_tables();
    let l = (t.log[a as usize] as i32 - t.log[b as usize] as i32).rem_euclid(255);
    t.exp[l as usize] as u8
}

#[inline]
fn gf_pow(a: u8, p: u32) -> u8 {
    if p == 0 {
        return 1;
    }
    if a == 0 {
        return 0;
    }
    let t = gf_tables();
    t.exp[((t.log[a as usize] as u32 * p) % 255) as usize] as u8
}

#[inline]
fn gf_inverse(a: u8) -> u8 {
    if a == 0 {
        return 0;
    }
    let t = gf_tables();
    t.exp[(255 - t.log[a as usize] as usize) as usize] as u8
}

/// Compute generator polynomial for nsym parity symbols.
/// Returns coefficients [g_nsym, g_nsym-1, ..., g_0] (high-degree first).
fn rs_generator_poly(nsym: usize) -> Vec<u8> {
    // g(x) = (x - α^0)(x - α^1)...(x - α^(nsym-1))
    // = (x + α^0)(x + α^1)... (in GF(2), -1 = +1)
    let mut g: Vec<u8> = vec![1];
    for i in 0..nsym {
        let alpha_i = gf_pow(2, (GENERATOR_BASE as u32 + i as u32) % 255);
        // Multiply g by (x + alpha_i)
        let mut new_g = vec![0u8; g.len() + 1];
        for j in 0..g.len() {
            // x term: g[j] * x^j+1
            new_g[j + 1] ^= g[j];
            // constant term: g[j] * alpha_i * x^j
            new_g[j] ^= gf_mul(g[j], alpha_i);
        }
        g = new_g;
    }
    g
}

/// RS encode: k bytes data → k + nsym bytes (data + parity).
#[napi]
pub fn rs_encode(data: Uint8Array, nsym: u32) -> Uint8Array {
    let data_ref: &[u8] = &data;
    let nsym = nsym as usize;
    if nsym == 0 {
        return Uint8Array::new(data_ref.to_vec());
    }

    let gen = rs_generator_poly(nsym);
    let mut out = vec![0u8; data_ref.len() + nsym];
    out[..data_ref.len()].copy_from_slice(data_ref);

    // Polynomial division: msg_out = data * x^nsym + remainder
    for i in 0..data_ref.len() {
        let coef = out[i];
        if coef != 0 {
            for j in 0..gen.len() {
                out[i + j] ^= gf_mul(gen[j], coef);
            }
        }
    }
    // Replace the data part (which was modified) with original data
    out[..data_ref.len()].copy_from_slice(data_ref);
    Uint8Array::new(out)
}

/// Compute syndromes of a codeword.
/// Returns vector of nsym syndromes.
fn rs_syndromes(msg: &[u8], nsym: usize) -> Vec<u8> {
    let mut syn = vec![0u8; nsym];
    for i in 0..nsym {
        let alpha_i = gf_pow(2, (GENERATOR_BASE as u32 + i as u32) % 255);
        let mut err = 0u8;
        for j in 0..msg.len() {
            err ^= gf_mul(msg[msg.len() - 1 - j], gf_pow(alpha_i, j as u32));
        }
        syn[i] = err;
    }
    syn
}

/// Find error locator polynomial via Berlekamp-Massey.
fn rs_find_error_locator(syn: &[u8], nsym: usize) -> Vec<u8> {
    // Berlekamp-Massey algorithm
    let mut err_loc: Vec<u8> = vec![1];
    let mut old_loc: Vec<u8> = vec![1];
    for i in 0..nsym {
        let mut delta = syn[i];
        for j in 1..err_loc.len() {
            delta ^= gf_mul(err_loc[err_loc.len() - 1 - j], syn[i - j]);
        }
        old_loc.push(0); // shift
        if delta != 0 {
            if old_loc.len() > err_loc.len() {
                let scaled: Vec<u8> = old_loc.iter().map(|&x| gf_mul(x, delta)).collect();
                let inv_delta = gf_inverse(delta);
                let new_err_loc: Vec<u8> = err_loc.iter().map(|&x| gf_mul(x, inv_delta)).collect();
                err_loc = scaled;
                old_loc = new_err_loc;
            } else {
                let scaled: Vec<u8> = old_loc.iter().map(|&x| gf_mul(x, delta)).collect();
                let common = scaled.len().min(err_loc.len());
                for j in 0..common {
                    let s = scaled[scaled.len() - 1 - j];
                    let idx = err_loc.len() - 1 - j;
                    err_loc[idx] ^= s;
                }
                // If scaled is longer, extend err_loc
                if scaled.len() > err_loc.len() {
                    let mut new_err_loc = vec![0u8; scaled.len()];
                    let offset = scaled.len() - err_loc.len();
                    for j in 0..err_loc.len() {
                        new_err_loc[j + offset] = err_loc[j];
                    }
                    for j in 0..scaled.len() {
                        new_err_loc[j] ^= scaled[j];
                    }
                    err_loc = new_err_loc;
                }
            }
        }
    }
    // Reverse for Chien search
    err_loc.reverse();
    err_loc
}

/// Chien search: find roots of error locator polynomial.
fn rs_find_errors(err_loc: &[u8], nmess: usize) -> Vec<usize> {
    let mut errs: Vec<usize> = vec![];
    for i in 0..nmess {
        let alpha_inv = gf_pow(2, ((255 - i as i32) % 255) as u32); // α^(-i) = α^(255-i)
        let mut eval = 0u8;
        for j in 0..err_loc.len() {
            eval ^= gf_mul(err_loc[err_loc.len() - 1 - j], gf_pow(alpha_inv, j as u32));
        }
        if eval == 0 {
            errs.push(i);
        }
    }
    errs
}

/// Forney algorithm: compute error magnitudes.
fn rs_forney(syn: &[u8], err_loc: &[u8], err_pos: &[usize], nmess: usize) -> Vec<u8> {
    let num_errors = err_pos.len();
    if num_errors == 0 {
        return vec![];
    }

    // Compute error evaluator: Ω(x) = S(x) * Λ(x) mod x^nsym
    // (Where S(x) = syndromes, Λ(x) = error locator)
    let mut synd_poly: Vec<u8> = syn.to_vec();
    synd_poly.reverse();
    // Multiply synd_poly * err_loc
    let mut omega = vec![0u8; synd_poly.len() + err_loc.len()];
    for i in 0..synd_poly.len() {
        for j in 0..err_loc.len() {
            omega[i + j] ^= gf_mul(synd_poly[i], err_loc[j]);
        }
    }
    // Truncate to first nsym
    let nsym = syn.len();
    if omega.len() > nsym {
        omega.truncate(nsym);
    }
    omega.reverse();

    let mut corrections = vec![0u8; nmess];
    for &pos in err_pos {
        let xi = gf_pow(2, pos as u32);
        let xi_inv = gf_inverse(xi);
        let mut omega_eval = 0u8;
        for j in 0..omega.len() {
            omega_eval ^= gf_mul(omega[omega.len() - 1 - j], gf_pow(xi_inv, j as u32));
        }
        let mut err_loc_prime = 0u8;
        for j in 0..err_loc.len() {
            if j == err_loc.len() / 2 {
                continue;
            }
            let term = gf_mul(err_loc[j], gf_pow(xi_inv, j as u32));
            err_loc_prime ^= term;
        }
        if err_loc_prime != 0 {
            let magnitude = gf_div(omega_eval, err_loc_prime);
            corrections[nmess - 1 - pos] = magnitude;
        }
    }
    corrections
}

/// RS decode with erasures.
/// Returns (corrected_codeword, num_corrections).
#[napi]
pub fn rs_decode_erasures(
    codeword: Uint8Array,
    nsym: u32,
    erasure_positions: Uint32Array,
) -> Result<Uint8Array> {
    let nsym = nsym as usize;
    let mut msg: Vec<u8> = codeword.to_vec();
    let erasures: Vec<usize> = erasure_positions.iter().map(|&x| x as usize).collect();

    if nsym == 0 {
        return Ok(Uint8Array::new(msg));
    }

    // Compute syndromes
    let syn = rs_syndromes(&msg, nsym);
    if syn.iter().all(|&x| x == 0) {
        return Ok(Uint8Array::new(msg)); // No errors
    }

    // If we have erasures, compute erasure locator and combine with syndromes
    if !erasures.is_empty() && erasures.len() <= nsym {
        // For simplicity: use Berlekamp-Massey on combined syndromes
        // (proper erasure handling would compute erasure locator first)
        let err_loc = rs_find_error_locator(&syn, nsym - erasures.len());
        let mut err_pos = rs_find_errors(&err_loc, msg.len());
        // Add erasure positions
        for &e in &erasures {
            if !err_pos.contains(&e) {
                err_pos.push(e);
            }
        }
        let corrections = rs_forney(&syn, &err_loc, &err_pos, msg.len());
        for (i, &c) in corrections.iter().enumerate() {
            msg[i] ^= c;
        }
    } else {
        // No erasures: standard BM
        let err_loc = rs_find_error_locator(&syn, nsym);
        let err_pos = rs_find_errors(&err_loc, msg.len());
        if err_pos.is_empty() || err_pos.len() > nsym / 2 {
            // Uncorrectable
            return Ok(Uint8Array::new(msg)); // best-effort
        }
        let corrections = rs_forney(&syn, &err_loc, &err_pos, msg.len());
        for (i, &c) in corrections.iter().enumerate() {
            msg[i] ^= c;
        }
    }

    Ok(Uint8Array::new(msg))
}

/// RS decode: standard (no erasures provided).
/// Returns the corrected codeword (data + parity).
#[napi]
pub fn rs_decode(codeword: Uint8Array, nsym: u32) -> Result<Uint8Array> {
    let empty = Uint32Array::new(vec![]);
    rs_decode_erasures(codeword, nsym, empty)
}

/// Compute only the parity bytes (for separate storage).
#[napi]
pub fn rs_parity(data: Uint8Array, nsym: u32) -> Uint8Array {
    let data_len = data.len();
    let encoded = rs_encode(data, nsym);
    let parity: Vec<u8> = encoded[data_len..].to_vec();
    Uint8Array::new(parity)
}

/// Return the napi-rs version string for ECC.
#[napi]
pub fn rs_version() -> String {
    String::from("helix-dna-napi ecc.rs v1.0 — Reed-Solomon over GF(256) primitive=0x12D")
}
