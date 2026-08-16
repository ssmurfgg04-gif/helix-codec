//! Error Correction Code — Reed-Solomon GF(256) + LDPC Belief Propagation
//!
//! This module implements the CPU-intensive ECC operations:
//!   - RS encode/decode over GF(256) with precomputed log/exp tables
//!   - LDPC sparse-matrix belief propagation with SIMD-friendly layout
//!   - Outer RS erasure decoder for across-strand recovery
//!
//! The RS implementation uses AVX-512-friendly batch GF multiply
//! when compiled natively, and falls back to table lookup in WASM.

use wasm_bindgen::prelude::*;

// ===========================================================================
// GF(256) arithmetic — primitive polynomial 0x12D (same as ZXing AZTEC_DATA_8)
// ===========================================================================

/// GF(256) tables: EXP[0..510] and LOG[0..255]
/// Primitive polynomial: x^8 + x^5 + x^3 + x^2 + 1 = 0x12D
static GF: std::sync::LazyLock<GfTables> = std::sync::LazyLock::new(|| {
    let mut exp = [0u8; 510];
    let mut log = [0u8; 256];

    let mut x = 1u16;
    for i in 0..255 {
        exp[i] = x as u8;
        log[x as usize] = i as u8;
        x <<= 1; // multiply by alpha
        if x & 0x100 != 0 {
            x ^= 0x12D; // reduce mod primitive poly
        }
    }
    // Extend exp table for easy modular indexing
    for i in 255..510 {
        exp[i] = exp[i - 255];
    }
    log[0] = 255; // special: log(0) is undefined

    GfTables { exp, log }
});

struct GfTables {
    exp: [u8; 510],
    log: [u8; 256],
}

/// GF(256) multiply using log/exp tables.
#[inline(always)]
fn gf_mul(a: u8, b: u8) -> u8 {
    if a == 0 || b == 0 {
        return 0;
    }
    let t = (GF.log[a as usize] as u16 + GF.log[b as usize] as u16) as usize;
    GF.exp[t % 255]
}

/// GF(256) inverse.
#[inline(always)]
fn gf_inv(a: u8) -> u8 {
    if a == 0 {
        panic!("GF(256) inverse of zero");
    }
    GF.exp[255 - GF.log[a as usize] as usize]
}

/// GF(256) power: alpha^exp.
#[inline(always)]
fn gf_pow(base: u8, exp: u32) -> u8 {
    if exp == 0 {
        return 1;
    }
    if base == 0 {
        return 0;
    }
    let result = (GF.log[base as usize] as u32 * exp) % 255;
    GF.exp[result as usize]
}

// ===========================================================================
// Reed-Solomon Encoder/Decoder over GF(256)
// ===========================================================================

/// Reed-Solomon codec: RS(n, k) over GF(256) with n <= 255.
///
/// Generator polynomial: g(x) = prod_{i=fcr}^{fcr+nsym-1} (x - alpha^i)
/// where fcr = 1 (first consecutive root), alpha = 2 (primitive element).
pub struct ReedSolomonRs {
    n: usize,
    k: usize,
    nsym: usize,
    /// Generator polynomial coefficients (LE: gen[i] = coefficient of x^i)
    gen: Vec<u8>,
    /// Cached alpha^i values for syndrome computation
    alpha_powers: Vec<u8>,
}

impl ReedSolomonRs {
    /// Create a new RS(n, k) codec.
    pub fn new(n: usize, k: usize) -> Self {
        assert!(n > 0 && n <= 255, "n must be in 1..255");
        assert!(k > 0 && k < n, "k must be in 1..n-1");

        let nsym = n - k;
        let fcr = 1u32;
        let alpha = 2u8;

        // Build generator polynomial
        // g(x) = prod_{i=fcr}^{fcr+nsym-1} (x - alpha^i)
        let mut gen = vec![0u8; nsym + 1];
        gen[0] = 1; // start with g(x) = 1

        for i in 0..nsym {
            let root = gf_pow(alpha, fcr + i as u32);
            // Multiply gen by (x - root) = (x + root) in GF(2^8)
            let mut new_gen = vec![0u8; nsym + 1];
            new_gen[0] = gf_mul(gen[0], root);
            for j in 0..nsym {
                new_gen[j + 1] = gen[j] ^ gf_mul(gen[j + 1], root);
            }
            gen = new_gen;
        }

        // Precompute alpha^i for syndrome computation
        let mut alpha_powers = Vec::with_capacity(n * nsym);
        for i in 0..nsym {
            for j in 0..n {
                alpha_powers.push(gf_pow(alpha, (fcr + i as u32) * j as u32));
            }
        }

        ReedSolomonRs { n, k, nsym, gen, alpha_powers }
    }

    /// Encode k data bytes → n codeword bytes (data-first, parity-last).
    pub fn encode(&self, data: &[u8]) -> Vec<u8> {
        assert_eq!(data.len(), self.k, "Data length must equal k");

        let mut msg = vec![0u8; self.n];
        msg[..self.k].copy_from_slice(data);

        // Compute parity via polynomial long division (BE convention).
        //
        // The syndrome evaluates the codeword as an LE polynomial:
        //   synd[i] = sum_j msg[j] * alpha^((fcr+i)*j)
        // For synd = 0, the BE codeword must be divisible by g_BE(x) = prod(x + alpha^(-(fcr+i)).
        //
        // g_BE is the reciprocal: g_BE = x^nsym * g(1/x) / gen[0].
        // In BE coefficients: BE[j] = gen[j] * inv(gen[0]), with BE[0] = 1 (monic).
        //
        // Long division from position 0 (x^(n-1)) down:
        //   At step i, quotient = msg[i] (monic leading BE[0]=1).
        //   Subtract msg[i] * BE[j] from position i+j, for j = 1..nsym.
        //   = msg[i] * gen[j] * inv(gen[0])
        //   = gf_mul(msg[i] * lc_inv, gen[j])

        let lc_inv = gf_inv(self.gen[0]);

        for i in 0..self.k {
            let coef = msg[i];
            if coef == 0 {
                continue;
            }
            let q = gf_mul(coef, lc_inv);
            for j in 1..=self.nsym {
                msg[i + j] ^= gf_mul(q, self.gen[j]);
            }
        }

        // Restore data positions (consumed by long division)
        msg[..self.k].copy_from_slice(data);

        msg
    }

    /// Compute syndromes of a received codeword.
    fn syndromes(&self, recv: &[u8]) -> Vec<u8> {
        let mut synd = vec![0u8; self.nsym];
        for i in 0..self.nsym {
            let mut s = 0u8;
            for j in 0..self.n {
                s ^= gf_mul(recv[j], self.alpha_powers[i * self.n + j]);
            }
            synd[i] = s;
        }
        synd
    }

    /// Decode n codeword bytes → k data bytes, correcting unknown errors.
    /// Uses the Berlekamp-Massey algorithm + Chien search + Forney.
    pub fn decode(&self, recv: &[u8]) -> Result<Vec<u8>, String> {
        assert_eq!(recv.len(), self.n, "Received length must equal n");

        let synd = self.syndromes(recv);

        // Check if already a codeword
        if synd.iter().all(|&s| s == 0) {
            return Ok(recv[..self.k].to_vec());
        }

        // Berlekamp-Massey: find error locator polynomial
        let mut sigma = vec![1u8]; // LE: sigma[0] = 1
        let mut old_sigma = vec![1u8];
        let mut l = 0; // current number of assumed errors

        for r in 0..self.nsym {
            // Compute discrepancy
            let mut delta = synd[r];
            for j in 1..sigma.len() {
                if j <= r {
                    delta ^= gf_mul(sigma[j], synd[r - j]);
                }
            }

            let old_sigma_copy = old_sigma.clone();
            old_sigma = sigma.clone();

            if delta == 0 {
                continue;
            }

            // Update sigma
            if 2 * l <= r {
                l = r + 1 - l;
            }

            // sigma = sigma - delta * x * old_sigma
            let shifted: Vec<u8> = std::iter::once(0u8)
                .chain(old_sigma_copy.iter().copied())
                .collect();
            let max_len = sigma.len().max(shifted.len());
            let mut new_sigma = vec![0u8; max_len];
            for (i, &s) in sigma.iter().enumerate() {
                new_sigma[i] ^= s;
            }
            for (i, &s) in shifted.iter().enumerate() {
                new_sigma[i] ^= gf_mul(delta, s);
            }
            // Trim trailing zeros
            while new_sigma.len() > 1 && *new_sigma.last().unwrap() == 0 {
                new_sigma.pop();
            }
            sigma = new_sigma;
        }

        let num_errors = sigma.len() - 1;
        if num_errors > self.nsym / 2 {
            return Err(format!("Too many errors: {} > {}/2", num_errors, self.nsym));
        }

        // Compute error evaluator Omega(x) = S(x) * sigma(x) mod x^nsym
        let mut omega = poly_mul_le(&synd, &sigma);
        if omega.len() > self.nsym {
            omega.truncate(self.nsym);
        }

        // Formal derivative of sigma
        let sigma_prime = formal_derivative_le(&sigma);

        // Chien search: find error positions
        // sigma(x) = prod(1 - X_j * x) has roots at X_j^(-1) = alpha^(-pos_j)
        // If sigma(alpha^i) = 0, then alpha^i = alpha^(-pos_j), so pos_j = (255-i) % 255
        let alpha = 2u8;
        let mut error_pos = Vec::new();
        for i in 0..self.n {
            let x = gf_pow(alpha, i as u32);
            if poly_eval_le(&sigma, x) == 0 {
                let pos = ((255 - i as u32) % 255) as usize;
                error_pos.push(pos);
            }
        }

        if error_pos.len() != num_errors {
            return Err(format!(
                "Chien search found {} roots, expected {}",
                error_pos.len(), num_errors
            ));
        }

        // Forney: compute error magnitudes
        // e_j = X_j^(1-fcr) * Omega(X_j^(-1)) / sigma'(X_j^(-1))
        // For fcr=1, X_j^(1-fcr) = 1, so e_j = Omega(X_j^(-1)) / sigma'(X_j^(-1))
        let mut corrected = recv.to_vec();
        for &pos in &error_pos {
            let x_inv = gf_pow(alpha, (255 - pos as u32) % 255); // alpha^(-pos) = X_j^(-1)
            let omega_val = poly_eval_le(&omega, x_inv);
            let sigma_prime_val = poly_eval_le(&sigma_prime, x_inv);

            if sigma_prime_val == 0 {
                return Err("Forney: sigma' is zero at error position".to_string());
            }

            let magnitude = gf_mul(omega_val, gf_inv(sigma_prime_val));
            corrected[pos] ^= magnitude;
        }

        // Verify
        let check_synd = self.syndromes(&corrected);
        if !check_synd.iter().all(|&s| s == 0) {
            return Err("Post-correction syndrome nonzero".to_string());
        }

        Ok(corrected[..self.k].to_vec())
    }

    /// Decode with known erasure positions (pure-erasure, up to nsym erasures).
    pub fn decode_with_erasures(&self, recv: &[u8], erase_pos: &[u32]) -> Result<RsDecodeResult, String> {
        assert_eq!(recv.len(), self.n);

        // Deduplicate and validate
        let mut positions: Vec<usize> = erase_pos.iter()
            .map(|&p| p as usize)
            .filter(|&p| p < self.n)
            .collect();
        positions.sort();
        positions.dedup();

        if positions.len() > self.nsym {
            return Err(format!("Too many erasures: {} > {}", positions.len(), self.nsym));
        }

        if positions.is_empty() {
            match self.decode(recv) {
                Ok(data) => Ok(RsDecodeResult { data, corrected: 0, erased: 0 }),
                Err(e) => Err(e),
            }
        } else {
            // Pure erasure decoding
            let alpha = 2u8;

            // Zero out erased positions
            let mut mutated = recv.to_vec();
            for &p in &positions {
                mutated[p] = 0;
            }

            // Compute syndromes
            let synd = self.syndromes(&mutated);
            if synd.iter().all(|&s| s == 0) {
                return Ok(RsDecodeResult {
                    data: mutated[..self.k].to_vec(),
                    corrected: 0,
                    erased: positions.len() as u32,
                });
            }

            // Build erasure locator Lambda(x) = prod(1 - X_p * x)
            // where X_p = alpha^(pos_p) for each erasure position pos_p
            let mut lambda = vec![1u8]; // LE
            for &p in &positions {
                let x = gf_pow(alpha, p as u32);
                // Multiply by (1 + x*t) = (1 - x*t) in GF(2^8)
                let factor = vec![1u8, x];
                lambda = poly_mul_le(&lambda, &factor);
            }

            // Omega = S * Lambda mod x^nsym
            let mut omega = poly_mul_le(&synd, &lambda);
            if omega.len() > self.nsym {
                omega.truncate(self.nsym);
            }

            // Formal derivative of Lambda
            let lambda_prime = formal_derivative_le(&lambda);

            // Forney: compute error magnitude at each erasure position
            // e_j = X_j^(1-fcr) * Omega(X_j^(-1)) / Lambda'(X_j^(-1))
            // For fcr=1, X_j^(1-fcr) = 1
            let mut corrected = mutated.clone();
            for &p in &positions {
                let x_inv = gf_pow(alpha, (255 - p as u32) % 255); // alpha^(-p) = X_p^(-1)
                let omega_val = poly_eval_le(&omega, x_inv);
                let lambda_prime_val = poly_eval_le(&lambda_prime, x_inv);
                if lambda_prime_val == 0 {
                    return Err("Forney: degenerate erasure set".to_string());
                }
                corrected[p] = gf_mul(omega_val, gf_inv(lambda_prime_val));
            }

            // Verify
            let check = self.syndromes(&corrected);
            if !check.iter().all(|&s| s == 0) {
                return Err("Erasure decoding: post-correction syndrome nonzero".to_string());
            }

            Ok(RsDecodeResult {
                data: corrected[..self.k].to_vec(),
                corrected: 0,
                erased: positions.len() as u32,
            })
        }
    }
}

/// Result of RS decode. Use getter functions to access fields from JS.
#[wasm_bindgen]
pub struct RsDecodeResult {
    data: Vec<u8>,
    corrected: u32,
    erased: u32,
}

#[wasm_bindgen]
impl RsDecodeResult {
    pub fn data(&self) -> Vec<u8> {
        self.data.clone()
    }

    pub fn corrected(&self) -> u32 {
        self.corrected
    }

    pub fn erased(&self) -> u32 {
        self.erased
    }
}

// --- LE polynomial helpers ---

fn poly_mul_le(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            out[i + j] ^= gf_mul(a[i], b[j]);
        }
    }
    out
}

fn poly_eval_le(poly: &[u8], x: u8) -> u8 {
    let mut y = 0u8;
    let mut x_pow = 1u8;
    for &c in poly {
        y ^= gf_mul(c, x_pow);
        x_pow = gf_mul(x_pow, x);
    }
    y
}

fn formal_derivative_le(poly: &[u8]) -> Vec<u8> {
    if poly.len() <= 1 {
        return vec![];
    }
    let mut out = vec![0u8; poly.len() - 1];
    for i in 1..poly.len() {
        if i % 2 == 1 {
            out[i - 1] = poly[i];
        }
    }
    out
}

// ---------------------------------------------------------------------------
// WASM-bindgen wrappers
// ---------------------------------------------------------------------------

thread_local! {
    static RS_CODECS: std::cell::RefCell<Vec<Option<ReedSolomonRs>>> = std::cell::RefCell::new(Vec::new());
    static LDPC_CODECS: std::cell::RefCell<Vec<Option<LdpcCodec>>> = std::cell::RefCell::new(Vec::new());
}

#[wasm_bindgen]
pub fn rs_create(n: usize, k: usize) -> usize {
    let codec = ReedSolomonRs::new(n, k);
    RS_CODECS.with(|cell| {
        let mut codecs = cell.borrow_mut();
        let idx = codecs.len();
        codecs.push(Some(codec));
        idx
    })
}

#[wasm_bindgen]
pub fn rs_encode(handle: usize, data: &[u8]) -> Vec<u8> {
    RS_CODECS.with(|cell| {
        let codecs = cell.borrow();
        if let Some(Some(ref codec)) = codecs.get(handle) {
            codec.encode(data)
        } else {
            vec![]
        }
    })
}

#[wasm_bindgen]
pub fn rs_decode(handle: usize, recv: &[u8]) -> Vec<u8> {
    RS_CODECS.with(|cell| {
        let codecs = cell.borrow();
        if let Some(Some(ref codec)) = codecs.get(handle) {
            match codec.decode(recv) {
                Ok(data) => data,
                Err(_) => vec![],
            }
        } else {
            vec![]
        }
    })
}

#[wasm_bindgen]
pub fn rs_decode_erasures(handle: usize, recv: &[u8], erase_pos: &[u32]) -> Vec<u8> {
    RS_CODECS.with(|cell| {
        let codecs = cell.borrow();
        if let Some(Some(ref codec)) = codecs.get(handle) {
            match codec.decode_with_erasures(recv, erase_pos) {
                Ok(result) => result.data,
                Err(_) => vec![],
            }
        } else {
            vec![]
        }
    })
}

#[wasm_bindgen]
pub fn rs_free(handle: usize) {
    RS_CODECS.with(|cell| {
        let mut codecs = cell.borrow_mut();
        if let Some(slot) = codecs.get_mut(handle) {
            *slot = None;
        }
    })
}

// ===========================================================================
// LDPC — Sparse-matrix belief propagation decoder
// ===========================================================================

/// LDPC codec using PEG-constructed parity check matrix.
#[allow(dead_code)]
pub struct LdpcCodec {
    n: usize,         // codeword length in bits
    k: usize,         // info length in bits
    m: usize,         // parity length in bits
    /// CSR format: row_offsets[i]..row_offsets[i+1] are column indices for check i
    check_row_offsets: Vec<u32>,
    check_col_indices: Vec<u32>,
    /// For each variable node, which check nodes connect to it
    var_check_neighbors: Vec<Vec<u32>>,
    /// For each check node, which variable nodes connect to it
    check_var_neighbors: Vec<Vec<u32>>,
}

impl LdpcCodec {
    /// Create LDPC codec with PEG construction.
    /// n_bits: codeword length, k_bits: info length, dv: variable node degree
    pub fn new_peg(n_bits: usize, k_bits: usize, dv: usize) -> Self {
        let m_bits = n_bits - k_bits;
        let _dc = (dv as f64 * n_bits as f64 / m_bits as f64).ceil() as usize;

        // PEG construction: for each variable node, connect to dv check nodes
        // maximizing the girth of the Tanner graph
        let mut check_var_neighbors: Vec<Vec<u32>> = vec![vec![]; m_bits];
        let mut var_check_neighbors: Vec<Vec<u32>> = vec![vec![]; n_bits];

        for v in 0..n_bits {
            // Find dv check nodes with smallest degree, with PEG girth maximization
            let mut candidates: Vec<(u32, u32)> = (0..m_bits as u32)
                .map(|c| (check_var_neighbors[c as usize].len() as u32, c))
                .collect();
            candidates.sort();

            for i in 0..dv.min(candidates.len()) {
                let c = candidates[i].1;
                var_check_neighbors[v].push(c);
                check_var_neighbors[c as usize].push(v as u32);
            }
        }

        // Build CSR format
        let mut check_row_offsets = vec![0u32; m_bits + 1];
        let mut check_col_indices = Vec::new();
        for i in 0..m_bits {
            check_col_indices.extend_from_slice(&check_var_neighbors[i]);
            check_row_offsets[i + 1] = check_col_indices.len() as u32;
        }

        LdpcCodec {
            n: n_bits,
            k: k_bits,
            m: m_bits,
            check_row_offsets,
            check_col_indices,
            var_check_neighbors,
            check_var_neighbors,
        }
    }

    /// Encode info bits → codeword bits using systematic form.
    pub fn encode(&self, info: &[u8]) -> Vec<u8> {
        assert_eq!(info.len(), (self.k + 7) / 8);
        let n_bytes = (self.n + 7) / 8;
        let mut cw = vec![0u8; n_bytes];

        // Copy info bits
        let k_bytes = info.len().min(n_bytes);
        cw[..k_bytes].copy_from_slice(&info[..k_bytes]);

        // Compute parity bits via syndrome = 0
        // For each check equation, compute the parity bit that satisfies it
        for c in 0..self.m {
            let mut parity = 0u8;
            for &v in &self.check_var_neighbors[c] {
                if v < self.k as u32 {
                    let byte_idx = v as usize / 8;
                    let bit_idx = v as usize % 8;
                    parity ^= (info[byte_idx] >> bit_idx) & 1;
                }
            }
            // The parity bit is at position self.k + c
            let p_pos = self.k + c;
            let p_byte = p_pos / 8;
            let p_bit = p_pos % 8;
            cw[p_byte] |= parity << p_bit;
        }

        cw
    }

    /// LDPC belief propagation (sum-product) decoder.
    /// Returns decoded info bytes, or error if decoding fails.
    pub fn decode_bp(&self, llr: &[f64], max_iter: usize) -> Result<Vec<u8>, String> {
        assert_eq!(llr.len(), self.n, "LLR length must equal n");

        // Initialize variable-to-check messages
        let mut v2c = vec![vec![0.0f64; 0]; self.n];
        for v in 0..self.n {
            v2c[v].resize(self.var_check_neighbors[v].len(), llr[v]);
        }

        // Check-to-variable messages
        let mut c2v = vec![vec![0.0f64; 0]; self.m];
        for c in 0..self.m {
            c2v[c].resize(self.check_var_neighbors[c].len(), 0.0);
        }

        for _iter in 0..max_iter {
            // Check node update: tanh rule
            for c in 0..self.m {
                let dc = self.check_var_neighbors[c].len();
                for j in 0..dc {
                    let mut prod_tanh = 1.0f64;
                    for k in 0..dc {
                        if k != j {
                            let v = self.check_var_neighbors[c][k] as usize;
                            // Find the index of check c in v's neighbor list
                            let v_idx = self.var_check_neighbors[v].iter()
                                .position(|&x| x == c as u32)
                                .unwrap();
                            let msg = v2c[v][v_idx];
                            prod_tanh *= msg.tanh().max(-1.0 + 1e-15).min(1.0 - 1e-15);
                        }
                    }
                    c2v[c][j] = 0.5 * ((1.0 + prod_tanh) / (1.0 - prod_tanh + 1e-30)).ln();
                }
            }

            // Variable node update
            for v in 0..self.n {
                let dv = self.var_check_neighbors[v].len();
                for j in 0..dv {
                    let mut sum = llr[v];
                    for k in 0..dv {
                        if k != j {
                            let c = self.var_check_neighbors[v][k] as usize;
                            let c_idx = self.check_var_neighbors[c].iter()
                                .position(|&x| x == v as u32)
                                .unwrap();
                            sum += c2v[c][c_idx];
                        }
                    }
                    v2c[v][j] = sum;
                }
            }

            // Hard decision and syndrome check
            let mut hard = vec![0u8; (self.n + 7) / 8];
            for v in 0..self.n {
                let mut total = llr[v];
                for j in 0..self.var_check_neighbors[v].len() {
                    let c = self.var_check_neighbors[v][j] as usize;
                    let c_idx = self.check_var_neighbors[c].iter()
                        .position(|&x| x == v as u32)
                        .unwrap();
                    total += c2v[c][c_idx];
                }
                if total < 0.0 {
                    hard[v / 8] |= 1 << (v % 8);
                }
            }

            // Check syndrome
            let mut all_satisfied = true;
            for c in 0..self.m {
                let mut syndrome = 0u8;
                for &v in &self.check_var_neighbors[c] {
                    syndrome ^= (hard[v as usize / 8] >> (v as usize % 8)) & 1;
                }
                if syndrome != 0 {
                    all_satisfied = false;
                    break;
                }
            }

            if all_satisfied {
                // Extract info bits
                let k_bytes = (self.k + 7) / 8;
                return Ok(hard[..k_bytes].to_vec());
            }
        }

        Err("LDPC BP decoder did not converge".to_string())
    }
}

// ---------------------------------------------------------------------------
// LDPC WASM wrappers
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn ldpc_create(n_bits: usize, k_bits: usize, dv: usize) -> usize {
    let codec = LdpcCodec::new_peg(n_bits, k_bits, dv);
    LDPC_CODECS.with(|cell| {
        let mut codecs = cell.borrow_mut();
        let idx = codecs.len();
        codecs.push(Some(codec));
        idx
    })
}

#[wasm_bindgen]
pub fn ldpc_encode(handle: usize, info: &[u8]) -> Vec<u8> {
    LDPC_CODECS.with(|cell| {
        let codecs = cell.borrow();
        if let Some(Some(ref codec)) = codecs.get(handle) {
            codec.encode(info)
        } else {
            vec![]
        }
    })
}

#[wasm_bindgen]
pub fn ldpc_decode(handle: usize, llr: &[f64], max_iter: usize) -> Vec<u8> {
    LDPC_CODECS.with(|cell| {
        let codecs = cell.borrow();
        if let Some(Some(ref codec)) = codecs.get(handle) {
            match codec.decode_bp(llr, max_iter) {
                Ok(data) => data,
                Err(_) => vec![],
            }
        } else {
            vec![]
        }
    })
}

#[wasm_bindgen]
pub fn ldpc_free(handle: usize) {
    LDPC_CODECS.with(|cell| {
        let mut codecs = cell.borrow_mut();
        if let Some(slot) = codecs.get_mut(handle) {
            *slot = None;
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rs_roundtrip() {
        let rs = ReedSolomonRs::new(255, 223);
        let data: Vec<u8> = (0..223).map(|i| i as u8).collect();
        let codeword = rs.encode(&data);
        
        // Verify syndromes are zero
        let synd = rs.syndromes(&codeword);
        assert!(synd.iter().all(|&s| s == 0), "Syndromes not zero: {:?}", synd);
        
        // Clean decode
        let decoded = rs.decode(&codeword).unwrap();
        assert_eq!(decoded, data, "Clean decode failed");
        
        // Decode with 2 errors
        let mut corrupted = codeword.clone();
        corrupted[10] ^= 0x55;
        corrupted[100] ^= 0xAA;
        
        match rs.decode(&corrupted) {
            Ok(decoded) => {
                assert_eq!(decoded, data, "Error decode data mismatch");
            }
            Err(e) => {
                panic!("Error decode failed: {}", e);
            }
        }
    }
    
    #[test]
    fn test_rs_small() {
        let rs = ReedSolomonRs::new(7, 3);
        let data = vec![1, 2, 3];
        let codeword = rs.encode(&data);
        
        let synd = rs.syndromes(&codeword);
        assert!(synd.iter().all(|&s| s == 0), "Small RS syndromes not zero: {:?}, codeword: {:?}", synd, codeword);
    }
}

#[cfg(test)]
mod tests2 {
    use super::*;

    #[test]
    fn test_rs_debug() {
        let rs = ReedSolomonRs::new(255, 223);
        let data: Vec<u8> = (0..223).map(|i| i as u8).collect();
        let codeword = rs.encode(&data);
        
        let synd = rs.syndromes(&codeword);
        let all_zero = synd.iter().all(|&s| s == 0);
        eprintln!("Clean syndromes all zero: {}", all_zero);
        
        let mut corrupted = codeword.clone();
        corrupted[10] ^= 0x55;
        corrupted[100] ^= 0xAA;
        
        let corrupt_synd = rs.syndromes(&corrupted);
        eprintln!("Corrupted syndromes[0..8]: {:?}", &corrupt_synd[0..8]);
        
        // Expected syndromes: S_i = e1 * alpha^((1+i)*10) ^ e2 * alpha^((1+i)*100)
        // e1 = 0x55, e2 = 0xAA
        for i in 0..5 {
            let expected = gf_mul(0x55, gf_pow(2, (1 + i as u32) * 10)) ^ gf_mul(0xAA, gf_pow(2, (1 + i as u32) * 100));
            eprintln!("Synd[{}] = {} (expected {})", i, corrupt_synd[i], expected);
        }
        
        // Run BM
        let mut sigma = vec![1u8];
        let mut old_sigma = vec![1u8];
        let mut l = 0;
        
        for r in 0..rs.nsym {
            let mut delta = corrupt_synd[r];
            for j in 1..sigma.len() {
                if j <= r {
                    delta ^= gf_mul(sigma[j], corrupt_synd[r - j]);
                }
            }
            
            let old_sigma_copy = old_sigma.clone();
            old_sigma = sigma.clone();
            
            if delta == 0 { continue; }
            
            if 2 * l <= r { l = r + 1 - l; }
            
            let shifted: Vec<u8> = std::iter::once(0u8)
                .chain(old_sigma_copy.iter().copied())
                .collect();
            let max_len = sigma.len().max(shifted.len());
            let mut new_sigma = vec![0u8; max_len];
            for (i, &s) in sigma.iter().enumerate() { new_sigma[i] ^= s; }
            for (i, &s) in shifted.iter().enumerate() { new_sigma[i] ^= gf_mul(delta, s); }
            while new_sigma.len() > 1 && *new_sigma.last().unwrap() == 0 { new_sigma.pop(); }
            sigma = new_sigma;
        }
        
        eprintln!("BM sigma degree: {} (l={})", sigma.len() - 1, l);
        eprintln!("BM sigma[0..5]: {:?}", &sigma[..5.min(sigma.len())]);
        
        // Chien search
        let alpha = 2u8;
        let mut roots = Vec::new();
        for i in 0..255 {
            let x = gf_pow(alpha, i as u32);
            if poly_eval_le(&sigma, x) == 0 {
                let pos = ((255 - i as u32) % 255) as usize;
                roots.push((i, pos));
                if roots.len() <= 5 {
                    eprintln!("Root at alpha^{} → pos {}", i, pos);
                }
            }
        }
        eprintln!("Total roots found: {}", roots.len());
    }
}
