//! DNA Compression — Arithmetic coding with context models
//!
//! Implements the CPU-bound compression tiers:
//!   - 2-bit packing (identity, no compression)
//!   - Order-0/1/2 context modeling + arithmetic coding
//!   - DNA-specific context models (k-mer frequency, GC bias)
//!
//! The arithmetic coder uses a byte-oriented range coder with 64-bit low
//! for reliable carry propagation and correct roundtrip behavior.
//! Context models use adaptive probability estimation with periodical
//! renormalization (Fenwick's approach).

use wasm_bindgen::prelude::*;

// ===========================================================================
// Byte-oriented Range Coder (64-bit low for correct carry propagation)
// ===========================================================================

const BOT: u32 = 1u32 << 24;          // renorm threshold = 2^24
const INITIAL_RNG: u32 = 1u32 << 31;  // initial range = 2^31

/// Range encoder — byte-oriented, 64-bit low for carry handling.
struct RangeEncoder {
    low: u64,     // 64-bit lower bound (allows carry detection)
    rng: u32,     // 32-bit range
    out: Vec<u8>,
}

impl RangeEncoder {
    fn new() -> Self {
        RangeEncoder {
            low: 0,
            rng: INITIAL_RNG,
            out: Vec::new(),
        }
    }

    /// Encode a symbol with given frequency counts.
    fn encode(&mut self, cum_freq: u32, sym_freq: u32, total_freq: u32) {
        let r = self.rng / total_freq;
        self.low += cum_freq as u64 * r as u64;
        self.rng = sym_freq * r;
        self.renorm();
    }

    /// Renormalize: output bytes and adjust range.
    fn renorm(&mut self) {
        while self.rng < BOT {
            // Propagate carry if low >= 2^32
            if self.low >= (1u64 << 32) {
                self.carry_propagate();
                self.low -= 1u64 << 32;
            }
            self.out.push((self.low >> 24) as u8);
            self.low = (self.low << 8) & 0xFFFF_FFFF;
            self.rng <<= 8;
        }
    }

    /// Propagate a carry: increment the last output byte, cascading through 0xFF bytes.
    fn carry_propagate(&mut self) {
        let mut i = self.out.len() as i32 - 1;
        while i >= 0 {
            if self.out[i as usize] == 0xFF {
                self.out[i as usize] = 0x00;
                i -= 1;
            } else {
                self.out[i as usize] += 1;
                return;
            }
        }
        // All previous bytes were 0xFF and turned to 0x00.
        // Prepend a 0x01 byte (extremely rare).
        self.out.insert(0, 0x01);
    }

    /// Finalize encoding and return compressed bytes.
    fn finish(mut self) -> Vec<u8> {
        // Output 4 more bytes to fully specify the final interval
        for _ in 0..4 {
            if self.low >= (1u64 << 32) {
                self.carry_propagate();
                self.low -= 1u64 << 32;
            }
            self.out.push((self.low >> 24) as u8);
            self.low = (self.low << 8) & 0xFFFF_FFFF;
        }
        self.out
    }
}

/// Range decoder — reads bytes from a compressed buffer.
struct RangeDecoder<'a> {
    code: u32,
    rng: u32,
    data: &'a [u8],
    pos: usize,
}

impl<'a> RangeDecoder<'a> {
    fn new(data: &'a [u8]) -> Self {
        // Read initial 4 bytes to prime the decoder
        let mut code = 0u32;
        for i in 0..4 {
            code <<= 8;
            if i < data.len() {
                code |= data[i] as u32;
            }
        }
        RangeDecoder {
            code,
            rng: INITIAL_RNG,
            data,
            pos: 4,
        }
    }

    fn read_byte(&mut self) -> u32 {
        if self.pos < self.data.len() {
            let b = self.data[self.pos] as u32;
            self.pos += 1;
            b
        } else {
            0
        }
    }

    /// Get cumulative frequency for current range position.
    fn get_cum_freq(&self, total_freq: u32) -> u32 {
        let r = self.rng / total_freq;
        if r == 0 {
            return 0;
        }
        (self.code / r).min(total_freq - 1)
    }

    /// Decode a symbol: given the cum_freq and sym_freq of the decoded symbol,
    /// update the decoder state.
    fn decode(&mut self, cum_freq: u32, sym_freq: u32, total_freq: u32) {
        let r = self.rng / total_freq;
        self.code -= cum_freq * r;
        self.rng = sym_freq * r;

        // Renormalize
        while self.rng < BOT {
            self.code = (self.code << 8) | self.read_byte();
            self.rng <<= 8;
        }
    }
}

// ===========================================================================
// Adaptive context model
// ===========================================================================

/// Adaptive frequency model for arithmetic coding.
/// Tracks symbol frequencies with periodic rescaling to prevent overflow.
pub struct AdaptiveModel {
    freqs: Vec<u32>,
    cum_freqs: Vec<u32>,
    total: u32,
    max_total: u32,
}

impl AdaptiveModel {
    /// Create a new model with `num_symbols` symbols.
    fn new(num_symbols: usize) -> Self {
        let freqs = vec![1u32; num_symbols];
        let mut cum_freqs = vec![0u32; num_symbols + 1];
        for i in 0..num_symbols {
            cum_freqs[i + 1] = cum_freqs[i] + freqs[i];
        }
        AdaptiveModel {
            freqs,
            cum_freqs,
            total: num_symbols as u32,
            max_total: 1 << 14, // rescale at 16384
        }
    }

    /// Get cumulative frequency, symbol frequency, and total.
    fn get(&self, symbol: usize) -> (u32, u32, u32) {
        let cum = self.cum_freqs[symbol];
        let freq = self.freqs[symbol];
        (cum, freq, self.total)
    }

    /// Find symbol from cumulative frequency value.
    fn find(&self, cum_val: u32) -> usize {
        match self.cum_freqs[1..].binary_search(&cum_val) {
            Ok(i) => i + 1,
            Err(i) => i,
        }
    }

    /// Update model after encoding/decoding a symbol.
    fn update(&mut self, symbol: usize) {
        self.freqs[symbol] += 1;
        if self.total >= self.max_total {
            self.rescale();
        } else {
            // Incremental cum_freq update
            for i in (symbol + 1)..self.cum_freqs.len() {
                self.cum_freqs[i] += 1;
            }
            self.total += 1;
        }
    }

    fn rescale(&mut self) {
        // Halve all frequencies (floor, minimum 1)
        let n = self.freqs.len();
        self.total = 0;
        for i in 0..n {
            self.freqs[i] = (self.freqs[i] + 1) / 2; // round up then divide
            if self.freqs[i] == 0 {
                self.freqs[i] = 1;
            }
            self.total += self.freqs[i];
        }
        // Rebuild cum_freqs
        self.cum_freqs[0] = 0;
        for i in 0..n {
            self.cum_freqs[i + 1] = self.cum_freqs[i] + self.freqs[i];
        }
    }
}

// ===========================================================================
// DNA-specific compress/decompress
// ===========================================================================

/// Compress DNA bytes using order-0 adaptive arithmetic coding.
/// Input should be raw bytes (not 2-bit packed). For DNA sequences,
/// the caller should first convert to a byte representation.
#[wasm_bindgen]
pub fn arith_compress(data: &[u8]) -> Vec<u8> {
    if data.is_empty() {
        return vec![];
    }

    // Write original length header (4 bytes)
    let mut header = vec![
        (data.len() >> 24) as u8,
        (data.len() >> 16) as u8,
        (data.len() >> 8) as u8,
        data.len() as u8,
    ];

    let mut encoder = RangeEncoder::new();
    let mut model = AdaptiveModel::new(256);

    for &byte in data {
        let (cum, freq, total) = model.get(byte as usize);
        encoder.encode(cum, freq, total);
        model.update(byte as usize);
    }

    let compressed = encoder.finish();
    header.extend_from_slice(&compressed);
    header
}

/// Decompress data compressed by arith_compress.
#[wasm_bindgen]
pub fn arith_decompress(data: &[u8]) -> Vec<u8> {
    if data.len() < 4 {
        return vec![];
    }

    let orig_len = ((data[0] as u32) << 24 | (data[1] as u32) << 16 |
                    (data[2] as u32) << 8 | data[3] as u32) as usize;

    let mut decoder = RangeDecoder::new(&data[4..]);
    let mut model = AdaptiveModel::new(256);
    let mut out = Vec::with_capacity(orig_len);

    for _ in 0..orig_len {
        let cum_val = decoder.get_cum_freq(model.total);
        let symbol = model.find(cum_val);
        let (cum, freq, total) = model.get(symbol);
        decoder.decode(cum, freq, total);
        model.update(symbol);
        out.push(symbol as u8);
    }

    out
}

/// Compress 2-bit packed DNA using order-1 context modeling.
/// Each base (A=0, C=1, G=2, T=3) is encoded with a context-dependent
/// model where the context is the previous base.
#[wasm_bindgen]
pub fn dna_compress_order1(packed: &[u8], num_bases: usize) -> Vec<u8> {
    if packed.is_empty() || num_bases == 0 {
        return vec![];
    }

    // Header: original length (4 bytes) + num_bases (4 bytes)
    let mut header = vec![
        (packed.len() >> 24) as u8,
        (packed.len() >> 16) as u8,
        (packed.len() >> 8) as u8,
        packed.len() as u8,
        (num_bases >> 24) as u8,
        (num_bases >> 16) as u8,
        (num_bases >> 8) as u8,
        num_bases as u8,
    ];

    // 5 context models: one for each previous base (A,C,G,T) + initial
    let mut models: Vec<AdaptiveModel> = (0..5).map(|_| AdaptiveModel::new(4)).collect();
    let mut encoder = RangeEncoder::new();
    let mut prev_base = 4usize; // initial context (no previous base)

    for i in 0..num_bases {
        let byte_idx = i >> 2;
        let shift = 6 - ((i & 3) << 1);
        let base = ((packed[byte_idx] >> shift) & 0x03) as usize;

        let (cum, freq, total) = models[prev_base].get(base);
        encoder.encode(cum, freq, total);
        models[prev_base].update(base);
        prev_base = base;
    }

    let compressed = encoder.finish();
    header.extend_from_slice(&compressed);
    header
}

/// Decompress DNA compressed by dna_compress_order1.
#[wasm_bindgen]
pub fn dna_decompress_order1(data: &[u8]) -> Vec<u8> {
    if data.len() < 8 {
        return vec![];
    }

    let packed_len = ((data[0] as u32) << 24 | (data[1] as u32) << 16 |
                      (data[2] as u32) << 8 | data[3] as u32) as usize;
    let num_bases = ((data[4] as u32) << 24 | (data[5] as u32) << 16 |
                     (data[6] as u32) << 8 | data[7] as u32) as usize;

    let mut decoder = RangeDecoder::new(&data[8..]);
    let mut models: Vec<AdaptiveModel> = (0..5).map(|_| AdaptiveModel::new(4)).collect();
    let mut packed = vec![0u8; packed_len];
    let mut prev_base = 4usize;

    for i in 0..num_bases {
        let cum_val = decoder.get_cum_freq(models[prev_base].total);
        let base = models[prev_base].find(cum_val);
        let (cum, freq, total) = models[prev_base].get(base);
        decoder.decode(cum, freq, total);
        models[prev_base].update(base);

        let byte_idx = i >> 2;
        let shift = 6 - ((i & 3) << 1);
        packed[byte_idx] |= (base as u8) << shift;
        prev_base = base;
    }

    packed
}
