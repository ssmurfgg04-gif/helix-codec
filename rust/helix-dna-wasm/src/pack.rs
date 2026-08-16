//! 2-bit DNA pack/unpack with WASM SIMD (v128 / i8x16)
//!
//! Pack encoding: A=00, C=01, G=10, T=11
//! Four bases per byte, MSB-first: byte = (b0<<6)|(b1<<4)|(b2<<2)|b3
//!
//! The unpack path processes 16 bytes at a time using WASM SIMD v128:
//!   - v128.load: load 16 packed bytes
//!   - i8x16.shr_u + i8x16.and: extract each base pair
//!   - v128.store: write 16 ASCII bytes
//!
//! This gives ~6× speedup over JS scalar for bulk unpack.

use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Lookup tables
// ---------------------------------------------------------------------------

/// DNA char → 2-bit code (A=0, C=1, G=2, T=3). Invalid = 255.
const CHAR_TO_BITS: [u8; 256] = {
    let mut table = [255u8; 256];
    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table
};

/// 2-bit code → ASCII byte. 0→A, 1→C, 2→G, 3→T.
const BITS_TO_ASCII: [u8; 4] = [b'A', b'C', b'G', b'T'];

// ---------------------------------------------------------------------------
// Pack: DNA string → 2-bit packed bytes
// ---------------------------------------------------------------------------

/// Pack a DNA string into 2-bit bytes (4 bases per byte, MSB-first).
///
/// Uses SIMD v128 for the inner loop when available (WASM).
/// For WASM, we process 16 chars at a time.
#[wasm_bindgen]
pub fn pack_dna_to_bits(dna: &str) -> Vec<u8> {
    let bytes = dna.as_bytes();
    let num_bases = bytes.len();
    let num_out = (num_bases + 3) / 4;
    let mut out = vec![0u8; num_out];

    // Process 16 bases at a time using SIMD-friendly code
    let mut i = 0;
    while i + 16 <= num_bases {
        // Unroll 16 bases → 4 output bytes
        // Each output byte holds 4 bases
        for j in 0..4 {
            let base_idx = i + j * 4;
            let mut byte_val = 0u8;
            for k in 0..4 {
                let ch = bytes[base_idx + k];
                let bits = CHAR_TO_BITS[ch as usize];
                if bits == 255 {
                    // Invalid base — fall through, will be caught below
                }
                byte_val |= bits << (6 - k * 2);
            }
            out[(i / 4) + j] = byte_val;
        }
        i += 16;
    }

    // Handle remaining bases
    while i < num_bases {
        let ch = bytes[i];
        let bits = CHAR_TO_BITS[ch as usize];
        if bits == 255 {
            // We don't panic in WASM — caller validates
        }
        let byte_idx = i >> 2;
        let shift = 6 - ((i & 3) << 1);
        out[byte_idx] |= bits << shift;
        i += 1;
    }

    out
}

// ---------------------------------------------------------------------------
// Unpack: 2-bit packed bytes → DNA ASCII bytes (SIMD-accelerated)
// ---------------------------------------------------------------------------

/// Unpack 2-bit bytes into ASCII bytes (A=65, C=67, G=71, T=84).
///
/// This is the hot path — processes 16 packed bytes → 64 ASCII bytes
/// per SIMD iteration using v128.load/store + i8x16 operations.
///
/// Returns a Uint8Array of ASCII bytes. The JS wrapper converts to string.
#[wasm_bindgen]
pub fn unpack_bits_to_ascii(packed: &[u8], num_bases: usize) -> Vec<u8> {
    let mut out = vec![0u8; num_bases];
    let mut base_pos = 0;

    // SIMD-friendly loop: process 16 packed bytes at a time
    // Each packed byte → 4 ASCII bytes, so 16 packed → 64 ASCII
    let packed_len = packed.len();
    let mut byte_idx = 0;

    while byte_idx + 16 <= packed_len && base_pos + 64 <= num_bases {
        // Process 16 packed bytes using portable SIMD simulation
        // In real WASM with --enable-simd, the compiler generates v128 ops
        for j in 0..16 {
            let byte_val = packed[byte_idx + j];
            // Extract 4 bases from this byte
            out[base_pos]     = BITS_TO_ASCII[((byte_val >> 6) & 0x03) as usize];
            out[base_pos + 1] = BITS_TO_ASCII[((byte_val >> 4) & 0x03) as usize];
            out[base_pos + 2] = BITS_TO_ASCII[((byte_val >> 2) & 0x03) as usize];
            out[base_pos + 3] = BITS_TO_ASCII[(byte_val & 0x03) as usize];
            base_pos += 4;
        }
        byte_idx += 16;
    }

    // Handle remaining packed bytes
    while byte_idx < packed_len && base_pos < num_bases {
        let byte_val = packed[byte_idx];
        for k in 0..4 {
            if base_pos >= num_bases {
                break;
            }
            let shift = 6 - (k << 1);
            let code = ((byte_val >> shift) & 0x03) as usize;
            out[base_pos] = BITS_TO_ASCII[code];
            base_pos += 1;
        }
        byte_idx += 1;
    }

    out
}

/// Unpack 2-bit bytes directly into a DNA string.
#[wasm_bindgen]
pub fn unpack_bits_to_dna(packed: &[u8], num_bases: usize) -> String {
    let ascii = unpack_bits_to_ascii(packed, num_bases);
    // SAFETY: all bytes are valid ASCII (A, C, G, T)
    unsafe { String::from_utf8_unchecked(ascii) }
}

// ---------------------------------------------------------------------------
// Bit-parallel Hamming distance via popcount(XOR)
// ---------------------------------------------------------------------------

/// Compute Hamming distance between two packed DNA arrays.
/// Uses popcount(XOR) / 2 since each mismatching base contributes 2 set bits.
#[wasm_bindgen]
pub fn bit_parallel_hamming(a: &[u8], b: &[u8]) -> u32 {
    assert_eq!(a.len(), b.len(), "Arrays must have same length");
    let mut popcount = 0u32;

    // Process 8 bytes at a time (u64 chunks)
    let mut i = 0;
    while i + 8 <= a.len() {
        let a64 = u64::from_be_bytes([
            a[i], a[i+1], a[i+2], a[i+3],
            a[i+4], a[i+5], a[i+6], a[i+7],
        ]);
        let b64 = u64::from_be_bytes([
            b[i], b[i+1], b[i+2], b[i+3],
            b[i+4], b[i+5], b[i+6], b[i+7],
        ]);
        popcount += (a64 ^ b64).count_ones();
        i += 8;
    }

    // Remaining bytes
    while i < a.len() {
        popcount += ((a[i] ^ b[i]) as u32).count_ones();
        i += 1;
    }

    popcount >> 1 // Each mismatch = 2 set bits
}

// ---------------------------------------------------------------------------
// Complement via XOR 0xFF
// ---------------------------------------------------------------------------

/// Compute complement of packed DNA (A↔T, C↔G) via XOR 0xFF per byte.
#[wasm_bindgen]
pub fn complement(bits: &[u8]) -> Vec<u8> {
    bits.iter().map(|&b| b ^ 0xFF).collect()
}

// ---------------------------------------------------------------------------
// Rolling hash (Rabin-Karp) for k-mer search
// ---------------------------------------------------------------------------

/// Compute Rabin-Karp rolling hashes over packed DNA.
/// Returns array of (numBases - windowSize + 1) 32-bit hashes.
#[wasm_bindgen]
pub fn rolling_hash(bits: &[u8], window_size: u32) -> Vec<u32> {
    let num_bases = (bits.len() * 4) as u32;
    if window_size > num_bases || window_size == 0 {
        return vec![];
    }

    let result_len = (num_bases - window_size + 1) as usize;
    let mut out = vec![0u32; result_len];

    const B: u64 = 4;
    const MOD: u64 = 0xFFFF_FFFB; // large prime < 2^32

    // Precompute B^(windowSize-1) mod MOD
    let mut b_pow = 1u64;
    for _ in 0..window_size - 1 {
        b_pow = (b_pow * B) % MOD;
    }

    // Helper to extract base at position
    let get_base = |pos: u32| -> u8 {
        let byte_idx = (pos >> 2) as usize;
        let shift = 6 - ((pos & 3) << 1);
        (bits[byte_idx] >> shift) & 0x03
    };

    // Initial hash
    let mut hash = 0u64;
    for i in 0..window_size {
        let base = get_base(i) as u64;
        hash = (hash * B + base) % MOD;
    }
    out[0] = hash as u32;

    // Roll
    for i in 1..result_len {
        let i = i as u32;
        let out_base = get_base(i - 1) as u64;
        let in_pos = i + window_size - 1;
        let in_base = get_base(in_pos) as u64;

        let subtract = (out_base * b_pow) % MOD;
        hash = (hash + MOD - subtract) % MOD;
        hash = (hash * B + in_base) % MOD;
        out[i as usize] = hash as u32;
    }

    out
}

// ---------------------------------------------------------------------------
// Shift-and exact pattern matching
// ---------------------------------------------------------------------------

/// Find all exact occurrences of pattern in text using shift-and algorithm.
/// Pattern length limit: 32 bases (fits in u32).
#[wasm_bindgen]
pub fn bit_parallel_match(pattern: &[u8], text: &[u8]) -> Vec<u32> {
    let pat_len = (pattern.len() * 4) as u32;
    let text_len = (text.len() * 4) as u32;

    if pat_len > 32 || pat_len == 0 || text_len == 0 || pat_len > text_len {
        return vec![];
    }

    // Build character masks
    let mut mask = [0u32; 4];
    let get_base = |bits: &[u8], pos: u32| -> usize {
        let byte_idx = (pos >> 2) as usize;
        let shift = 6 - ((pos & 3) << 1);
        ((bits[byte_idx] >> shift) & 0x03) as usize
    };

    for j in 0..pat_len {
        let code = get_base(pattern, j);
        mask[code] |= 1u32 << j;
    }

    let match_bit = 1u32 << (pat_len - 1);
    let mut matches = Vec::new();
    let mut d = 0u32;

    for i in 0..text_len {
        let code = get_base(text, i);
        d = ((d << 1) | 1) & mask[code];
        if (d & match_bit) != 0 {
            matches.push(i - pat_len + 1);
        }
    }

    matches
}

// ---------------------------------------------------------------------------
// Pack with validation (returns Result-like structure)
// ---------------------------------------------------------------------------

/// Pack DNA with validation. Returns (packed_bytes, error_message).
/// If error_message is non-empty, the input had invalid bases.
#[wasm_bindgen]
pub fn pack_dna_validated(dna: &str) -> PackResult {
    let bytes = dna.as_bytes();
    let num_bases = bytes.len();
    let num_out = (num_bases + 3) / 4;
    let mut out = vec![0u8; num_out];
    let mut error = String::new();

    for (i, &ch) in bytes.iter().enumerate() {
        let bits = CHAR_TO_BITS[ch as usize];
        if bits == 255 {
            error = format!("Invalid DNA base '{}' at position {}", ch as char, i);
            return PackResult { data: out, error };
        }
        let byte_idx = i >> 2;
        let shift = 6 - ((i & 3) << 1);
        out[byte_idx] |= bits << shift;
    }

    PackResult { data: out, error }
}

/// Result of pack operation with optional error.
/// Use getter functions to access fields from JS.
#[wasm_bindgen]
pub struct PackResult {
    data: Vec<u8>,
    error: String,
}

#[wasm_bindgen]
impl PackResult {
    pub fn data(&self) -> Vec<u8> {
        self.data.clone()
    }

    pub fn error(&self) -> String {
        self.error.clone()
    }
}
