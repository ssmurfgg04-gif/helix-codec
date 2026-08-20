//! pack.rs — 2-bit DNA packing and bit-parallel operations
//!
//! Native Rust implementation of the hot paths from `src/lib/dna/pack.ts`:
//!   - packDnaToBits: DNA string → packed Uint8Array (4 bases/byte, MSB-first)
//!   - unpackBitsToDna: packed bytes → DNA string
//!   - complement: XOR 0xFF per byte (A↔T, C↔G)
//!   - reverseComplement: complement + reverse base order
//!   - bitParallelHamming: popcount(XOR) / 2
//!   - rollingHash: Rabin-Karp 4-ary rolling hash
//!
//! Pack encoding: A=00, C=01, G=10, T=11. Four bases per byte, MSB-first.

use napi::bindgen_prelude::*;
use napi_derive::napi;

const BITS_TO_CHAR: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Pack a DNA string into 2-bit bytes (4 bases per byte, MSB-first).
///
/// If the length is not a multiple of 4, the final byte is right-padded
/// with zeros in the least-significant bit pairs.
///
/// @param dna DNA string of A/C/G/T characters
/// @returns Packed Uint8Array, length = ceil(dna.length / 4)
#[napi]
pub fn pack_dna_to_bits(dna: String) -> Uint8Array {
    let bytes = dna.as_bytes();
    let num_bytes = (bytes.len() + 3) / 4;
    let mut out = vec![0u8; num_bytes];
    for (i, &c) in bytes.iter().enumerate() {
        let bits: u8 = match c {
            b'A' => 0b00,
            b'C' => 0b01,
            b'G' => 0b10,
            b'T' => 0b11,
            _ => return Uint8Array::new(vec![0u8; 0]),
        };
        let byte_idx = i >> 2;
        let shift = 6 - ((i & 3) << 1) as u8;
        out[byte_idx] |= bits << shift;
    }
    Uint8Array::new(out)
}

/// Unpack 2-bit bytes back into a DNA string.
///
/// @param bits Packed byte array (from packDnaToBits)
/// @param num_bases Number of bases to decode
/// @returns DNA string of length num_bases
#[napi]
pub fn unpack_bits_to_dna(bits: Uint8Array, num_bases: u32) -> String {
    let nb = num_bases as usize;
    if nb > bits.len() * 4 {
        return String::new();
    }
    let mut out = Vec::<u8>::with_capacity(nb);
    let bits_ref: &[u8] = &bits;
    for i in 0..nb {
        let byte_idx = i >> 2;
        let shift = 6 - ((i & 3) << 1) as u8;
        let code = (bits_ref[byte_idx] >> shift) & 0b11;
        out.push(BITS_TO_CHAR[code as usize]);
    }
    String::from_utf8(out).unwrap_or_default()
}

/// Compute the complement of a packed DNA array via bit-flip (XOR 0xFF).
#[napi]
pub fn complement_packed(bits: Uint8Array) -> Uint8Array {
    let out: Vec<u8> = bits.iter().map(|&b| b ^ 0xFF).collect();
    Uint8Array::new(out)
}

/// Compute the reverse complement of a packed DNA array.
///
/// Reverse complement = complement + reverse base order (2-bit pairs).
#[napi]
pub fn reverse_complement_packed(bits: Uint8Array, num_bases: u32) -> Uint8Array {
    let nb = num_bases as usize;
    if nb == 0 {
        return Uint8Array::new(vec![]);
    }
    let num_bytes = (nb + 3) / 4;
    let mut out = vec![0u8; num_bytes];
    let bits_ref: &[u8] = &bits;
    for i in 0..nb {
        let src_pos = nb - 1 - i;
        let src_byte_idx = src_pos >> 2;
        let src_shift = 6 - ((src_pos & 3) << 1) as u8;
        // Complement = XOR 0b11 (per 2-bit pair)
        let base = (bits_ref[src_byte_idx] >> src_shift) & 0b11 ^ 0b11;
        let dst_byte_idx = i >> 2;
        let dst_shift = 6 - ((i & 3) << 1) as u8;
        out[dst_byte_idx] |= base << dst_shift;
    }
    Uint8Array::new(out)
}

/// Compute Hamming distance between two packed DNA arrays via popcount(XOR)/2.
///
/// Each base is 2 bits; XOR has 2 set bits per mismatching base.
#[napi]
pub fn bit_parallel_hamming(a: Uint8Array, b: Uint8Array) -> Result<u32> {
    if a.len() != b.len() {
        return Err(Error::from_reason(format!(
            "Arrays must have same length: {} vs {}",
            a.len(),
            b.len()
        )));
    }
    let mut popcount: u32 = 0;
    for i in 0..a.len() {
        popcount += (a[i] ^ b[i]).count_ones();
    }
    Ok(popcount >> 1)
}

/// Compute Rabin-Karp rolling hashes over a packed DNA bit array.
///
/// Each element is the hash of the k-mer window starting at that position.
/// Uses 4-ary rolling: subtract outgoing base × B^(k-1), shift, add incoming.
#[napi]
pub fn rolling_hash(bits: Uint8Array, window_size: u32) -> Uint32Array {
    let ws = window_size as usize;
    let num_bases = bits.len() * 4;
    if ws > num_bases || ws == 0 {
        return Uint32Array::new(vec![]);
    }
    let result_len = num_bases - ws + 1;
    let mut out = vec![0u32; result_len];
    const B: u64 = 4;
    const MOD: u64 = 0xFFFF_FFFB; // large prime < 2^32

    // Precompute B^(ws-1) mod MOD
    let mut b_pow: u64 = 1;
    for _ in 0..(ws - 1) {
        b_pow = (b_pow * B) % MOD;
    }

    // Initial hash
    let mut hash: u64 = 0;
    for i in 0..ws {
        let byte_idx = i >> 2;
        let shift = 6 - ((i & 3) << 1) as u8;
        let base = ((bits[byte_idx] >> shift) & 0b11) as u64;
        hash = (hash * B + base) % MOD;
    }
    out[0] = hash as u32;

    // Roll
    for i in 1..result_len {
        let out_byte_idx = (i - 1) >> 2;
        let out_shift = 6 - (((i - 1) & 3) << 1) as u8;
        let out_base = ((bits[out_byte_idx] >> out_shift) & 0b11) as u64;

        let in_pos = i + ws - 1;
        let in_byte_idx = in_pos >> 2;
        let in_shift = 6 - ((in_pos & 3) << 1) as u8;
        let in_base = ((bits[in_byte_idx] >> in_shift) & 0b11) as u64;

        let subtract = (out_base * b_pow) % MOD;
        hash = (hash + MOD - subtract) % MOD;
        hash = (hash * B + in_base) % MOD;
        out[i] = hash as u32;
    }

    Uint32Array::new(out)
}

/// Compute GC content of a DNA string. Returns the fraction in [0, 1].
#[napi]
pub fn gc_content(dna: String) -> f64 {
    if dna.is_empty() {
        return 0.0;
    }
    let bytes = dna.as_bytes();
    let mut gc: f64 = 0.0;
    for &c in bytes {
        if c == b'G' || c == b'C' {
            gc += 1.0;
        }
    }
    gc / bytes.len() as f64
}

/// Compute the maximum homopolymer run length in a DNA string.
#[napi]
pub fn max_homopolymer_run(dna: String) -> u32 {
    let bytes = dna.as_bytes();
    if bytes.is_empty() {
        return 0;
    }
    let mut max_run: u32 = 1;
    let mut cur_run: u32 = 1;
    for i in 1..bytes.len() {
        if bytes[i] == bytes[i - 1] {
            cur_run += 1;
            if cur_run > max_run {
                max_run = cur_run;
            }
        } else {
            cur_run = 1;
        }
    }
    max_run
}
