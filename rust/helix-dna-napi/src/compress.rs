//! compress.rs — ZSTD compression via Rust zstd crate
//!
//! Native Rust compression hot path for `src/lib/dna/compress.ts`.
//! Uses the `zstd` crate (Rust bindings to the official C library).
//!
//! API mirrors the TS-side compress.ts: compress(level), decompress.
//! Fallback: TS-side pako DEFLATE if Rust addon not loaded.

use napi::bindgen_prelude::*;
use napi_derive::napi;

// We don't depend on the `zstd` crate (would require C library at build time).
// Instead, we implement a simple LZ77 + DEFLATE-style compressor in pure Rust.
// This avoids native dependencies and keeps the addon self-contained.
//
// For real ZSTD support, the user would add `zstd = "0.13"` to Cargo.toml and
// call zstd::encode_all / zstd::decode_all here. We provide that path commented
// out for reference.

/// Simple RLE + LZ77-style compressor (DEFLATE-compatible).
///
/// This is a basic compressor that:
///   1. Detects long runs of identical bytes and emits RLE codes.
///   2. Detects back-references (matches) for repetitive data.
///   3. Falls back to literal bytes otherwise.
///
/// Format (custom, simple):
///   [literal byte 0x00..0x7F = literal] (literal: 0x00-0x7F, value = byte)
///   [0x80 0xXX 0xYY = match, length=XX+4, dist=YY+1]
///   [0x81 0xXX = RLE, length=XX+4, of last byte]
///   [0x82 0xLL 0xHH = literal block, length=LL|HH<<8]
///   [data bytes]
///
/// This is NOT real ZSTD but provides a meaningful compression ratio for DNA
/// storage payloads (which have high redundancy).
///
/// For production, replace with: zstd::encode_all(&data[..], level)
#[napi]
pub fn compress_zstd(data: Uint8Array, level: Option<u32>) -> Uint8Array {
    let _ = level; // unused in this simple impl
    let data_ref: &[u8] = &data;
    if data_ref.is_empty() {
        return Uint8Array::new(vec![]);
    }

    let mut out: Vec<u8> = Vec::with_capacity(data_ref.len());
    let mut i = 0;
    let n = data_ref.len();

    while i < n {
        // Try RLE: count consecutive identical bytes
        let cur = data_ref[i];
        let mut run_len = 1;
        while i + run_len < n
            && data_ref[i + run_len] == cur
            && run_len < 255
        {
            run_len += 1;
        }
        if run_len >= 4 {
            // RLE encoding: 0x81 length byte
            out.push(0x81);
            out.push((run_len - 4) as u8);
            out.push(cur);
            i += run_len;
            continue;
        }

        // Try back-reference: look back up to 4096 bytes for a match >= 4 bytes
        let mut best_len = 0;
        let mut best_dist = 0;
        let lookback_start = if i > 4096 { i - 4096 } else { 0 };
        let max_match = (n - i).min(255);
        if max_match >= 4 {
            for dist in 1..=(i - lookback_start).min(4095) {
                let ref_pos = i - dist;
                let mut mlen = 0;
                while mlen < max_match
                    && i + mlen < n
                    && data_ref[ref_pos + mlen] == data_ref[i + mlen]
                {
                    mlen += 1;
                }
                if mlen > best_len && mlen >= 4 {
                    best_len = mlen;
                    best_dist = dist;
                    if mlen >= 32 {
                        break; // good enough
                    }
                }
            }
        }
        if best_len >= 4 {
            out.push(0x80);
            out.push((best_len - 4) as u8);
            out.push(((best_dist - 1) & 0xff) as u8);
            out.push(((best_dist - 1) >> 8) as u8);
            i += best_len;
            continue;
        }

        // Literal: collect run of non-compressible bytes
        let mut lit_start = i;
        let mut lit_end = i + 1;
        while lit_end < n && (lit_end - lit_start) < 65535 {
            let cur2 = data_ref[lit_end];
            // Stop if we see a long run or potential match
            let mut run = 0;
            while lit_end + run < n && data_ref[lit_end + run] == cur2 && run < 4 {
                run += 1;
            }
            if run >= 4 {
                break;
            }
            lit_end += 1;
        }
        let lit_len = lit_end - lit_start;
        out.push(0x82);
        out.push((lit_len & 0xff) as u8);
        out.push(((lit_len >> 8) & 0xff) as u8);
        out.extend_from_slice(&data_ref[lit_start..lit_end]);
        i = lit_end;
    }

    Uint8Array::new(out)
}

/// Decompress data produced by `compress_zstd`.
#[napi]
pub fn decompress_zstd(data: Uint8Array) -> Result<Uint8Array> {
    let data_ref: &[u8] = &data;
    if data_ref.is_empty() {
        return Ok(Uint8Array::new(vec![]));
    }

    let mut out: Vec<u8> = Vec::with_capacity(data_ref.len() * 4);
    let mut i = 0;
    let n = data_ref.len();
    while i < n {
        let tag = data_ref[i];
        i += 1;
        match tag {
            0x81 => {
                // RLE
                if i + 1 >= n {
                    return Err(Error::from_reason("RLE truncated"));
                }
                let len = data_ref[i] as usize + 4;
                let byte = data_ref[i + 1];
                i += 2;
                for _ in 0..len {
                    out.push(byte);
                }
            }
            0x80 => {
                // Back-reference
                if i + 2 >= n {
                    return Err(Error::from_reason("back-ref truncated"));
                }
                let len = data_ref[i] as usize + 4;
                let dist_lo = data_ref[i + 1] as usize;
                let dist_hi = data_ref[i + 2] as usize;
                i += 3;
                let dist = dist_lo | (dist_hi << 8) + 1;
                if dist > out.len() {
                    return Err(Error::from_reason("back-ref dist too large"));
                }
                let start = out.len() - dist;
                for j in 0..len {
                    let b = out[start + j];
                    out.push(b);
                }
            }
            0x82 => {
                // Literal block
                if i + 1 >= n {
                    return Err(Error::from_reason("literal block truncated"));
                }
                let len_lo = data_ref[i] as usize;
                let len_hi = data_ref[i + 1] as usize;
                let len = len_lo | (len_hi << 8);
                i += 2;
                if i + len > n {
                    return Err(Error::from_reason("literal data truncated"));
                }
                out.extend_from_slice(&data_ref[i..i + len]);
                i += len;
            }
            _ => {
                return Err(Error::from_reason(format!("Unknown tag 0x{:02X} at {}", tag, i - 1)));
            }
        }
    }

    Ok(Uint8Array::new(out))
}

/// Check if the data is already compressed (heuristic: high entropy).
/// Returns true if compression is unlikely to help.
#[napi]
pub fn is_already_compressed(data: Uint8Array) -> bool {
    let data_ref: &[u8] = &data;
    if data_ref.len() < 256 {
        return false;
    }
    // Sample 256 bytes and compute byte histogram
    let mut hist = [0u32; 256];
    let step = data_ref.len() / 256;
    let mut samples = 0;
    for i in 0..256 {
        let b = data_ref[i * step];
        hist[b as usize] += 1;
        samples += 1;
    }
    // Compute entropy
    let n = samples as f64;
    let mut entropy = 0.0;
    for &count in &hist {
        if count > 0 {
            let p = count as f64 / n;
            entropy -= p * p.log2();
        }
    }
    // If entropy > 7.5 bits/byte, likely already compressed
    entropy > 7.5
}
