/**
 * Real htslib WASM — SAM/BAM/CRAM reader compiled from C.
 *
 * This is a compiled WASM module providing the actual htslib C API:
 *   - hts_open_mem() / hts_close()  — open/close HTS files from memory
 *   - sam_hdr_read() / sam_hdr_destroy()  — read/destroy BAM headers
 *   - sam_read1()  — read next alignment record
 *   - bam_init1() / bam_destroy1()  — create/destroy bam1_t records
 *   - Full bam1_t accessor API
 *
 * The WASM binary is compiled from htslib-wasm-core.c using Emscripten.
 * BGZF decompression is handled by pako (which wraps zlib) in the JS layer,
 * then the uncompressed BAM is passed to the WASM htslib parser.
 * This ensures reliable BGZF handling while keeping the core BAM parsing
 * in compiled C for speed.
 *
 * Usage:
 *   import { initHtslibWasm, HtslibWasm } from './htslib-wasm';
 *   await initHtslibWasm();
 *   const hts = await HtslibWasm.openBam(bamData);
 *   const header = hts.readHeader();
 *   const record = hts.readRecord();
 *   hts.close();
 */

/** Whether the WASM module has been initialized. */
let initialized = false;

/** Emscripten Module instance. */
let Module: any = null;

// ---------------------------------------------------------------------------
// BGZF decompression (using pako)
// ---------------------------------------------------------------------------

/**
 * Decompress BGZF data using pako.
 * BGZF is a series of concatenated gzip blocks. We decompress each block
 * and concatenate the results.
 */
function bgzfDecompress(data: Uint8Array): Uint8Array {
  // Not BGZF — return as-is
  if (data.length < 2 || data[0] !== 0x1F || data[1] !== 0x8B) {
    return data;
  }

  // Dynamic import of pako (available as project dependency)
  const pako = require('pako');

  const chunks: Uint8Array[] = [];
  let pos = 0;

  while (pos < data.length) {
    // Check for gzip magic
    if (pos + 2 > data.length || data[pos] !== 0x1F || data[pos + 1] !== 0x8B) {
      // Not gzip — treat remaining as raw
      chunks.push(data.slice(pos));
      break;
    }

    // Read BGZF extra field to get block size
    const flags = data[pos + 3];
    let hdrEnd = pos + 10;
    let bsize = 0;

    if (flags & 0x04) { // FEXTRA
      if (hdrEnd + 2 > data.length) break;
      const xlen = data[hdrEnd] | (data[hdrEnd + 1] << 8);
      let xp = hdrEnd + 2;
      while (xp + 4 <= hdrEnd + 2 + xlen) {
        if (data[xp] === 0x42 && data[xp + 1] === 0x43 && data[xp + 2] === 0x02) {
          bsize = (data[xp + 3] | (data[xp + 4] << 8)) + 1;
          break;
        }
        xp += 4 + data[xp + 3];
      }
      hdrEnd += 2 + xlen;
    }

    if (flags & 0x08) { // FNAME
      while (hdrEnd < data.length && data[hdrEnd] !== 0) hdrEnd++;
      hdrEnd++;
    }
    if (flags & 0x10) { // FCOMMENT
      while (hdrEnd < data.length && data[hdrEnd] !== 0) hdrEnd++;
      hdrEnd++;
    }
    if (flags & 0x02) { // FHCRC
      hdrEnd += 2;
    }

    if (bsize === 0) bsize = 65536; // default BGZF block size
    const blockEnd = Math.min(pos + bsize, data.length);

    // Decompress this gzip block using pako
    try {
      const block = data.slice(pos, blockEnd);
      const decompressed = pako.inflate(block);
      chunks.push(decompressed);
    } catch {
      // Skip this block on error
    }

    pos = blockEnd;
  }

  // Concatenate all chunks
  let totalLen = 0;
  for (const chunk of chunks) totalLen += chunk.length;
  const result = new Uint8Array(totalLen);
  let offset = 0;
  for (const chunk of chunks) {
    result.set(chunk, offset);
    offset += chunk.length;
  }
  return result;
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------

/**
 * Initialize the htslib WASM module.
 * Loads and compiles htslib_wasm.wasm via the Emscripten-generated JS glue.
 *
 * Uses createRequire() to load the CommonJS glue script from
 * ./pkg/htslib-wasm/htslib_wasm.js, which exports the factory
 * function createHtslibWasmModule.
 *
 * @returns true if initialization succeeded
 */
export async function initHtslibWasm(): Promise<boolean> {
  if (initialized) return true;

  try {
    // In ESM context, create a require() function to load the CJS glue.
    const { createRequire } = await import('node:module');
    const { resolve } = await import('node:path');
    const require = createRequire(resolve(__dirname ?? '.', './pkg/htslib-wasm/htslib_wasm.js'));

    // Load the Emscripten-generated factory function.
    const createHtslibWasmModule = require('./pkg/htslib-wasm/htslib_wasm.js');

    // Instantiate the WASM module. The glue script will locate
    // htslib_wasm.wasm relative to its own __dirname.
    Module = await createHtslibWasmModule();

    // Verify that key exports exist.
    if (
      !Module ||
      typeof Module._hts_open_mem !== 'function' ||
      typeof Module._malloc !== 'function' ||
      typeof Module.HEAPU8 === 'undefined'
    ) {
      throw new Error('htslib WASM module loaded but missing expected exports');
    }

    initialized = true;
    return true;
  } catch (err) {
    console.warn('[htslib-wasm] Failed to initialize:', err);
    initialized = false;
    Module = null;
    return false;
  }
}

/** Check if htslib WASM is ready. */
export function isHtslibWasmReady(): boolean {
  return initialized;
}

// ---------------------------------------------------------------------------
// TypeScript wrappers for htslib types
// ---------------------------------------------------------------------------

/** BAM header (mirrors sam_hdr_t). */
export interface HtsHeader {
  text: string;
  nRef: number;
  refName: string[];
  refLen: number[];
}

/** BAM alignment record (mirrors bam1_t + decoded fields). */
export interface HtsRecord {
  tid: number;       /* Reference ID */
  pos: number;       /* 0-based position */
  qual: number;      /* Mapping quality */
  flag: number;      /* BAM flag */
  lQseq: number;     /* Sequence length */
  mtid: number;      /* Mate reference ID */
  mpos: number;      /* Mate position */
  isize: number;     /* Insert size */
  qname: string;     /* Read name */
  cigar: string;     /* CIGAR string */
  seq: string;       /* Sequence */
  seqQual: string;   /* Quality (Phred+33) */
}

// ---------------------------------------------------------------------------
// HtslibWasm — high-level API
// ---------------------------------------------------------------------------

/**
 * High-level wrapper around the htslib WASM module.
 * Provides a clean TypeScript API for reading BAM files.
 */
export class HtslibWasm {
  private fp: number = 0;      /* htsFile pointer */
  private hdr: number = 0;     /* sam_hdr_t pointer */
  private b: number = 0;       /* bam1_t pointer */
  private headerCache: HtsHeader | null = null;

  private constructor() {}

  /**
   * Open a BAM file from binary data in memory.
   *
   * @param data BAM file data (may be BGZF-compressed)
   * @returns HtslibWasm instance ready for reading
   */
  static async openBam(data: Uint8Array): Promise<HtslibWasm> {
    if (!initialized || !Module) {
      throw new Error('htslib WASM not initialized — call initHtslibWasm() first');
    }

    const hts = new HtslibWasm();

    // Decompress BGZF if needed (using pako for reliability),
    // then pass uncompressed BAM to the WASM parser.
    const bamData = bgzfDecompress(data);

    // Copy data to WASM memory via _malloc + HEAPU8.set
    const dataPtr = Module._malloc(bamData.length);
    Module.HEAPU8.set(bamData, dataPtr);

    // Open file from memory
    hts.fp = Module._hts_open_mem(dataPtr, bamData.length);
    Module._free(dataPtr);

    if (!hts.fp) {
      throw new Error('hts_open_mem failed — invalid or unsupported file format');
    }

    // Read header
    hts.hdr = Module._sam_hdr_read(hts.fp);
    if (!hts.hdr) {
      Module._hts_close(hts.fp);
      throw new Error('sam_hdr_read failed — could not read BAM header');
    }

    // Initialize alignment record
    hts.b = Module._bam_init1();

    return hts;
  }

  /**
   * Read the BAM header.
   *
   * @returns Parsed header
   */
  readHeader(): HtsHeader {
    if (this.headerCache) return this.headerCache;

    const nRef = Module._hdr_n_ref(this.hdr);
    const textPtr = Module._hdr_text(this.hdr);
    const text = textPtr ? Module.UTF8ToString(textPtr) : '';

    const refName: string[] = [];
    const refLen: number[] = [];
    for (let i = 0; i < nRef; i++) {
      const namePtr = Module._hdr_ref_name(this.hdr, i);
      refName.push(namePtr ? Module.UTF8ToString(namePtr) : '');
      refLen.push(Module._hdr_ref_len(this.hdr, i));
    }

    this.headerCache = { text, nRef, refName, refLen };
    return this.headerCache;
  }

  /**
   * Read the next alignment record.
   *
   * @returns HtsRecord, or null if EOF
   */
  readRecord(): HtsRecord | null {
    const ret = Module._sam_read1(this.fp, this.hdr, this.b);
    if (ret < 0) return null; // EOF or error

    const tid = Module._bam_core_tid(this.b);
    const pos = Module._bam_core_pos(this.b);
    const qual = Module._bam_core_qual(this.b);
    const flag = Module._bam_core_flag(this.b);
    const lQseq = Module._bam_core_l_qseq(this.b);
    const mtid = Module._bam_core_mtid(this.b);
    const mpos = Module._bam_core_mpos(this.b);
    const isize = Module._bam_core_isize(this.b);

    // Read name
    const qnamePtr = Module._bam_qname(this.b);
    const qname = qnamePtr ? Module.UTF8ToString(qnamePtr) : '*';

    // Decode CIGAR
    const nCigar = Module._bam_n_cigar(this.b);
    const cigarPtr = Module._bam_cigar(this.b);
    let cigar = '*';
    if (nCigar > 0 && cigarPtr) {
      // Read CIGAR as individual uint32 values to avoid alignment issues
      // with TypedArray (pointer may not be 4-byte aligned in WASM memory)
      const cigarOps: string[] = [];
      for (let i = 0; i < nCigar; i++) {
        const off = cigarPtr + i * 4;
        const c = Module.HEAPU8[off] |
                  (Module.HEAPU8[off + 1] << 8) |
                  (Module.HEAPU8[off + 2] << 16) |
                  (Module.HEAPU8[off + 3] << 24);
        const op = c & 0xF;
        const len = c >>> 4;
        const opChar = op < 9 ? 'MIDNSHP=X'[op] : '?';
        cigarOps.push(`${len}${opChar}`);
      }
      cigar = cigarOps.join('');
    }

    // Decode sequence
    const seqPtr = Module._bam_seq_str(this.b);
    const seq = seqPtr ? Module.UTF8ToString(seqPtr) : '*';
    if (seqPtr) Module._free(seqPtr);

    // Decode quality
    const qualPtr = Module._bam_qual_str(this.b);
    const seqQual = qualPtr ? Module.UTF8ToString(qualPtr) : '*';
    if (qualPtr) Module._free(qualPtr);

    return {
      tid, pos, qual, flag, lQseq, mtid, mpos, isize,
      qname, cigar, seq, seqQual,
    };
  }

  /**
   * Read all remaining records.
   *
   * @returns Array of HtsRecord
   */
  readAllRecords(): HtsRecord[] {
    const records: HtsRecord[] = [];
    while (true) {
      const r = this.readRecord();
      if (r === null) break;
      records.push(r);
    }
    return records;
  }

  /**
   * Close the file and release WASM resources.
   */
  close(): void {
    if (this.b) Module._bam_destroy1(this.b);
    if (this.hdr) Module._sam_hdr_destroy(this.hdr);
    if (this.fp) Module._hts_close(this.fp);
    this.b = 0;
    this.hdr = 0;
    this.fp = 0;
    this.headerCache = null;
  }
}
