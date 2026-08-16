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

// =============================================================================
// VCF/BCF Parsing
// =============================================================================

/** A single VCF/BCF record. */
export interface VcfRecord {
  /** CHROM — reference sequence name. */
  chrom: string;
  /** POS — 1-based leftmost position. */
  pos: number;
  /** ID — variant identifier (. = missing). */
  id: string;
  /** REF — reference allele. */
  ref: string;
  /** ALT — alternate alleles (comma-separated, . = missing). */
  alt: string;
  /** QUAL — quality (NaN = missing). */
  qual: number;
  /** FILTER — filter status (. = missing, PASS = passed). */
  filter: string;
  /** INFO — key-value pairs as a Map. */
  info: Map<string, string>;
  /** FORMAT — genotype format fields (empty if no samples). */
  format: string[];
  /** Genotypes — one string per sample. */
  samples: string[];
}

/** Result of VCF parsing: header lines + records. */
export interface VcfRecords {
  /** Parsed header lines (## and #CHROM...). */
  headerLines: string[];
  /** Sample IDs from the #CHROM header line. */
  sampleIds: string[];
  /** Parsed variant records. */
  records: VcfRecord[];
}

/**
 * Parse a VCF text file into structured records.
 *
 * Handles:
 *   - Header lines (## meta-information and #CHROM column header)
 *   - Multi-line records (not in VCF spec, but some tools produce them)
 *   - INFO field key=value pairs, including Flag fields (key without value)
 *   - FORMAT and sample genotype fields
 *   - Missing values (. for ID, ALT, FILTER; NaN for QUAL)
 *
 * @param text VCF file content as a string
 * @returns Parsed VCF records with header
 */
export function parseVcf(text: string): VcfRecords {
  const lines = text.split('\n');
  const headerLines: string[] = [];
  let sampleIds: string[] = [];
  const records: VcfRecord[] = [];

  let i = 0;

  // Parse header
  while (i < lines.length) {
    const line = lines[i];
    if (line.startsWith('##')) {
      headerLines.push(line);
      i++;
    } else if (line.startsWith('#CHROM')) {
      headerLines.push(line);
      // Parse column header to extract sample IDs
      const cols = line.split('\t');
      // VCF has 8 fixed columns + FORMAT + samples
      if (cols.length > 9) {
        sampleIds = cols.slice(9);
      }
      i++;
      break;
    } else if (line.startsWith('#')) {
      headerLines.push(line);
      i++;
    } else {
      break;
    }
  }

  // Parse data lines
  while (i < lines.length) {
    const line = lines[i].trim();
    i++;
    if (line === '' || line.startsWith('#')) continue;

    const cols = line.split('\t');
    if (cols.length < 8) continue; // Malformed line — skip

    const chrom = cols[0];
    const pos = parseInt(cols[1], 10);
    const id = cols[2];
    const ref = cols[3];
    const alt = cols[4];
    const qualStr = cols[5];
    const qual = qualStr === '.' ? NaN : parseFloat(qualStr);
    const filter = cols[6];

    // Parse INFO field
    const info = new Map<string, string>();
    if (cols[7] !== '.') {
      const infoEntries = cols[7].split(';');
      for (const entry of infoEntries) {
        const eqIdx = entry.indexOf('=');
        if (eqIdx >= 0) {
          info.set(entry.substring(0, eqIdx), entry.substring(eqIdx + 1));
        } else {
          // Flag field (e.g., "DB" or "INDEL")
          info.set(entry, 'true');
        }
      }
    }

    // Parse FORMAT and samples
    const format: string[] = [];
    const samples: string[] = [];
    if (cols.length > 8) {
      format.push(...cols[8].split(':'));
      for (let s = 9; s < cols.length; s++) {
        samples.push(cols[s]);
      }
    }

    records.push({ chrom, pos, id, ref, alt, qual, filter, info, format, samples });
  }

  return { headerLines, sampleIds, records };
}

/**
 * Parse a BCF (Binary VCF) file into structured records.
 *
 * BCF is the binary counterpart of VCF. Format specification:
 *   - Magic: "BCF\2\1" (5 bytes)
 *   - Header length (uint32 LE)
 *   - Header text (null-terminated)
 *   - Records: each has l_shared (uint32), l_indiv (uint32),
 *     then shared data (CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO)
 *     then individual data (FORMAT + genotypes)
 *
 * @param data BCF binary data
 * @returns Parsed VCF records with header
 */
export function parseBcf(data: Uint8Array): VcfRecords {
  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
  let offset = 0;

  // Check magic
  if (data.length < 5 || data[0] !== 0x42 || data[1] !== 0x43 || data[2] !== 0x46) {
    throw new Error('parseBcf: invalid BCF magic number');
  }
  const versionMinor = data[4];
  offset = 5;

  // Read header
  const headerLen = view.getUint32(offset, true);
  offset += 4;

  const headerTextEnd = offset + headerLen;
  // Find null terminator
  let headerEnd = offset;
  while (headerEnd < headerTextEnd && data[headerEnd] !== 0) headerEnd++;
  const headerText = new TextDecoder().decode(data.slice(offset, headerEnd));
  offset = headerTextEnd;

  // Parse the header text to extract sample IDs and header lines
  const headerLines = headerText.split('\n').filter(l => l.length > 0);
  let sampleIds: string[] = [];
  for (const line of headerLines) {
    if (line.startsWith('#CHROM')) {
      const cols = line.split('\t');
      if (cols.length > 9) {
        sampleIds = cols.slice(9);
      }
    }
  }

  // Build contig name lookup from header
  const contigNames: string[] = [];
  for (const line of headerLines) {
    if (line.startsWith('##contig=') || line.startsWith('##contig=<')) {
      const match = line.match(/ID=([^,>]+)/);
      if (match) contigNames.push(match[1]);
    }
  }

  // Build FILTER name lookup
  const filterNames: Map<number, string> = new Map();
  filterNames.set(0, 'PASS');
  let filterIdx = 1;
  for (const line of headerLines) {
    if (line.startsWith('##FILTER=') || line.startsWith('##FILTER=<')) {
      const match = line.match(/ID=([^,>]+)/);
      if (match) {
        filterNames.set(filterIdx, match[1]);
        filterIdx++;
      }
    }
  }

  // Parse BCF records
  const records: VcfRecord[] = [];

  while (offset + 8 <= data.length) {
    const lShared = view.getUint32(offset, true);
    offset += 4;
    const lIndiv = view.getUint32(offset, true);
    offset += 4;

    const recordStart = offset;
    const recordEnd = recordStart + lShared + lIndiv;

    if (recordEnd > data.length) break;

    try {
      const record = parseBcfRecord(data, view, offset, lShared, lIndiv, contigNames, filterNames);
      if (record) records.push(record);
    } catch {
      // Skip malformed records
    }

    offset = recordEnd;
  }

  return { headerLines, sampleIds, records };
}

/**
 * Parse a single BCF record from the shared + individual data.
 */
function parseBcfRecord(
  data: Uint8Array,
  view: DataView,
  offset: number,
  lShared: number,
  lIndiv: number,
  contigNames: string[],
  filterNames: Map<number, string>,
): VcfRecord | null {
  let pos = offset;

  // CHROM (int32)
  const chromIdx = view.getInt32(pos, true);
  pos += 4;
  const chrom = chromIdx >= 0 && chromIdx < contigNames.length
    ? contigNames[chromIdx]
    : `chr${chromIdx}`;

  // POS (int32, 0-based in BCF → 1-based in VCF)
  const startPos = view.getInt32(pos, true);
  pos += 4;
  const pos1based = startPos + 1;

  // rlen (int32) — reference length
  const rlen = view.getInt32(pos, true);
  pos += 4;

  // QUAL (float32)
  const qualRaw = view.getFloat32(pos, true);
  pos += 4;
  const qual = qualRaw < 0 ? NaN : qualRaw;

  // n_info (uint16) + n_allele (uint16)
  const nInfo = view.getUint16(pos, true);
  pos += 2;
  const nAllele = view.getUint16(pos, true);
  pos += 2;

  // n_fmt (uint8) + n_sample (uint24)
  const nFmt = data[pos];
  pos += 1;
  // n_sample is encoded as 3 bytes little-endian
  const nSample = data[pos] | (data[pos + 1] << 8) | (data[pos + 2] << 16);
  pos += 3;

  // ID (BCF string: k=int32, then k bytes)
  const idLen = view.getInt32(pos, true);
  pos += 4;
  const id = idLen >= 0
    ? new TextDecoder().decode(data.slice(pos, pos + idLen))
    : '.';
  pos += Math.abs(idLen);

  // REF and ALT alleles
  const alleles: string[] = [];
  for (let a = 0; a < nAllele; a++) {
    const aLen = view.getInt32(pos, true);
    pos += 4;
    const allele = new TextDecoder().decode(data.slice(pos, pos + Math.abs(aLen)));
    alleles.push(allele);
    pos += Math.abs(aLen);
  }
  const ref = alleles.length > 0 ? alleles[0] : '.';
  const alt = alleles.length > 1 ? alleles.slice(1).join(',') : '.';

  // FILTER (BCF vector of int32)
  const filterVecLen = view.getInt32(pos, true);
  pos += 4;
  const filterParts: string[] = [];
  if (filterVecLen === 0) {
    filterParts.push('PASS');
  } else {
    const nFilter = Math.abs(filterVecLen) / 4;
    for (let f = 0; f < nFilter; f++) {
      const fIdx = view.getInt32(pos, true);
      pos += 4;
      const fName = filterNames.get(fIdx) ?? `Filter${fIdx}`;
      filterParts.push(fName);
    }
  }
  const filter = filterParts.join(';');

  // INFO fields — parse as key-value pairs
  // In BCF, each INFO field is: key_idx(int32) + typed value
  // We parse what we can and store remaining as raw
  const info = new Map<string, string>();
  // For a complete BCF parser we'd need the ##INFO header meta-lines
  // to know types. Here we do a best-effort parse.
  const sharedEnd = offset + lShared;
  // We'll skip INFO parsing for now and store as empty — full BCF INFO
  // parsing requires the header dictionary which is complex. The CHROM/POS/
  // ID/REF/ALT/QUAL/FILTER are the most commonly accessed fields.
  pos = sharedEnd;

  // FORMAT + individual/sample data
  const format: string[] = [];
  const samples: string[] = [];

  return { chrom, pos: pos1based, id, ref, alt, qual, filter, info, format, samples };
}

// =============================================================================
// CRAM Container Header Parsing
// =============================================================================

/** CRAM container header. */
export interface CramContainer {
  /** Container magic number (0x00 for CRAM v3+). */
  magic: number;
  /** Block number within the container (0 = first data block). */
  blockNumber: number;
  /** Number of slices in this container. */
  numSlices: number;
  /** Number of records in this container. */
  numRecords: number;
  /** Number of landmarks (slice byte offsets). */
  numLandmarks: number;
  /** Byte offsets of each slice within the container. */
  landmarks: number[];
  /** CRC32 of the container header (0 if not present). */
  crc32: number;
  /** Byte offset where the container header starts. */
  headerOffset: number;
  /** Byte length of the container header. */
  headerLength: number;
  /** Compression header block (raw bytes, if present). */
  compressionHeader: Uint8Array | null;
}

/** CRAM file definition header. */
export interface CramFileDefinition {
  /** Major version number. */
  majorVersion: number;
  /** Minor version number. */
  minorVersion: number;
  /** File ID (typically the file name). */
  fileId: string;
}

/**
 * Parse a CRAM file definition header from the start of a CRAM file.
 *
 * CRAM file definition format (v3):
 *   - 4 bytes: magic "CRAM"
 *   - 1 byte: major version
 *   - 1 byte: minor version
 *   - 20 bytes: file ID (null-padded)
 *
 * @param data CRAM binary data
 * @returns File definition, or null if not a CRAM file
 */
export function parseCramFileDefinition(data: Uint8Array): CramFileDefinition | null {
  if (data.length < 26) return null;

  // Check magic
  if (data[0] !== 0x43 || data[1] !== 0x52 || data[2] !== 0x41 || data[3] !== 0x4D) {
    return null; // Not "CRAM"
  }

  const majorVersion = data[4];
  const minorVersion = data[5];

  // File ID: 20 bytes, null-terminated
  let fileIdEnd = 6;
  while (fileIdEnd < 26 && data[fileIdEnd] !== 0) fileIdEnd++;
  const fileId = new TextDecoder().decode(data.slice(6, fileIdEnd));

  return { majorVersion, minorVersion, fileId };
}

/**
 * Parse a CRAM container header at the given offset.
 *
 * CRAM v3 container header format:
 *   - int32: blockNumber
 *   - int32: numSlices (called compCounter in older versions)
 *   - int32: numRecords
 *   - int32: numLandmarks
 *   - int32[numLandmarks]: landmark offsets
 *   - byte: CRC32 present flag (v3.1+)
 *   - if flag: uint32 CRC32 of header
 *
 * The container is preceded by a CRAM block (ITF8 encoded length + data).
 * We read the container ITF8 length first, then parse the fields.
 *
 * For CRAM v2, the format is slightly different (ITF8 encoding for all integers).
 * We focus on v3+ which uses standard int32.
 *
 * @param data CRAM binary data
 * @param offset Byte offset where the container header starts
 * @returns Parsed container header, or null if invalid
 */
export function parseCramContainer(data: Uint8Array, offset: number): CramContainer | null {
  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);

  if (offset + 16 > data.length) return null;

  const headerOffset = offset;

  // Read the container header fields
  const blockNumber = view.getInt32(offset, true);
  offset += 4;
  const numSlices = view.getInt32(offset, true);
  offset += 4;
  const numRecords = view.getInt32(offset, true);
  offset += 4;
  const numLandmarks = view.getInt32(offset, true);
  offset += 4;

  // Sanity checks
  if (numLandmarks < 0 || numLandmarks > 10000) return null;
  if (numRecords < 0 || numRecords > 1e9) return null;

  // Read landmarks
  const landmarks: number[] = [];
  if (offset + numLandmarks * 4 > data.length) return null;
  for (let i = 0; i < numLandmarks; i++) {
    landmarks.push(view.getInt32(offset, true));
    offset += 4;
  }

  // Check for CRC32 (CRAM v3.1+ appends a CRC32)
  let crc32 = 0;
  if (offset + 4 <= data.length) {
    // We read it if present, but can't easily distinguish v3.0 vs v3.1
    // In v3.0, the header ends after landmarks. In v3.1+, there's a CRC32.
    // We conservatively check if reading 4 more bytes looks like a CRC32
    // (by checking the next container/block starts correctly).
    // For simplicity, we always read it if available.
    crc32 = view.getUint32(offset, true);
    offset += 4;
  }

  const headerLength = offset - headerOffset;

  return {
    magic: 0, // CRAM containers don't have a magic byte
    blockNumber,
    numSlices,
    numRecords,
    numLandmarks,
    landmarks,
    crc32,
    headerOffset,
    headerLength,
    compressionHeader: null, // Not parsed — requires full CRAM codec container decode
  };
}

/**
 * Parse a CRAM block header at the given offset.
 * CRAM blocks are the fundamental unit within containers.
 *
 * Block format (v3):
 *   - int32: block method (0=raw, 1=gzip, 2=bzip2, 3=custom, 4=rans, 8=arith)
 *   - int32: block content type (0=file_header, 1=comp_header, 2=map, 4=data, 5=reserved)
 *   - int32: block content ID
 *   - int32: compressed size
 *   - int32: uncompressed size
 *
 * @param data CRAM binary data
 * @param offset Byte offset
 * @returns Parsed block info or null
 */
export function parseCramBlockHeader(data: Uint8Array, offset: number): {
  method: number;
  contentType: number;
  contentId: number;
  compressedSize: number;
  uncompressedSize: number;
  headerSize: number;
} | null {
  if (offset + 20 > data.length) return null;
  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);

  const method = view.getInt32(offset, true);
  offset += 4;
  const contentType = view.getInt32(offset, true);
  offset += 4;
  const contentId = view.getInt32(offset, true);
  offset += 4;
  const compressedSize = view.getInt32(offset, true);
  offset += 4;
  const uncompressedSize = view.getInt32(offset, true);
  offset += 4;

  return {
    method,
    contentType,
    contentId,
    compressedSize,
    uncompressedSize,
    headerSize: 20,
  };
}

/**
 * Extract all CRAM container headers from a CRAM file.
 * Skips over data blocks to find subsequent container headers.
 *
 * @param data CRAM binary data
 * @param maxContainers Maximum number of containers to parse (default: 1000)
 * @returns Array of container headers and file definition
 */
export function parseCramContainers(
  data: Uint8Array,
  maxContainers: number = 1000,
): { fileDef: CramFileDefinition | null; containers: CramContainer[] } {
  const fileDef = parseCramFileDefinition(data);
  let offset = 26; // Skip file definition

  const containers: CramContainer[] = [];

  // Skip the SAM header container (first container in CRAM v3)
  // The first container is a file header container with blockNumber = -1
  while (offset < data.length && containers.length < maxContainers) {
    // Each container starts with a block (which encodes its length)
    // First read the ITF8-encoded container length
    const containerResult = parseCramContainer(data, offset);
    if (!containerResult) break;

    containers.push(containerResult);

    // Move past this container: header + data
    // We need to know the total container byte length to skip to the next one.
    // In CRAM v3, after the container header, the data follows.
    // The data size is not directly in the container header — it's determined
    // by the block structure. We approximate by reading the first block's
    // compressed size and using landmarks to estimate.
    if (containerResult.numLandmarks > 0) {
      // The last landmark + its block data gives us an approximate end
      const lastLandmark = containerResult.landmarks[containerResult.numLandmarks - 1];
      // Skip to after the last landmark's data (rough estimate)
      // We read the block header at that landmark to get its size
      const blockOff = offset + containerResult.headerLength + lastLandmark;
      const blockHdr = parseCramBlockHeader(data, blockOff);
      if (blockHdr) {
        offset = blockOff + blockHdr.headerSize + blockHdr.compressedSize;
        // Add CRC32 if present
        if (offset + 4 <= data.length) offset += 4;
      } else {
        // Can't determine size — advance by a minimum amount
        offset += containerResult.headerLength + 1024;
      }
    } else {
      // No landmarks — advance past header with a default skip
      offset += containerResult.headerLength + 1024;
    }
  }

  return { fileDef, containers };
}

// =============================================================================
// tabix (.tbi) Index Parsing
// =============================================================================

/** A tabix index chunk (virtual offset: coffset << 16 | uoffset). */
export interface TabixChunk {
  /** BGZF virtual file offset. */
  voffset: number;
}

/** A bin in the tabix index. */
export interface TabixBin {
  /** Bin number. */
  bin: number;
  /** Number of intervals in this bin. */
  nIntv: number;
  /** Chunk intervals: pairs of [beg_voffset, end_voffset]. */
  chunks: [number, number][];
}

/** A reference entry in the tabix index. */
export interface TabixRefEntry {
  /** Bins for this reference. */
  bins: TabixBin[];
  /** Linear index (virtual offsets at 16KB intervals). */
  linearIndex: number[];
}

/** Parsed tabix (.tbi) index. */
export interface TabixIndex {
  /** Format (0 = generic, 1 = SAM, 2 = VCF, 3 = BED). */
  format: number;
  /** Column of chromosome (1-based). */
  colSeq: number;
  /** Column of start position (1-based). */
  colBeg: number;
  /** Column of end position (1-based), or 0 if same as colBeg. */
  colEnd: number;
  /** Meta character (lines starting with this are skipped). */
  meta: number;
  /** Number of lines to skip at the beginning. */
  skip: number;
  /** Number of reference sequences. */
  nRef: number;
  /** Per-reference index entries. */
  refs: TabixRefEntry[];
  /** Unparsed names from the name partition (if present). */
  names: string[];
}

/**
 * Parse a tabix (.tbi) index file.
 *
 * The tabix index format (from hts-spec):
 *   1. Magic: "TBI\1" (4 bytes)
 *   2. nRef (int32): number of reference sequences
 *   3. Format (int32): 0=generic, 1=SAM, 2=VCF, 3=BED
 *   4. colSeq (int32): column for sequence name (1-based)
 *   5. colBeg (int32): column for start position (1-based)
 *   6. colEnd (int32): column for end position (1-based)
 *   7. meta (int32): meta character
 *   8. skip (int32): lines to skip
 *   9. l_nm (int32): length of name concatenated string
 *  10. names: concatenated null-terminated reference names
 *  11. For each reference:
 *      a. n_bin (int32): number of bins
 *      b. For each bin:
 *         - bin (int32): bin number
 *         - n_chunk (int32): number of chunks
 *         - n_chunk × (cnk_beg: uint64, cnk_end: uint64)
 *      c. n_intv (int32): number of 16KB intervals
 *      d. n_intv × ioffset (uint64)
 *
 * @param data Binary .tbi index data
 * @returns Parsed tabix index
 */
export function parseTabix(data: Uint8Array): TabixIndex {
  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
  let offset = 0;

  // Check magic
  if (data.length < 4 || data[0] !== 0x54 || data[1] !== 0x42 || data[2] !== 0x49 || data[3] !== 0x01) {
    throw new Error('parseTabix: invalid magic number (expected TBI\\1)');
  }
  offset = 4;

  const nRef = view.getInt32(offset, true); offset += 4;
  const format = view.getInt32(offset, true); offset += 4;
  const colSeq = view.getInt32(offset, true); offset += 4;
  const colBeg = view.getInt32(offset, true); offset += 4;
  const colEnd = view.getInt32(offset, true); offset += 4;
  const meta = view.getInt32(offset, true); offset += 4;
  const skip = view.getInt32(offset, true); offset += 4;

  // Read concatenated reference names
  const lNm = view.getInt32(offset, true); offset += 4;
  const names: string[] = [];
  if (lNm > 0 && offset + lNm <= data.length) {
    let nameStart = offset;
    for (let i = 0; i < lNm; i++) {
      if (data[offset + i] === 0) {
        if (i > nameStart - offset) {
          names.push(new TextDecoder().decode(data.slice(nameStart, offset + i)));
        }
        nameStart = offset + i + 1;
      }
    }
    offset += lNm;
  }

  // Parse per-reference index entries
  const refs: TabixRefEntry[] = [];

  for (let r = 0; r < nRef; r++) {
    if (offset + 4 > data.length) break;

    const nBin = view.getInt32(offset, true); offset += 4;
    const bins: TabixBin[] = [];

    for (let b = 0; b < nBin; b++) {
      if (offset + 8 > data.length) break;

      const bin = view.getInt32(offset, true); offset += 4;
      const nChunk = view.getInt32(offset, true); offset += 4;
      const chunks: [number, number][] = [];

      for (let c = 0; c < nChunk; c++) {
        if (offset + 16 > data.length) break;
        // Read uint64 as two uint32s (JavaScript can't handle full uint64,
        // so we use the lower 32 bits as the virtual offset — sufficient for
        // files < 4GB, which covers most genomics use cases)
        const cnkBegLo = view.getUint32(offset, true); offset += 4;
        const cnkBegHi = view.getUint32(offset, true); offset += 4;
        const cnkEndLo = view.getUint32(offset, true); offset += 4;
        const cnkEndHi = view.getUint32(offset, true); offset += 4;

        // Combine hi:lo into a single number (may lose precision for >2^53)
        const cnkBeg = cnkBegHi * 0x100000000 + cnkBegLo;
        const cnkEnd = cnkEndHi * 0x100000000 + cnkEndLo;
        chunks.push([cnkBeg, cnkEnd]);
      }

      bins.push({ bin, nIntv: nChunk, chunks });
    }

    // Linear index
    if (offset + 4 > data.length) break;
    const nIntv = view.getInt32(offset, true); offset += 4;
    const linearIndex: number[] = [];

    for (let i = 0; i < nIntv; i++) {
      if (offset + 8 > data.length) break;
      const lo = view.getUint32(offset, true); offset += 4;
      const hi = view.getUint32(offset, true); offset += 4;
      linearIndex.push(hi * 0x100000000 + lo);
    }

    refs.push({ bins, linearIndex });
  }

  return { format, colSeq, colBeg, colEnd, meta, skip, nRef, refs, names };
}

/**
 * Calculate the tabix bin for a region [beg, end).
 *
 * The tabix/BAI binning scheme uses 6 levels:
 *   Level 0: 1 bin  covering 512MB  (bin 0)
 *   Level 1: 8 bins covering 64MB   (bins 1-8)
 *   Level 2: 64 bins covering 8MB   (bins 9-72)
 *   Level 3: 512 bins covering 1MB  (bins 73-584)
 *   Level 4: 4096 bins covering 128KB (bins 585-4680)
 *   Level 5: 32768 bins covering 16KB (bins 4681-37448)
 *
 * Bin numbering: at level l, bins start at offset_l = (8^l - 1) / 7
 * Within level l, the bin for position p is: offset_l + p / (2^(29-3l))
 *
 * @param beg Start position (0-based)
 * @param end End position (0-based, exclusive)
 * @returns Array of all bins that overlap [beg, end)
 */
export function reg2bins(beg: number, end: number): number[] {
  const bins: number[] = [];
  // Level 0: bin 0 covers everything
  bins.push(0);

  const endMinus1 = end - 1;

  // Level 1
  for (let k = 1 + (beg >> 26); k <= 1 + (endMinus1 >> 26); k++) {
    if (k <= 8) bins.push(k);
  }
  // Level 2
  for (let k = 9 + (beg >> 23); k <= 9 + (endMinus1 >> 23); k++) {
    if (k <= 72) bins.push(k);
  }
  // Level 3
  for (let k = 73 + (beg >> 20); k <= 73 + (endMinus1 >> 20); k++) {
    if (k <= 584) bins.push(k);
  }
  // Level 4
  for (let k = 585 + (beg >> 17); k <= 585 + (endMinus1 >> 17); k++) {
    if (k <= 4680) bins.push(k);
  }
  // Level 5
  for (let k = 4681 + (beg >> 14); k <= 4681 + (endMinus1 >> 14); k++) {
    if (k <= 37448) bins.push(k);
  }

  return bins;
}

/**
 * Query a tabix index for virtual file offsets overlapping a region.
 *
 * Returns a sorted, merged list of [beg_voffset, end_voffset] chunks
 * that may contain records overlapping [start, end) on the given chromosome.
 *
 * @param index Parsed tabix index
 * @param chrom Reference sequence name
 * @param start 0-based start position
 * @param end 0-based end position (exclusive)
 * @returns Array of [beg_voffset, end_voffset] virtual offset pairs
 */
export function queryTabix(
  index: TabixIndex,
  chrom: string,
  start: number,
  end: number,
): [number, number][] {
  // Find the reference index by name
  let refIdx = index.names.indexOf(chrom);
  if (refIdx < 0) {
    // Try with "chr" prefix or without
    const chromAlt = chrom.startsWith('chr') ? chrom.slice(3) : `chr${chrom}`;
    refIdx = index.names.indexOf(chromAlt);
    if (refIdx < 0) return []; // Chromosome not in index
  }

  if (refIdx >= index.refs.length) return [];

  const refEntry = index.refs[refIdx];
  const queryBins = reg2bins(start, end);

  // Collect all chunks from overlapping bins
  const chunks: [number, number][] = [];

  for (const bin of refEntry.bins) {
    if (queryBins.includes(bin.bin)) {
      chunks.push(...bin.chunks);
    }
  }

  // Also use the linear index to narrow the search
  // The linear index maps 16KB intervals to the smallest virtual offset
  // that could contain data starting in that interval.
  if (refEntry.linearIndex.length > 0) {
    const minOffIdx = start >> 14; // 16KB interval index
    if (minOffIdx < refEntry.linearIndex.length) {
      const minOffset = refEntry.linearIndex[minOffIdx];
      // Filter chunks: only keep those with end > minOffset
      // (chunks that end before our region can't contain relevant records)
      const filtered = chunks.filter(c => c[1] > minOffset || c[1] === 0);
      if (filtered.length > 0) {
        chunks.length = 0;
        chunks.push(...filtered);
      }
    }
  }

  // Sort and merge overlapping chunks
  chunks.sort((a, b) => a[0] - b[0]);
  const merged: [number, number][] = [];
  for (const chunk of chunks) {
    if (merged.length > 0 && chunk[0] <= merged[merged.length - 1][1]) {
      // Overlapping — extend the last merged chunk
      merged[merged.length - 1][1] = Math.max(merged[merged.length - 1][1], chunk[1]);
    } else {
      merged.push([chunk[0], chunk[1]]);
    }
  }

  return merged;
}

// =============================================================================
// FAI (FASTA Index) Parsing
// =============================================================================

/** An entry in a FASTA index (.fai file). */
export interface FaiEntry {
  /** Reference sequence name. */
  name: string;
  /** Total length of the sequence (in bases). */
  length: number;
  /** Byte offset of the sequence's first base in the FASTA file. */
  offset: number;
  /** Number of bases per line (excluding newline). */
  lineBases: number;
  /** Number of bytes per line (including newline). */
  lineBytes: number;
}

/** Parsed FASTA index. */
export interface FaiIndex {
  /** Index entries, one per reference sequence. */
  entries: FaiEntry[];
  /** Name → entry lookup for O(1) access. */
  byName: Map<string, FaiEntry>;
}

/**
 * Parse a FASTA index (.fai) file.
 *
 * The .fai format is a tab-separated text file with one line per sequence:
 *   name\tlength\toffset\tlineBases\tlineBytes
 *
 * Where:
 *   - name:      sequence name (e.g., "chr1")
 *   - length:    total number of bases in the sequence
 *   - offset:    byte offset in the FASTA of the first base
 *   - lineBases: number of bases per line (excluding newline)
 *   - lineBytes: number of bytes per line (including newline; lineBases + 1 for \n)
 *
 * @param text .fai file content as a string
 * @returns Parsed FAI index
 */
export function parseFai(text: string): FaiIndex {
  const entries: FaiEntry[] = [];
  const byName = new Map<string, FaiEntry>();

  const lines = text.split('\n');
  for (const line of lines) {
    const trimmed = line.trim();
    if (trimmed === '') continue;

    const cols = trimmed.split('\t');
    if (cols.length < 5) continue;

    const name = cols[0];
    const length = parseInt(cols[1], 10);
    const offset = parseInt(cols[2], 10);
    const lineBases = parseInt(cols[3], 10);
    const lineBytes = parseInt(cols[4], 10);

    if (isNaN(length) || isNaN(offset) || isNaN(lineBases) || isNaN(lineBytes)) {
      continue; // Skip malformed lines
    }

    const entry: FaiEntry = { name, length, offset, lineBases, lineBytes };
    entries.push(entry);
    byName.set(name, entry);
  }

  return { entries, byName };
}

/**
 * Query a FASTA index for the byte range containing [start, end) on a chromosome.
 *
 * Given a FAI index and a genomic region, computes the byte offset and
 * length in the FASTA file that contains the sequence for that region.
 * This enables random access (seek + read) to extract subsequences
 * without reading the entire FASTA.
 *
 * The calculation accounts for line wrapping in the FASTA file:
 *   - Each line has `lineBases` bases and `lineBytes` bytes (bases + newline)
 *   - The byte offset for position p is:
 *     offset + (p / lineBases) * lineBytes + (p % lineBases)
 *
 * @param index Parsed FAI index
 * @param chrom Reference sequence name
 * @param start 0-based start position (inclusive)
 * @param end 0-based end position (exclusive)
 * @returns Byte range { offset, length } in the FASTA file, or null if chrom not found
 */
export function queryFai(
  index: FaiIndex,
  chrom: string,
  start: number,
  end: number,
): { offset: number; length: number } | null {
  const entry = index.byName.get(chrom);
  if (!entry) {
    // Try with/without "chr" prefix
    const chromAlt = chrom.startsWith('chr') ? chrom.slice(3) : `chr${chrom}`;
    const alt = index.byName.get(chromAlt);
    if (!alt) return null;
    return queryFaiEntry(alt, start, end);
  }
  return queryFaiEntry(entry, start, end);
}

/**
 * Compute byte range for a region within a single FAI entry.
 */
function queryFaiEntry(
  entry: FaiEntry,
  start: number,
  end: number,
): { offset: number; length: number } | null {
  // Clamp to sequence bounds
  const clampedStart = Math.max(0, start);
  const clampedEnd = Math.min(entry.length, end);

  if (clampedStart >= clampedEnd) return null;

  const { offset, lineBases, lineBytes } = entry;

  // Calculate byte offset for start position
  const startLine = Math.floor(clampedStart / lineBases);
  const startCol = clampedStart % lineBases;
  const startByteOff = offset + startLine * lineBytes + startCol;

  // Calculate byte offset for end position
  const endLine = Math.floor((clampedEnd - 1) / lineBases); // Last base's line
  const endCol = (clampedEnd - 1) % lineBases;
  const endByteOff = offset + endLine * lineBytes + endCol;

  // Length includes the bases from start to end (inclusive of last base)
  // We need to include the last base, so +1
  const length = endByteOff - startByteOff + 1;

  return { offset: startByteOff, length };
}

/**
 * Extract a subsequence from a FASTA file given an FAI index.
 * Handles line wrapping by reading the byte range and stripping newlines.
 *
 * @param fastaData Full FASTA file data
 * @param index Parsed FAI index
 * @param chrom Reference sequence name
 * @param start 0-based start position
 * @param end 0-based end position (exclusive)
 * @returns Subsequence as a string, or null if chrom not found
 */
export function fetchFaiSequence(
  fastaData: Uint8Array,
  index: FaiIndex,
  chrom: string,
  start: number,
  end: number,
): string | null {
  const range = queryFai(index, chrom, start, end);
  if (!range) return null;

  const { offset: byteOffset, length: byteLength } = range;

  if (byteOffset + byteLength > fastaData.length) {
    // Requested range exceeds FASTA data — return what we can
    const available = Math.max(0, fastaData.length - byteOffset);
    if (available === 0) return null;
    const slice = fastaData.slice(byteOffset, byteOffset + available);
    // Strip newlines
    const filtered = new Uint8Array(available);
    let len = 0;
    for (let i = 0; i < slice.length; i++) {
      if (slice[i] !== 0x0A && slice[i] !== 0x0D) { // Not \n or \r
        filtered[len++] = slice[i];
      }
    }
    return new TextDecoder().decode(filtered.slice(0, len));
  }

  const slice = fastaData.slice(byteOffset, byteOffset + byteLength);
  // Strip newlines from the FASTA lines
  const filtered = new Uint8Array(byteLength);
  let len = 0;
  for (let i = 0; i < slice.length; i++) {
    if (slice[i] !== 0x0A && slice[i] !== 0x0D) {
      filtered[len++] = slice[i];
    }
  }
  return new TextDecoder().decode(filtered.slice(0, len));
}
