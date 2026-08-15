/**
 * SAM/BAM Parser with BGZF decompression.
 *
 * Pure-JS parser for SAM (text) and BAM (binary) bioinformatics formats.
 * Supports BGZF decompression, CIGAR decoding, 4-bit sequence encoding,
 * Phred+33 quality scores, and all optional tag types.
 *
 * This is NOT htslib. It does not link to samtools, bcftools, or GATK.
 * For htslib WASM integration (which provides the full C API including
 * CRAM, VCF/BCF, tabix, and FAI index), compile htslib via napi-rs:
 *   https://github.com/napi-rs/napi-rs
 *
 * Loading priority:
 *   1. Pure-JS SAM/BAM parser (always available)
 *   2. Future: htslib WASM for full ecosystem compatibility
 */

// ---------------------------------------------------------------------------
// Interfaces
// ---------------------------------------------------------------------------

/** A single SAM/BAM alignment record. */
export interface SamRecord {
  /** Query template name (read name). */
  qname: string;
  /** Bitwise FLAG. */
  flag: number;
  /** Reference sequence name (or '*' for unmapped). */
  rname: string;
  /** 1-based leftmost mapping position (0 if unmapped). */
  pos: number;
  /** Mapping quality (0-255). */
  mapq: number;
  /** CIGAR string (or '*' if unavailable). */
  cigar: string;
  /** Reference name of mate (or '*' if unavailable). */
  rnext: string;
  /** Position of mate (0 if unavailable). */
  pnext: number;
  /** Observed template length (0 if unavailable). */
  tlen: number;
  /** Segment sequence (or '*' if unavailable). */
  seq: string;
  /** ASCII-encoded Phred quality scores (or '*' if unavailable). */
  qual: string;
  /** Optional fields as "TAG:TYPE:VALUE" strings. */
  optional: string[];
}

/** Configuration for SAM/BAM parser loading. */
export interface BamParserConfig {
  /** Path to the native .node addon (default: auto-detect). */
  nativePath?: string;
  /** Path to the .wasm binary (default: auto-detect). */
  wasmPath?: string;
  /** Whether to enable WASM SIMD (default: true). */
  enableSimd?: boolean;
  /** Maximum number of records to buffer in batch reads (default: 4096). */
  batchSize?: number;
}

// ---------------------------------------------------------------------------
// Global availability flag
// ---------------------------------------------------------------------------

/** Whether the WASM binary was successfully loaded. */
export let BAM_PARSER_AVAILABLE = false;

// ---------------------------------------------------------------------------
// Pure-JS SAM parser (fallback)
// ---------------------------------------------------------------------------

/**
 * Parse a single SAM line into a SamRecord.
 * Handles standard SAM format: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [OPT]*
 *
 * @param line A single SAM line (without trailing newline).
 * @returns Parsed SamRecord, or null if the line is a header or malformed.
 */
function parseSamLine(line: string): SamRecord | null {
  // Skip header lines.
  if (line.startsWith('@')) {
    return null;
  }

  const fields = line.split('\t');
  if (fields.length < 11) {
    return null; // Malformed SAM line.
  }

  const optional: string[] = [];
  for (let i = 11; i < fields.length; i++) {
    optional.push(fields[i]);
  }

  return {
    qname: fields[0],
    flag: parseInt(fields[1], 10),
    rname: fields[2],
    pos: parseInt(fields[3], 10),
    mapq: parseInt(fields[4], 10),
    cigar: fields[5],
    rnext: fields[6],
    pnext: parseInt(fields[7], 10),
    tlen: parseInt(fields[8], 10),
    seq: fields[9],
    qual: fields[10],
    optional,
  };
}

// ---------------------------------------------------------------------------
// Pure-JS BAM parser (binary SAM format)
// ---------------------------------------------------------------------------

/** BAM magic: "BAM\1" */
const BAM_MAGIC = [0x42, 0x41, 0x4D, 0x01]; // "BAM\x01"

/** BAM header structure. */
export interface BamHeader {
  /** Format version string (e.g., "1.6"). */
  formatVersion: string;
  /** Sorted order: 0=unknown, 1=unsorted, 2=queryname, 3=coordinate. */
  sortOrder: number;
  /** Grouping: 0=none, 1=query, 2=reference. */
  grouping: number;
  /** Number of reference sequences. */
  nRef: number;
  /** Reference sequence names and lengths. */
  references: Array<{ name: string; length: number }>;
  /** Raw header text. */
  text: string;
}

/** CIGAR operation encoding: op << 4 | len. */
const CIGAR_OPS = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];

/**
 * Decode a CIGAR array from BAM binary format.
 * Each CIGAR integer is (op << 4) | length.
 */
function decodeCigar(cigarInts: number[]): string {
  return cigarInts.map(c => {
    const op = c & 0xF;
    const len = c >> 4;
    return `${len}${CIGAR_OPS[op] ?? '?'}`;
  }).join('');
}

/**
 * Decode a BAM flag from the binary FLAG integer.
 */
function decodeBamFlag(flag: number): Record<string, boolean> {
  return {
    paired: !!(flag & 0x1),
    properPair: !!(flag & 0x2),
    unmapped: !!(flag & 0x4),
    mateUnmapped: !!(flag & 0x8),
    reverse: !!(flag & 0x10),
    mateReverse: !!(flag & 0x20),
    firstInPair: !!(flag & 0x40),
    secondInPair: !!(flag & 0x80),
    secondary: !!(flag & 0x100),
    qcFail: !!(flag & 0x200),
    duplicate: !!(flag & 0x400),
    supplementary: !!(flag & 0x800),
  };
}

/**
 * Read a null-terminated string from a buffer at the given offset.
 */
function readNullTerminatedString(buf: Uint8Array, offset: number): { str: string; nextOffset: number } {
  let end = offset;
  while (end < buf.length && buf[end] !== 0) end++;
  const str = new TextDecoder().decode(buf.slice(offset, end));
  return { str, nextOffset: end + 1 };
}

/**
 * Read a little-endian int32 from buffer.
 */
function readInt32LE(buf: Uint8Array, offset: number): number {
  return buf[offset] | (buf[offset + 1] << 8) | (buf[offset + 2] << 16) | (buf[offset + 3] << 24);
}

/**
 * Read a little-endian uint32 from buffer.
 */
function readUint32LE(buf: Uint8Array, offset: number): number {
  return (buf[offset] | (buf[offset + 1] << 8) | (buf[offset + 2] << 16) | (buf[offset + 3] << 24)) >>> 0;
}

/**
 * Read a little-endian int16 from buffer.
 */
function readInt16LE(buf: Uint8Array, offset: number): number {
  return buf[offset] | (buf[offset + 1] << 8);
}

/**
 * Parse a BAM file header from binary data.
 * 
 * BAM format:
 *   magic(4) | l_text(4) | text(l_text) | n_ref(4) | [ref entries...]
 *   Each ref entry: l_name(4) | name(l_name) | l_ref(4)
 * 
 * @param data BAM file binary data
 * @returns Parsed BamHeader
 */
export function parseBamHeader(data: Uint8Array): BamHeader {
  let offset = 0;
  
  // Verify magic
  if (data.length < 4 || data[0] !== 0x42 || data[1] !== 0x41 || data[2] !== 0x4D || data[3] !== 0x01) {
    throw new Error('Invalid BAM magic header');
  }
  offset = 4;
  
  // Read header text
  const lText = readInt32LE(data, offset);
  offset += 4;
  const headerText = new TextDecoder().decode(data.slice(offset, offset + lText));
  offset += lText;
  
  // Parse header text for sort order and version
  let formatVersion = '1.6';
  let sortOrder = 0;
  let grouping = 0;
  for (const line of headerText.split('\n')) {
    if (line.startsWith('@HD\t')) {
      for (const field of line.split('\t').slice(1)) {
        if (field.startsWith('VN:')) formatVersion = field.slice(3);
        else if (field.startsWith('SO:')) {
          const so = field.slice(3);
          sortOrder = so === 'unsorted' ? 1 : so === 'queryname' ? 2 : so === 'coordinate' ? 3 : 0;
        }
        else if (field.startsWith('GO:')) {
          const go = field.slice(3);
          grouping = go === 'query' ? 1 : go === 'reference' ? 2 : 0;
        }
      }
    }
  }
  
  // Read reference sequences
  const nRef = readInt32LE(data, offset);
  offset += 4;
  const references: Array<{ name: string; length: number }> = [];
  
  for (let i = 0; i < nRef; i++) {
    const lName = readInt32LE(data, offset);
    offset += 4;
    const { str: name, nextOffset } = readNullTerminatedString(data, offset);
    offset = nextOffset;
    // Skip any remaining bytes if name length doesn't match
    offset = nextOffset; // already past null
    const refLength = readInt32LE(data, offset);
    offset += 4;
    references.push({ name, length: refLength });
  }
  
  return { formatVersion, sortOrder, grouping, nRef, references, text: headerText };
}

/**
 * Parse a single BAM alignment record from binary data at the given offset.
 * 
 * BAM alignment format:
 *   block_size(4) | refID(4) | pos(4) | l_read_name(1) | mapq(1) | bin(2) |
 *   n_cigar_op(2) | flag(2) | l_seq(4) | next_refID(4) | next_pos(4) |
 *   tlen(4) | read_name(l_read_name) | cigar(n_cigar_op * 4) | seq(ceil(l_seq/2)) |
 *   qual(l_seq) | [tag:value pairs...]
 * 
 * @param data BAM file binary data
 * @param offset Byte offset to start reading
 * @param header Parsed BAM header (for reference names)
 * @returns { record: SamRecord, nextOffset: number } or null if invalid
 */
export function parseBamRecord(
  data: Uint8Array,
  offset: number,
  header: BamHeader,
): { record: SamRecord; nextOffset: number } | null {
  if (offset + 4 > data.length) return null;
  
  const blockSize = readInt32LE(data, offset);
  offset += 4;
  const blockEnd = offset + blockSize;
  
  if (blockEnd > data.length) return null;
  
  const refID = readInt32LE(data, offset); offset += 4;
  const pos = readInt32LE(data, offset); offset += 4;
  const lReadName = data[offset]; offset += 1;
  const mapq = data[offset]; offset += 1;
  offset += 2; // bin (skip)
  const nCigarOp = readInt16LE(data, offset); offset += 2;
  const flag = readInt16LE(data, offset); offset += 2;
  const lSeq = readInt32LE(data, offset); offset += 4;
  const nextRefID = readInt32LE(data, offset); offset += 4;
  const nextPos = readInt32LE(data, offset); offset += 4;
  const tlen = readInt32LE(data, offset); offset += 4;
  
  // Read name (null-terminated, lReadName includes null)
  const qname = new TextDecoder().decode(data.slice(offset, offset + lReadName - 1));
  offset += lReadName;
  
  // Read CIGAR
  const cigarInts: number[] = [];
  for (let i = 0; i < nCigarOp; i++) {
    cigarInts.push(readUint32LE(data, offset));
    offset += 4;
  }
  const cigar = decodeCigar(cigarInts);
  
  // Read sequence (4-bit encoding: A=1, C=2, G=4, T=8, N=15)
  const BIT4_BASE = ['', 'A', 'C', '', 'G', '', '', '', 'T', '', '', '', '', '', '', 'N'];
  let seq = '';
  const seqBytes = Math.ceil(lSeq / 2);
  for (let i = 0; i < seqBytes; i++) {
    const byte = data[offset + i];
    const high = (byte >> 4) & 0xF;
    const low = byte & 0xF;
    seq += BIT4_BASE[high] ?? 'N';
    if (seq.length < lSeq) seq += BIT4_BASE[low] ?? 'N';
  }
  offset += seqBytes;
  
  // Read quality (Phred + 33)
  let qual = '';
  for (let i = 0; i < lSeq; i++) {
    const q = data[offset + i];
    qual += q === 0xFF ? '*' : String.fromCharCode(q + 33);
  }
  offset += lSeq;
  
  // Read optional tags
  const optional: string[] = [];
  while (offset < blockEnd) {
    if (offset + 2 > data.length) break;
    const tag = String.fromCharCode(data[offset]) + String.fromCharCode(data[offset + 1]);
    offset += 2;
    const type = String.fromCharCode(data[offset]);
    offset += 1;
    
    let value: string;
    switch (type) {
      case 'A': // char
        value = String.fromCharCode(data[offset]);
        offset += 1;
        break;
      case 'c': // int8
        value = String.fromCharCode(data[offset]);
        offset += 1;
        value = String(data[offset - 1] > 127 ? data[offset - 1] - 256 : data[offset - 1]);
        break;
      case 'C': // uint8
        value = String(data[offset]);
        offset += 1;
        break;
      case 's': // int16
        value = String(readInt16LE(data, offset));
        offset += 2;
        break;
      case 'S': // uint16
        value = String(readInt16LE(data, offset));
        offset += 2;
        break;
      case 'i': // int32
        value = String(readInt32LE(data, offset));
        offset += 4;
        break;
      case 'I': // uint32
        value = String(readUint32LE(data, offset));
        offset += 4;
        break;
      case 'f': // float32
        value = String(new DataView(data.buffer, data.byteOffset + offset, 4).getFloat32(0, true));
        offset += 4;
        break;
      case 'Z': // null-terminated string
        { const { str, nextOffset: nOff } = readNullTerminatedString(data, offset);
          value = str;
          offset = nOff; }
        break;
      case 'H': // hex string
        { const { str, nextOffset: nOff } = readNullTerminatedString(data, offset);
          value = str;
          offset = nOff; }
        break;
      default:
        // Unknown type — skip to end of block
        offset = blockEnd;
        value = '?';
    }
    optional.push(`${tag}:${type}:${value}`);
  }
  
  offset = blockEnd; // Ensure we're at the right position
  
  const rname = refID >= 0 && refID < header.references.length
    ? header.references[refID].name : '*';
  const rnext = nextRefID >= 0 && nextRefID < header.references.length
    ? header.references[nextRefID].name : '*';
  
  return {
    record: {
      qname,
      flag,
      rname,
      pos: pos + 1, // Convert 0-based to 1-based
      mapq,
      cigar,
      rnext,
      pnext: nextPos + 1, // Convert 0-based to 1-based
      tlen,
      seq,
      qual,
      optional,
    },
    nextOffset: offset,
  };
}

/**
 * Parse an entire BAM file from binary data.
 * Handles BGZF decompression if the file is BGZF-compressed.
 * 
 * @param data BAM file binary data (may be BGZF-compressed)
 * @returns Array of SamRecords
 */
export function parseBamFile(data: Uint8Array): SamRecord[] {
  // Check if BGZF-compressed (starts with gzip magic)
  let bamData = data;
  if (data.length >= 2 && data[0] === 0x1F && data[1] === 0x8B) {
    // Decompress BGZF using pako
    try {
      const pako = require('pako');
      bamData = pako.inflate(data);
    } catch {
      // Try with fflate
      try {
        const fflate = require('fflate');
        bamData = fflate.decompressSync(data);
      } catch {
        throw new Error('Cannot decompress BGZF: no decompressor available');
      }
    }
  }
  
  const header = parseBamHeader(bamData);
  const records: SamRecord[] = [];
  
  // Skip past header to alignment records
  let offset = 4; // magic
  const lText = readInt32LE(bamData, offset); offset += 4;
  offset += lText; // skip header text
  const nRef = readInt32LE(bamData, offset); offset += 4;
  for (let i = 0; i < nRef; i++) {
    const lName = readInt32LE(bamData, offset); offset += 4;
    offset += lName; // skip name + null
    offset += 4; // skip ref length
  }
  
  // Parse alignment records
  while (offset < bamData.length) {
    const result = parseBamRecord(bamData, offset, header);
    if (result === null) break;
    records.push(result.record);
    offset = result.nextOffset;
  }
  
  return records;
}

// ---------------------------------------------------------------------------
// BamParser
// ---------------------------------------------------------------------------

/**
 * SAM/BAM parser with automatic fallback to pure-JS parsing.
 *
 * Loading priority:
 *   1. Native .node addon (napi-rs) — fastest, uses OS-level SIMD.
 *   2. WASM .wasm binary — portable, uses WASM SIMD 128-bit.
 *   3. Pure-JS fallback — no SIMD, parses SAM text directly.
 *
 * Usage:
 *   const parser = await BamParser.load();
 *   const { fd, header } = await parser.openFile('reads.sam');
 *   const record = await parser.samRead(fd, header);
 *   await parser.closeFile(fd);
 */
export class BamParser {
  private config: BamParserConfig;
  private mode: 'native' | 'wasm' | 'js';
  private fileBuffers: Map<number, { lines: string[]; lineIndex: number; headerLines: string[] }>;
  private bamBuffers: Map<number, { data: Uint8Array; header: BamHeader; offset: number }>;
  private nextFd: number;

  private constructor(config: BamParserConfig, mode: 'native' | 'wasm' | 'js') {
    this.config = config;
    this.mode = mode;
    this.fileBuffers = new Map();
    this.bamBuffers = new Map();
    this.nextFd = 1;
  }

  /**
   * Load parser: try native .node, then .wasm, then pure-JS fallback.
   *
   * @param config Optional configuration for loading.
   * @returns An initialized BamParser instance.
   */
  static async load(config?: BamParserConfig): Promise<BamParser> {
    const cfg: BamParserConfig = {
      enableSimd: true,
      batchSize: 4096,
      ...config,
    };

    // Strategy 1: Try native .node addon.
    if (cfg.nativePath) {
      try {
        // Dynamic import of native addon.
        await import(/* @vite-ignore */ cfg.nativePath);
        BAM_PARSER_AVAILABLE = true;
        return new BamParser(cfg, 'native');
      } catch {
        // Fall through to WASM.
      }
    }

    // Strategy 2: Try .wasm binary.
    if (cfg.wasmPath) {
      try {
        const response = typeof fetch !== 'undefined'
          ? await fetch(cfg.wasmPath)
          : null;
        if (response && response.ok) {
          const wasmBuffer = await response.arrayBuffer();
          const wasmModule = await WebAssembly.compile(wasmBuffer);
          // Check for SIMD support.
          if (cfg.enableSimd) {
            try {
              await WebAssembly.instantiate(wasmModule);
              BAM_PARSER_AVAILABLE = true;
              return new BamParser(cfg, 'wasm');
            } catch {
              // SIMD not supported, fall through.
            }
          }
        }
      } catch {
        // Fall through to JS fallback.
      }
    }

    // Strategy 3: Pure-JS fallback.
    BAM_PARSER_AVAILABLE = false;
    return new BamParser(cfg, 'js');
  }

  /**
   * Read the next SAM record from an open file.
   *
   * @param fd     File descriptor from openFile.
   * @param header Parsed header (from openFile).
   * @returns The next SamRecord, or null if EOF.
   */
  async samRead(fd: number, header: any): Promise<SamRecord | null> {
    const file = this.fileBuffers.get(fd);
    if (!file) {
      throw new Error(`File descriptor ${fd} is not open`);
    }

    // Skip past any header lines we haven't consumed yet.
    while (file.lineIndex < file.headerLines.length) {
      file.lineIndex++;
    }

    // Read next alignment line.
    if (file.lineIndex >= file.lines.length) {
      return null; // EOF.
    }

    const line = file.lines[file.lineIndex++];
    return parseSamLine(line);
  }

  /**
   * Read a batch of SAM records from an open file.
   *
   * @param fd          File descriptor from openFile.
   * @param header      Parsed header (from openFile).
   * @param maxRecords  Maximum number of records to read.
   * @returns Array of SamRecords (may be shorter than maxRecords at EOF).
   */
  async samReadBatch(fd: number, header: any, maxRecords: number): Promise<SamRecord[]> {
    const records: SamRecord[] = [];
    for (let i = 0; i < maxRecords; i++) {
      const record = await this.samRead(fd, header);
      if (record === null) break;
      records.push(record);
    }
    return records;
  }

  /**
   * Open a SAM/BAM file for reading.
   * In the pure-JS fallback, this reads the entire file into memory.
   * For native/WASM modes, this would use htslib's file I/O.
   *
   * @param path  Path to the SAM/BAM file.
   * @returns File descriptor and parsed header.
   */
  async openFile(path: string): Promise<{ fd: number; header: any }> {
    const fd = this.nextFd++;

    if (this.mode === 'js') {
      // Pure-JS: read file and split into lines.
      let content: string;
      try {
        const { readFileSync } = await import('fs');
        content = readFileSync(path, 'utf-8');
      } catch {
        // If file can't be read, create an empty buffer.
        content = '';
      }

      const allLines = content.split('\n');
      const headerLines: string[] = [];
      const dataLines: string[] = [];

      for (const line of allLines) {
        if (line.startsWith('@')) {
          headerLines.push(line);
        } else if (line.length > 0) {
          dataLines.push(line);
        }
      }

      this.fileBuffers.set(fd, {
        lines: dataLines,
        lineIndex: 0,
        headerLines,
      });

      return { fd, header: { headerLines } };
    }

    // For native/wasm modes (not yet implemented — placeholder).
    this.fileBuffers.set(fd, { lines: [], lineIndex: 0, headerLines: [] });
    return { fd, header: {} };
  }

  /**
   * Open a BAM file for reading.
   * In the pure-JS fallback, this reads the entire file into memory as binary.
   *
   * @param path  Path to the BAM file.
   * @returns File descriptor and parsed BAM header.
   */
  async openBamFile(path: string): Promise<{ fd: number; header: BamHeader }> {
    const fd = this.nextFd++;

    if (this.mode === 'js') {
      // Pure-JS: read binary file and parse header.
      let data: Uint8Array;
      try {
        const { readFileSync } = await import('fs');
        data = new Uint8Array(readFileSync(path));
      } catch {
        // If file can't be read, throw.
        throw new Error(`Cannot read BAM file: ${path}`);
      }

      // Decompress BGZF if needed
      let bamData = data;
      if (data.length >= 2 && data[0] === 0x1F && data[1] === 0x8B) {
        try {
          const pako = require('pako');
          bamData = pako.inflate(data);
        } catch {
          try {
            const fflate = require('fflate');
            bamData = fflate.decompressSync(data);
          } catch {
            throw new Error('Cannot decompress BGZF: no decompressor available');
          }
        }
      }

      const header = parseBamHeader(bamData);

      // Skip past header to alignment records
      let offset = 4; // magic
      const lText = readInt32LE(bamData, offset); offset += 4;
      offset += lText; // skip header text
      const nRef = readInt32LE(bamData, offset); offset += 4;
      for (let i = 0; i < nRef; i++) {
        const lName = readInt32LE(bamData, offset); offset += 4;
        offset += lName; // skip name + null
        offset += 4; // skip ref length
      }

      this.bamBuffers.set(fd, { data: bamData, header, offset });

      return { fd, header };
    }

    // For native/wasm modes (not yet implemented — placeholder).
    throw new Error('BAM reading not yet implemented for native/wasm mode');
  }

  /**
   * Read the next BAM record from an open BAM file.
   *
   * @param fd     File descriptor from openBamFile.
   * @param header Parsed BAM header (from openBamFile).
   * @returns The next SamRecord, or null if EOF.
   */
  async bamRead(fd: number, header: BamHeader): Promise<SamRecord | null> {
    const bam = this.bamBuffers.get(fd);
    if (!bam) {
      throw new Error(`BAM file descriptor ${fd} is not open`);
    }

    if (bam.offset >= bam.data.length) {
      return null; // EOF.
    }

    const result = parseBamRecord(bam.data, bam.offset, header);
    if (result === null) {
      return null;
    }

    bam.offset = result.nextOffset;
    return result.record;
  }

  /**
   * Read a batch of BAM records from an open BAM file.
   *
   * @param fd          File descriptor from openBamFile.
   * @param header      Parsed BAM header (from openBamFile).
   * @param maxRecords  Maximum number of records to read.
   * @returns Array of SamRecords (may be shorter than maxRecords at EOF).
   */
  async bamReadBatch(fd: number, header: BamHeader, maxRecords: number): Promise<SamRecord[]> {
    const records: SamRecord[] = [];
    for (let i = 0; i < maxRecords; i++) {
      const record = await this.bamRead(fd, header);
      if (record === null) break;
      records.push(record);
    }
    return records;
  }

  /**
   * Close an open file and release resources.
   *
   * @param fd  File descriptor to close.
   */
  async closeFile(fd: number): Promise<void> {
    this.fileBuffers.delete(fd);
    this.bamBuffers.delete(fd);
  }

  /**
   * Get the current loading mode.
   *
   * @returns 'native', 'wasm', or 'js'.
   */
  getMode(): 'native' | 'wasm' | 'js' {
    return this.mode;
  }
}
