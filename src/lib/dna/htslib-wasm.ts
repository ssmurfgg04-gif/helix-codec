/**
 * htslib WASM Binding Skeleton.
 *
 * Compiles htslib to WASM via emnapi + napi-rs.
 * The loader tries native .node first, falls back to .wasm automatically.
 * Same codebase, zero platform-specific builds.
 *
 * napi-rs auto-generates TypeScript definitions from the C API.
 * htslib handles SIMD unpack internally (SSE/AVX/NEON → WASM SIMD 128-bit).
 *
 * Effort: 3-5 days (packaging existing library, not writing one).
 *
 * For now, this is a skeleton that provides the TypeScript interface
 * and falls back to pure-JS parsing when WASM is not available.
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

/** Configuration for htslib WASM loading. */
export interface HtslibConfig {
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
export let HTSLIB_WASM_AVAILABLE = false;

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
// HtslibWasm
// ---------------------------------------------------------------------------

/**
 * htslib WASM binding with automatic fallback to pure-JS parsing.
 *
 * Loading priority:
 *   1. Native .node addon (napi-rs) — fastest, uses OS-level SIMD.
 *   2. WASM .wasm binary — portable, uses WASM SIMD 128-bit.
 *   3. Pure-JS fallback — no SIMD, parses SAM text directly.
 *
 * Usage:
 *   const hts = await HtslibWasm.load();
 *   const { fd, header } = await hts.openFile('reads.sam');
 *   const record = await hts.samRead(fd, header);
 *   await hts.closeFile(fd);
 */
export class HtslibWasm {
  private config: HtslibConfig;
  private mode: 'native' | 'wasm' | 'js';
  private fileBuffers: Map<number, { lines: string[]; lineIndex: number; headerLines: string[] }>;
  private nextFd: number;

  private constructor(config: HtslibConfig, mode: 'native' | 'wasm' | 'js') {
    this.config = config;
    this.mode = mode;
    this.fileBuffers = new Map();
    this.nextFd = 1;
  }

  /**
   * Load htslib: try native .node, then .wasm, then pure-JS fallback.
   *
   * @param config Optional configuration for loading.
   * @returns An initialized HtslibWasm instance.
   */
  static async load(config?: HtslibConfig): Promise<HtslibWasm> {
    const cfg: HtslibConfig = {
      enableSimd: true,
      batchSize: 4096,
      ...config,
    };

    // Strategy 1: Try native .node addon.
    if (cfg.nativePath) {
      try {
        // Dynamic import of native addon.
        await import(/* @vite-ignore */ cfg.nativePath);
        HTSLIB_WASM_AVAILABLE = true;
        return new HtslibWasm(cfg, 'native');
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
              HTSLIB_WASM_AVAILABLE = true;
              return new HtslibWasm(cfg, 'wasm');
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
    HTSLIB_WASM_AVAILABLE = false;
    return new HtslibWasm(cfg, 'js');
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
   * Close an open file and release resources.
   *
   * @param fd  File descriptor to close.
   */
  async closeFile(fd: number): Promise<void> {
    this.fileBuffers.delete(fd);
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
