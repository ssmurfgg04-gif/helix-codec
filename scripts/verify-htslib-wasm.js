#!/usr/bin/env node
/**
 * Verify that the htslib WASM module works correctly.
 *
 * This script:
 *   1. Loads the htslib WASM module via createRequire (CJS glue)
 *   2. Creates a minimal BAM in memory
 *   3. Opens it with _hts_open_mem
 *   4. Reads the header
 *   5. Verifies the reference count and name
 *   6. Reads alignment records
 *   7. Prints "htslib WASM: REAL ✓" on success
 */

import fs from 'node:fs';
import path from 'node:path';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// ---- Step 1: Load the htslib WASM module ----

const wasmDir = path.join(__dirname, '..', 'src', 'lib', 'dna', 'pkg', 'htslib-wasm');
const require = createRequire(path.join(wasmDir, 'htslib_wasm.js'));

let Module;
try {
  const createHtslibWasmModule = require('./htslib_wasm.js');
  Module = await createHtslibWasmModule();
  console.log('[1] WASM module loaded successfully');
} catch (err) {
  console.error('[1] FAILED to load WASM module:', err.message);
  process.exit(1);
}

// Verify key exports exist
const requiredExports = [
  '_hts_open_mem', '_sam_hdr_read', '_sam_read1',
  '_bam_init1', '_bam_destroy1', '_hts_close',
  '_sam_hdr_destroy', '_hdr_n_ref', '_hdr_ref_name',
  '_hdr_ref_len', '_hdr_text', '_bam_core_tid', '_bam_core_pos',
  '_bam_core_qual', '_bam_core_flag', '_bam_core_l_qseq',
  '_bam_qname', '_bam_cigar', '_bam_n_cigar',
  '_bam_seq_str', '_bam_qual_str', '_malloc', '_free',
];

for (const name of requiredExports) {
  if (typeof Module[name] !== 'function') {
    console.error(`[1] Missing export: ${name}`);
    process.exit(1);
  }
}
console.log(`[1] All ${requiredExports.length} expected exports present`);
console.log('[1] HEAPU8 available:', !!Module.HEAPU8);

// ---- Step 2: Create minimal BAM in memory ----

// We'll reuse the generated file, or build one inline
const bamPath = path.join(__dirname, 'test-data', 'minimal.bam');
let bamData;

if (fs.existsSync(bamPath)) {
  bamData = new Uint8Array(fs.readFileSync(bamPath));
  console.log(`[2] Loaded test BAM from ${bamPath} (${bamData.length} bytes)`);
} else {
  // Build a minimal BAM inline
  console.log('[2] test-data/minimal.bam not found, building inline...');

  const headerText = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n';
  const headerBytes = Buffer.from(headerText, 'utf-8');

  // Build reference entry
  const refName = Buffer.from('chr1\0', 'utf-8');

  // Build alignment record
  const qnameBytes = Buffer.from('read1\0', 'utf-8');
  const cigarBuf = Buffer.alloc(4);
  cigarBuf.writeUInt32LE((10 << 4) | 0); // 10M

  // 4-bit seq encoding: ACGTACGTAC
  const seqEnc = Buffer.from([0x12, 0x48, 0x12, 0x48, 0x12]); // A=1,C=2,G=4,T=8
  const qualEnc = Buffer.alloc(10, 30); // Phred 30

  const varData = Buffer.concat([qnameBytes, cigarBuf, seqEnc, qualEnc]);
  const blockSize = 32 + varData.length;

  const record = Buffer.alloc(4 + blockSize);
  let off = 0;
  record.writeInt32LE(blockSize, off); off += 4;
  record.writeInt32LE(0, off); off += 4;      // refID = 0
  record.writeInt32LE(0, off); off += 4;      // pos = 0
  record[off++] = qnameBytes.length;           // l_read_name
  record[off++] = 30;                          // mapq
  record.writeUInt16LE(4681, off); off += 2;  // bin
  record.writeUInt16LE(1, off); off += 2;     // n_cigar_op
  record.writeUInt16LE(0, off); off += 2;     // flag
  record.writeInt32LE(10, off); off += 4;     // l_seq
  record.writeInt32LE(-1, off); off += 4;     // next_refID
  record.writeInt32LE(-1, off); off += 4;     // next_pos
  record.writeInt32LE(0, off); off += 4;      // tlen
  varData.copy(record, off);

  // Assemble full BAM
  const parts = [
    Buffer.from([0x42, 0x41, 0x4D, 0x01]),  // BAM\1
    Buffer.alloc(4),                           // l_text placeholder
    headerBytes,
    Buffer.alloc(4),                           // n_ref placeholder
  ];
  const ltBuf = Buffer.alloc(4);
  ltBuf.writeInt32LE(headerBytes.length);
  parts[1] = ltBuf;

  const nrBuf = Buffer.alloc(4);
  nrBuf.writeInt32LE(1);
  parts[3] = nrBuf;

  const lnBuf = Buffer.alloc(4);
  lnBuf.writeInt32LE(refName.length);
  const rlBuf = Buffer.alloc(4);
  rlBuf.writeInt32LE(1000);

  parts.push(lnBuf, refName, rlBuf, record);

  bamData = new Uint8Array(Buffer.concat(parts));
  console.log(`[2] Built inline BAM (${bamData.length} bytes)`);
}

// ---- Step 3: Open BAM with _hts_open_mem ----

const dataPtr = Module._malloc(bamData.length);
Module.HEAPU8.set(bamData, dataPtr);
console.log(`[3] Copied ${bamData.length} bytes to WASM memory at ptr=${dataPtr}`);

const fp = Module._hts_open_mem(dataPtr, bamData.length);
Module._free(dataPtr);

if (!fp) {
  console.error('[3] _hts_open_mem returned NULL — invalid BAM format');
  process.exit(1);
}
console.log(`[3] _hts_open_mem succeeded: fp=${fp}`);

// ---- Step 4: Read the header ----

const hdr = Module._sam_hdr_read(fp);
if (!hdr) {
  console.error('[4] _sam_hdr_read returned NULL — could not read header');
  Module._hts_close(fp);
  process.exit(1);
}
console.log(`[4] _sam_hdr_read succeeded: hdr=${hdr}`);

// ---- Step 5: Verify reference count and name ----

const nRef = Module._hdr_n_ref(hdr);
console.log(`[5] Reference count: ${nRef}`);

if (nRef !== 1) {
  console.error(`[5] Expected 1 reference, got ${nRef}`);
  Module._sam_hdr_destroy(hdr);
  Module._hts_close(fp);
  process.exit(1);
}

const namePtr = Module._hdr_ref_name(hdr, 0);
const refNameStr = namePtr ? Module.UTF8ToString(namePtr) : '';
const refLen = Module._hdr_ref_len(hdr, 0);
console.log(`[5] Reference 0: name="${refNameStr}", length=${refLen}`);

if (refNameStr !== 'chr1') {
  console.error(`[5] Expected ref name "chr1", got "${refNameStr}"`);
  Module._sam_hdr_destroy(hdr);
  Module._hts_close(fp);
  process.exit(1);
}

if (refLen !== 1000) {
  console.error(`[5] Expected ref length 1000, got ${refLen}`);
  Module._sam_hdr_destroy(hdr);
  Module._hts_close(fp);
  process.exit(1);
}

const textPtr = Module._hdr_text(hdr);
const headerText = textPtr ? Module.UTF8ToString(textPtr) : '';
console.log(`[5] Header text: ${JSON.stringify(headerText)}`);

// ---- Step 6: Read alignment records ----

const b = Module._bam_init1();
if (!b) {
  console.error('[6] _bam_init1 returned NULL');
  Module._sam_hdr_destroy(hdr);
  Module._hts_close(fp);
  process.exit(1);
}
console.log(`[6] _bam_init1 succeeded: b=${b}`);

let recordCount = 0;
while (true) {
  const ret = Module._sam_read1(fp, hdr, b);
  if (ret < 0) break;

  recordCount++;
  const tid = Module._bam_core_tid(b);
  const pos = Module._bam_core_pos(b);
  const qual = Module._bam_core_qual(b);
  const flag = Module._bam_core_flag(b);
  const lQseq = Module._bam_core_l_qseq(b);

  const qnamePtr = Module._bam_qname(b);
  const qname = qnamePtr ? Module.UTF8ToString(qnamePtr) : '*';

  const nCigar = Module._bam_n_cigar(b);
  const cigarPtr = Module._bam_cigar(b);
  let cigar = '*';
  if (nCigar > 0 && cigarPtr) {
    // Read CIGAR as individual uint32 values to avoid alignment issues
    const ops = [];
    for (let i = 0; i < nCigar; i++) {
      const off = cigarPtr + i * 4;
      const c = Module.HEAPU8[off] |
                (Module.HEAPU8[off + 1] << 8) |
                (Module.HEAPU8[off + 2] << 16) |
                (Module.HEAPU8[off + 3] << 24);
      const op = c & 0xF;
      const len = c >>> 4;
      ops.push(`${len}${'MIDNSHP=X'[op] || '?'}`);
    }
    cigar = ops.join('');
  }

  const seqPtr = Module._bam_seq_str(b);
  const seqStr = seqPtr ? Module.UTF8ToString(seqPtr) : '*';
  if (seqPtr) Module._free(seqPtr);

  const qualPtr = Module._bam_qual_str(b);
  const qualStr = qualPtr ? Module.UTF8ToString(qualPtr) : '*';
  if (qualPtr) Module._free(qualPtr);

  console.log(`[6] Record ${recordCount}: ${qname} flag=${flag} tid=${tid} pos=${pos} mapq=${qual} cigar=${cigar} seq=${seqStr} qual=${qualStr}`);
}

console.log(`[6] Total records read: ${recordCount}`);

if (recordCount !== 1) {
  console.error(`[6] Expected 1 record, got ${recordCount}`);
  Module._bam_destroy1(b);
  Module._sam_hdr_destroy(hdr);
  Module._hts_close(fp);
  process.exit(1);
}

// ---- Cleanup ----

Module._bam_destroy1(b);
Module._sam_hdr_destroy(hdr);
Module._hts_close(fp);
console.log('[7] Cleanup complete');

// ---- Success ----

console.log('');
console.log('htslib WASM: REAL ✓');
