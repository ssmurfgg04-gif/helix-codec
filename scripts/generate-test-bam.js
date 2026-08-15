#!/usr/bin/env node
/**
 * Generate a valid uncompressed BAM file for testing.
 *
 * BAM binary format (uncompressed):
 *   magic(4) | l_text(4) | text(l_text) | n_ref(4) | [ref entries...]
 *   [alignment records...]
 *
 * Each ref entry: l_name(4) | name(l_name, NUL-terminated) | l_ref(4)
 *
 * Each alignment record:
 *   block_size(4) | refID(4) | pos(4) | l_read_name(1) | mapq(1) | bin(2) |
 *   n_cigar_op(2) | flag(2) | l_seq(4) | next_refID(4) | next_pos(4) |
 *   tlen(4) | read_name(l_read_name) | cigar(n_cigar_op*4) |
 *   seq(ceil(l_seq/2)) | qual(l_seq) | [aux...]
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// ---- helpers ----

function writeInt32LE(buf, offset, val) {
  buf[offset]     = val & 0xFF;
  buf[offset + 1] = (val >> 8) & 0xFF;
  buf[offset + 2] = (val >> 16) & 0xFF;
  buf[offset + 3] = (val >> 24) & 0xFF;
  return offset + 4;
}

function writeUint16LE(buf, offset, val) {
  buf[offset]     = val & 0xFF;
  buf[offset + 1] = (val >> 8) & 0xFF;
  return offset + 2;
}

// ---- header ----

const headerText = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n';
const headerBytes = Buffer.from(headerText, 'utf-8');
const l_text = headerBytes.length;

// Reference sequences
const refs = [{ name: 'chr1', length: 1000 }];

// ---- alignment record ----

const qname = 'read1';
const flag = 0;         // mapped
const refID = 0;        // chr1
const pos = 0;          // 0-based
const mapq = 30;
const nCigarOp = 1;     // 10M
const cigarOp = (10 << 4) | 0; // 10M: len=10, op=M(0)
const seq = 'ACGTACGTAC';
const lSeq = seq.length;
const qualVal = 30;     // Phred 30 for all bases
const nextRefID = -1;
const nextPos = -1;
const tlen = 0;

// Compute BAI bin for [beg=0, end=10)
// hts_reg2bin: level 5 offset = ((1<<15)-1)/7 = 4681; bin = 4681 + (beg>>14)
const bin = 4681 + (0 >> 14);

// Encode read name (NUL-terminated)
const qnameBytes = Buffer.from(qname + '\0', 'utf-8');
const lReadName = qnameBytes.length;

// Encode sequence in 4-bit BAM encoding
// A=1, C=2, G=4, T=8, N=15
const baseMap = { A: 1, C: 2, G: 4, T: 8, N: 15 };
const seqEncodedLen = Math.ceil(lSeq / 2);
const seqEncoded = Buffer.alloc(seqEncodedLen);
for (let i = 0; i < lSeq; i += 2) {
  const hi = baseMap[seq[i]] || 15;
  const lo = (i + 1 < lSeq) ? (baseMap[seq[i + 1]] || 15) : 0;
  seqEncoded[i / 2] = (hi << 4) | lo;
}

// Encode quality (raw Phred values, not +33)
const qualEncoded = Buffer.alloc(lSeq, qualVal);

// Variable-length data
const varData = Buffer.concat([
  qnameBytes,                               // read name
  (() => { const b = Buffer.alloc(4); b.writeUInt32LE(cigarOp, 0); return b; })(),  // CIGAR
  seqEncoded,                               // sequence
  qualEncoded,                              // quality
]);

// block_size = 32 (fixed fields) + varData.length
const blockSize = 32 + varData.length;

// ---- assemble the BAM file ----

// Calculate total size
const refNameBytes = Buffer.from(refs[0].name + '\0', 'utf-8');
const lName = refNameBytes.length;

const totalSize =
  4 +                     // magic
  4 + l_text +            // header text
  4 +                     // n_ref
  (4 + lName + 4) +       // ref entries
  4 + blockSize;          // alignment record (block_size prefix + data)

const buf = Buffer.alloc(totalSize);
let off = 0;

// Magic: "BAM\1"
buf[off++] = 0x42; // B
buf[off++] = 0x41; // A
buf[off++] = 0x4D; // M
buf[off++] = 0x01; // \1

// Header text
off = writeInt32LE(buf, off, l_text);
headerBytes.copy(buf, off);
off += l_text;

// Number of reference sequences
off = writeInt32LE(buf, off, refs.length);

// Reference entries
for (const ref of refs) {
  const rn = Buffer.from(ref.name + '\0', 'utf-8');
  off = writeInt32LE(buf, off, rn.length);
  rn.copy(buf, off);
  off += rn.length;
  off = writeInt32LE(buf, off, ref.length);
}

// ---- alignment record ----

// block_size
off = writeInt32LE(buf, off, blockSize);

// Fixed fields (32 bytes)
off = writeInt32LE(buf, off, refID);      // refID
off = writeInt32LE(buf, off, pos);         // pos
buf[off++] = lReadName;                    // l_read_name
buf[off++] = mapq;                         // mapq
off = writeUint16LE(buf, off, bin);        // bin
off = writeUint16LE(buf, off, nCigarOp);   // n_cigar_op
off = writeUint16LE(buf, off, flag);       // flag
off = writeInt32LE(buf, off, lSeq);        // l_seq
off = writeInt32LE(buf, off, nextRefID);   // next_refID
off = writeInt32LE(buf, off, nextPos);     // next_pos
off = writeInt32LE(buf, off, tlen);        // tlen

// Variable-length data
varData.copy(buf, off);
off += varData.length;

// Write the file
const outPath = path.join(__dirname, 'test-data', 'minimal.bam');
fs.writeFileSync(outPath, buf.slice(0, off));
console.log(`Created ${outPath} (${off} bytes)`);
console.log(`  Header text: ${JSON.stringify(headerText)}`);
console.log(`  1 alignment: ${qname} mapped to chr1:0, seq=${seq}, qual=${qualVal}`);
