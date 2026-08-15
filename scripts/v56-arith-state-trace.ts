#!/usr/bin/env bun
/**
 * v56 — Direct state tracing for 2 bytes of 0xFF.
 * Patches the ArithmeticDecoder and ArithmeticEncoder to log state.
 */

// Inline patched classes (copy of markov-arithmetic.ts with logging)

const BASES = ["A", "C", "G", "T"] as const;
type Base = typeof BASES[number];

const NUM_STATE_BITS = 24;
const STATE_MASK = 0xFFFFFF;
const HALF_STATE = 0x800000;
const QUARTER_STATE = 0x400000;

class BitInputStream {
  private bits: number[] = [];
  private idx = 0;
  constructor(data: Uint8Array) {
    for (let i = 0; i < data.length; i++)
      for (let b = 7; b >= 0; b--)
        this.bits.push((data[i] >> b) & 1);
  }
  readBit(): number { return this.idx < this.bits.length ? this.bits[this.idx++] : 0; }
  get position() { return this.idx; }
}

class BitOutputStream {
  private bits: number[] = [];
  writeBit(bit: number) { this.bits.push(bit & 1); }
  toBytes(n: number): Uint8Array {
    const out = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      let byte = 0;
      for (let b = 0; b < 8; b++) byte = (byte << 1) | (this.bits[i * 8 + b] ?? 0);
      out[i] = byte;
    }
    return out;
  }
  get length() { return this.bits.length; }
  get bitsArr() { return this.bits; }
}

// Patched ArithmeticDecoder (encoder role) — logs state
class ArithmeticDecoder {
  private low = 0;
  private high = STATE_MASK;
  private code = 0;
  private input: BitInputStream;
  constructor(input: BitInputStream) {
    this.input = input;
    for (let i = 0; i < NUM_STATE_BITS; i++)
      this.code = (((this.code << 1) | this.input.readBit()) & STATE_MASK) >>> 0;
  }
  read(numSymbols: number, pos: number): number {
    const range = (this.high - this.low + 1) >>> 0;
    const symbolWidth = Math.floor(range / numSymbols);
    const valueOffset = ((this.code - this.low) >>> 0);
    let symbolIdx = Math.floor(valueOffset / symbolWidth);
    if (symbolIdx >= numSymbols) symbolIdx = numSymbols - 1;
    this.high = (this.low + symbolWidth * (symbolIdx + 1) - 1) >>> 0;
    this.low = (this.low + symbolWidth * symbolIdx) >>> 0;
    while (true) {
      if (this.high < HALF_STATE) {}
      else if (this.low >= HALF_STATE) {
        this.code = (this.code - HALF_STATE) >>> 0;
        this.low = (this.low - HALF_STATE) >>> 0;
        this.high = (this.high - HALF_STATE) >>> 0;
      } else if (this.low >= QUARTER_STATE && this.high < (QUARTER_STATE * 3)) {
        this.code = (this.code - QUARTER_STATE) >>> 0;
        this.low = (this.low - QUARTER_STATE) >>> 0;
        this.high = (this.high - QUARTER_STATE) >>> 0;
      } else break;
      this.low = ((this.low << 1) & STATE_MASK) >>> 0;
      this.high = (((this.high << 1) | 1) & STATE_MASK) >>> 0;
      this.code = (((this.code << 1) | this.input.readBit()) & STATE_MASK) >>> 0;
    }
    console.log(`  ENC pos=${pos} sym=${symbolIdx} n=${numSymbols} low=0x${this.low.toString(16).padStart(6,'0')} high=0x${this.high.toString(16).padStart(6,'0')} code=0x${this.code.toString(16).padStart(6,'0')} bitpos=${this.input.position}`);
    return symbolIdx;
  }
}

// Patched ArithmeticEncoder (decoder role) — logs state
class ArithmeticEncoder {
  private low = 0;
  private high = STATE_MASK;
  private pendingBits = 0;
  private output: BitOutputStream;
  constructor(output: BitOutputStream) { this.output = output; }
  write(symbolIdx: number, numSymbols: number, pos: number) {
    const range = (this.high - this.low + 1) >>> 0;
    const symbolWidth = Math.floor(range / numSymbols);
    this.high = (this.low + symbolWidth * (symbolIdx + 1) - 1) >>> 0;
    this.low = (this.low + symbolWidth * symbolIdx) >>> 0;
    let bitsOutput = 0;
    while (true) {
      if (this.high < HALF_STATE) {
        this.outputBit(0); bitsOutput++;
      } else if (this.low >= HALF_STATE) {
        this.outputBit(1); bitsOutput++;
        this.low = (this.low - HALF_STATE) >>> 0;
        this.high = (this.high - HALF_STATE) >>> 0;
      } else if (this.low >= QUARTER_STATE && this.high < (QUARTER_STATE * 3)) {
        this.pendingBits++;
      } else break;
      this.low = ((this.low << 1) & STATE_MASK) >>> 0;
      this.high = (((this.high << 1) | 1) & STATE_MASK) >>> 0;
    }
    console.log(`  DEC pos=${pos} sym=${symbolIdx} n=${numSymbols} low=0x${this.low.toString(16).padStart(6,'0')} high=0x${this.high.toString(16).padStart(6,'0')} pend=${this.pendingBits} outBits=${bitsOutput} totalBits=${this.output.length}`);
  }
  private outputBit(bit: number) {
    this.output.writeBit(bit);
    while (this.pendingBits > 0) { this.output.writeBit(1 - bit); this.pendingBits--; }
  }
  finish() {
    this.pendingBits++;
    if (this.low < QUARTER_STATE) this.outputBit(0);
    else this.outputBit(1);
  }
}

// Test: 2 bytes of 0xFF
const data = new Uint8Array([0xFF, 0xFF]);
const numBases = 19;

console.log("=== ENCODER (ArithmeticDecoder) ===");
const encInput = new BitInputStream(data);
const encoder = new ArithmeticDecoder(encInput);
const symbols: number[] = [];
const allowedCounts: number[] = [];
let prev: Base = "A";
let runLength = 0;

for (let pos = 0; pos < numBases; pos++) {
  let allowed: Base[];
  if (runLength >= 3) allowed = BASES.filter(b => b !== prev);
  else allowed = [...BASES];
  const sym = encoder.read(allowed.length, pos);
  symbols.push(sym);
  allowedCounts.push(allowed.length);
  const base = allowed[Math.min(sym, allowed.length - 1)];
  if (base === prev) runLength++;
  else { runLength = 1; prev = base; }
}

console.log("\n=== DECODER (ArithmeticEncoder) ===");
const decOutput = new BitOutputStream();
const decoder = new ArithmeticEncoder(decOutput);

prev = "A";
runLength = 0;
for (let pos = 0; pos < numBases; pos++) {
  const allowed = runLength >= 3 ? BASES.filter(b => b !== prev) : [...BASES];
  decoder.write(symbols[pos], allowed.length, pos);
  const base = allowed[Math.min(symbols[pos], allowed.length - 1)];
  if (base === prev) runLength++;
  else { runLength = 1; prev = base; }
}
decoder.finish();

const result = decOutput.toBytes(2);
console.log(`\n=== RESULT ===`);
console.log(`Expected: FF FF`);
console.log(`Got:      ${result[0].toString(16).padStart(2,'0')} ${result[1].toString(16).padStart(2,'0')}`);
console.log(`Total output bits: ${decOutput.length}`);
console.log(`Output bits: ${decOutput.bitsArr.join("")}`);

// Show the expected bits
const expBits: number[] = [];
for (let i = 0; i < 2; i++) for (let b = 7; b >= 0; b--) expBits.push((data[i] >> b) & 1);
console.log(`Expected bits: ${expBits.join("")}`);
