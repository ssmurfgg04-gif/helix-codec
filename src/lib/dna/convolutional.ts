/**
 * Convolutional Inner Code with Viterbi Decoder
 *
 * A rate-1/2 convolutional code with memory 2, decoded via Viterbi algorithm.
 * Designed for indel tolerance in DNA storage (HEDGES approach).
 *
 * The convolutional encoder processes the input bit stream and produces 2 output
 * bits per input bit. The Viterbi decoder finds the maximum-likelihood input
 * sequence given the received (noisy) output bits.
 *
 * For DNA storage, the key advantage over block codes is that convolutional
 * codes can recover from insertions and deletions (with modifications), while
 * block codes (RS, LDPC) only handle substitutions.
 *
 * HEDGES (Press 2020) uses a convolutional code with greedy tree search to
 * handle indels. This implementation uses the standard Viterbi algorithm for
 * substitutions, and can be extended with synchronization markers for indels.
 *
 * Code parameters:
 *   - Rate: 1/2 (1 input bit → 2 output bits)
 *   - Memory: 2 (4 states)
 *   - Generator polynomials: G1 = 7 (111), G2 = 5 (101)
 *   - Constraint length: K = 3
 *   - Free distance: 5 (corrects 2 errors per constraint length)
 *
 * For DNA storage, we use this as an OPTIONAL inner code alongside LDPC:
 *   Payload → LDPC (substitution) → Convolutional (indel) → DNA mapping
 *
 * Reference:
 *   - Press et al. (2020). "HEDGES error-correcting code for DNA storage."
 *     PNAS 117:31.
 *   - Viterbi (1967). "Error bounds for convolutional codes and an
 *     asymptotically optimum decoding algorithm." IEEE Trans. IT 13:2.
 */

export interface ConvolutionalConfig {
  /** Rate (1/n outputs per input bit). Currently only 1/2 supported. */
  rate: number;
  /** Memory (constraint length - 1). Default 2. */
  memory: number;
  /** Generator polynomials (octal). Default [7, 5] for rate 1/2. */
  generators: number[];
}

export const DEFAULT_CONV_CONFIG: ConvolutionalConfig = {
  rate: 2,
  memory: 2,
  generators: [7, 5], // G1=111, G2=101
};

/**
 * Convolutional encoder (rate 1/2, memory 2).
 *
 * State machine:
 *   State = last `memory` input bits (2 bits = 4 states)
 *   Input: 1 bit
 *   Output: 2 bits (G1, G2)
 *   Next state: shift in the input bit
 *
 * Generator polynomials (for memory 2):
 *   G1 = 111 (octal 7): output = bit[0] XOR bit[1] XOR bit[2]
 *   G2 = 101 (octal 5): output = bit[0] XOR bit[2]
 *   where bit[0] = current input, bit[1] = previous, bit[2] = before that
 */
export class ConvolutionalCode {
  readonly rate: number;
  readonly memory: number;
  readonly numStates: number;
  readonly generators: number[];

  // Precomputed transition table: [state][input] → [output, nextState]
  private transitions: { output: number; nextState: number }[][];

  constructor(cfg: ConvolutionalConfig = DEFAULT_CONV_CONFIG) {
    this.rate = cfg.rate;
    this.memory = cfg.memory;
    this.numStates = 1 << cfg.memory;
    this.generators = cfg.generators;
    this.transitions = this.buildTransitionTable();
  }

  private buildTransitionTable(): { output: number; nextState: number }[][] {
    const table: { output: number; nextState: number }[][] = [];
    for (let state = 0; state < this.numStates; state++) {
      table[state] = [];
      for (let input = 0; input < 2; input++) {
        // Build the register: [input, state bits...]
        const reg = (input << this.memory) | state;
        let output = 0;
        for (let g = 0; g < this.generators.length; g++) {
          const gen = this.generators[g];
          let bit = 0;
          for (let b = 0; b < this.memory + 1; b++) {
            if ((gen >> b) & 1) bit ^= (reg >> b) & 1;
          }
          output = (output << 1) | bit;
        }
        const nextState = (reg >> 1) & (this.numStates - 1);
        table[state][input] = { output, nextState };
      }
    }
    return table;
  }

  /**
   * Encode a bit array using the convolutional code.
   *
   * The output includes a zero-tail of length `memory` to flush the encoder
   * back to state 0, which lets the Viterbi decoder traceback from state 0.
   *
   * @param inputBits Input bits (length = k)
   * @returns Output bits (length = (k + memory) * rate), including tail
   */
  encode(inputBits: Uint8Array): Uint8Array {
    // BUGFIX v52: original allocated inputBits.length * rate but the zero-tail
    // loop writes memory additional rate-bit blocks → out-of-bounds writes.
    // Allocate room for the tail.
    const totalSteps = inputBits.length + this.memory;
    const outputBits = new Uint8Array(totalSteps * this.rate);
    let state = 0;

    for (let i = 0; i < inputBits.length; i++) {
      const input = inputBits[i] & 1;
      const { output, nextState } = this.transitions[state][input];
      for (let r = 0; r < this.rate; r++) {
        outputBits[i * this.rate + r] = (output >> (this.rate - 1 - r)) & 1;
      }
      state = nextState;
    }

    // Flush with zero tail (memory bits) — drives encoder back to state 0
    for (let i = 0; i < this.memory; i++) {
      const { output, nextState } = this.transitions[state][0];
      for (let r = 0; r < this.rate; r++) {
        outputBits[(inputBits.length + i) * this.rate + r] = (output >> (this.rate - 1 - r)) & 1;
      }
      state = nextState;
    }

    return outputBits;
  }

  /**
   * Decode using Viterbi algorithm (maximum-likelihood sequence decoder).
   *
   * BUGFIX v52: The original traceback tried to recover the predecessor state
   * by scanning transitions[prev][input].nextState === state, but multiple
   * `prev` values can lead to the same (state, input) pair (e.g., for memory=2,
   * rate=1/2: transitions[0][0].nextState === transitions[1][0].nextState === 0).
   * This caused the traceback to always pick the first match, producing
   * all-zero output. Fix: store the predecessor state explicitly in the path.
   *
   * @param receivedBits Received bits (length = k * rate, with tail)
   * @param numInfoBits Number of original information bits (excluding tail)
   * @returns Decoded information bits
   */
  decode(receivedBits: Uint8Array, numInfoBits: number): Uint8Array {
    const totalBits = numInfoBits + this.memory; // including tail
    const decoded = new Uint8Array(numInfoBits);

    // Viterbi DP tables
    // pathMetric[state] = accumulated metric for the best path ending at `state`
    let pathMetric = new Float64Array(this.numStates).fill(Infinity);
    pathMetric[0] = 0; // start at state 0

    // paths[step] is a pair of Int8Arrays:
    //   pathsInput[step][state] = the input bit that led to this state at this step
    //   pathsPrev[step][state]  = the previous state (before the transition)
    const pathsInput: Int8Array[] = [];
    const pathsPrev: Int8Array[] = [];

    for (let step = 0; step < totalBits; step++) {
      const newMetric = new Float64Array(this.numStates).fill(Infinity);
      const stepInput = new Int8Array(this.numStates).fill(-1);
      const stepPrev = new Int8Array(this.numStates).fill(-1);

      for (let state = 0; state < this.numStates; state++) {
        if (pathMetric[state] === Infinity) continue;

        for (let input = 0; input < 2; input++) {
          const { output, nextState } = this.transitions[state][input];
          // Compute Hamming distance between expected output and received bits
          let dist = 0;
          for (let r = 0; r < this.rate; r++) {
            const expectedBit = (output >> (this.rate - 1 - r)) & 1;
            const receivedBit = receivedBits[step * this.rate + r] ?? 0;
            if (expectedBit !== receivedBit) dist++;
          }
          const metric = pathMetric[state] + dist;
          if (metric < newMetric[nextState]) {
            newMetric[nextState] = metric;
            stepInput[nextState] = input;
            stepPrev[nextState] = state; // record predecessor explicitly
          }
        }
      }

      pathMetric = newMetric;
      pathsInput.push(stepInput);
      pathsPrev.push(stepPrev);
    }

    // Traceback from state 0 (due to zero tail)
    let state = 0;
    for (let step = totalBits - 1; step >= 0; step--) {
      const input = pathsInput[step][state];
      const prev = pathsPrev[step][state];
      if (step < numInfoBits && input >= 0) {
        decoded[step] = input;
      }
      if (prev >= 0) {
        state = prev;
      }
    }

    return decoded;
  }
}

/**
 * Convert bytes to bits and back.
 */
export function bytesToBits(data: Uint8Array): Uint8Array {
  const bits = new Uint8Array(data.length * 8);
  for (let i = 0; i < data.length; i++) {
    for (let b = 0; b < 8; b++) {
      bits[i * 8 + b] = (data[i] >> (7 - b)) & 1;
    }
  }
  return bits;
}

export function bitsToBytes(bits: Uint8Array): Uint8Array {
  const out = new Uint8Array(Math.floor(bits.length / 8));
  for (let i = 0; i < out.length; i++) {
    let byte = 0;
    for (let b = 0; b < 8; b++) {
      byte |= bits[i * 8 + b] << (7 - b);
    }
    out[i] = byte;
  }
  return out;
}

/**
 * v52: Byte-oriented wrapper around the convolutional code.
 *
 * Used by the HEDGES-style pipeline:
 *
 *   ENCODE:  data → LDPC → conv-encode(bits) → DNA mapping
 *   DECODE:  DNA → bits → conv-Viterbi-decode → LDPC → data
 *
 * The convolutional code is the INNER code (closer to the channel) and
 * provides indel tolerance via Viterbi maximum-likelihood sequence estimation.
 * LDPC is the OUTER code and corrects residual substitutions.
 *
 * Rate: 1/2 (1 input bit → 2 output bits). Each input byte (8 bits) becomes
 * 16 bits = 2 bytes of conv-encoded output, plus a small tail (memory*rate
 * bits = 4 bits = 0.5 bytes, padded up to 1 byte).
 *
 * For a payload of N bytes, the conv-encoded output is:
 *   ceil((N*8 + memory) * rate / 8) bytes ≈ 2N + 1 bytes
 *
 * The decoder recovers the original N bytes from the 2N+1 byte stream.
 */
export class ConvolutionalInnerCode {
  private conv: ConvolutionalCode;
  readonly inputBytes: number;
  readonly outputBytes: number;

  constructor(inputBytes: number, cfg: ConvolutionalConfig = DEFAULT_CONV_CONFIG) {
    this.conv = new ConvolutionalCode(cfg);
    this.inputBytes = inputBytes;
    // Output bits = (inputBits + memory) * rate, padded to byte boundary
    const inputBits = inputBytes * 8;
    const outputBits = (inputBits + this.conv.memory) * this.conv.rate;
    this.outputBytes = Math.ceil(outputBits / 8);
  }

  /**
   * Encode N input bytes to ~2N+1 output bytes.
   * The output is padded with zero bits to the next byte boundary.
   */
  encode(data: Uint8Array): Uint8Array {
    if (data.length !== this.inputBytes) {
      throw new Error(`ConvolutionalInnerCode.encode: expected ${this.inputBytes} bytes, got ${data.length}`);
    }
    const inputBits = bytesToBits(data);
    const outputBits = this.conv.encode(inputBits);
    // Pad to byte boundary (zero-padding at the end)
    const out = new Uint8Array(this.outputBytes);
    for (let i = 0; i < outputBits.length; i++) {
      out[i >> 3] |= outputBits[i] << (7 - (i & 7));
    }
    return out;
  }

  /**
   * Decode ~2N+1 received bytes back to N bytes via Viterbi MLSE.
   *
   * The decoder treats the channel as a BSC (binary symmetric channel) —
   * substitutions only. For indel tolerance, the HMM-based Viterbi
   * preprocessor (viterbi-preprocess.ts) must run BEFORE this decoder
   * to convert indels into substitutions.
   */
  decode(received: Uint8Array): { decoded: Uint8Array; corrected: number } {
    if (received.length < this.outputBytes) {
      // Pad with zeros if short (defensive — shouldn't happen in practice)
      const padded = new Uint8Array(this.outputBytes);
      padded.set(received);
      received = padded;
    }
    // Unpack bits (only the first outputBits bits are meaningful)
    const totalBits = (this.inputBytes * 8 + this.conv.memory) * this.conv.rate;
    const receivedBits = new Uint8Array(totalBits);
    for (let i = 0; i < totalBits; i++) {
      receivedBits[i] = (received[i >> 3] >> (7 - (i & 7))) & 1;
    }
    const decodedBits = this.conv.decode(receivedBits, this.inputBytes * 8);
    return { decoded: bitsToBytes(decodedBits), corrected: 0 };
  }
}
