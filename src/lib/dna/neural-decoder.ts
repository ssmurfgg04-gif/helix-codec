/**
 * Neural Decoder Interface — ONNX-runtime-web adapter
 *
 * Provides an interface for ML-based decoders (Neural Polar Decoders, Astra
 * learned-BP, AGNES GNN seed-chaining). Loads pre-trained models via
 * ONNX-runtime-web (WASM) and runs inference.
 *
 * This is the adapter layer — the actual models need to be trained separately
 * (PyTorch/TensorFlow) and exported to ONNX format. The interface supports:
 *   - Loading ONNX models
 *   - Running inference on input tensors
 *   - Fallback to algorithmic decoder if model unavailable
 *
 * Reference:
 *   - Aharoni & Pfister (2025). "Neural Polar Decoders for DNA Storage."
 *     arXiv:2506.17076. github.com/zivaharoni/neural-polar-decoders-dna-storage
 *   - Maan et al. (2025). "Astra: Learning BP on Tanner Graphs."
 *     arXiv:2408.07038.
 *   - ONNX Runtime: github.com/microsoft/onnxruntime
 */

export interface NeuralModelConfig {
  /** Path or URL to the ONNX model file. */
  modelPath: string;
  /** Input tensor name. */
  inputName: string;
  /** Output tensor name. */
  outputName: string;
  /** Input shape (e.g., [1, 252] for 252-bit codewords). */
  inputShape: number[];
  /** Output shape. */
  outputShape: number[];
}

export interface NeuralDecodeResult {
  /** Decoded output (bits or bases). */
  output: Float32Array;
  /** Inference time in ms. */
  inferenceMs: number;
  /** Whether a model was used (false = fallback). */
  usedModel: boolean;
  /** Confidence scores per output position (0-1). */
  confidence?: Float32Array;
}

/**
 * Neural Decoder — loads ONNX models and runs inference.
 *
 * Falls back to an algorithmic decoder if:
 *   - ONNX runtime is not available
 *   - Model file is not found
 *   - Inference fails
 */
export class NeuralDecoder {
  private config: NeuralModelConfig;
  private session: any = null;
  private available = false;

  constructor(config: NeuralModelConfig) {
    this.config = config;
    this.init();
  }

  private async init(): Promise<void> {
    try {
      // Try to load ONNX runtime (dynamic import to avoid hard dependency)
      const ort = await import("onnxruntime-web").catch(() => null);
      if (!ort) {
        console.warn("ONNX runtime not available — neural decoder will use fallback");
        return;
      }

      this.session = await ort.InferenceSession.create(this.config.modelPath);
      this.available = true;
    } catch (e) {
      console.warn("Failed to load neural model:", (e as Error).message);
      this.available = false;
    }
  }

  /**
   * Run inference on input data.
   * Falls back to algorithmic decoder if model unavailable.
   */
  async decode(
    input: Float32Array,
    fallback: (input: Float32Array) => Float32Array,
  ): Promise<NeuralDecodeResult> {
    if (!this.available || !this.session) {
      const t0 = Date.now();
      const output = fallback(input);
      return {
        output,
        inferenceMs: Date.now() - t0,
        usedModel: false,
      };
    }

    try {
      const t0 = Date.now();
      const tensor = new (await import("onnxruntime-web")).Tensor(
        "float32",
        input,
        this.config.inputShape,
      );
      const results = await this.session.run({ [this.config.inputName]: tensor });
      const output = results[this.config.outputName].data as Float32Array;

      // Compute confidence from output softmax (if applicable)
      const confidence = this.computeConfidence(output);

      return {
        output,
        inferenceMs: Date.now() - t0,
        usedModel: true,
        confidence,
      };
    } catch (e) {
      console.warn("Neural inference failed, using fallback:", (e as Error).message);
      const output = fallback(input);
      return { output, inferenceMs: 0, usedModel: false };
    }
  }

  private computeConfidence(output: Float32Array): Float32Array {
    // Simple softmax confidence: max(p) per position
    const numPositions = output.length / 4; // assuming 4 classes (A,C,G,T)
    const confidence = new Float32Array(numPositions);
    for (let i = 0; i < numPositions; i++) {
      let maxVal = -Infinity;
      for (let j = 0; j < 4; j++) {
        if (output[i * 4 + j] > maxVal) maxVal = output[i * 4 + j];
      }
      // Softmax: exp(x) / sum(exp(x))
      let sum = 0;
      for (let j = 0; j < 4; j++) sum += Math.exp(output[i * 4 + j] - maxVal);
      confidence[i] = Math.exp(maxVal - maxVal) / sum; // = 1/sum... wait, this is wrong
      // Actually: confidence = max(softmax) = exp(max) / sum(exp)
      confidence[i] = 1 / sum; // since exp(max - max) = 1
    }
    return confidence;
  }

  /** Check if the neural model is loaded and available. */
  isAvailable(): boolean {
    return this.available;
  }
}

/**
 * Stub for Neural Polar Decoder (Aharoni 2025).
 * In production, load the trained ONNX model from:
 *   github.com/zivaharoni/neural-polar-decoders-dna-storage
 */
export function createNeuralPolarDecoder(): NeuralDecoder {
  return new NeuralDecoder({
    modelPath: "/models/neural_polar_dna.onnx",
    inputName: "llr",
    outputName: "info_bits",
    inputShape: [1, 252],
    outputShape: [1, 208],
  });
}

/**
 * Stub for Astra learned-BP decoder (Maan 2025).
 */
export function createAstraDecoder(): NeuralDecoder {
  return new NeuralDecoder({
    modelPath: "/models/astra_bp.onnx",
    inputName: "llr",
    outputName: "codeword",
    inputShape: [1, 252],
    outputShape: [1, 252],
  });
}

/**
 * Stub for AGNES GNN seed-chaining (RawHash3, arXiv:2510.16013).
 * Would replace the k-mer prefilter in the HMM alignment pipeline.
 */
export function createAGNESChainer(): NeuralDecoder {
  return new NeuralDecoder({
    modelPath: "/models/agnes_seed_chaining.onnx",
    inputName: "reads",
    outputName: "chains",
    inputShape: [1, 200],
    outputShape: [1, 15],
  });
}
