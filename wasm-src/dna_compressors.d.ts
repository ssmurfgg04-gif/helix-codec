// DNA Compressors WASM — TypeScript declarations
export interface DnaCompressorsModule {
  _dna_compress(algo: number, input: number, input_len: number, output: number, output_cap: number): number;
  _dna_decompress(algo: number, input: number, input_len: number, output: number, output_cap: number): number;
  _dna_compressor_name(algo: number): number;
  _dna_compressor_count(): number;
  _malloc(size: number): number;
  _free(ptr: number): void;
  HEAPU8: Uint8Array;
  UTF8ToString(ptr: number): string;
}

export default function createModule(): Promise<DnaCompressorsModule>;
