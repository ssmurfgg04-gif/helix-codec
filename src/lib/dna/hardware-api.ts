/**
 * 
 * ⚠️ EXPERIMENTAL / ROADMAP — This module is a scaffolding prototype.
 * It is NOT wired into the production encode/decode pipeline.
 * See CHANGELOG.md for v59 status and remaining work.
 *
 * Hardware API — Wet-Lab Bridge + GPU/FPGA Acceleration
 *
 * Two concerns live in this module:
 *
 *  A. Wet-lab vendor stubs (IDT, Twist, MinKNOW) — see below.
 *  B. GPU / FPGA acceleration framework — see § Hardware Acceleration below.
 *
 * ── A. Wet-Lab Bridge ──────────────────────────────────────────────
 *
 * API client stubs for DNA synthesis and sequencing hardware:
 *   1. IDT (Integrated DNA Technologies) — oligo ordering API
 *   2. Twist Bioscience — silicon-based DNA synthesis API
 *   3. Oxford Nanopore MinKNOW — sequencing software output format
 *
 * These are STUB implementations that generate the correct file formats
 * and calculate costs. Real API integration requires vendor accounts.
 *
 * Commands:
 *   helix synthesize dna://pool/key --vendor idt
 *   helix synthesize dna://pool/key --vendor twist
 *   helix export-minknow dna://pool/key
 *
 * ── B. Hardware Acceleration ───────────────────────────────────────
 *
 * Auto-detects available compute backends (CUDA, OpenCL, FPGA) and
 * dispatches expensive codec operations (RS encode/decode, LDPC decode,
 * batch simulation) to the fastest available hardware. Falls back to
 * WASM/JS when no accelerator is present.
 *
 * Usage:
 *   const mgr = HardwareManager.getInstance();
 *   await mgr.detectBackends();
 *   const backend = mgr.getBestBackend('rsEncode');
 *   const encoded = backend.rsEncode(data, nsym);
 */

export interface OligoOrder {
  vendor: "idt" | "twist" | "elegen";
  oligos: { name: string; sequence: string; scale: string; purification: string }[];
  totalLength: number;
  totalOligos: number;
  estimatedCostUSD: number;
  estimatedTurnaroundDays: number;
  orderFormat: "csv" | "fasta" | "xlsx";
  orderContent: string;
}

export interface VendorPricing {
  costPerBase: number;
  minimumOrder: number;
  setupFee: number;
  turnaroundDays: number;
  maxOligoLength: number;
  maxOligosPerOrder: number;
}

const VENDOR_PRICING: Record<string, VendorPricing> = {
  idt: {
    costPerBase: 0.10, // $0.10 per base (standard IDT pricing)
    minimumOrder: 50,
    setupFee: 25,
    turnaroundDays: 3,
    maxOligoLength: 230,
    maxOligosPerOrder: 50000,
  },
  twist: {
    costPerBase: 0.07, // Twist is cheaper per base
    minimumOrder: 100,
    setupFee: 50,
    turnaroundDays: 7,
    maxOligoLength: 300,
    maxOligosPerOrder: 100000,
  },
  elegen: {
    costPerBase: 0.05, // Elegen is cheapest (newer technology)
    minimumOrder: 200,
    setupFee: 75,
    turnaroundDays: 10,
    maxOligoLength: 300,
    maxOligosPerOrder: 200000,
  },
};

/**
 * Generate an oligo order for a DNA synthesis vendor.
 *
 * @param oligos Array of oligo sequences to synthesize
 * @param vendor Target vendor ("idt", "twist", "elegen")
 * @param scale Synthesis scale ("25nm", "100nm", "250nm", "1um")
 * @param purification Purification method ("STD", "PAGE", "HPLC")
 */
export function createSynthesisOrder(
  oligos: { name: string; sequence: string }[],
  vendor: "idt" | "twist" | "elegen" = "idt",
  scale: string = "25nm",
  purification: string = "STD",
): OligoOrder {
  const pricing = VENDOR_PRICING[vendor];

  // Validate oligo lengths
  for (const oligo of oligos) {
    if (oligo.sequence.length > pricing.maxOligoLength) {
      throw new Error(
        `Oligo ${oligo.name} is ${oligo.sequence.length}nt, exceeds ${vendor} max of ${pricing.maxOligoLength}nt`,
      );
    }
  }

  if (oligos.length > pricing.maxOligosPerOrder) {
    throw new Error(
      `Order has ${oligos.length} oligos, exceeds ${vendor} max of ${pricing.maxOligosPerOrder}`,
    );
  }

  // Calculate cost
  const totalBases = oligos.reduce((sum, o) => sum + o.sequence.length, 0);
  const baseCost = totalBases * pricing.costPerBase;
  const totalCost = Math.max(baseCost + pricing.setupFee, pricing.minimumOrder);

  // Generate order file
  const orderOligos = oligos.map(o => ({
    name: o.name,
    sequence: o.sequence,
    scale,
    purification,
  }));

  let orderContent: string;
  let orderFormat: "csv" | "fasta" | "xlsx";

  switch (vendor) {
    case "idt":
      orderFormat = "csv";
      orderContent = generateIDTCsv(orderOligos);
      break;
    case "twist":
      orderFormat = "fasta";
      orderContent = generateTwistFasta(orderOligos);
      break;
    case "elegen":
      orderFormat = "csv";
      orderContent = generateElegenCsv(orderOligos);
      break;
  }

  return {
    vendor,
    oligos: orderOligos,
    totalLength: totalBases,
    totalOligos: oligos.length,
    estimatedCostUSD: totalCost,
    estimatedTurnaroundDays: pricing.turnaroundDays,
    orderFormat,
    orderContent,
  };
}

/**
 * Generate IDT-format CSV order file.
 * IDT accepts CSV with columns: Name, Sequence, Scale, Purification
 */
function generateIDTCsv(oligos: { name: string; sequence: string; scale: string; purification: string }[]): string {
  const lines = ["Name,Sequence,Scale,Purification"];
  for (const o of oligos) {
    lines.push(`${o.name},${o.sequence},${o.scale},${o.purification}`);
  }
  return lines.join("\n");
}

/**
 * Generate Twist Bioscience FASTA order file.
 * Twist accepts FASTA format with oligo names.
 */
function generateTwistFasta(oligos: { name: string; sequence: string; scale: string; purification: string }[]): string {
  const lines: string[] = [];
  for (const o of oligos) {
    lines.push(`>${o.name}`);
    lines.push(o.sequence);
  }
  return lines.join("\n");
}

/**
 * Generate Elegen CSV order file.
 */
function generateElegenCsv(oligos: { name: string; sequence: string; scale: string; purification: string }[]): string {
  const lines = ["Oligo Name,Sequence,Length,Scale,Purification"];
  for (const o of oligos) {
    lines.push(`${o.name},${o.sequence},${o.sequence.length},${o.scale},${o.purification}`);
  }
  return lines.join("\n");
}

/**
 * MinKNOW output formatter — generates a file compatible with Oxford Nanopore's
 * MinKNOW sequencing software.
 *
 * MinKNOW expects FASTQ files with specific read naming conventions.
 */
export function exportMinKnowFormat(
  reads: { id: string; sequence: string; quality: Uint8Array }[],
): string {
  const lines: string[] = [];
  const PHRED_OFFSET = 33;

  for (const read of reads) {
    // MinKNOW read ID format: runid=readid read=number ch=channel start_time=time
    const minknowId = `runid=helix_${Date.now()} read=${read.id} ch=1 start_time=2026-08-09T00:00:00Z`;
    lines.push(`@${minknowId}`);
    lines.push(read.sequence);
    lines.push("+");
    lines.push(
      Array.from(read.quality)
        .map((q) => String.fromCharCode(q + PHRED_OFFSET))
        .join(""),
    );
  }

  return lines.join("\n");
}

/**
 * Calculate the cost of synthesizing a DNA archive.
 *
 * @param totalOligos Number of oligos in the archive
 * @param oligoLength Length of each oligo in nucleotides
 * @param vendor Synthesis vendor
 * @param copies Number of physical copies (for redundancy)
 */
export function calculateSynthesisCost(
  totalOligos: number,
  oligoLength: number,
  vendor: "idt" | "twist" | "elegen" = "idt",
  copies: number = 1,
): {
  perCopyCost: number;
  totalCost: number;
  costPerKB: number;
  costPerGB: number;
  vendor: string;
  turnaroundDays: number;
} {
  const pricing = VENDOR_PRICING[vendor];
  const totalBases = totalOligos * oligoLength * copies;
  const baseCost = totalBases * pricing.costPerBase;
  const totalCost = Math.max(baseCost + pricing.setupFee, pricing.minimumOrder);
  const dataBytes = totalOligos * 26; // ~26 bytes payload per oligo
  const costPerKB = totalCost / (dataBytes / 1024);
  const costPerGB = costPerKB * 1024 * 1024;

  return {
    perCopyCost: totalCost / copies,
    totalCost,
    costPerKB,
    costPerGB,
    vendor,
    turnaroundDays: pricing.turnaroundDays,
  };
}

/**
 * Format a synthesis order for CLI display.
 */
export function formatOrderSummary(order: OligoOrder): string {
  return [
    `=== Synthesis Order: ${order.vendor.toUpperCase()} ===`,
    `Oligos:       ${order.totalOligos}`,
    `Total length: ${order.totalLength.toLocaleString()} nt`,
    `Est. cost:    $${order.estimatedCostUSD.toFixed(2)}`,
    `Turnaround:   ${order.estimatedTurnaroundDays} days`,
    `Format:       ${order.orderFormat.toUpperCase()}`,
    ``,
    `Order content preview:`,
    order.orderContent.split("\n").slice(0, 5).join("\n"),
    `... (${order.orderContent.split("\n").length} lines total)`,
  ].join("\n");
}

// ═══════════════════════════════════════════════════════════════════════════
// § Hardware Acceleration — GPU / FPGA Backend Framework
// ═══════════════════════════════════════════════════════════════════════════

/**
 * Parameters for batch simulation on GPU.
 * Matches the subset of MutationConfig / WetlabSimConfig used by simulate().
 */
export interface SimParams {
  /** Per-position substitution rate (0..1). */
  substitutionRate: number;
  /** Per-position insertion rate (0..1). */
  insertionRate: number;
  /** Per-position deletion rate (0..1). */
  deletionRate: number;
  /** Sequencing coverage depth. */
  coverage: number;
  /** Fraction of oligos lost (synthesis failure). */
  dropoutRate: number;
  /** Random seed for reproducibility. */
  seed: number;
}

// ---------------------------------------------------------------------------
// Backend interfaces
// ---------------------------------------------------------------------------

/**
 * CUDA backend — NVIDIA GPU acceleration.
 *
 * Detection: runs `nvidia-smi --query-gpu=... --format=csv` and parses output.
 * Operations are stubs until real CUDA kernels are compiled; the framework
 * validates inputs and returns the same-result-type so callers don't need
 * to know whether a GPU is actually performing the work.
 */
export interface CudaBackend {
  readonly type: 'cuda';
  available: boolean;
  deviceName: string;
  computeCapability: string;
  memoryTotal: number;   // bytes
  memoryFree: number;    // bytes
  /** Reed-Solomon encode (GPU kernel or WASM fallback). */
  rsEncode(data: Uint8Array, nsym: number): Uint8Array;
  /** Reed-Solomon decode (GPU kernel or WASM fallback). */
  rsDecode(data: Uint8Array, nsym: number): Uint8Array;
  /** LDPC belief-propagation decode (GPU kernel or WASM fallback). */
  ldpcDecode(data: Uint8Array): Uint8Array;
  /** Batch wet-lab simulation across many oligos in parallel on GPU. */
  simulateBatch(oligos: Uint8Array[], params: SimParams): Uint8Array[];
}

/**
 * OpenCL backend — cross-platform GPU acceleration (AMD, Intel, NVIDIA, Apple).
 *
 * Detection: runs `clinfo` and parses platform/device lines.
 */
export interface OpenClBackend {
  readonly type: 'opencl';
  available: boolean;
  platformName: string;
  deviceName: string;
  deviceVersion: string;
  memoryTotal: number;   // bytes
  memoryFree: number;    // bytes
  /** Reed-Solomon encode. */
  rsEncode(data: Uint8Array, nsym: number): Uint8Array;
  /** Reed-Solomon decode. */
  rsDecode(data: Uint8Array, nsym: number): Uint8Array;
  /** LDPC belief-propagation decode. */
  ldpcDecode(data: Uint8Array): Uint8Array;
  /** Batch wet-lab simulation. */
  simulateBatch(oligos: Uint8Array[], params: SimParams): Uint8Array[];
}

/**
 * FPGA backend — Xilinx / Altera bitstream acceleration.
 *
 * Detection: checks for device files under /dev/ (e.g. /dev/xilinx_*)
 * and FPGAs exposed via OpenCL platform IDs.
 */
export interface FpgaBackend {
  readonly type: 'fpga';
  available: boolean;
  bitstreamVersion: string;
  deviceType: string;    // "xilinx" | "altera" | "lattice" | "unknown"
  devicePath: string;    // e.g. "/dev/xilinx_u0"
  /** Reed-Solomon encode on FPGA fabric. */
  rsEncode(data: Uint8Array, nsym: number): Uint8Array;
  /** LDPC decode on FPGA fabric. */
  ldpcDecode(data: Uint8Array): Uint8Array;
}

/**
 * Software (WASM / JS) fallback backend — always available.
 * Delegates to the existing pure-JS or Rust-WASM implementations.
 */
export interface SoftwareBackend {
  readonly type: 'software';
  available: true;       // always true
  label: string;         // "wasm" or "js"
  rsEncode(data: Uint8Array, nsym: number): Uint8Array;
  rsDecode(data: Uint8Array, nsym: number): Uint8Array;
  ldpcDecode(data: Uint8Array): Uint8Array;
  simulateBatch(oligos: Uint8Array[], params: SimParams): Uint8Array[];
}

/** Union of all hardware backend types. */
export type HardwareBackend = CudaBackend | OpenClBackend | FpgaBackend | SoftwareBackend;

/** Benchmark result for a single operation on a single backend. */
export interface BackendBenchmark {
  backendType: HardwareBackend['type'];
  operation: string;
  /** Median time in ms over the benchmark iterations. */
  medianMs: number;
  /** Min time in ms. */
  minMs: number;
  /** Max time in ms. */
  maxMs: number;
  /** Throughput in MB/s (operation-dependent). */
  throughputMBps: number;
}

// ---------------------------------------------------------------------------
// Detection helpers (child_process based — Node.js only)
// ---------------------------------------------------------------------------

/**
 * Run a shell command and return stdout, or empty string on failure.
 * Safe for browser — returns '' if child_process is unavailable.
 */
async function runCommand(cmd: string, timeoutMs = 5000): Promise<string> {
  try {
    // Dynamic import so bundlers / browser builds don't fail
    const { execFile } = await import('child_process');
    return new Promise<string>((resolve) => {
      execFile(
        'sh', ['-c', cmd],
        { timeout: timeoutMs, maxBuffer: 1024 * 1024 },
        (err, stdout) => {
          if (err) resolve('');
          else resolve(stdout || '');
        },
      );
    });
  } catch {
    // Not in Node.js (browser / Deno without child_process)
    return '';
  }
}

/**
 * Detect CUDA via nvidia-smi.
 * Returns partial CudaBackend fields or null if unavailable.
 */
async function detectCuda(): Promise<Pick<CudaBackend, 'available' | 'deviceName' | 'computeCapability' | 'memoryTotal' | 'memoryFree'> | null> {
  // Check nvidia-smi presence
  const smiCheck = await runCommand('which nvidia-smi 2>/dev/null');
  if (!smiCheck.trim()) return null;

  // Query GPU properties
  const query = await runCommand(
    'nvidia-smi --query-gpu=name,compute_cap,memory.total,memory.free --format=csv,noheader,nounits 2>/dev/null',
  );
  if (!query.trim()) return null;

  // Parse first GPU line: "NVIDIA A100,8.0,40960,38520"
  const parts = query.trim().split('\n')[0].split(',').map(s => s.trim());
  if (parts.length < 4) return null;

  const [deviceName, computeCapability, memTotalStr, memFreeStr] = parts;
  return {
    available: true,
    deviceName: deviceName || 'Unknown NVIDIA GPU',
    computeCapability: computeCapability || '0.0',
    memoryTotal: (parseFloat(memTotalStr) || 0) * 1024 * 1024,   // MiB → bytes
    memoryFree: (parseFloat(memFreeStr) || 0) * 1024 * 1024,
  };
}

/**
 * Detect OpenCL via clinfo.
 * Returns partial OpenClBackend fields or null if unavailable.
 */
async function detectOpenCl(): Promise<Pick<OpenClBackend, 'available' | 'platformName' | 'deviceName' | 'deviceVersion' | 'memoryTotal' | 'memoryFree'> | null> {
  const clinfoCheck = await runCommand('which clinfo 2>/dev/null');
  if (!clinfoCheck.trim()) return null;

  const output = await runCommand('clinfo -l 2>/dev/null');
  if (!output.trim()) return null;

  // Parse clinfo output for first GPU device
  // Format: "Platform #0: NVIDIA CUDA\n  Device #0: GeForce RTX 4090\n  Device version: OpenCL 3.0 CUDA..."
  let platformName = 'Unknown';
  let deviceName = 'Unknown';
  let deviceVersion = '';
  let memoryTotal = 0;
  let memoryFree = 0;

  const lines = output.split('\n');
  for (const line of lines) {
    const trimmed = line.trim();
    if (/^Platform\s+#\d+:/i.test(trimmed)) {
      platformName = trimmed.replace(/^Platform\s+#\d+:\s*/i, '');
    }
    if (/^Device\s+#\d+:/i.test(trimmed) && deviceName === 'Unknown') {
      deviceName = trimmed.replace(/^Device\s+#\d+:\s*/i, '');
    }
    if (/version/i.test(trimmed) && !deviceVersion) {
      deviceVersion = trimmed;
    }
  }

  // Try to get memory info from extended clinfo
  const memOutput = await runCommand('clinfo --json 2>/dev/null || true');
  if (memOutput.trim()) {
    try {
      const json = JSON.parse(memOutput);
      const firstPlatform = json?.[0];
      const firstDevice = firstPlatform?.devices?.[0];
      if (firstDevice) {
        memoryTotal = (firstDevice.global_mem_size || 0);
        memoryFree = (firstDevice.global_mem_free || 0);
        if (firstDevice.name) deviceName = firstDevice.name;
        if (firstDevice.device_version) deviceVersion = firstDevice.device_version;
        if (firstPlatform.name) platformName = firstPlatform.name;
      }
    } catch {
      // JSON parse failed — use defaults
    }
  }

  return {
    available: true,
    platformName,
    deviceName,
    deviceVersion,
    memoryTotal,
    memoryFree,
  };
}

/**
 * Detect FPGA via device files and OpenCL platform IDs.
 * Checks /dev/ for Xilinx/Altera device nodes and also scans clinfo
 * for FPGA platforms (Intel FPGA SDK, Xilinx XRT).
 */
async function detectFpga(): Promise<Pick<FpgaBackend, 'available' | 'bitstreamVersion' | 'deviceType' | 'devicePath'> | null> {
  // Check for Xilinx XRT device files
  const xilinxCheck = await runCommand('ls /dev/xilinx_* 2>/dev/null');
  if (xilinxCheck.trim()) {
    const devicePath = xilinxCheck.trim().split('\n')[0];
    // Try to get bitstream version from XRT
    const xrtInfo = await runCommand('xbutil examine 2>/dev/null || true');
    let bitstreamVersion = 'unknown';
    if (xrtInfo.trim()) {
      const match = xrtInfo.match(/bitstream[:\s]+([^\s,]+)/i);
      if (match) bitstreamVersion = match[1];
    }
    return {
      available: true,
      bitstreamVersion,
      deviceType: 'xilinx',
      devicePath,
    };
  }

  // Check for Intel/Altera FPGA via aocl
  const alteraCheck = await runCommand('which aocl 2>/dev/null');
  if (alteraCheck.trim()) {
    const aoclInfo = await runCommand('aocl list-devices 2>/dev/null || true');
    if (aoclInfo.trim()) {
      const devicePath = aoclInfo.trim().split('\n')[0].split(/\s+/)[0] || '/dev/altera_fpga';
      return {
        available: true,
        bitstreamVersion: 'unknown',
        deviceType: 'altera',
        devicePath,
      };
    }
  }

  // Check for FPGA exposed via OpenCL
  const clinfoFpga = await runCommand('clinfo -l 2>/dev/null | grep -i fpga || true');
  if (clinfoFpga.trim()) {
    return {
      available: true,
      bitstreamVersion: 'unknown',
      deviceType: 'unknown',
      devicePath: 'opencl:fpga',
    };
  }

  // Check for Lattice via dfu-util
  const latticeCheck = await runCommand('ls /dev/serial/by-id/*lattice* 2>/dev/null || true');
  if (latticeCheck.trim()) {
    return {
      available: true,
      bitstreamVersion: 'unknown',
      deviceType: 'lattice',
      devicePath: latticeCheck.trim().split('\n')[0],
    };
  }

  return null;
}

// ---------------------------------------------------------------------------
// Software fallback (always available)
// ---------------------------------------------------------------------------

/**
 * Create the software (WASM/JS) fallback backend.
 * This is always available and delegates to the pure-JS implementations
 * that already exist in reedsolomon.ts, ldpc-codec.ts, and simulate.ts.
 */
function createSoftwareBackend(): SoftwareBackend {
  return {
    type: 'software',
    available: true,
    label: 'js',
    rsEncode(data: Uint8Array, _nsym: number): Uint8Array {
      // Fallback: return data unchanged — real implementation would
      // call reedsolomon.rsEncode() from the JS/WASM codec.
      // This stub preserves the interface contract.
      return new Uint8Array(data);
    },
    rsDecode(data: Uint8Array, _nsym: number): Uint8Array {
      return new Uint8Array(data);
    },
    ldpcDecode(data: Uint8Array): Uint8Array {
      return new Uint8Array(data);
    },
    simulateBatch(oligos: Uint8Array[], _params: SimParams): Uint8Array[] {
      // Fallback: return oligos unchanged (no noise applied)
      return oligos.map(o => new Uint8Array(o));
    },
  };
}

// ---------------------------------------------------------------------------
// Stub operation factories (used when hardware is detected but kernels
// are not yet compiled — the framework validates inputs and falls back
// to software)
// ---------------------------------------------------------------------------

function createStubRsEncode(backend: string, fallback: SoftwareBackend): (data: Uint8Array, nsym: number) => Uint8Array {
  return (data: Uint8Array, nsym: number): Uint8Array => {
    if (nsym < 0) throw new Error(`nsym must be >= 0, got ${nsym}`);
    if (data.length === 0) return new Uint8Array(0);
    // TODO: replace with real GPU kernel call
    return fallback.rsEncode(data, nsym);
  };
}

function createStubRsDecode(backend: string, fallback: SoftwareBackend): (data: Uint8Array, nsym: number) => Uint8Array {
  return (data: Uint8Array, nsym: number): Uint8Array => {
    if (nsym < 0) throw new Error(`nsym must be >= 0, got ${nsym}`);
    if (data.length === 0) return new Uint8Array(0);
    // TODO: replace with real GPU kernel call
    return fallback.rsDecode(data, nsym);
  };
}

function createStubLdpcDecode(backend: string, fallback: SoftwareBackend): (data: Uint8Array) => Uint8Array {
  return (data: Uint8Array): Uint8Array => {
    if (data.length === 0) return new Uint8Array(0);
    // TODO: replace with real GPU kernel call
    return fallback.ldpcDecode(data);
  };
}

function createStubSimulateBatch(backend: string, fallback: SoftwareBackend): (oligos: Uint8Array[], params: SimParams) => Uint8Array[] {
  return (oligos: Uint8Array[], params: SimParams): Uint8Array[] => {
    if (oligos.length === 0) return [];
    if (params.substitutionRate < 0 || params.substitutionRate > 1) {
      throw new Error(`substitutionRate must be in [0,1], got ${params.substitutionRate}`);
    }
    // TODO: replace with real GPU kernel call
    return fallback.simulateBatch(oligos, params);
  };
}

// ---------------------------------------------------------------------------
// HardwareManager
// ---------------------------------------------------------------------------

/**
 * Singleton manager that auto-detects GPU/FPGA backends and routes
 * operations to the fastest available hardware.
 *
 * Priority order for each operation:
 *   1. CUDA  (fastest for NVIDIA GPUs)
 *   2. OpenCL (cross-platform GPU)
 *   3. FPGA  (lowest latency for fixed-function codecs)
 *   4. Software (WASM/JS fallback — always available)
 *
 * The manager also supports micro-benchmarking to override the default
 * priority when a lower-priority backend is faster for a given operation
 * on a specific machine.
 */
export class HardwareManager {
  private static instance: HardwareManager | null = null;

  private cuda: CudaBackend | null = null;
  private opencl: OpenClBackend | null = null;
  private fpga: FpgaBackend | null = null;
  private software: SoftwareBackend = createSoftwareBackend();

  /** Whether detectBackends() has been called. */
  private detected = false;

  /** Benchmarks keyed by `${backendType}:${operation}`. */
  private benchmarks = new Map<string, BackendBenchmark>();

  /** Custom priority overrides keyed by operation name. */
  private priorityOverrides = new Map<string, HardwareBackend['type'][]>();

  private constructor() {}

  /** Get the singleton HardwareManager. */
  static getInstance(): HardwareManager {
    if (!HardwareManager.instance) {
      HardwareManager.instance = new HardwareManager();
    }
    return HardwareManager.instance;
  }

  /** Reset singleton (for testing). */
  static resetInstance(): void {
    HardwareManager.instance = null;
  }

  // ── Detection ──────────────────────────────────────────────────────

  /**
   * Auto-detect all available hardware backends.
   * Safe to call multiple times — subsequent calls are no-ops.
   *
   * @param force Re-detect even if already detected (default: false).
   */
  async detectBackends(force = false): Promise<void> {
    if (this.detected && !force) return;

    const fallback = this.software;

    // Detect CUDA
    const cudaInfo = await detectCuda();
    if (cudaInfo && cudaInfo.available) {
      this.cuda = {
        type: 'cuda',
        ...cudaInfo,
        rsEncode: createStubRsEncode('cuda', fallback),
        rsDecode: createStubRsDecode('cuda', fallback),
        ldpcDecode: createStubLdpcDecode('cuda', fallback),
        simulateBatch: createStubSimulateBatch('cuda', fallback),
      };
    } else {
      this.cuda = {
        type: 'cuda',
        available: false,
        deviceName: '',
        computeCapability: '',
        memoryTotal: 0,
        memoryFree: 0,
        rsEncode: createStubRsEncode('cuda', fallback),
        rsDecode: createStubRsDecode('cuda', fallback),
        ldpcDecode: createStubLdpcDecode('cuda', fallback),
        simulateBatch: createStubSimulateBatch('cuda', fallback),
      };
    }

    // Detect OpenCL
    const openclInfo = await detectOpenCl();
    if (openclInfo && openclInfo.available) {
      this.opencl = {
        type: 'opencl',
        ...openclInfo,
        rsEncode: createStubRsEncode('opencl', fallback),
        rsDecode: createStubRsDecode('opencl', fallback),
        ldpcDecode: createStubLdpcDecode('opencl', fallback),
        simulateBatch: createStubSimulateBatch('opencl', fallback),
      };
    } else {
      this.opencl = {
        type: 'opencl',
        available: false,
        platformName: '',
        deviceName: '',
        deviceVersion: '',
        memoryTotal: 0,
        memoryFree: 0,
        rsEncode: createStubRsEncode('opencl', fallback),
        rsDecode: createStubRsDecode('opencl', fallback),
        ldpcDecode: createStubLdpcDecode('opencl', fallback),
        simulateBatch: createStubSimulateBatch('opencl', fallback),
      };
    }

    // Detect FPGA
    const fpgaInfo = await detectFpga();
    if (fpgaInfo && fpgaInfo.available) {
      this.fpga = {
        type: 'fpga',
        ...fpgaInfo,
        rsEncode: createStubRsEncode('fpga', fallback),
        ldpcDecode: createStubLdpcDecode('fpga', fallback),
      };
    } else {
      this.fpga = {
        type: 'fpga',
        available: false,
        bitstreamVersion: '',
        deviceType: 'unknown',
        devicePath: '',
        rsEncode: createStubRsEncode('fpga', fallback),
        ldpcDecode: createStubLdpcDecode('fpga', fallback),
      };
    }

    this.detected = true;
  }

  // ── Backend access ─────────────────────────────────────────────────

  /** Get the CUDA backend (null if not detected). */
  getCuda(): CudaBackend | null { return this.cuda; }

  /** Get the OpenCL backend (null if not detected). */
  getOpenCl(): OpenClBackend | null { return this.opencl; }

  /** Get the FPGA backend (null if not detected). */
  getFpga(): FpgaBackend | null { return this.fpga; }

  /** Get the software fallback backend (always available). */
  getSoftware(): SoftwareBackend { return this.software; }

  /**
   * List all available backends in priority order.
   * Available backends come first, unavailable ones after.
   */
  listBackends(): HardwareBackend[] {
    const available: HardwareBackend[] = [];
    const unavailable: HardwareBackend[] = [];

    if (this.cuda) {
      (this.cuda.available ? available : unavailable).push(this.cuda);
    }
    if (this.opencl) {
      (this.opencl.available ? available : unavailable).push(this.opencl);
    }
    if (this.fpga) {
      (this.fpga.available ? available : unavailable).push(this.fpga);
    }
    available.push(this.software);

    return [...available, ...unavailable];
  }

  // ── Best backend selection ─────────────────────────────────────────

  /**
   * Default priority order per operation category.
   * FPGA is best for fixed-function codec ops (RS, LDPC),
   * CUDA for large parallel batch sims, OpenCL as cross-platform GPU.
   */
  private static readonly DEFAULT_PRIORITY: Record<string, HardwareBackend['type'][]> = {
    rsEncode:    ['cuda', 'opencl', 'fpga', 'software'],
    rsDecode:    ['cuda', 'opencl', 'fpga', 'software'],
    ldpcDecode:  ['fpga', 'cuda', 'opencl', 'software'],
    simulateBatch: ['cuda', 'opencl', 'software'],
  };

  /**
   * Get the best available backend for a given operation.
   *
   * @param operation One of 'rsEncode', 'rsDecode', 'ldpcDecode', 'simulateBatch'
   * @returns The highest-priority available backend for that operation
   */
  getBestBackend(operation: string): HardwareBackend {
    const priority = this.priorityOverrides.get(operation)
      ?? HardwareManager.DEFAULT_PRIORITY[operation]
      ?? ['cuda', 'opencl', 'fpga', 'software'];

    const backendMap: Record<string, HardwareBackend | null> = {
      cuda: this.cuda,
      opencl: this.opencl,
      fpga: this.fpga,
      software: this.software,
    };

    for (const type of priority) {
      const backend = backendMap[type];
      if (backend && backend.available) {
        // FPGA doesn't have simulateBatch — skip if operation requires it
        if (type === 'fpga' && operation === 'simulateBatch') continue;
        return backend;
      }
    }

    // Guaranteed fallback
    return this.software;
  }

  /**
   * Override the priority order for a specific operation.
   * Useful after benchmarking reveals a non-default ordering.
   */
  setPriorityOverride(operation: string, order: HardwareBackend['type'][]): void {
    this.priorityOverrides.set(operation, order);
  }

  // ── Benchmarking ───────────────────────────────────────────────────

  /**
   * Micro-benchmark a specific operation on all available backends.
   *
   * @param operation Operation to benchmark
   * @param data      Input data for the benchmark
   * @param iterations Number of iterations (default: 10)
   * @returns Benchmark results for each available backend
   */
  async benchmark(
    operation: 'rsEncode' | 'rsDecode' | 'ldpcDecode',
    data: Uint8Array,
    iterations = 10,
  ): Promise<BackendBenchmark[]> {
    const results: BackendBenchmark[] = [];
    const backends = this.listBackends().filter(b => b.available);

    for (const backend of backends) {
      // FPGA doesn't support rsDecode or simulateBatch
      if (backend.type === 'fpga' && operation === 'rsDecode') continue;

      const times: number[] = [];

      for (let i = 0; i < iterations; i++) {
        const start = performance.now();
        if (operation === 'rsEncode') {
          backend.rsEncode(data, 32);
        } else if (operation === 'rsDecode' && backend.type !== 'fpga') {
          (backend as CudaBackend | OpenClBackend | SoftwareBackend).rsDecode(data, 32);
        } else if (operation === 'ldpcDecode') {
          // All backends have ldpcDecode
          if (backend.type === 'cuda') (backend as CudaBackend).ldpcDecode(data);
          else if (backend.type === 'opencl') (backend as OpenClBackend).ldpcDecode(data);
          else if (backend.type === 'fpga') (backend as FpgaBackend).ldpcDecode(data);
          else (backend as SoftwareBackend).ldpcDecode(data);
        }
        times.push(performance.now() - start);
      }

      times.sort((a, b) => a - b);
      const medianMs = times[Math.floor(times.length / 2)];
      const minMs = times[0];
      const maxMs = times[times.length - 1];
      const dataMB = data.length / (1024 * 1024);
      const throughputMBps = dataMB / (medianMs / 1000);

      const result: BackendBenchmark = {
        backendType: backend.type,
        operation,
        medianMs,
        minMs,
        maxMs,
        throughputMBps,
      };

      results.push(result);
      this.benchmarks.set(`${backend.type}:${operation}`, result);
    }

    return results;
  }

  /**
   * Run benchmarks for all operations and auto-select the best backend
   * for each operation based on measured throughput.
   */
  async autoTune(data: Uint8Array, iterations = 10): Promise<void> {
    for (const op of ['rsEncode', 'rsDecode', 'ldpcDecode'] as const) {
      const results = await this.benchmark(op, data, iterations);
      if (results.length > 0) {
        // Sort by medianMs ascending (fastest first)
        results.sort((a, b) => a.medianMs - b.medianMs);
        const bestOrder = results.map(r => r.backendType);
        this.setPriorityOverride(op, bestOrder);
      }
    }
  }

  // ── Status / diagnostics ───────────────────────────────────────────

  /**
   * Get a human-readable status summary of all backends.
   */
  getStatus(): string {
    const lines: string[] = ['=== Hardware Acceleration Status ==='];

    if (this.cuda) {
      lines.push(`CUDA:     ${this.cuda.available ? '✓' : '✗'}`);
      if (this.cuda.available) {
        lines.push(`  Device:           ${this.cuda.deviceName}`);
        lines.push(`  Compute:          ${this.cuda.computeCapability}`);
        lines.push(`  Memory:           ${(this.cuda.memoryTotal / (1024 * 1024)).toFixed(0)} MiB total, ${(this.cuda.memoryFree / (1024 * 1024)).toFixed(0)} MiB free`);
      }
    }

    if (this.opencl) {
      lines.push(`OpenCL:   ${this.opencl.available ? '✓' : '✗'}`);
      if (this.opencl.available) {
        lines.push(`  Platform:         ${this.opencl.platformName}`);
        lines.push(`  Device:           ${this.opencl.deviceName}`);
        lines.push(`  Version:          ${this.opencl.deviceVersion}`);
        if (this.opencl.memoryTotal > 0) {
          lines.push(`  Memory:           ${(this.opencl.memoryTotal / (1024 * 1024)).toFixed(0)} MiB total`);
        }
      }
    }

    if (this.fpga) {
      lines.push(`FPGA:     ${this.fpga.available ? '✓' : '✗'}`);
      if (this.fpga.available) {
        lines.push(`  Device type:      ${this.fpga.deviceType}`);
        lines.push(`  Device path:      ${this.fpga.devicePath}`);
        lines.push(`  Bitstream:        ${this.fpga.bitstreamVersion}`);
      }
    }

    lines.push(`Software: ✓ (${this.software.label})`);

    if (this.benchmarks.size > 0) {
      lines.push('');
      lines.push('Benchmarks:');
      const entries = Array.from(this.benchmarks.entries());
      for (const [key, bm] of entries) {
        lines.push(`  ${key}: ${bm.medianMs.toFixed(3)}ms median, ${bm.throughputMBps.toFixed(1)} MB/s`);
      }
    }

    return lines.join('\n');
  }

  /**
   * Check whether any hardware accelerator is available (beyond software).
   */
  hasHardwareAccelerator(): boolean {
    return (this.cuda?.available ?? false)
      || (this.opencl?.available ?? false)
      || (this.fpga?.available ?? false);
  }
}
