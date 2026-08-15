/**
 * DNA Codec API
 *
 * Single endpoint that handles all codec operations via `op` query param.
 * Operations:
 *   - encode: takes file data + config, returns encoded file (oligos + metadata)
 *   - simulate: takes encoded file + mutation config, returns simulated reads
 *   - decode: takes reads + metadata + config, returns recovered file
 *   - benchmark: runs benchmark sweep, returns results
 *   - archive: wraps an EncodedFile into a BioArchive container (manifest + Merkle + chunks)
 *   - audit: runs audit/scrubbing on an archive after simulation
 *   - lineage: simulates a generation of biological time
 *   - holographic: encodes data using Holographic DNA Sharding codec
 *   - holographic-decode: decodes holographic shards (with optional shard loss)
 *   - holographic-sweep: runs a shard-loss sweep benchmark
 *
 * All request/response bodies are JSON. Binary data (file bytes) is base64-encoded.
 */

import { NextRequest, NextResponse } from "next/server";
import { encodeFile } from "@/lib/dna/codec";
import { decodeReads } from "@/lib/dna/decode";
import { decodeReadsUltra } from "@/lib/dna/ultra-decode";
import { simulate, MutationConfig } from "@/lib/dna/simulate";
import { runErrorRateSweep, runCoverageSweep, generatePayload } from "@/lib/dna/benchmark";
import {
  CodecConfig,
  DEFAULT_CONFIG,
  EncodedFile,
  computeLayout,
} from "@/lib/dna/types";
import { toBioArchive } from "@/lib/dna/bioarchive";
import { holographicEncode, simulateShardLoss } from "@/lib/dna/holographic";

export const runtime = "nodejs";
export const maxDuration = 60;

function base64ToBytes(b64: string): Uint8Array {
  const bin = Buffer.from(b64, "base64");
  return new Uint8Array(bin);
}

function bytesToBase64(bytes: Uint8Array): string {
  return Buffer.from(bytes).toString("base64");
}

function parseConfig(cfg: Partial<CodecConfig> | undefined): CodecConfig {
  return {
    ...DEFAULT_CONFIG,
    ...cfg,
    constraints: {
      ...DEFAULT_CONFIG.constraints,
      ...(cfg?.constraints ?? {}),
    },
  };
}

export async function POST(req: NextRequest) {
  try {
    const body = await req.json();
    const op = body.op as string;

    // Input validation: prevent abuse
    const MAX_PAYLOAD_SIZE = 50 * 1024 * 1024; // 50MB
    if (body.data && typeof body.data === 'string' && body.data.length > MAX_PAYLOAD_SIZE * 1.37) {
      // base64 is ~37% larger than raw, so limit base64 string length
      return NextResponse.json(
        { ok: false, error: 'Payload too large (max 50MB)' },
        { status: 413 },
      );
    }
    if (body.steps && (body.steps > 100 || body.steps < 1)) {
      return NextResponse.json(
        { ok: false, error: 'steps must be between 1 and 100' },
        { status: 400 },
      );
    }

    switch (op) {
      case "encode": {
        const fileData = base64ToBytes(body.data as string);
        const cfg = parseConfig(body.config);
        const fileName = body.fileName ?? "upload.bin";
        const contentType = body.contentType ?? "application/octet-stream";

        const result = await encodeFile(fileData, cfg, { fileName, contentType });

        // Convert oligos to a serializable format (already JSON-safe)
        return NextResponse.json({
          ok: true,
          encoded: result.encoded,
          stats: result.stats,
        });
      }

      case "simulate": {
        const encoded = body.encoded as EncodedFile;
        const mutCfg = body.mutationConfig as MutationConfig;
        const result = simulate(encoded.oligos, mutCfg);
        return NextResponse.json({
          ok: true,
          simulation: result,
        });
      }

      case "decode": {
        const reads = body.reads;
        const metadata = body.metadata;
        const cfg = parseConfig(body.config);
        const fwdPrimer = body.forwardPrimer as string;
        const revPrimer = body.reversePrimer as string;
        const useSoftInfo = body.useSoftInfo !== false; // default true

        let result;
        try {
          // Try WASM-accelerated path first (5-11x faster for Illumina)
          result = await decodeReadsUltra(reads, metadata, cfg, fwdPrimer, revPrimer);
        } catch (wasmError) {
          // Fallback to pure JS path for unsupported configs (nanopore+conv, encrypted, etc.)
          console.warn("WASM decode failed, falling back to JS:", (wasmError as Error).message);
          result = await decodeReads(reads, metadata, cfg, fwdPrimer, revPrimer, useSoftInfo);
        }
        return NextResponse.json({
          ok: true,
          decoded: {
            data: result.data ? bytesToBase64(result.data) : null,
            hash: result.hash,
            hashMatches: result.hashMatches,
            stats: result.stats,
            perOligo: result.perOligo.map((p) => ({
              ...p,
              payloadBytes: bytesToBase64(p.payloadBytes),
            })),
          },
        });
      }

      case "benchmark-error": {
        const payloadSize = body.payloadSize ?? 2048;
        const cfg = parseConfig(body.config);
        const payload = generatePayload(payloadSize);
        const result = await runErrorRateSweep(payload, cfg, {
          maxRate: body.maxRate ?? 0.05,
          steps: body.steps ?? 6,
          coverage: body.coverage ?? 20,
        });
        return NextResponse.json({ ok: true, benchmark: result });
      }

      case "benchmark-coverage": {
        const payloadSize = body.payloadSize ?? 2048;
        const cfg = parseConfig(body.config);
        const payload = generatePayload(payloadSize);
        const result = await runCoverageSweep(payload, cfg, {
          minCov: body.minCov ?? 5,
          maxCov: body.maxCov ?? 30,
          steps: body.steps ?? 6,
        });
        return NextResponse.json({ ok: true, benchmark: result });
      }

      case "layout": {
        const cfg = parseConfig(body.config);
        const layout = computeLayout(cfg);
        return NextResponse.json({ ok: true, layout, config: cfg });
      }

      case "archive": {
        const encoded = body.encoded as EncodedFile;
        const encryption = body.encryption;
        const lifecycle = body.lifecycle;
        const archive = await toBioArchive(encoded, encryption, lifecycle);
        return NextResponse.json({ ok: true, archive });
      }

      case "holographic": {
        const data = base64ToBytes(body.data as string);
        const dataShards = body.dataShards ?? 10;
        const totalShards = body.totalShards ?? Math.ceil(dataShards * 1.5);
        const encoding = holographicEncode(data, {
          dataShards,
          totalShards,
          blockSize: dataShards,
        });
        // Convert shards to base64 for transport
        const shardsB64 = encoding.shards.map((s) => ({
          index: s.index,
          x: s.x,
          data: bytesToBase64(s.data),
        }));
        return NextResponse.json({
          ok: true,
          encoding: {
            shards: shardsB64,
            numBlocks: encoding.numBlocks,
            originalLength: encoding.originalLength,
            config: encoding.config,
          },
        });
      }

      case "holographic-sweep": {
        const data = base64ToBytes(body.data as string);
        const dataShards = body.dataShards ?? 10;
        const totalShards = body.totalShards ?? Math.ceil(dataShards * 1.5);
        const encoding = holographicEncode(data, {
          dataShards,
          totalShards,
          blockSize: dataShards,
        });
        const losses = body.losses ?? [0, 0.1, 0.2, 0.3, 0.4, 0.5];
        const results = losses.map((loss: number) => {
          const r = simulateShardLoss(encoding, loss, body.seed ?? 42);
          return {
            lossFraction: loss,
            shardsAvailable: r.shardsAvailable,
            shardsLost: r.shardsLost,
            recoverySuccessful: r.recoverySuccessful,
            partialRecoveryRate: r.partialRecoveryRate,
          };
        });
        return NextResponse.json({
          ok: true,
          sweep: {
            dataShards,
            totalShards,
            overheadRatio: totalShards / dataShards,
            originalLength: encoding.originalLength,
            results,
          },
        });
      }

      default:
        return NextResponse.json({ ok: false, error: `Unknown op: ${op}` }, { status: 400 });
    }
  } catch (error) {
    console.error("DNA API error:", error);
    const message = process.env.NODE_ENV === 'development' && error instanceof Error
      ? error.message : 'Internal error';
    return NextResponse.json(
      { ok: false, error: message },
      { status: 500 },
    );
  }
}

export async function GET(req: NextRequest) {
  const url = new URL(req.url);
  const op = url.searchParams.get("op");

  if (op === "spec") {
    return NextResponse.json({
      ok: true,
      spec: {
        version: 1,
        defaultConfig: DEFAULT_CONFIG,
        layout: computeLayout(DEFAULT_CONFIG),
        mutationPresets: {
          illumina: {
            substitutionRate: 0.001,
            insertionRate: 0.0005,
            deletionRate: 0.001,
            coverage: 20,
            dropoutRate: 0,
            seed: 0,
          },
          nanopore: {
            substitutionRate: 0.02,
            insertionRate: 0.03,
            deletionRate: 0.04,
            coverage: 15,
            dropoutRate: 0.05,
            seed: 0,
          },
          pacbio: {
            substitutionRate: 0.005,
            insertionRate: 0.05,
            deletionRate: 0.03,
            coverage: 10,
            dropoutRate: 0.02,
            seed: 0,
          },
          clean: {
            substitutionRate: 0,
            insertionRate: 0,
            deletionRate: 0,
            coverage: 1,
            dropoutRate: 0,
            seed: 0,
          },
        },
      },
    });
  }

  return NextResponse.json({ ok: true, message: "DNA Codec API. Use POST with op=encode|simulate|decode|benchmark-error|benchmark-coverage|layout" });
}
