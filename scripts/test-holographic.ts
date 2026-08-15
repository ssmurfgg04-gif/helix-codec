// Test the Holographic DNA Sharding Codec.
import {
  holographicEncode,
  holographicDecode,
  simulateShardLoss,
  HolographicConfig,
} from "../src/lib/dna/holographic";

function assert(cond: boolean, msg: string) {
  if (!cond) {
    console.error("FAIL:", msg);
    process.exit(1);
  } else {
    console.log("PASS:", msg);
  }
}

function equalArrays(a: Uint8Array, b: Uint8Array): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

console.log("\n=== Holographic Codec Tests ===\n");

// Test 1: Basic encode/decode (no shard loss)
{
  const data = new Uint8Array(100);
  for (let i = 0; i < 100; i++) data[i] = (i * 7 + 13) & 0xff;
  const config: HolographicConfig = { dataShards: 10, totalShards: 15, blockSize: 10 };
  const encoding = holographicEncode(data, config);
  console.log(`Encoded ${data.length} bytes into ${encoding.shards.length} shards (${encoding.numBlocks} blocks)`);

  const recovered = holographicDecode(encoding.shards, encoding);
  assert(equalArrays(recovered, data), "Basic encode/decode (no loss) matches");
}

// Test 2: 1/3 shard loss (5 of 15 shards lost), should still recover
{
  const data = new Uint8Array(100);
  for (let i = 0; i < 100; i++) data[i] = (i * 7 + 13) & 0xff;
  const config: HolographicConfig = { dataShards: 10, totalShards: 15, blockSize: 10 };
  const encoding = holographicEncode(data, config);

  const result = simulateShardLoss(encoding, 1 / 3, 42);
  console.log(`\n1/3 shard loss: ${result.shardsLost} lost, ${result.shardsAvailable} available`);
  assert(result.recoverySuccessful, "1/3 shard loss recovery succeeds");
  assert(result.partialRecoveryRate === 1.0, "1/3 shard loss: 100% recovery");
}

// Test 3: Exactly K shards (10 of 15), should still recover
{
  const data = new Uint8Array(100);
  for (let i = 0; i < 100; i++) data[i] = (i * 7 + 13) & 0xff;
  const config: HolographicConfig = { dataShards: 10, totalShards: 15, blockSize: 10 };
  const encoding = holographicEncode(data, config);

  // Use only first 10 shards
  const available = encoding.shards.slice(0, 10);
  const recovered = holographicDecode(available, encoding);
  assert(equalArrays(recovered, data), "Exactly K shards recovery matches");
}

// Test 4: Less than K shards (9 of 15), should fail
{
  const data = new Uint8Array(100);
  for (let i = 0; i < 100; i++) data[i] = (i * 7 + 13) & 0xff;
  const config: HolographicConfig = { dataShards: 10, totalShards: 15, blockSize: 10 };
  const encoding = holographicEncode(data, config);

  try {
    const available = encoding.shards.slice(0, 9);
    holographicDecode(available, encoding);
    assert(false, "Should fail with < K shards");
  } catch (e) {
    assert(true, `Fails correctly with < K shards: ${(e as Error).message}`);
  }
}

// Test 5: Larger data, multiple blocks
{
  const data = new Uint8Array(1000);
  for (let i = 0; i < 1000; i++) data[i] = (i * 31 + 17) & 0xff;
  const config: HolographicConfig = { dataShards: 20, totalShards: 30, blockSize: 20 };
  const encoding = holographicEncode(data, config);
  console.log(`\nLarger test: ${data.length} bytes -> ${encoding.shards.length} shards, ${encoding.numBlocks} blocks`);

  const result = simulateShardLoss(encoding, 0.3, 42); // lose 30%
  console.log(`30% loss: ${result.shardsLost} lost, recovered=${result.recoverySuccessful}`);
  assert(result.recoverySuccessful, "Larger data 30% loss recovery succeeds");
}

// Test 6: Compare overhead: holographic 1.5x vs RS 2x
{
  console.log("\n--- Overhead comparison ---");
  const data = new Uint8Array(500);
  for (let i = 0; i < 500; i++) data[i] = (i * 31 + 17) & 0xff;

  // Holographic: K=10, N=15 (1.5x overhead)
  const holoConfig: HolographicConfig = { dataShards: 10, totalShards: 15, blockSize: 10 };
  const holoEnc = holographicEncode(data, holoConfig);
  console.log(`Holographic: ${holoEnc.shards.length} shards for ${data.length} bytes (overhead ${holoEnc.shards.length / holoConfig.dataShards}x)`);

  // Test recovery at various loss rates
  for (const loss of [0.0, 0.1, 0.2, 0.3, 0.4]) {
    const result = simulateShardLoss(holoEnc, loss, 42);
    console.log(`  Loss ${loss * 100}%: ${result.recoverySuccessful ? "RECOVERED" : "FAILED"} (${result.shardsAvailable}/${holoEnc.shards.length} shards)`);
  }
}

console.log("\nAll holographic codec tests passed.");
