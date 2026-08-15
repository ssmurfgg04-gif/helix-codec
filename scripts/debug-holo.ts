import { holographicEncode, holographicDecode } from "../src/lib/dna/holographic";

const data = new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
const config = { dataShards: 10, totalShards: 15, blockSize: 10 };
const enc = holographicEncode(data, config);
console.log("Shards:");
for (const s of enc.shards) {
  console.log(`  x=${s.x}: ${Array.from(s.data).join(",")}`);
}
const recovered = holographicDecode(enc.shards, enc);
console.log("Original: ", Array.from(data).join(","));
console.log("Recovered:", Array.from(recovered).join(","));
console.log("Match:", JSON.stringify(Array.from(data)) === JSON.stringify(Array.from(recovered)));
