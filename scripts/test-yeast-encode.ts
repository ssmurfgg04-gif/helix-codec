import { encodeFile } from "../src/lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

const content = fs.readFileSync(path.join(__dirname, "..", "datasets", "medium", "yeast.fa"), "utf-8");
const seq = content.split("\n").filter(l => !l.startsWith(">") && l.trim()).join("");
const data = new Uint8Array(Buffer.from(seq, "utf-8"));
console.log(`Yeast data: ${data.length} bytes`);
console.log(`Encoding at ${new Date().toISOString()}`);

const t0 = Date.now();
const result = await encodeFile(data, V51_DEFAULT_CONFIG, { fileName: "yeast.fa", contentType: "application/octet-stream" });
const encMs = Date.now() - t0;

console.log(`Encode done: ${encMs}ms, ${result.stats.oligoCount} oligos, density=${result.stats.netDensityBitsPerNt.toFixed(3)}, retries=${result.stats.screeningRetries}`);

// Check constraints
let gcV = 0, hpV = 0, gcMin = 1, gcMax = 0, maxHp = 0;
for (const o of result.encoded.oligos) {
  if (o.gc < 0.4 || o.gc > 0.6) gcV++;
  if (o.maxHomopolymer > 3) hpV++;
  if (o.gc < gcMin) gcMin = o.gc;
  if (o.gc > gcMax) gcMax = o.gc;
  if (o.maxHomopolymer > maxHp) maxHp = o.maxHomopolymer;
}
console.log(`Constraints: GC=[${gcMin.toFixed(3)},${gcMax.toFixed(3)}], maxHp=${maxHp}, gcViol=${gcV}, hpViol=${hpV}`);
