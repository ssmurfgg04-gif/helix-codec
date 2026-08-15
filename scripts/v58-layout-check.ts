import { ULTIMATE_V55_DENSITY_CONFIG } from "../src/lib/dna/presets";
import { computeLayoutAuto } from "../src/lib/dna/types";

const layout = computeLayoutAuto(ULTIMATE_V55_DENSITY_CONFIG);
console.log("config:", JSON.stringify({
  oligoLength: ULTIMATE_V55_DENSITY_CONFIG.oligoLength,
  mappingMode: ULTIMATE_V55_DENSITY_CONFIG.mappingMode,
  innerParityBytes: ULTIMATE_V55_DENSITY_CONFIG.innerParityBytes,
}, null, 2));
console.log("layout:", JSON.stringify(layout));
const innerN = layout.addressBytes + layout.payloadBytes + layout.innerParityBytes;
console.log("innerN:", innerN, "totalInnerBytes:", layout.totalInnerBytes);

// Check what capacity the encoder sees
const innerDnaLen = layout.totalInnerBytes * 4;
const ARITH_CAPACITY_RATE = 1.95;
const blockSize = Math.floor(innerDnaLen / 2);
const bTotal = Math.max(2, Math.floor((blockSize * ARITH_CAPACITY_RATE) / 8));
const bData = bTotal - 1;
const numBlocks = Math.floor(innerDnaLen / blockSize);
const cap = numBlocks * bData;
console.log(`blockSize=${blockSize}, bTotal=${bTotal}, bData=${bData}, numBlocks=${numBlocks}, capacity=${cap}`);
console.log(`innerN=${innerN}, capacity=${cap}, fits=${cap >= innerN}`);
