/**
 * Deep decode diagnostic: trace the full decode pipeline
 */
import { encodeFile } from "../src/lib/dna/codec";
import { decodeReads } from "../src/lib/dna/decode";
import { simulate, PRESET_CLEAN } from "../src/lib/dna/simulate";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { bytesToSplitConstrainedDna, splitConstrainedDnaToBytesWithErasure } from "../src/lib/dna/constrained-mapping";
import { whitenAddress, unwhitenAddress, homopolymerSafeDnaToAddress } from "../src/lib/dna/mapping";
import { crc16Bytes, verifyCrc16 } from "../src/lib/dna/crc16";

async function main() {
  const data = new Uint8Array(Buffer.from("Hello World!"));
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  
  const enc = await encodeFile(data, cfg, { fileName: "test", contentType: "text/plain" });
  console.log("Encoded:", enc.stats.oligoCount, "oligos");
  
  const primerLen = cfg.primerLength;
  const totalInnerBytes = enc.encoded.metadata.innerRS.n + 2;
  const addressBytes = 4;
  
  // Manual per-oligo decode
  let successCount = 0;
  let failCount = 0;
  
  for (let i = 0; i < enc.encoded.oligos.length; i++) {
    const oligo = enc.encoded.oligos[i];
    const innerDna = oligo.sequence.slice(primerLen, oligo.sequence.length - primerLen);
    
    // Step 1: Address extraction
    const addressDna = innerDna.slice(0, 16);
    const addrBytes = homopolymerSafeDnaToAddress(addressDna);
    const unwhitened = unwhitenAddress(addrBytes);
    const index = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
    const seed = unwhitened[3];
    
    // Step 2: Full constrained decode
    try {
      const result = splitConstrainedDnaToBytesWithErasure(innerDna, 3, addressBytes, totalInnerBytes);
      const innerBlock = result.data;
      const erasures = result.erasures.filter(e => e).length;
      
      // Verify CRC
      const crcOk = verifyCrc16(innerBlock);
      
      // Check address in decoded block
      const decodedAddr = innerBlock.slice(0, 4);
      const addrMatch = addrBytes.every((b, j) => b === decodedAddr[j]);
      
      if (i < 3) {
        console.log(`Oligo ${i}: index=${index} seed=${seed} crcOk=${crcOk} addrMatch=${addrMatch} erasures=${erasures} blockLen=${innerBlock.length}`);
      }
      
      if (crcOk) successCount++;
      else failCount++;
    } catch (e: any) {
      console.log(`Oligo ${i}: DECODE ERROR: ${e.message}`);
      failCount++;
    }
  }
  
  console.log(`\nPer-oligo decode: ${successCount} OK, ${failCount} FAIL`);
  
  // Now try full pipeline
  const sim = simulate(enc.encoded.oligos, { ...PRESET_CLEAN, coverage: 30, simulator: "basic" });
  const dec = await decodeReads(sim.reads, enc.encoded.metadata, cfg, enc.encoded.forwardPrimer, enc.encoded.reversePrimer);
  console.log(`Full pipeline: dataLen=${dec.data?.length} hashOk=${dec.hashMatches}`);
}

main().catch(e => { console.error(e); process.exit(1); });
