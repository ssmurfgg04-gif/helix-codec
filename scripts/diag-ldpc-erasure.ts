/**
 * Test LDPC decode with erasures
 */
import { encodeFile } from "../src/lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { splitConstrainedDnaToBytesWithErasure } from "../src/lib/dna/constrained-mapping";
import { LDPCInnerCode } from "../src/lib/dna/ldpc-codec";
import { verifyCrc16, crc16Bytes } from "../src/lib/dna/crc16";
import { unwhitenAddress, xorWithSeed } from "../src/lib/dna/mapping";

async function main() {
  const data = new Uint8Array(Buffer.from("Hello World!"));
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const enc = await encodeFile(data, cfg, { fileName: "test", contentType: "text/plain" });
  
  const primerLen = cfg.primerLength;
  const innerN = enc.encoded.metadata.innerRS.n;
  const innerK = enc.encoded.metadata.innerRS.k;
  const totalInnerBytes = innerN + 2;
  const addressBytes = 4;
  
  console.log("innerRS:", JSON.stringify(enc.encoded.metadata.innerRS));
  console.log("innerN:", innerN, "innerK:", innerK, "totalInner:", totalInnerBytes);
  
  // Create LDPC decoder
  const ldpc = new LDPCInnerCode({ n: innerN, k: innerK });
  
  for (let i = 0; i < enc.encoded.oligos.length; i++) {
    const oligo = enc.encoded.oligos[i];
    const innerDna = oligo.sequence.slice(primerLen, oligo.sequence.length - primerLen);
    
    // Constrained decode
    const result = splitConstrainedDnaToBytesWithErasure(innerDna, 3, 4, totalInnerBytes);
    const innerBlock = result.data;
    const erasures = result.erasures;
    
    // Extract the LDPC codeword (innerN bytes, before CRC)
    const codeword = innerBlock.slice(0, innerN);
    const crc = innerBlock.slice(innerN, innerN + 2);
    
    // Check CRC before LDPC
    const crcOk = verifyCrc16(innerBlock);
    console.log(`\nOligo ${i}: CRC before LDPC: ${crcOk ? 'OK' : 'FAIL'}`);
    
    // Try LDPC decode with erasures
    const erasePositions: number[] = [];
    for (let j = 0; j < erasures.length && j < innerN * 8; j++) {
      if (erasures[j]) erasePositions.push(j);
    }
    console.log(`  Erasure positions (${erasePositions.length}):`, erasePositions.slice(0, 10).join(','), '...');
    
    try {
      const decoded = ldpc.decode(codeword, erasePositions);
      // Rebuild inner block with decoded codeword
      const newBlock = new Uint8Array(totalInnerBytes);
      newBlock.set(decoded, 0);
      newBlock.set(crc16Bytes(decoded), innerN);
      
      const crcAfterLdpc = verifyCrc16(newBlock);
      console.log(`  LDPC decode OK, CRC after: ${crcAfterLdpc ? 'OK' : 'FAIL'}`);
      
      // Check if decoded data matches original
      const addr = decoded.slice(0, 4);
      const unwhitened = unwhitenAddress(addr);
      const index = (unwhitened[0] << 16) | (unwhitened[1] << 8) | unwhitened[2];
      const seed = unwhitened[3];
      console.log(`  Index: ${index} (expected: ${oligo.index}), Seed: ${seed}`);
    } catch (e: any) {
      console.log(`  LDPC decode FAILED: ${e.message}`);
    }
  }
}

main().catch(e => { console.error(e); process.exit(1); });
