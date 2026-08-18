/**
 * Check erasure positions from constrained decode
 */
import { encodeFile } from "../src/lib/dna/codec";
import { V51_DEFAULT_CONFIG } from "../src/lib/dna/presets";
import { splitConstrainedDnaToBytesWithErasure } from "../src/lib/dna/constrained-mapping";

async function main() {
  const data = new Uint8Array(Buffer.from("Hello World!"));
  const cfg = { ...V51_DEFAULT_CONFIG, maxRetries: 10 };
  const enc = await encodeFile(data, cfg, { fileName: "test", contentType: "text/plain" });
  
  const primerLen = cfg.primerLength;
  const totalInnerBytes = enc.encoded.metadata.innerRS.n + 2;
  
  for (let i = 0; i < enc.encoded.oligos.length; i++) {
    const oligo = enc.encoded.oligos[i];
    const innerDna = oligo.sequence.slice(primerLen, oligo.sequence.length - primerLen);
    
    const result = splitConstrainedDnaToBytesWithErasure(innerDna, 3, 4, totalInnerBytes);
    const erasures = result.erasures;
    
    // Count erasures in address vs payload region
    const addrErasureBits = erasures.slice(0, 32).filter(e => e).length; // 4 bytes * 8 bits
    const payloadErasureBits = erasures.slice(32).filter(e => e).length;
    
    console.log(`Oligo ${i}: total=${erasures.filter(e=>e).length} addr_bits=${addrErasureBits} payload_bits=${payloadErasureBits} dataLen=${result.data.length}`);
  }
}

main().catch(e => { console.error(e); process.exit(1); });
