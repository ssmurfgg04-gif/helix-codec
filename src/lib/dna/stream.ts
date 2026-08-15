/** Streaming Encode/Decode — O(chunkSize) memory */
export function* createChunkIterator(data:Uint8Array,chunkSize:number=1048576):Generator<Uint8Array> { for (let o=0;o<data.length;o+=chunkSize) yield data.slice(o,Math.min(o+chunkSize,data.length)); }
export async function* streamEncode(data:Uint8Array,encodeChunk:(c:Uint8Array)=>Promise<string[]>,chunkSize:number=1048576):AsyncGenerator<string[]> { for (const chunk of createChunkIterator(data,chunkSize)) yield await encodeChunk(chunk); }
export async function* streamDecode(oligos:string[],decodeBatch:(b:string[])=>Promise<Uint8Array>,batchSize:number=512):AsyncGenerator<Uint8Array> { for (let i=0;i<oligos.length;i+=batchSize) yield await decodeBatch(oligos.slice(i,i+batchSize)); }
