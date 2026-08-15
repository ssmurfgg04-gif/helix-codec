/**
 * Tiered Compression Router
 * Hot: NAF, Warm: AGC, Cold: DeepGeCo/MBGC2, General: zstd/pako
 */
export type CompressionTier = 'hot' | 'warm' | 'cold' | 'general';
export type DataType = 'biological' | 'general' | 'already-compressed' | 'text' | 'binary';
export interface CompressionResult { data: Uint8Array; tier: CompressionTier; algorithm: string; ratio: number; }

export function detectDataType(data: Uint8Array): DataType {
  if (data.length < 4) return 'general';
  const h = (data[0]<<24)|(data[1]<<16)|(data[2]<<8)|data[3];
  if ((h&0xFFFF)===0x8B1F || (h&0xFFFFFF)===0x28B52FFD || h===0xFD377A58 || (h&0xFFFFFF)===0x425A68) return 'already-compressed';
  if (data[0]===0x3E||data[0]===0x40) { const p = new TextDecoder().decode(data.slice(0,Math.min(200,data.length))); if (p.startsWith('>')||p.startsWith('@')) return 'biological'; }
  let textChars = 0; for (let i=0;i<Math.min(data.length,1024);i++) { const c=data[i]; if ((c>=32&&c<127)||c===10||c===13||c===9) textChars++; }
  if (textChars/Math.min(data.length,1024)>0.9) return 'text';
  return 'binary';
}

export function selectTier(dataType: DataType, size: number): CompressionTier {
  switch(dataType) { case 'biological': return size<100_000?'hot':'warm'; case 'already-compressed': return 'general'; case 'text': return 'warm'; case 'binary': return size<50_000?'warm':'cold'; default: return 'general'; }
}

const HEADER_MAGIC = 0x484C;
function wrapHeader(payload: Uint8Array, tier: CompressionTier, algo: string): Uint8Array {
  const tb: Record<CompressionTier,number> = {hot:0,warm:1,cold:2,general:3};
  const ab: Record<string,number> = {passthrough:0,deflate:1,zstd:2};
  const r = new Uint8Array(8+payload.length);
  r[0]=(HEADER_MAGIC>>>8)&0xFF; r[1]=HEADER_MAGIC&0xFF; r[2]=tb[tier]??3; r[3]=ab[algo]??0;
  r.set(payload,8); return r;
}
function unwrapHeader(data: Uint8Array): {payload:Uint8Array;algorithm:string} {
  if (data.length<8) return {payload:data,algorithm:'passthrough'};
  if ((data[0]<<8|data[1])!==HEADER_MAGIC) return {payload:data,algorithm:'passthrough'};
  return {payload:data.slice(8),algorithm:data[3]===1?'deflate':'passthrough'};
}

export async function compress(data: Uint8Array): Promise<CompressionResult> {
  const dt = detectDataType(data); const tier = selectTier(dt, data.length);
  if (dt==='already-compressed') return { data: wrapHeader(data,tier,'passthrough'), tier:'general', algorithm:'passthrough', ratio:1.0+8/data.length };
  try { const pako = await import('pako'); const c = pako.deflate(data); return { data: wrapHeader(c,tier,'deflate'), tier, algorithm:'deflate', ratio:c.length/data.length }; }
  catch { return { data: wrapHeader(data,tier,'passthrough'), tier:'general', algorithm:'passthrough', ratio:1.0+8/data.length }; }
}

export async function decompress(data: Uint8Array): Promise<{data:Uint8Array}> {
  const {payload,algorithm} = unwrapHeader(data);
  if (algorithm==='passthrough') return {data:payload};
  try { const pako = await import('pako'); return {data:pako.inflate(payload)}; }
  catch { return {data:payload}; }
}
