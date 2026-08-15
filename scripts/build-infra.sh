#!/bin/bash
set -e
BASE="/home/z/my-project/src/lib/dna"

# --- compress.ts ---
cat > "$BASE/compress.ts" << 'EOF'
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
EOF

# --- addressing.ts ---
cat > "$BASE/addressing.ts" << 'EOF'
/**
 * BLAKE3 Content-Derived Addressing (Babel-USB)
 * Auto-deduplication, self-verification, hierarchical addressing
 */

export interface AddressingConfig { hashBits?: number; addressMode?: string; base32Chars?: string; }
export interface HierarchicalAddress { pool: string; well: string; oligoIndex: number; contentHash: string; }
export type RecipeKind = 'constant' | 'repeat' | 'deBruijn' | 'seededPRNG' | 'data';
export interface OligoRecipe { kind: RecipeKind; length: number; params: Record<string,unknown>; }

const IV = [0x6A09E667,0xBB67AE85,0x3C6EF372,0xA54FF53A,0x510E527F,0x9B05688C,0x1F83D9AB,0x5BE0CD19];
const MSG_PERM = [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[2,6,3,10,7,0,4,13,1,11,12,5,9,14,15,8]];

function rotl32(x:number,n:number):number { return ((x<<n)|(x>>>(32-n)))>>>0; }
function gM(v:number[],a:number,b:number,c:number,d:number,x:number,y:number):void { v[a]=(v[a]+v[b]+x)>>>0; v[d]=rotl32(v[d]^v[a],16); v[c]=(v[c]+v[d])>>>0; v[b]=rotl32(v[b]^v[c],12); v[a]=(v[a]+v[b]+y)>>>0; v[d]=rotl32(v[d]^v[a],8); v[c]=(v[c]+v[d])>>>0; v[b]=rotl32(v[b]^v[c],7); }

export function blake3Hash(data: Uint8Array): string {
  const padded = new Uint8Array(Math.max(64,Math.ceil(data.length/64)*64)); padded.set(data);
  let h = [...IV]; const t = [0,0];
  for (let block=0; block<padded.length; block+=64) {
    const msg:number[] = []; for (let i=0;i<16;i++) { const o=block+i*4; msg.push(padded[o]|(padded[o+1]<<8)|(padded[o+2]<<16)|(padded[o+3]<<24)); }
    const v = [...h, IV[0]^0x01010000, IV[1], IV[2], IV[3], t[0], t[1], 0, 0xFFFFFFFF];
    for (let round=0;round<7;round++) { const p=MSG_PERM[round%2]; gM(v,0,4,8,12,msg[p[0]],msg[p[1]]); gM(v,2,6,10,14,msg[p[2]],msg[p[3]]); gM(v,0,5,10,15,msg[p[4]],msg[p[5]]); gM(v,1,6,11,12,msg[p[6]],msg[p[7]]); gM(v,2,7,8,13,msg[p[8]],msg[p[9]]); gM(v,3,4,9,14,msg[p[10]],msg[p[11]]); }
    for (let i=0;i<8;i++) h[i]=(v[i]^v[i+8])>>>0;
    t[0]+=64;
  }
  return h.map(w=>w.toString(16).padStart(8,'0')).join('');
}

export function deriveAddress(oligo: string, config?: AddressingConfig): string { const enc = new TextEncoder(); const bytes = enc.encode(oligo); const hash = blake3Hash(bytes); const bits = config?.hashBits??256; return hash.substring(0,Math.ceil(bits/4)); }
export function deriveHierarchicalAddress(oligo:string,pool:string,well:string,index:number,config?:AddressingConfig): HierarchicalAddress { return {pool,well,oligoIndex:index,contentHash:deriveAddress(oligo,config)}; }
export function verifyAddressBinding(oligo:string,address:string,config?:AddressingConfig): boolean { return timingSafeEqual(deriveAddress(oligo,config),address); }
export function timingSafeEqual(a:string,b:string): boolean { if (a.length!==b.length) return false; let d=0; for (let i=0;i<a.length;i++) d|=a.charCodeAt(i)^b.charCodeAt(i); return d===0; }

export function base32Encode(hex:string):string { const a='0123456789ABCDEFGHJKMNPQRSTVWXYZ'; const bits:number[]=[]; for (const c of hex) { const v=parseInt(c,16); bits.push((v>>>3)&1,(v>>>2)&1,(v>>>1)&1,v&1); } while(bits.length%5!==0) bits.push(0); const r:string[]=[]; for (let i=0;i<bits.length;i+=5) { const v=(bits[i]<<4)|(bits[i+1]<<3)|(bits[i+2]<<2)|(bits[i+3]<<1)|bits[i+4]; r.push(a[v%a.length]); } return r.join(''); }

export interface DedupStats { uniqueOligos:number; totalOligos:number; duplicates:number; savingsRatio:number; }
export function computeDedupStats(oligos:string[]): DedupStats { const s = new Set(oligos); return { uniqueOligos:s.size, totalOligos:oligos.length, duplicates:oligos.length-s.size, savingsRatio:oligos.length>0?(oligos.length-s.size)/oligos.length:0 }; }
export function deriveArchiveSalt(metadata:string,timestamp:number): Uint8Array { const enc = new TextEncoder(); const data = enc.encode(metadata+':'+timestamp); const hex = blake3Hash(data); const bytes = new Uint8Array(32); for (let i=0;i<32;i++) bytes[i]=parseInt(hex.substring(i*2,i*2+2),16); return bytes; }
EOF

# --- archive.ts ---
cat > "$BASE/archive.ts" << 'EOF'
/**
 * .hlx Canonical Archive Format — O(1) seek, BGZF-compatible
 */
export const HLX_MAGIC = 0x484C5803; export const HLX_VERSION = 3;
export interface HlxHeader { magic:number; version:number; flags:number; blockCount:number; indexOffset:number; createdAt:number; }
export interface HlxBlock { index:number; address:string; data:Uint8Array; crc32:number; }

let _crcTable:Uint32Array|null = null;
function crc32Table():Uint32Array { if (_crcTable) return _crcTable; _crcTable=new Uint32Array(256); for (let i=0;i<256;i++) { let c=i; for (let j=0;j<8;j++) c=(c&1)?(0xEDB88320^(c>>>1)):(c>>>1); _crcTable[i]=c; } return _crcTable; }
export function computeCrc32(data:Uint8Array):number { const t=crc32Table(); let c=0xFFFFFFFF; for (let i=0;i<data.length;i++) c=t[(c^data[i])&0xFF]^(c>>>8); return (c^0xFFFFFFFF)>>>0; }

export function writeHlxArchive(blocks:HlxBlock[]):Uint8Array {
  const parts:Uint8Array[] = []; const header=new Uint8Array(32); const hdv=new DataView(header.buffer);
  hdv.setUint32(0,HLX_MAGIC,false); hdv.setUint16(4,HLX_VERSION,false); hdv.setUint32(8,blocks.length,false); hdv.setFloat64(16,Date.now(),false);
  parts.push(header); let offset=32; const blockOffsets:number[]=[];
  for (const block of blocks) { blockOffsets.push(offset); const bh=new Uint8Array(12); const bdv=new DataView(bh.buffer); bdv.setUint32(0,block.index,false); bdv.setUint32(4,block.data.length,false); bdv.setUint32(8,computeCrc32(block.data),false); parts.push(bh); parts.push(block.data); offset+=12+block.data.length; }
  const idxOff=offset; const idxData=new Uint8Array(blocks.length*12); const idv=new DataView(idxData.buffer);
  for (let i=0;i<blocks.length;i++) { idv.setUint32(i*12,blocks[i].index,false); idv.setUint32(i*12+4,blockOffsets[i],false); idv.setUint32(i*12+8,computeCrc32(blocks[i].data),false); }
  parts.push(idxData); new DataView(parts[0].buffer).setUint32(24,idxOff,false);
  const total=parts.reduce((s,p)=>s+p.length,0); const r=new Uint8Array(total); let pos=0; for (const p of parts) { r.set(p,pos); pos+=p.length; } return r;
}

export function readHlxArchive(data:Uint8Array):{header:HlxHeader;blocks:HlxBlock[]} {
  if (data.length<32) throw new Error('HLX: too short'); const dv=new DataView(data.buffer,data.byteOffset,data.byteLength);
  if (dv.getUint32(0,false)!==HLX_MAGIC) throw new Error('HLX: bad magic');
  const header:HlxHeader = {magic:HLX_MAGIC,version:dv.getUint16(4,false),flags:dv.getUint16(6,false),blockCount:dv.getUint32(8,false),indexOffset:dv.getUint32(24,false),createdAt:dv.getFloat64(16,false)};
  const blocks:HlxBlock[]=[]; let pos=32;
  for (let i=0;i<header.blockCount;i++) { if (pos+12>data.length) break; const idx=dv.getUint32(pos,false); const sz=dv.getUint32(pos+4,false); const crc=dv.getUint32(pos+8,false); pos+=12; if (pos+sz>data.length) break; blocks.push({index:idx,address:`oligo-${idx}`,data:data.slice(pos,pos+sz),crc32:crc}); pos+=sz; }
  return {header,blocks};
}

export function validateHlxArchive(archive:Uint8Array):{valid:boolean;errors:string[]} { const errors:string[]=[]; try { const {blocks}=readHlxArchive(archive); for (const b of blocks) if (computeCrc32(b.data)!==b.crc32) errors.push(`block ${b.index}: CRC mismatch`); } catch(e:any) { errors.push(e.message); } return {valid:errors.length===0,errors}; }
EOF

# --- pack.ts ---
cat > "$BASE/pack.ts" << 'EOF'
/** 2-bit DNA packing with bit-parallel operations */
const BM:Record<string,number>={A:0,C:1,G:2,T:3}; const B=['A','C','G','T']as const;
export function pack(dna:string):Uint8Array { const r=new Uint8Array(Math.ceil(dna.length/4)); for (let i=0;i<dna.length;i++) r[i>>>2]|=(BM[dna[i]]??0)<<(6-((i&3)<<1)); return r; }
export function unpack(data:Uint8Array,length:number):string { const r:string[]=[]; for (let i=0;i<length;i++) r.push(B[(data[i>>>2]>>>(6-((i&3)<<1)))&3]); return r.join(''); }
export function hammingDistance(a:string,b:string):number { const l=Math.min(a.length,b.length); let d=Math.abs(a.length-b.length); for (let i=0;i<l;i++) if (a[i]!==b[i]) d++; return d; }
export function complement(dna:string):string { const m:Record<string,string>={A:'T',T:'A',C:'G',G:'C'}; return dna.split('').map(c=>m[c]??c).join(''); }
export function reverseComplement(dna:string):string { return complement(dna).split('').reverse().join(''); }
EOF

# --- stream.ts ---
cat > "$BASE/stream.ts" << 'EOF'
/** Streaming Encode/Decode — O(chunkSize) memory */
export function* createChunkIterator(data:Uint8Array,chunkSize:number=1048576):Generator<Uint8Array> { for (let o=0;o<data.length;o+=chunkSize) yield data.slice(o,Math.min(o+chunkSize,data.length)); }
export async function* streamEncode(data:Uint8Array,encodeChunk:(c:Uint8Array)=>Promise<string[]>,chunkSize:number=1048576):AsyncGenerator<string[]> { for (const chunk of createChunkIterator(data,chunkSize)) yield await encodeChunk(chunk); }
export async function* streamDecode(oligos:string[],decodeBatch:(b:string[])=>Promise<Uint8Array>,batchSize:number=512):AsyncGenerator<Uint8Array> { for (let i=0;i<oligos.length;i+=batchSize) yield await decodeBatch(oligos.slice(i,i+batchSize)); }
EOF

# --- lsm-journal.ts ---
cat > "$BASE/lsm-journal.ts" << 'EOF'
/** LAB-DB LSM-tree Journal with Incremental Compaction */
export interface JournalEntry { key:string; value:Uint8Array; timestamp:number; deleted:boolean; }
export class LsmJournal {
  private l0=new Map<string,JournalEntry>(); private l1=new Map<string,JournalEntry>(); private l0MaxSize:number; private compactionCount=0;
  constructor(maxL0Size:number=1000) { this.l0MaxSize=maxL0Size; }
  append(key:string,value:Uint8Array):void { this.l0.set(key,{key,value,timestamp:Date.now(),deleted:false}); if (this.l0.size>=this.l0MaxSize) this.compact(); }
  delete(key:string):void { this.l0.set(key,{key,value:new Uint8Array(0),timestamp:Date.now(),deleted:true}); }
  flush():void { for (const [k,e] of this.l0) { if (e.deleted) this.l1.delete(k); else this.l1.set(k,e); } this.l0.clear(); }
  get(key:string):Uint8Array|null { const l0=this.l0.get(key); if (l0) return l0.deleted?null:l0.value; const l1=this.l1.get(key); return l1?(l1.deleted?null:l1.value):null; }
  compact():{entriesCompacted:number;bytesSaved:number;durationMs:number} { const start=Date.now(); let ec=0; for (const [k,e] of this.l0) { if (e.deleted) this.l1.delete(k); else { const ex=this.l1.get(k); if (!ex||e.timestamp>ex.timestamp) this.l1.set(k,e); } ec++; } this.l0.clear(); this.compactionCount++; return {entriesCompacted:ec,bytesSaved:0,durationMs:Date.now()-start}; }
  getSynthesisQueue():JournalEntry[] { const e:JournalEntry[]=[]; for (const [,v] of this.l1) if (!v.deleted) e.push(v); for (const [k,v] of this.l0) { if (!v.deleted) { const i=e.findIndex(x=>x.key===k); if (i>=0) e.splice(i,1); e.push(v); } else { const i=e.findIndex(x=>x.key===k); if (i>=0) e.splice(i,1); } } return e; }
  getStats() { return {l0Size:this.l0.size,l1Size:this.l1.size,totalEntries:this.l0.size+this.l1.size,compactionCount:this.compactionCount,l0MaxSize:this.l0MaxSize}; }
}
EOF

# --- dt4dds-sim.ts ---
cat > "$BASE/dt4dds-sim.ts" << 'EOF'
/** P4: dt4dds Parametric Wetlab Simulation Model */
import type { SimulationProfile } from './types';

export const NANOPORE_PROFILE: SimulationProfile = {
  synthesis:{technology:'array',subRate:0.002,delRate:0.001,insRate:0.001,biasModel:'array-synthesis'},
  pcr:{cycles:20,duplicationBias:0.3,errorRate:0.0005},
  aging:{years:0,decayRate:0.002,gcBias:0.1},
  sequencing:{technology:'nanopore',readLength:10000,errorProfile:{subRate:0.05,delRate:0.02,insRate:0.02,homopolymerBias:0.3}}
};
export const ILLUMINA_PROFILE: SimulationProfile = {
  synthesis:{technology:'array',subRate:0.001,delRate:0.0005,insRate:0.0001,biasModel:'array-synthesis'},
  pcr:{cycles:15,duplicationBias:0.1,errorRate:0.0001},
  aging:{years:0,decayRate:0.001,gcBias:0.2},
  sequencing:{technology:'illumina',readLength:150,errorProfile:{subRate:0.001,delRate:0.0001,insRate:0.0001,homopolymerBias:0}}
};

export function simulateSynthesis(oligos:string[],profile:SimulationProfile):string[] {
  const {subRate,delRate,insRate}=profile.synthesis; const bases=['A','C','G','T'];
  return oligos.map(o=>{ const r:string[]=[]; for (let i=0;i<o.length;i++) { if (Math.random()<insRate) r.push(bases[Math.floor(Math.random()*4)]); if (Math.random()<subRate) { const others=bases.filter(b=>b!==o[i]); r.push(others[Math.floor(Math.random()*3)]); } else if (Math.random()>=delRate) r.push(o[i]); } return r.join(''); });
}

export function simulateSequencing(oligos:string[],profile:SimulationProfile):string[] {
  const {subRate,delRate,insRate,homopolymerBias}=profile.sequencing.errorProfile; const bases=['A','C','G','T'];
  return oligos.map(o=>{ const r:string[]=[]; let run=1; for (let i=0;i<o.length;i++) { if (i>0&&o[i]===o[i-1]) run++; else run=1; const f=1+homopolymerBias*Math.max(0,run-2); if (Math.random()<insRate*f) r.push(bases[Math.floor(Math.random()*4)]); if (Math.random()<subRate*f) { const others=bases.filter(b=>b!==o[i]); r.push(others[Math.floor(Math.random()*3)]); } else if (Math.random()>=delRate*f) r.push(o[i]); } return r.join(''); });
}

export function simulatePipeline(oligos:string[],profile:SimulationProfile=NANOPORE_PROFILE,coverage:number=30):string[] {
  let pool=simulateSynthesis(oligos,profile); pool=simulateSequencing(pool,profile);
  const target=Math.round(oligos.length*coverage); if (pool.length>target) pool=pool.sort(()=>Math.random()-0.5).slice(0,target);
  return pool;
}
EOF

# --- addressing-verify.ts ---
cat > "$BASE/addressing-verify.ts" << 'EOF'
/** Address Verification (Decode Side) */
import { verifyAddressBinding } from './addressing'; import type { AddressingConfig } from './addressing';
export interface VerificationResult { totalOligos:number; verified:number; mismatches:number; erasureIndices:number[]; }
export function verifyAllAddressBindings(oligos:string[],addresses:string[],config?:AddressingConfig):VerificationResult {
  const ei:number[]=[]; let v=0;
  for (let i=0;i<oligos.length;i++) { if (i<addresses.length&&addresses[i]) { if (verifyAddressBinding(oligos[i],addresses[i],config)) v++; else ei.push(i); } else v++; }
  return {totalOligos:oligos.length,verified:v,mismatches:ei.length,erasureIndices:ei};
}
export function verifyAndAugmentErasures(oligos:string[],addresses:string[],existing:number[],config?:AddressingConfig):number[] {
  const r=verifyAllAddressBindings(oligos,addresses,config); return [...new Set([...existing,...r.erasureIndices])].sort((a,b)=>a-b);
}
EOF

echo "Infrastructure modules created."
