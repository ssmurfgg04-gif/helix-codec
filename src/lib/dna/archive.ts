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
