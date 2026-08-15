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
