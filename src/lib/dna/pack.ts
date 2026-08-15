/** 2-bit DNA packing with bit-parallel operations */
const BM:Record<string,number>={A:0,C:1,G:2,T:3}; const B=['A','C','G','T']as const;
export function pack(dna:string):Uint8Array { const r=new Uint8Array(Math.ceil(dna.length/4)); for (let i=0;i<dna.length;i++) r[i>>>2]|=(BM[dna[i]]??0)<<(6-((i&3)<<1)); return r; }
export function unpack(data:Uint8Array,length:number):string { const r:string[]=[]; for (let i=0;i<length;i++) r.push(B[(data[i>>>2]>>>(6-((i&3)<<1)))&3]); return r.join(''); }
export function hammingDistance(a:string,b:string):number { const l=Math.min(a.length,b.length); let d=Math.abs(a.length-b.length); for (let i=0;i<l;i++) if (a[i]!==b[i]) d++; return d; }
export function complement(dna:string):string { const m:Record<string,string>={A:'T',T:'A',C:'G',G:'C'}; return dna.split('').map(c=>m[c]??c).join(''); }
export function reverseComplement(dna:string):string { return complement(dna).split('').reverse().join(''); }
