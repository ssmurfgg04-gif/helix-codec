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
