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
