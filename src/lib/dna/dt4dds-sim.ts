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
