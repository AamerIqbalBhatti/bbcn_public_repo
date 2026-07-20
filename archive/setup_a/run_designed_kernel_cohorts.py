"""
run_designed_kernel_cohorts.py
Design-small / test-big, per patient, on the REPAIRED network.
  DESIGN: seed the 9-node switch from each patient, run the REAL forward-stabilization
          kernel designer (forward_stab_kernel_design.design_patient_kernel) restricted
          to druggable nodes -> the per-patient drug nomination (the CDE).
  TEST  : apply that designed kernel on the full 136-node repaired delayed net
          (honest_bbcn, death-engaged accumulator), STATIC (hold) and DYNAMIC (withdraw).
This threads all three obstacles: repaired rules (both sides), PHLPP (both sides),
accumulator (full-net test phase).
"""
import sys, time, random
from collections import Counter
sys.path.insert(0,'setup_a'); sys.path.insert(0,'setup_b/code')
import pandas as pd
from bbcn import honest_bbcn as HB
import forward_stab_kernel_design as FSK

NODES=HB.NODES; lags=HB.make_lags(); APOP=FSK.APOP_TARGET
SW_STATE=['AKT1','PTEN','MTOR','PDPK1','PIK3CA','MDM2','TP53','FOXO3']  # PHLPP derived
COH={'TCGA':'setup_a/data/binarized/tcga_brca_1082x135.csv',
     'METABRIC':'setup_a/data/binarized/metabric_1980x135.csv',
     'ISPY2':'setup_a/data/binarized/ispy2_988x135.csv'}

def seed_switch(b):
    U=HB.patient_U(b)
    x0={nd:int(b.get(nd,0)) for nd in SW_STATE}
    x0['PHLPP']=int(x0['FOXO3'] and not U)             # repaired PHLPP init
    I=dict(SRC=int(b.get('SRC',0)),RHEB=int(b.get('RHEB',0)),IGF1R=int(b.get('IGF1R',0)),
           RTK_up=int(b.get('GRB2',0) or b.get('IRS1',0)),
           CDKN2A=int(b.get('CDKN2A',0)),E2F1=int(b.get('E2F1',0)),
           ATM=int(b.get('ATM',0)),ATR=int(b.get('ATR',0)))
    return x0,I,U

def apop(s): return s['CASP3']==1 and s['AKT1']==0
def fullnet(b,U,clamp,T1,T2,theta=10,cmax=20):
    bus={n:int(b.get(n,0)) for n in NODES}
    for k,v in clamp.items(): bus[k]=v
    c=0;commit=0;hist=[bus]
    for t in range(T1+T2):
        active=(clamp if t<T1 else {})
        nb=HB.step_honest(HB._view(hist,t,lags),U)
        c=min(c+1,cmax) if HB.commit_signal(nb) else max(c-1,0); commit=commit or (c>=theta)
        if commit:
            for bk in HB.BRAKES: nb[bk]=0
        for k,v in active.items(): nb[k]=v
        hist.append(nb)
    return hist
def static_reached(b,U,clamp): return apop(fullnet(b,U,clamp,260,0)[-1])
def dynamic_durable(b,U,clamp):
    h=fullnet(b,U,clamp,120,180); return (sum(apop(s) for s in h[-40:])/40)>0.5

random.seed(11); N=200; t0=time.time()
print("REPAIRED net | designed druggable kernel per patient (forward-stab CDE) | static+dynamic test")
print(f"{'cohort':9}| N | already | resistant | kernel | designed-clamp full-net  | nominated genes")
print(f"{'':9}|   |  apop   |  (need k) | found  |  static / dynamic        |")
print("-"*104)
agg={}
for name,path in COH.items():
    df=pd.read_csv(path,index_col=0); idx=random.sample(range(len(df)),N)
    already=res=found=0; stat=dyn=0; genes=Counter()
    for i in idx:
        b={n:int(df.iloc[i].get(n,0)) for n in NODES}
        x0,I,U=seed_switch(b)
        result,kernels=FSK.design_patient_kernel(x0,I,U,method='stabilize')
        if not kernels and result=='APOPTOTIC':
            already+=1; continue
        res+=1
        if not kernels: continue          # resistant but no druggable kernel found
        found+=1
        clamp={nd:APOP[nd] for bx in kernels for nd in kernels[bx]}
        for nd in clamp: genes[nd]+=1
        stat+=static_reached(b,U,clamp); dyn+=dynamic_durable(b,U,clamp)
    p=lambda x,d:f"{round(100*x/d)}%" if d else "  -"
    nm=", ".join(f"{g}:{genes[g]}" for g,_ in genes.most_common())
    print(f"{name:9}|{N:3}| {p(already,N):>6} | {p(res,N):>8} | {p(found,res):>5} | "
          f"{p(stat,found):>6} / {p(dyn,found):<8}      | {nm}")
    agg[name]=dict(already=already,res=res,found=found,stat=stat,dyn=dyn,genes=dict(genes),N=N)
print(f"\n[{time.time()-t0:.0f}s, N={N}/cohort]  static=hold the designed clamp; dynamic=withdraw then check durable")
import json; open('/home/claude/designed_kernel_agg.json','w').write(json.dumps(agg))
