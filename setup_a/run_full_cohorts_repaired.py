"""
run_full_cohorts_repaired.py — FULL N, all three cohorts, on the repaired delayed
network with the death-engaged accumulator and U derived per patient. No sampling.
"""
import sys, time; sys.path.insert(0,'.')
import pandas as pd
from bbcn import honest_bbcn as HB
NODES=HB.NODES; lags=HB.make_lags()
COH={'TCGA':'data/binarized/tcga_brca_1082x135.csv',
     'METABRIC':'data/binarized/metabric_1980x135.csv',
     'ISPY2':'data/binarized/ispy2_988x135.csv'}
def apop(s): return s['CASP3']==1 and s['AKT1']==0

def course(b,U,hold,T1,T2=180,theta=10,cmax=20):
    bus={n:int(b.get(n,0)) for n in NODES}
    for k,v in (hold or {}).items(): bus[k]=v
    c=0;commit=0;history=[bus]
    for t in range(T1+T2):
        active=(hold if t<T1 else {})
        nb=HB.step_honest(HB._view(history,t,lags),U)
        c=min(c+1,cmax) if HB.commit_signal(nb) else max(c-1,0); commit=commit or (1 if c>=theta else 0)
        if commit:
            for bk in HB.BRAKES: nb[bk]=0
        for k,v in (active or {}).items(): nb[k]=v
        history.append(nb)
    return (sum(apop(s) for s in history[-40:])/40)>0.5

print(f"{'cohort':9} |   N  | U=1  | drug flips | genotoxic flips | drug-resist==U")
print("-"*78)
T0=time.time()
for name,path in COH.items():
    df=pd.read_csv(path, index_col=0); N=len(df)
    u=drug=gen=match=0
    for i in range(N):
        b={n:int(df.iloc[i].get(n,0)) for n in NODES}
        U=HB.patient_U(b); u+=U
        df_flip=course(b,U,{'AKT1':0},120); drug+=df_flip
        gen+=course(b,U,{'ATM':1,'ATR':1},120)
        if (U==1 and not df_flip) or (U==0 and df_flip): match+=1
    pc=lambda x:f"{round(100*x/N)}%"
    print(f"{name:9} | {N:4} | {pc(u):>4} | {pc(drug):>9}  | {pc(gen):>13}   | {match}/{N} ({pc(match)})")
print(f"\n[full run, {time.time()-T0:.0f}s]  network: repaired delayed BBCN, death-engaged accumulator, U per patient")
