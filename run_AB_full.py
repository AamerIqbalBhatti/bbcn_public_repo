import sys, time, json
from collections import Counter
sys.path.insert(0,'setup_a'); sys.path.insert(0,'setup_b/code')
import pandas as pd
from bbcn import honest_bbcn as HB
import forward_stab_kernel_design as FSK
from bbcn_switch import FAST, MID, SLOW, rules as sw_rules, label as sw_label
NODES=HB.NODES; lags=HB.make_lags(); APOP=FSK.APOP_TARGET
SW_STATE=['AKT1','PTEN','MTOR','PDPK1','PIK3CA','MDM2','TP53','FOXO3']
COH={'TCGA':'setup_a/data/binarized/tcga_brca_1082x135.csv',
     'METABRIC':'setup_a/data/binarized/metabric_1980x135.csv',
     'ISPY2':'setup_a/data/binarized/ispy2_988x135.csv'}
def seed_switch(b):
    U=HB.patient_U(b); x0={nd:int(b.get(nd,0)) for nd in SW_STATE}; x0['PHLPP']=int(x0['FOXO3'] and not U)
    I=dict(SRC=int(b.get('SRC',0)),RHEB=int(b.get('RHEB',0)),IGF1R=int(b.get('IGF1R',0)),
           RTK_up=int(b.get('GRB2',0) or b.get('IRS1',0)),CDKN2A=int(b.get('CDKN2A',0)),
           E2F1=int(b.get('E2F1',0)),ATM=int(b.get('ATM',0)),ATR=int(b.get('ATR',0)))
    return x0,I,U
def sw_course(x0,I,U,clamp,hold,total=200):
    s=dict(x0)
    for nd,v in clamp.items(): s[nd]=v
    tail=[]
    for t in range(1,total+1):
        snap=dict(s); nf=sw_rules(snap,I,U)
        for nd in FAST: s[nd]=nf[nd]
        if t%5==0: s['MDM2']=nf['MDM2']
        if t%25==0:
            for nd in SLOW: s[nd]=nf[nd]
        if t<=hold:
            for nd,v in clamp.items(): s[nd]=v
        tail.append(sw_label(s))
    return tail
def B_static(x0,I,U,c): return sw_course(x0,I,U,c,200)[-1]=='APOPTOTIC'
def B_dyn(x0,I,U,c):
    t=sw_course(x0,I,U,c,120)[-40:]; return sum(l=='APOPTOTIC' for l in t)/len(t)>0.5
def apop(s): return s['CASP3']==1 and s['AKT1']==0
def fullnet(b,U,clamp,T1,T2,theta=10,cmax=20):
    bus={n:int(b.get(n,0)) for n in NODES}
    for k,v in clamp.items(): bus[k]=v
    c=0;commit=0;hist=[bus]
    for t in range(T1+T2):
        nb=HB.step_honest(HB._view(hist,t,lags),U)
        c=min(c+1,cmax) if HB.commit_signal(nb) else max(c-1,0); commit=commit or (c>=theta)
        if commit:
            for bk in HB.BRAKES: nb[bk]=0
        if t<T1:
            for k,v in clamp.items(): nb[k]=v
        hist.append(nb)
    return hist
def A_dyn(b,U,c):
    h=fullnet(b,U,c,120,180); return (sum(apop(s) for s in h[-40:])/40)>0.5

name=sys.argv[1]; start=int(sys.argv[2]); end=int(sys.argv[3])
df=pd.read_csv(COH[name],index_col=0); n=len(df)
end = n if end<0 else min(end,n)
t0=time.time(); res=found=0; Bs=Bd=Asingle=Acombo=0; genes=Counter()
for i in range(start,end):
    b={nd:int(df.iloc[i].get(nd,0)) for nd in NODES}
    x0,I,U=seed_switch(b)
    result,kernels=FSK.design_patient_kernel(x0,I,U,method='stabilize')
    if not kernels and result=='APOPTOTIC': continue
    res+=1
    if not kernels: continue
    found+=1
    clamp={nd:APOP[nd] for bx in kernels for nd in kernels[bx]}
    for nd in clamp: genes[nd]+=1
    Bs+=B_static(x0,I,U,clamp); Bd+=B_dyn(x0,I,U,clamp)
    Asingle+=A_dyn(b,U,clamp)
    Acombo +=A_dyn(b,U,{**clamp,'ATM':1,'ATR':1})   # designed drug + p53 engagement, then withdraw
out=dict(name=name,start=start,end=end,n=end-start,res=res,found=found,
         Bs=Bs,Bd=Bd,Asingle=Asingle,Acombo=Acombo,genes=dict(genes),secs=round(time.time()-t0))
json.dump(out,open(f'/home/claude/ab_{name}_{start}.json','w'))
print(out)
