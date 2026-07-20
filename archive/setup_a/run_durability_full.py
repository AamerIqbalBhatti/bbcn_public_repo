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
def seed(b):
    U=HB.patient_U(b); x0={nd:int(b.get(nd,0)) for nd in SW_STATE}; x0['PHLPP']=int(x0['FOXO3'] and not U)
    I=dict(SRC=int(b.get('SRC',0)),RHEB=int(b.get('RHEB',0)),IGF1R=int(b.get('IGF1R',0)),
           RTK_up=int(b.get('GRB2',0) or b.get('IRS1',0)),CDKN2A=int(b.get('CDKN2A',0)),
           E2F1=int(b.get('E2F1',0)),ATM=int(b.get('ATM',0)),ATR=int(b.get('ATR',0)))
    return x0,I,U
def sw_dyn(x0,I,U,clamp,hold=120,total=200):
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
    return sum(l=='APOPTOTIC' for l in tail[-40:])/40>0.5
def A_dyn(b,U,clamp,T1=120,T2=180,theta=10,cmax=20):
    bus={n:int(b.get(n,0)) for n in NODES}
    for k,v in clamp.items(): bus[k]=v
    hist=[bus]
    for t in range(T1+T2):
        nb=HB.step_honest(HB._view(hist,t,lags),U,cascade=True)   # Step 16: real cascade, no brake-zeroing
        if t<T1:
            for k,v in clamp.items(): nb[k]=v
        hist.append(nb)
    tail=hist[-40:]
    strict=sum(s['CASP3']==1 and s['AKT1']==0 for s in tail)/40>0.5
    commitd=sum(s['CASP3']==1 for s in tail)/40>0.5
    return strict, commitd

t0=time.time(); res={}
for name,path in COH.items():
    df=pd.read_csv(path,index_col=0); n=len(df)
    found=resist=strictc=commitc=swc=0
    for i in range(n):
        b={nd:int(df.iloc[i].get(nd,0)) for nd in NODES}
        x0,I,U=seed(b)
        r,kernels=FSK.design_patient_kernel(x0,I,U,method='stabilize')
        if not kernels and r=='APOPTOTIC': continue
        resist+=1
        if not kernels: continue
        found+=1
        clamp={nd:APOP[nd] for bx in kernels for nd in kernels[bx]}
        st,cm=A_dyn(b,U,clamp); strictc+=st; commitc+=cm
        swc+=sw_dyn(x0,I,U,clamp)
    res[name]=dict(n=n,resist=resist,found=found,strict=strictc,commit=commitc,switch=swc)
    p=lambda x,d: f"{round(100*x/d)}%" if d else "-"
    print(f"{name:9} N={n:4} | resistant {p(resist,n):>4} | kernel {p(found,resist):>4} | "
          f"SWITCH durable {p(swc,found):>4} | FULLNET strict {p(strictc,found):>4} | FULLNET commit {p(commitc,found):>4}")
json.dump(res,open('/home/claude/durability_full.json','w'))
print(f"[{time.time()-t0:.0f}s]")
