import sys, json, time, os
sys.path.insert(0,'.'); sys.path.insert(0,'../setup_b/code')
import pandas as pd
from bbcn import honest_bbcn as HB
import forward_stab_kernel_design as FSK
from bbcn_switch import FAST, MID, SLOW, rules as sw_rules, label as sw_label
NODES=HB.NODES; lags=HB.make_lags(); APOP=FSK.APOP_TARGET
SW_STATE=['AKT1','PTEN','MTOR','PDPK1','PIK3CA','MDM2','TP53','FOXO3']
COH={'TCGA':'data/binarized/tcga_brca_1082x135.csv',
     'METABRIC':'data/binarized/metabric_1980x135.csv',
     'ISPY2':'data/binarized/ispy2_988x135.csv'}
OUT='/home/claude/fulln_e40.json'
def seed(b):
    U=HB.patient_clamp_off(b); x0={nd:int(b.get(nd,0)) for nd in SW_STATE}; x0['PHLPP']=int(x0['FOXO3'] and not U)
    I=dict(SRC=int(b.get('SRC',0)),RHEB=int(b.get('RHEB',0)),IGF1R=int(b.get('IGF1R',0)),
           RTK_up=int(b.get('GRB2',0) or b.get('IRS1',0)),CDKN2A=int(b.get('CDKN2A',0)),
           E2F1=int(b.get('E2F1',0)),ATM=int(b.get('ATM',0)),ATR=int(b.get('ATR',0)))
    return x0,I,U
def sw_dyn(x0,I,U,clamp,hold=120,total=200):
    s=dict(x0)
    for nd,v in clamp.items(): s[nd]=v
    tail=[]
    for t in range(1,total+1):
        nf=sw_rules(dict(s),I,U)
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
    c=0;commit=0;hist=[bus]
    for t in range(T1+T2):
        nb=HB.step_honest(HB._view(hist,t,lags),U)
        c=min(c+1,cmax) if HB.commit_signal(nb) else max(c-1,0); commit=commit or (c>=theta)
        if commit:
            for bk in HB.BRAKES: nb[bk]=0
        if t<T1:
            for k,v in clamp.items(): nb[k]=v
        hist.append(nb)
    tail=hist[-40:]
    strict=sum(s['CASP3']==1 and s['AKT1']==0 for s in tail)/40>0.5
    commitd=sum(s['CASP3']==1 for s in tail)/40>0.5
    return strict, commitd
cohort=sys.argv[1]; a=int(sys.argv[2]); bb=int(sys.argv[3])
df=pd.read_csv(COH[cohort],index_col=0); bb=min(bb,len(df))
acc=json.load(open(OUT)) if os.path.exists(OUT) else {}
d=acc.get(cohort, dict(n=0,resist=0,found=0,strict=0,commit=0,switch=0,done=0))
t0=time.time()
for i in range(a,bb):
    b={nd:int(df.iloc[i].get(nd,0)) for nd in NODES}
    x0,I,U=seed(b)
    r,kernels=FSK.design_patient_kernel(x0,I,U,method='stabilize')
    d['n']+=1
    if not kernels and r=='APOPTOTIC': continue
    d['resist']+=1
    if not kernels: continue
    d['found']+=1
    clamp={nd:APOP[nd] for bx in kernels for nd in kernels[bx]}
    st,cm=A_dyn(b,U,clamp); d['strict']+=st; d['commit']+=cm
    d['switch']+=sw_dyn(x0,I,U,clamp)
d['done']=max(d['done'],bb)
acc[cohort]=d; json.dump(acc,open(OUT,'w'))
dt=time.time()-t0
print(f"{cohort} [{a}:{bb}] {bb-a} pts in {dt:.0f}s ({dt/max(1,bb-a)*1000:.0f} ms/pt) | n={d['n']} resist={d['resist']} found={d['found']} done_to={d['done']}")
