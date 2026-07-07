import sys, json, time, os
sys.path.insert(0,'.'); sys.path.insert(0,'../setup_b/code')
import pandas as pd
from bbcn import honest_bbcn as HB
import durability_clock as DC
import forward_stab_kernel_design as FSK
NODES=HB.NODES; APOP=FSK.APOP_TARGET
SW=['AKT1','PTEN','MTOR','PDPK1','PIK3CA','MDM2','TP53','FOXO3']
COH={'TCGA':'data/binarized/tcga_brca_1082x135.csv',
     'METABRIC':'data/binarized/metabric_1980x135.csv',
     'ISPY2':'data/binarized/ispy2_988x135.csv'}
OUT='/home/claude/fulln_frontload.json'
def seed(b):
    U=HB.patient_clamp_off(b); x0={nd:int(b.get(nd,0)) for nd in SW}; x0['PHLPP']=int(x0['FOXO3'] and not U)
    I=dict(SRC=int(b.get('SRC',0)),RHEB=int(b.get('RHEB',0)),IGF1R=int(b.get('IGF1R',0)),
           RTK_up=int(b.get('GRB2',0) or b.get('IRS1',0)),CDKN2A=int(b.get('CDKN2A',0)),
           E2F1=int(b.get('E2F1',0)),ATM=int(b.get('ATM',0)),ATR=int(b.get('ATR',0)))
    return x0,I,U
cohort=sys.argv[1]; a=int(sys.argv[2]); bb=int(sys.argv[3])
df=pd.read_csv(COH[cohort],index_col=0); bb=min(bb,len(df))
acc=json.load(open(OUT)) if os.path.exists(OUT) else {}
d=acc.get(cohort, dict(found=0,esc=0,front=0,done=0))
t0=time.time()
for i in range(a,bb):
    b={nd:int(df.iloc[i].get(nd,0)) for nd in NODES}
    x0,I,U=seed(b)
    r,kernels=FSK.design_patient_kernel(x0,I,U,method='stabilize')
    if not kernels: continue
    d['found']+=1
    d['esc']  +=DC.run_capped_on_clock(b,U,kernels,APOP,pulsed=True,front_loaded=False)['durable']
    d['front']+=DC.run_capped_on_clock(b,U,kernels,APOP,pulsed=True,front_loaded=True)['durable']
d['done']=max(d['done'],bb)
acc[cohort]=d; json.dump(acc,open(OUT,'w'))
print(f"{cohort} [{a}:{bb}] in {time.time()-t0:.0f}s | found={d['found']} done_to={d['done']}")
