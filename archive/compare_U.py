"""compare_U.py COHORT [START END]  -- E40 durability under BOTH U definitions, full N (or slice)."""
import sys, time
sys.path.insert(0, 'setup_a'); sys.path.insert(0, 'setup_b/code')
import pandas as pd
from bbcn import honest_bbcn as HB
import forward_stab_kernel_design as FSK
from bbcn_switch import FAST, MID, SLOW, rules as sw_rules, label as sw_label
NODES = HB.NODES; lags = HB.make_lags(); APOP = FSK.APOP_TARGET
SW_STATE = ['AKT1','PTEN','MTOR','PDPK1','PIK3CA','MDM2','TP53','FOXO3']
AKT_ACT = ['AKT1','MTOR','RPS6KB1','EIF4EBP1']; FOXO_NUC = ['FOXO3','FOXO1']
COH = {'TCGA':'setup_a/data/binarized/tcga_brca_1082x135.csv',
       'METABRIC':'setup_a/data/binarized/metabric_1980x135.csv',
       'ISPY2':'setup_a/data/binarized/ispy2_988x135.csv'}

def U_any(b):  return 1 if (any(int(b.get(k,0)) for k in AKT_ACT) and any(int(b.get(k,0)) for k in FOXO_NUC)) else 0
def U_mean(b):
    a = sum(int(b.get(k,0)) for k in AKT_ACT)/len(AKT_ACT); f = sum(int(b.get(k,0)) for k in FOXO_NUC)/len(FOXO_NUC)
    return 1 if (a>0.5 and f>0.5) else 0

def seedU(b, U):
    x0 = {nd:int(b.get(nd,0)) for nd in SW_STATE}; x0['PHLPP'] = int(x0['FOXO3'] and not U)
    I = dict(SRC=int(b.get('SRC',0)),RHEB=int(b.get('RHEB',0)),IGF1R=int(b.get('IGF1R',0)),
             RTK_up=int(b.get('GRB2',0) or b.get('IRS1',0)),CDKN2A=int(b.get('CDKN2A',0)),
             E2F1=int(b.get('E2F1',0)),ATM=int(b.get('ATM',0)),ATR=int(b.get('ATR',0)))
    return x0, I

def A_dyn(b, U, clamp, T1=120, T2=180, theta=10, cmax=20):
    bus = {n:int(b.get(n,0)) for n in NODES}
    for k,v in clamp.items(): bus[k]=v
    c=0; commit=0; hist=[bus]
    for t in range(T1+T2):
        nb = HB.step_honest(HB._view(hist,t,lags), U)
        c = min(c+1,cmax) if HB.commit_signal(nb) else max(c-1,0); commit = commit or (c>=theta)
        if commit:
            for bk in HB.BRAKES: nb[bk]=0
        if t<T1:
            for k,v in clamp.items(): nb[k]=v
        hist.append(nb)
    return sum(s['CASP3']==1 for s in hist[-40:])/40 > 0.5   # commit durability

def run(b, Ufn):
    U = Ufn(b); x0, I = seedU(b, U)
    r, kernels = FSK.design_patient_kernel(x0, I, U, method='stabilize')
    if not kernels and r == 'APOPTOTIC': return (0, 0)          # not resistant
    if not kernels: return (1, 0)                               # resistant, no kernel
    clamp = {nd:APOP[nd] for bx in kernels for nd in kernels[bx]}
    return (1, int(A_dyn(b, U, clamp)))                        # resistant, commit?

name = sys.argv[1]; df = pd.read_csv(COH[name], index_col=0); n = len(df)
s = int(sys.argv[2]) if len(sys.argv) > 2 else 0
e = int(sys.argv[3]) if len(sys.argv) > 3 else n
t0 = time.time()
res = {'any':[0,0,0], 'mean':[0,0,0]}   # [n, resistant, commit]
for i in range(s, e):
    b = {nd:int(df.iloc[i].get(nd,0)) for nd in NODES}
    for tag, Ufn in (('any',U_any), ('mean',U_mean)):
        rz, cm = run(b, Ufn); res[tag][0]+=1; res[tag][1]+=rz; res[tag][2]+=cm
for tag in ('any','mean'):
    N, rz, cm = res[tag]
    print(f"{name:9} {tag:5} | N={N:4} | resistant {rz:4} ({round(100*rz/N)}%) | "
          f"commit {cm:4} ({round(100*cm/rz) if rz else 0}% of resistant)")
print(f"[{name} {time.time()-t0:.0f}s]")
