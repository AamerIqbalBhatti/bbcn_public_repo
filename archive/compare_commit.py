"""Phase 5 step 16 -- DECISIVE COMPARISON (measurement only, no production edits).

For every resistant+kerneled patient, compute durability two ways:
  FLAG    : the abstract COMMIT flag (commit_signal/theta accumulator -> zero BRAKES)
  CASCADE : the real p53->PUMA/NOXA->SMAC->executioner cascade (step_honest cascade=True),
            with NO brake-zeroing.
Both for the idealized continuous hold (E40) and the clinical pulsed clock (E42).

Usage:  python compare_commit.py <TCGA|METABRIC|ISPY2> [start] [limit]
"""
import sys
sys.path.insert(0, 'setup_a'); sys.path.insert(0, 'setup_b/code')
import pandas as pd
from bbcn import honest_bbcn as HB
from bbcn import clinical_clock as CK
import forward_stab_kernel_design as FSK
import durability_clock as DC

NODES = HB.NODES; lags = HB.make_lags(); APOP = FSK.APOP_TARGET
SW_STATE = ['AKT1','PTEN','MTOR','PDPK1','PIK3CA','MDM2','TP53','FOXO3']
COH = {'TCGA':'setup_a/data/binarized/tcga_brca_1082x135.csv',
       'METABRIC':'setup_a/data/binarized/metabric_1980x135.csv',
       'ISPY2':'setup_a/data/binarized/ispy2_988x135.csv'}

def seed(b):
    U = HB.patient_clamp_off(b); x0 = {nd: int(b.get(nd,0)) for nd in SW_STATE}
    x0['PHLPP'] = int(x0['FOXO3'] and not U)
    I = dict(SRC=int(b.get('SRC',0)), RHEB=int(b.get('RHEB',0)), IGF1R=int(b.get('IGF1R',0)),
             RTK_up=int(b.get('GRB2',0) or b.get('IRS1',0)), CDKN2A=int(b.get('CDKN2A',0)),
             E2F1=int(b.get('E2F1',0)), ATM=int(b.get('ATM',0)), ATR=int(b.get('ATR',0)))
    return x0, I, U

# --- idealized continuous hold (mirror regression A_dyn exactly) ---
def ideal_flag(b, U, clamp, T1=120, T2=180, theta=10, cmax=20):
    bus = {n: int(b.get(n,0)) for n in NODES}
    for k, v in clamp.items(): bus[k] = v
    c = 0; commit = 0; hist = [bus]
    for t in range(T1+T2):
        nb = HB.step_honest(HB._view(hist, t, lags), U)
        c = min(c+1, cmax) if HB.commit_signal(nb) else max(c-1, 0); commit = commit or (c >= theta)
        if commit:
            for bk in HB.BRAKES: nb[bk] = 0
        if t < T1:
            for k, v in clamp.items(): nb[k] = v
        hist.append(nb)
    return sum(s['CASP3'] == 1 for s in hist[-40:]) / 40 > 0.5

def ideal_cascade(b, U, clamp, T1=120, T2=180):
    bus = {n: int(b.get(n,0)) for n in NODES}
    for k, v in clamp.items(): bus[k] = v
    hist = [bus]
    for t in range(T1+T2):
        nb = HB.step_honest(HB._view(hist, t, lags), U, cascade=True)   # real cascade, NO brake-zeroing
        if t < T1:
            for k, v in clamp.items(): nb[k] = v
        hist.append(nb)
    return sum(s['CASP3'] == 1 for s in hist[-40:]) / 40 > 0.5

# --- clinical pulsed clock ---
def clock_flag(b, U, clamp):
    return DC.A_dyn_clock(b, U, clamp, pulsed=True)        # current flag path

def clock_cascade(b, U, clamp):
    bus = {n: int(b.get(n,0)) for n in NODES}
    for k, v in clamp.items(): bus[k] = v
    hist = [bus]
    for t in range(CK.total_days()):
        day = t + 1; dosed = CK.is_on_day(day)
        nb = HB.step_honest(HB._view(hist, t, lags), U, cascade=True)
        if dosed:
            for k, v in clamp.items(): nb[k] = v
        hist.append(nb)
    tail = hist[-CK.observation_days():]
    return sum(s['CASP3'] == 1 for s in tail) / len(tail) > 0.5

def main():
    name = sys.argv[1]
    start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
    df = pd.read_csv(COH[name], index_col=0)
    end = min(start + limit, len(df))
    found = idf = idc = ckf = ckc = 0
    for i in range(start, end):
        b = {nd: int(df.iloc[i].get(nd, 0)) for nd in NODES}
        x0, I, U = seed(b)
        r, kernels = FSK.design_patient_kernel(x0, I, U, method='stabilize')
        if not kernels: continue
        found += 1
        clamp = {nd: APOP[nd] for bx in kernels for nd in kernels[bx]}
        idf += ideal_flag(b, U, clamp)
        idc += ideal_cascade(b, U, clamp)
        ckf += clock_flag(b, U, clamp)
        ckc += clock_cascade(b, U, clamp)
    print(f"{name} rows[{start}:{end}] kerneled={found}")
    print(f"  IDEAL   flag={idf} ({100*idf/found:.0f}%)  cascade={idc} ({100*idc/found:.0f}%)")
    print(f"  CLOCK   flag={ckf} ({100*ckf/found:.0f}%)  cascade={ckc} ({100*ckc/found:.0f}%)")
    print(f"RAW {name} {start} {end} {found} {idf} {idc} {ckf} {ckc}")

if __name__ == '__main__':
    main()
