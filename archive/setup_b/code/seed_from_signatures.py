#!/usr/bin/env python3
"""
seed_from_signatures.py  --  I-SPY2-grade multi-gene SIGNATURE front-end for the BBCN switch
=============================================================================================
WHAT THIS DOES
--------------
Instead of binarizing single genes at their median (the old, lossy front-end), this seeds the
switch's signaling nodes from MULTI-GENE, PATHWAY-ANCHORED activity signatures -- the same
class of signal-extraction the I-SPY2 team uses for their validated response predictors.
We compute a continuous robust-z score per signature, then THRESHOLD it (stay Boolean).

Ledger ref: this is the "level playing field with I-SPY2" front-end (session 2026-06-11).
The drug-overlap concordance (our kernels nominate PI3K/AKT/mTOR, overlapping I-SPY2's
MK-2206 arm) is ROBUST to this front-end. pCR prediction is NOT rescued by it (see ledger 79).

INPUTS (edit paths at bottom):
  - GSE194040 I-SPY2 expression matrix (.txt.gz, genes x samples, sample IDs = ResIDs)
  - ISPY2_response_table.csv  (patient_id, arm, pcr, hr, her2, mp)  [recovered from GEO series matrix]
  - bbcn_switch.py + forward_stab_kernel_design.py on the path
"""
import numpy as np, csv, sys, gzip
sys.path.insert(0, ".")  # expects bbcn_switch.py + forward_stab_kernel_design.py alongside
from bbcn_switch import NODES, simulate, label

# pull the loaders + design fns from the kernel-design module (everything before the demo marker)
exec(open("forward_stab_kernel_design.py").read().split("# run small first")[0])

# ----------------------------------------------------------------------------
# 1. SIGNATURE DEFINITIONS  (multi-gene, pathway-anchored; literature gene-sets)
#    For a paper, swap these for Wolf et al. 2022 Table S1 gene members verbatim.
# ----------------------------------------------------------------------------
SIGNATURES = {
  'AKT_MTOR_ACT': ['AKT1','AKT2','AKT3','MTOR','RPS6KB1','RPS6','EIF4EBP1','RPTOR','PDPK1'],
  'PTEN_loss'   : ['PTEN','INPP4B'],                 # low score => PTEN-pathway loss
  'FOXO_NUC'    : ['FOXO1','FOXO3','FOXO4'],
  'PROLIF'      : ['MKI67','CCNB1','CCNE1','PCNA','MCM2','AURKA','BUB1','CDC20','UBE2C','TOP2A','CCNB2','CDK1','FOXM1'],
  'DAMAGE'      : ['ATM','ATR','CHEK1','CHEK2','BRCA1','BRCA2','RAD51','MDC1','H2AFX','MRE11','RBBP8','PARP1'],
  'RTK_up'      : ['EGFR','ERBB2','ERBB3','IGF1R','GRB2','GRB7','INSR','IRS1'],
  'P53_ACT'     : ['TP53','CDKN1A','MDM2','BAX','BBC3','PMAIP1'],
}

def robust_z(v):
    """median/MAD robust z-score (the project default, ledger 73)."""
    m = np.nanmedian(v); mad = np.nanmedian(np.abs(v - m))
    sd = 1.4826 * mad if mad > 0 else (np.nanstd(v) or 1.0)
    z = (v - m) / sd; z[np.isnan(z)] = 0
    return z

def load_full_matrix(path, need_genes):
    """Load only the rows we need from the big gene x sample matrix."""
    expr = {}
    with gzip.open(path, "rt") as f:
        samples = f.readline().rstrip("\n").split("\t")[1:]
        for line in f:
            p = line.rstrip("\n").split("\t")
            g = p[0].strip().strip('"')
            if g in need_genes:
                expr[g] = np.array([float(x) if x not in ("", "NA") else np.nan for x in p[1:]])
    return samples, expr


def run(expr_path, response_csv):
    need = set(NODES)
    for v in SIGNATURES.values():
        need |= set(v)
    samples, expr = load_full_matrix(expr_path, need)
    n = len(samples); idx = {s: i for i, s in enumerate(samples)}
    print(f"loaded {n} samples, {len(expr)} signature/node genes")

    # 2. continuous signature scores, then THRESHOLD at cohort median (stay Boolean)
    def sig(name):
        g = [x for x in SIGNATURES[name] if x in expr]
        return np.mean([robust_z(expr[x]) for x in g], axis=0) if g else np.zeros(n)
    S = {k: sig(k) for k in SIGNATURES}
    ON = lambda x: (x > np.median(x)).astype(int)
    akt_on   = ON(S['AKT_MTOR_ACT'])
    pten_off = (S['PTEN_loss'] < np.median(S['PTEN_loss'])).astype(int)
    foxo_on  = ON(S['FOXO_NUC'])
    prolif_on= ON(S['PROLIF'])
    damage_on= ON(S['DAMAGE'])
    rtk_on   = ON(S['RTK_up'])
    p53_on   = ON(S['P53_ACT'])
    U_sig    = (akt_on & foxo_on)   # uncoupled self-latch condition (ledger 70)

    # 3. seed switch nodes from signatures, route each patient, design kernels
    def route(i):
        x0 = {nd: 0 for nd in NODES}
        x0['AKT1']=int(akt_on[i]); x0['MTOR']=int(akt_on[i]); x0['PDPK1']=int(akt_on[i]); x0['PIK3CA']=int(akt_on[i])
        x0['PTEN']=int(1-pten_off[i]); x0['FOXO3']=int(foxo_on[i]); x0['PHLPP']=int(foxo_on[i] and not U_sig[i])
        x0['TP53']=int(p53_on[i])
        I = dict(SRC=0, RHEB=int(akt_on[i]), IGF1R=int(rtk_on[i]), RTK_up=int(rtk_on[i]),
                 CDKN2A=int(prolif_on[i]), E2F1=int(prolif_on[i]), ATM=int(damage_on[i]), ATR=int(damage_on[i]))
        base = label(simulate(x0, I, int(U_sig[i])))
        if base == 'APOPTOTIC':
            return 'apoptotic', {}
        attr, kernels = design_patient_kernel(x0, I, int(U_sig[i]))
        return ('flippable' if attr == 'APOPTOTIC' else 'resistant'), kernels

    # 4. report: routing + 3-way pCR on MK-2206 arm + kernel drug nominations
    resp = {}
    with open(response_csv) as f:
        for r in csv.DictReader(f):
            resp[r['patient_id']] = r
    mk  = [p for p in resp if 'MK-2206' in resp[p]['arm'] and 'Trastuzumab' not in resp[p]['arm']
           and p in idx and resp[p]['pcr'] in ('0','1')]
    from collections import Counter, defaultdict
    allids = [p for p in resp if p in idx and resp[p]['pcr'] in ('0','1')]
    dist = Counter(route(idx[p])[0] for p in allids)
    print("routing:", {k: f"{dist[k]} ({100*dist[k]/len(allids):.0f}%)" for k in ['apoptotic','flippable','resistant']})
    g = defaultdict(lambda: [0,0])
    for p in mk:
        grp, _ = route(idx[p]); g[grp][0]+=int(resp[p]['pcr']); g[grp][1]+=1
    print("MK-2206 arm 3-way pCR:")
    for k in ['apoptotic','flippable','resistant']:
        pc,nn = g[k]; print(f"  {k:10s} n={nn:3d} pCR={pc} ({100*pc/max(nn,1):.0f}%)")
    kg = Counter()
    for p in allids:
        for v in route(idx[p])[1].values():
            for gene in ('PIK3CA','AKT1','PDPK1','MTOR','PTEN'):
                if gene in v: kg[gene]+=1
    print("kernel drug nominations:", dict(kg.most_common()))


if __name__ == "__main__":
    # EDIT THESE PATHS to wherever you keep the data:
    EXPR = "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988_txt.gz"
    RESP = "ISPY2_response_table.csv"
    run(EXPR, RESP)
