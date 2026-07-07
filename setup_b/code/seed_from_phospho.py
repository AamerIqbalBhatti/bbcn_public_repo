#!/usr/bin/env python3
"""
seed_from_phospho.py  --  PHOSPHO-PROTEIN (RPPA) activity front-end for the BBCN switch
========================================================================================
WHAT THIS DOES  (this is the "Bee" run, session 2026-06-11)
----------------------------------------------------------
Seeds the switch's signaling nodes from the ACTUAL phospho-protein activity state
(I-SPY2 RPPA data, GSE196093) -- the real signal we had been proxying with mRNA all along.
Every signaling node has a direct phospho readout (see PHOSPHO->NODE mapping below).

WHY IT MATTERS
--------------
This was the definitive test of the project's recurring "expression != signaling" question.
RESULT: phospho seeding did NOT predict pCR any better than mRNA. Both fail equally.
=> the bottleneck was NEVER signal type; it is ENDPOINT DISTANCE (pCR is a distal whole-tumor
   outcome -- immune/stroma/PK -- beyond the single-cell decision the switch models).
The DRUG-OVERLAP concordance (kernels nominate PI3K/AKT/mTOR, overlapping the MK-2206 arm)
remains robust here too. Ledger ref: section 80.

One intriguing PHOSPHO-ONLY lead (hypothesis, n tiny): apoptotic-primed patients respond to
MK-2206 (41% pCR) far better than to chemo-only control (12%) -- a +29pp treatment interaction.
Only visible with phospho seeding. n=29 vs 91, NOT a result -- needs power + a Fisher test.

PHOSPHO -> SWITCH NODE MAPPING (RPPA endpoint -> node):
  AKT1   <- AKT.S473, AKT.T308
  MTOR   <- mTOR.S2448, p70S6K.T389, S6RP.S240.S244, X4EBP1.S65   (activity + downstream output)
  FOXO3  <- FOXO1.S256, FOXO3a.S253, FOXO1.T24.FOXO3a.T32   (HIGH phospho => FOXO INACTIVE)
  GSK3B  <- GSK3aB.S21.S9         (HIGH => GSK3 inhibited => AKT active)
  TSC2   <- Tuberin.TSC2.Y1571    (HIGH => TSC2 inactivated => mTOR de-repressed)
  PTEN   <- PTEN.S380, PTEN.total
  TP53   <- p53.S15, p53.total ;  MDM2 <- MDM2.S166
  DAMAGE <- ATM.S1981, ATR.S428, CHK1.S345, CHK2.S33.S35
  apoptosis execution <- Caspase.3/7/9.cleaved, PARP.cleaved, BAD.S136
  U (uncoupled self-latch, ledger 70) = AKT-active AND FOXO-inactive(phospho)

INPUT: GSE196093_series_matrix.gz  (RPPA component; carries phospho values + arm + pcr)
"""
import numpy as np, sys, gzip
sys.path.insert(0, ".")
from bbcn_switch import NODES, simulate, label
exec(open("forward_stab_kernel_design.py").read().split("# run small first")[0])


def parse_rppa(path):
    """Parse the GSE196093 RPPA series matrix: phospho values + per-sample arm/pcr."""
    chars = {}; data = {}; intable = False
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!Sample_characteristics"):
                vals = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                key = vals[0].split(":")[0].strip()
                chars[key] = [v.split(":",1)[1].strip() if ":" in v else "" for v in vals]
            if line.startswith("!series_matrix_table_begin"): intable = True; continue
            if line.startswith("!series_matrix_table_end"):   intable = False; continue
            if intable:
                p = line.rstrip("\n").split("\t")
                name = p[0].strip().strip('"')
                if name != "ID_REF":
                    data[name] = np.array([float(x) if x not in ("","NA","null") else np.nan for x in p[1:]])
    return chars, data


def run(rppa_path):
    chars, data = parse_rppa(rppa_path)
    n = len(next(iter(data.values())))
    arm = chars.get("arm", [""]*n); pcr = chars.get("pcr", [""]*n)
    print(f"RPPA: {n} samples, {len(data)} phospho endpoints")

    def rz(v):
        m = np.nanmedian(v); mad = np.nanmedian(np.abs(v-m))
        sd = 1.4826*mad if mad>0 else (np.nanstd(v) or 1.0)
        z = (v-m)/sd; z[np.isnan(z)] = 0; return z
    def get(*names):
        arrs = [rz(data[k]) for k in names if k in data]
        return np.mean(arrs, axis=0) if arrs else np.zeros(n)

    akt_act   = get('AKT.S473','AKT.T308')
    mtor_act  = get('mTOR.S2448','p70S6K.T389','S6RP.S240.S244','X4EBP1.S65')
    foxo_phos = get('FOXO1.S256','FOXO1.T24.FOXO3a.T32','FOXO3a.S253')   # HIGH = FOXO inactive
    pten_act  = get('PTEN.S380','PTEN.total')
    p53_act   = get('p53.S15','p53.total')
    damage    = get('ATM.S1981','ATR.S428','CHK1.S345','CHK2.S33.S35')
    prolif    = get('Cyclin.B1.total','Cyclin.A2.total','RB.S780')

    ON = lambda x: (x > np.median(x)).astype(int)
    akt_on=ON(akt_act); mtor_on=ON(mtor_act); foxo_inact=ON(foxo_phos)
    pten_on=ON(pten_act); p53_on=ON(p53_act); dmg_on=ON(damage); prolif_on=ON(prolif)
    U = (akt_on & foxo_inact)

    def route(i):
        x0 = {nd: 0 for nd in NODES}
        x0['AKT1']=int(akt_on[i]); x0['MTOR']=int(mtor_on[i]); x0['PDPK1']=int(akt_on[i]); x0['PIK3CA']=int(akt_on[i])
        x0['PTEN']=int(pten_on[i]); x0['FOXO3']=int(1-foxo_inact[i]); x0['PHLPP']=int((1-foxo_inact[i]) and not U[i])
        x0['TP53']=int(p53_on[i])
        I = dict(SRC=0, RHEB=int(mtor_on[i]), IGF1R=int(akt_on[i]), RTK_up=int(akt_on[i]),
                 CDKN2A=int(prolif_on[i]), E2F1=int(prolif_on[i]), ATM=int(dmg_on[i]), ATR=int(dmg_on[i]))
        base = label(simulate(x0, I, int(U[i])))
        if base == 'APOPTOTIC':
            return 'apoptotic', {}
        attr, kernels = design_patient_kernel(x0, I, int(U[i]))
        return ('flippable' if attr=='APOPTOTIC' else 'resistant'), kernels

    from collections import Counter, defaultdict
    valid = [i for i in range(n) if pcr[i] in ('0','1')]
    mk    = [i for i in valid if 'MK-2206' in arm[i] and 'Trastuzumab' not in arm[i]]
    ctrl  = [i for i in valid if arm[i].strip()=='Paclitaxel']
    dist = Counter(route(i)[0] for i in valid)
    print("routing:", {k: f"{dist[k]} ({100*dist[k]/len(valid):.0f}%)" for k in ['apoptotic','flippable','resistant']})
    def split3(ids, lab):
        g = defaultdict(lambda: [0,0])
        for i in ids:
            grp,_ = route(i); g[grp][0]+=int(pcr[i]); g[grp][1]+=1
        print(f"[{lab}]")
        for k in ['apoptotic','flippable','resistant']:
            pc,nn=g[k]; print(f"  {k:10s} n={nn:3d} pCR={pc} ({100*pc/max(nn,1):.0f}%)")
    split3(mk, "MK-2206 arm");  split3(ctrl, "Control arm")
    kg = Counter()
    for i in valid:
        for v in route(i)[1].values():
            for g in ('PIK3CA','AKT1','PDPK1','MTOR','PTEN'):
                if g in v: kg[g]+=1
    print("kernel drug nominations:", dict(kg.most_common()))


if __name__ == "__main__":
    run("GSE196093_series_matrix.gz")   # EDIT path as needed
