"""
cohort_pipeline.py — The full safeguarded patient pipeline (ledger sec.66-70).

Stage 0  expression -> Boolean initial state (multi-gene activity seeds)
Stage 1  four redundant phenotype sensors (proliferation / apoptosis / resistance / invasion)
Stage 2  three-box multirate switch -> attractor (APOPTOTIC / SURVIVAL / mixed)
         with safeguards: damage-ACTIVITY gating, mutation correction, uncoupled-state flag U
Stage 3  minimal druggable kernel (single-input reachability)

Reproduces the cohort table in the manuscript (TCGA-BRCA, METABRIC, I-SPY2).
Expression / mutation files are NOT bundled (large, controlled-access); see README
for download instructions and point the loaders at your local copies.
"""

from __future__ import annotations
import gzip
import time
import numpy as np
from collections import Counter

from bbcn_switch import FAST, MID, SLOW, NODES, rules, simulate, label

# ---------------------------------------------------------------------------
# Gene sets (multi-gene activity signatures and sensors)
# ---------------------------------------------------------------------------
SEED_SIG = {
    "AKT1":  ["AKT1", "MTOR", "RPS6KB1", "RPS6", "EIF4EBP1", "PDPK1", "GSK3B"],
    "TP53":  ["CDKN1A", "GADD45A", "BBC3", "MDM2", "SESN1", "RRM2B"],
    "MTOR":  ["MTOR", "RPS6KB1", "EIF4EBP1"], "PDPK1": ["PDPK1", "AKT1"],
    "PIK3CA": ["PIK3CA", "AKT1", "PDPK1"], "PTEN": ["PTEN"], "MDM2": ["MDM2"],
    "PHLPP": ["PHLPP1", "PHLPP2", "FOXO3"], "FOXO3": ["FOXO3", "FOXO1"],
}
AKT_ACT = ["AKT1", "MTOR", "RPS6KB1", "RPS6", "EIF4EBP1"]   # for U flag
FOXO_NUC = ["FOXO3", "FOXO1"]                                # for U flag
DAMAGE_SIG = ["CHEK2", "CHEK1", "H2AFX", "MDC1", "RAD51", "BRCA1", "FANCD2",
              "RNF168", "TP53BP1", "GADD45A", "CDKN1A", "ATM", "ATR",
              "MRE11", "NBN", "RAD50"]
PROLIF = ["MKI67", "CCNB1", "CCNE1", "PCNA", "MCM2", "AURKA"]
APOP_PRO = ["BAX", "BAK1", "BBC3", "PMAIP1", "BCL2L11", "BID"]
APOP_ANTI = ["BCL2", "BCL2L1", "MCL1"]
RESIST = ["ABCB1", "ABCC1", "ABCG2"]
EXTRA = ["SRC", "RHEB", "IGF1R", "GRB2", "IRS1", "EGFR", "ERBB2", "CDKN2A", "E2F1"]

NEED = set(EXTRA + PROLIF + APOP_PRO + APOP_ANTI + DAMAGE_SIG + AKT_ACT + FOXO_NUC + RESIST)
for v in SEED_SIG.values():
    NEED |= set(v)

DRUG = ["RTK_up", "SRC", "RHEB", "IGF1R", "ATM", "ATR"]
DAMAGE_PCT = 0.80   # damage-activity gating threshold


# ---------------------------------------------------------------------------
# Loaders (cBioPortal / GEO formats)
# ---------------------------------------------------------------------------
def load_matrix(path, idcols=2, gzipped=False, header_idcols=None):
    """Load a gene x sample expression matrix; keep only NEED genes.

    header_idcols handles files whose header row has a different count of leading
    non-sample cells than data rows (the I-SPY2 GEO matrix: 988 sample IDs in the
    header, but data rows are gene + 988 values). Pass header_idcols=0, idcols=1.
    """
    hc = idcols if header_idcols is None else header_idcols
    opener = (lambda p: gzip.open(p, "rt")) if gzipped else open
    rows = {}
    with opener(path) as f:
        samples = f.readline().rstrip("\n").split("\t")[hc:]
        for line in f:
            p = line.rstrip("\n").split("\t")
            g = p[0].strip()
            if g in NEED:
                vals = [float(x) if x not in ("", "NA", "NaN") else np.nan
                        for x in p[idcols:]]
                if len(vals) != len(samples):
                    vals = (vals + [np.nan] * len(samples))[:len(samples)]
                rows[g] = np.array(vals)
    return samples, rows


def load_mutations(path):
    """Load a cBioPortal MAF; return (TP53_LOF_samples, PIK3CA_missense_samples)."""
    tp53, pik = set(), set()
    with open(path) as f:
        line = f.readline()
        while line.startswith("#"):
            line = f.readline()
        h = line.rstrip("\n").split("\t")
        gi, bi = h.index("Hugo_Symbol"), h.index("Tumor_Sample_Barcode")
        vi = h.index("Variant_Classification")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(gi, bi, vi):
                continue
            if p[gi] == "TP53" and p[vi] not in ("Silent", "3'UTR", "5'UTR", "Intron", "RNA"):
                tp53.add(p[bi])
            if p[gi] == "PIK3CA" and p[vi] == "Missense_Mutation":
                pik.add(p[bi])
    return tp53, pik


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------
def run_cohort(name, samples, rows, tp53_mut=None, pik_mut=None, idtrim=lambda x: x,
               method="mad"):
    """Run the full safeguarded pipeline over one cohort; print summary + checks.

    method : per-gene normalisation, one of 'mad' (default, robust), 'zscore', 'median'.
    """
    # guard: align samples AND data arrays to common width
    _w = min([len(v) for v in rows.values()]) if rows else len(samples)
    _w = min(_w, len(samples))
    samples = samples[:_w]
    for _g in list(rows.keys()):
        rows[_g] = rows[_g][:_w]
    n = _w
    t0 = time.time()

    def z(g):
        """Per-gene cohort normalisation. `method` selects centring/scaling:
          'mad'    : robust z = (x - median) / (1.4826 * MAD)   [DEFAULT, outlier-robust]
          'zscore' : classic z = (x - mean) / std
          'median' : hard split, +1 above median else -1 (discards magnitude)
        """
        if g not in rows:
            return np.zeros(n)
        v = rows[g]
        if method == "zscore":
            c = np.nanmean(v); s = np.nanstd(v) or 1.0
            r = (v - c) / s
        elif method == "median":
            c = np.nanmedian(v)
            r = np.where(v > c, 1.0, -1.0)
        else:  # 'mad' (default): robust z-score, preferred for skewed RNA-seq
            c = np.nanmedian(v)
            mad = np.nanmedian(np.abs(v - c))
            s = 1.4826 * mad or 1.0
            r = (v - c) / s
        r[np.isnan(r)] = 0.0
        return r

    Z = {g: z(g) for g in NEED}

    def seed(node):
        gs = [Z[g] for g in SEED_SIG[node] if g in rows]
        return (np.mean(gs, axis=0) > 0).astype(int) if gs else np.zeros(n, int)

    SEED = {nd: seed(nd) for nd in SEED_SIG}
    akt_act = np.mean([Z[g] for g in AKT_ACT if g in rows], axis=0)
    foxo_nuc = np.mean([Z[g] for g in FOXO_NUC if g in rows], axis=0)
    U = ((akt_act > 0) & (foxo_nuc > 0)).astype(int)               # uncoupled flag
    damage = np.mean([Z[g] for g in DAMAGE_SIG if g in rows], axis=0)
    dcut = np.quantile(damage, DAMAGE_PCT)
    prolif = np.mean([Z[g] for g in PROLIF if g in rows], axis=0)
    resist = (np.mean([Z[g] for g in RESIST if g in rows], axis=0)
              if any(g in rows for g in RESIST) else np.zeros(n))
    have_mut = tp53_mut is not None

    routing, kernels, targeted, records = [], [], 0, []
    for i in range(n):
        sid = idtrim(samples[i])
        dmg = int(damage[i] > dcut)
        rtk = int((Z["GRB2"][i] > 0) or (Z["IRS1"][i] > 0)) if "GRB2" in rows else 0
        if have_mut and sid in pik_mut:
            rtk = 1
        I = dict(SRC=int(Z["SRC"][i] > 0) if "SRC" in rows else 0,
                 RHEB=int(Z["RHEB"][i] > 0) if "RHEB" in rows else 0,
                 IGF1R=int(Z["IGF1R"][i] > 0) if "IGF1R" in rows else 0,
                 RTK_up=rtk,
                 CDKN2A=int(Z["CDKN2A"][i] > 0) if "CDKN2A" in rows else 0,
                 E2F1=int(Z["E2F1"][i] > 0) if "E2F1" in rows else 0,
                 ATM=dmg, ATR=dmg)
        x0 = {nd: int(SEED[nd][i]) for nd in NODES}
        if have_mut and sid in tp53_mut:      # mutation safeguard: p53 LOF inactive
            x0["TP53"] = 0
        lab = label(simulate(x0, I, int(U[i])))
        routing.append(lab)
        records.append((sid, lab, int(U[i])))
        if lab != "APOPTOTIC":                # minimal single-input kernel
            for k in DRUG:
                J = dict(I); J[k] = 0 if k in ("RTK_up", "SRC", "RHEB", "IGF1R") else 1
                if label(simulate(x0, J, int(U[i]))) == "APOPTOTIC":
                    kernels.append(k)
                    if k in ("RTK_up", "SRC", "RHEB", "IGF1R"):
                        targeted += 1
                    break
            else:
                kernels.append("none")

    c = Counter(routing); rt = np.array(routing)
    print(f"\n{name}: N={n} ({time.time()-t0:.1f}s) | uncoupled U=1 in {100*U.mean():.0f}%")
    print(f"  ROUTING: APOP {100*c['APOPTOTIC']/n:.0f}% | SURV {100*c['SURVIVAL']/n:.0f}% "
          f"| mixed {100*c['mixed']/n:.0f}%")
    sm, am = rt == "SURVIVAL", rt == "APOPTOTIC"
    surv_unc = res_surv = res_apop = None
    if sm.sum():
        surv_unc = float(U[sm].mean())
        res_surv = float(resist[sm].mean()); res_apop = float(resist[am].mean())
        print(f"  BIO-CHECK: SURV-routed uncoupled-rate {100*surv_unc:.0f}% vs cohort {100*U.mean():.0f}%")
        print(f"             resistance-sig SURV {res_surv:+.2f}z vs APOP {res_apop:+.2f}z")
    return dict(name=name, n=n, method=method, counts=dict(c), U=float(U.mean()),
                apop_pct=round(100*c['APOPTOTIC']/n), surv_pct=round(100*c['SURVIVAL']/n),
                mixed_pct=round(100*c['mixed']/n),
                surv_uncoupled_pct=(round(100*surv_unc) if surv_unc is not None else None),
                resist_surv=res_surv, resist_apop=res_apop, records=records)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Run the BBCN safeguarded cohort pipeline.")
    ap.add_argument("--expr", required=True, help="expression matrix (gene x sample)")
    ap.add_argument("--muts", default=None, help="cBioPortal MAF (optional)")
    ap.add_argument("--name", default="cohort")
    ap.add_argument("--gzip", action="store_true")
    ap.add_argument("--idtrim", type=int, default=0, help="truncate sample IDs to N chars (0=off)")
    ap.add_argument("--method", default="mad", choices=["mad", "zscore", "median"],
                    help="per-gene normalisation (default: mad = robust z, outlier-resistant)")
    a = ap.parse_args()
    samples, rows = load_matrix(a.expr, gzipped=a.gzip)
    tp53 = pik = None
    if a.muts:
        tp53, pik = load_mutations(a.muts)
    trim = (lambda x: x[:a.idtrim]) if a.idtrim else (lambda x: x)
    run_cohort(a.name, samples, rows, tp53, pik, idtrim=trim, method=a.method)
