#!/usr/bin/env python3
"""
make_tnbc_setupB.py -- receptor-low (TNBC-like) re-affirmation panel for Setup B.
Routes patients through the bistable switch on the shipped switch inputs (no
mutation safeguard, so percentages sit within ~2 points of the Table values),
then compares the full cohort with the receptor-low subset on routing and on the
headline survival-routed uncoupled enrichment. Subset is defined by low ESR1,
PGR, ERBB2 in the binarised Setup A data, joined by sample ID.
Output: preprint/figures/figB_tnbc_panel.png
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "setup_b", "code"))
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cohort_pipeline as CP

TEAL, GOLD, NAVY, BG, SLATE = '#1B6E8B', '#C9992B', '#0D1F3C', '#FAFAFA', '#4A6170'
TSV = {"TCGA": "tcga_brca_switch_inputs.tsv", "METABRIC": "metabric_switch_inputs.tsv", "ISPY2": "ispy2_switch_inputs.tsv"}
BIN = {"TCGA": "tcga_brca_1082x135.csv", "METABRIC": "metabric_1980x135.csv", "ISPY2": "ispy2_988x135.csv"}
LBL = {"TCGA": "TCGA", "METABRIC": "METABRIC", "ISPY2": "I-SPY2"}


def recneg_ids(cohort):
    B = pd.read_csv(os.path.join(HERE, "setup_a", "data", "binarized", BIN[cohort]))
    m = (B["ESR1"] == 0) & (B["PGR"] == 0) & (B["ERBB2"] == 0)
    return set(B.loc[m, "SAMPLE_ID"].astype(str)), set(B["SAMPLE_ID"].astype(str))


def stats(recs):
    n = len(recs)
    routes = [r[1] for r in recs]
    apop = round(100 * sum(x == "APOPTOTIC" for x in routes) / n)
    surv = round(100 * sum(x == "SURVIVAL" for x in routes) / n)
    mixed = round(100 * sum(x == "mixed" for x in routes) / n)
    surv_recs = [r for r in recs if r[1] == "SURVIVAL"]
    surv_unc = round(100 * np.mean([r[2] for r in surv_recs])) if surv_recs else 0
    return dict(n=n, apop=apop, surv=surv, mixed=mixed, surv_unc=surv_unc)


res = {}
for cohort in ["TCGA", "METABRIC", "ISPY2"]:
    s, r = CP.load_matrix(os.path.join(HERE, "setup_b", "data", "samples", TSV[cohort]), idcols=1, gzipped=False)
    out = CP.run_cohort(cohort, s, r, None, None, idtrim=lambda x: x, method="mad")
    recs = out["records"]                       # [(sid, route, U)]
    neg, allids = recneg_ids(cohort)
    full = recs
    sub = [rec for rec in recs if str(rec[0]) in neg]
    res[cohort] = {"full": stats(full), "tnbc": stats(sub), "tnbc_n": len(sub)}
    print(cohort, "full", res[cohort]["full"], "| tnbc", res[cohort]["tnbc"])

json.dump(res, open(os.path.join(HERE, "tnbc_setupB.json"), "w"), indent=1)

metrics = [("apop", "Apoptotic-routed"), ("surv", "Survival-routed"),
           ("surv_unc", "Survival-routed\nin uncoupled state")]
cohorts = ["TCGA", "METABRIC", "ISPY2"]
fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), facecolor=BG)
x = np.arange(len(cohorts)); w = 0.38
for ax, (key, title) in zip(axes, metrics):
    fv = [res[c]["full"][key] for c in cohorts]
    tv = [res[c]["tnbc"][key] for c in cohorts]
    b1 = ax.bar(x - w/2, fv, w, label="full cohort", color=TEAL, edgecolor=NAVY, linewidth=0.6)
    b2 = ax.bar(x + w/2, tv, w, label="TNBC-like (receptor-low)", color=GOLD, edgecolor=NAVY, linewidth=0.6)
    for b in list(b1) + list(b2):
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+1.5, f"{int(b.get_height())}",
                ha="center", va="bottom", fontsize=8, color=SLATE)
    ax.set_title(title, fontsize=11, fontweight="bold", color=NAVY)
    ax.set_xticks(x); ax.set_xticklabels([LBL[c] for c in cohorts], fontsize=9)
    ax.set_ylim(0, 105); ax.set_facecolor(BG)
    ax.spines[["top", "right"]].set_visible(False)
    if key == "apop":
        ax.set_ylabel("% of patients", fontsize=9, color=SLATE)
axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9, frameon=False)
fig.suptitle("Setup B re-affirmation on the receptor-low (TNBC-like) subset",
             fontsize=13, fontweight="bold", color=NAVY, y=1.02)
ns = ", ".join(f"{LBL[c]} n={res[c]['tnbc_n']}" for c in cohorts)
fig.text(0.5, -0.06, "Routing on shipped switch inputs (no mutation safeguard; within ~2 points of Table values).  "
         "TNBC-like joined by sample ID.  " + ns + ".",
         ha="center", fontsize=8, color=SLATE)
out = os.path.join(HERE, "preprint", "figures", "figB_tnbc_panel.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=BG)
print("WROTE", out)
