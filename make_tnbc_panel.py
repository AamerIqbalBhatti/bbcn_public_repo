#!/usr/bin/env python3
"""
make_tnbc_panel.py -- receptor-low (TNBC-like) re-affirmation panel for Setup A.
Defines the subset by low ESR1, PGR, ERBB2 in the shipped binarised data (a
transcriptomic proxy for triple-negative status; no external PAM50 labels), then
runs the SAME Setup A static controller (ranked) on the full cohort and the
subset, and plots the individual-vs-joint contrast. Output:
  preprint/figures/figA_tnbc_panel.png
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
SETUP_A = os.path.join(HERE, "setup_a")
sys.path.insert(0, SETUP_A)
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bbcn import controller as C, harness as H

TEAL, GOLD, RED, NAVY, LGREY, BG, SLATE = (
    '#1B6E8B', '#C9992B', '#C0392B', '#0D1F3C', '#ECF0F1', '#FAFAFA', '#4A6170')
BIN = {"TCGA": "tcga_brca_1082x135.csv", "METABRIC": "metabric_1980x135.csv", "ISPY2": "ispy2_988x135.csv"}
LBL = {"TCGA": "TCGA", "METABRIC": "METABRIC", "ISPY2": "I-SPY2"}
PHEN = ["Resistance_OFF", "Apoptosis_ON", "Proliferation_OFF"]
NODES = set(H.ALL_NODES)


def run_group(pts):
    n = len(pts)
    iso = {p: 0 for p in PHEN}; allthree = 0
    for init in pts:
        o = C.run_patient(init, mode="isolated", family="ideal", kernel_method="ranked")
        for p in PHEN: iso[p] += int(o[p])
    for init in pts:
        o = C.run_patient(init, mode="sequenced", family="compatible", kernel_method="ranked")
        allthree += int(o.get("all_three", False))
    return dict(n=n,
                apop=round(100*iso["Apoptosis_ON"]/n),
                prolif=round(100*iso["Proliferation_OFF"]/n),
                allthree=round(100*allthree/n))


res = {}
for cohort, fn in BIN.items():
    B = pd.read_csv(os.path.join(SETUP_A, "data", "binarized", fn))
    cols = [c for c in B.columns if c in NODES]
    recneg = (B["ESR1"] == 0) & (B["PGR"] == 0) & (B["ERBB2"] == 0)
    full = [{c: int(r[c]) for c in cols} for _, r in B.iterrows()]
    sub = [{c: int(r[c]) for c in cols} for _, r in B[recneg].iterrows()]
    res[cohort] = {"full": run_group(full), "tnbc": run_group(sub),
                   "tnbc_n": int(recneg.sum()), "tnbc_pct": round(100*recneg.mean())}
    print(cohort, res[cohort])

json.dump(res, open(os.path.join(HERE, "tnbc_panel.json"), "w"), indent=1)

# ---- plot: 3 metrics x (full vs TNBC-like), grouped by cohort ----
metrics = [("apop", "Apoptosis-ON\n(isolated)"), ("prolif", "Proliferation-OFF\n(isolated)"),
           ("allthree", "All three\n(sequenced)")]
cohorts = ["TCGA", "METABRIC", "ISPY2"]
fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), facecolor=BG)
x = np.arange(len(cohorts)); w = 0.38
for ax, (key, title) in zip(axes, metrics):
    full_v = [res[c]["full"][key] for c in cohorts]
    tnbc_v = [res[c]["tnbc"][key] for c in cohorts]
    b1 = ax.bar(x - w/2, full_v, w, label="full cohort", color=TEAL, edgecolor=NAVY, linewidth=0.6)
    b2 = ax.bar(x + w/2, tnbc_v, w, label="TNBC-like (receptor-low)", color=GOLD, edgecolor=NAVY, linewidth=0.6)
    for b in list(b1) + list(b2):
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+1.5, f"{int(b.get_height())}",
                ha="center", va="bottom", fontsize=8, color=SLATE)
    ax.set_title(title, fontsize=11, fontweight="bold", color=NAVY)
    ax.set_xticks(x); ax.set_xticklabels([LBL[c] for c in cohorts], fontsize=9)
    ax.set_ylim(0, 105); ax.set_facecolor(BG)
    ax.spines[["top", "right"]].set_visible(False)
    if key == "apop":
        ax.set_ylabel("% of patients reaching a fixed point", fontsize=9, color=SLATE)
axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9, frameon=False)
fig.suptitle("Setup A re-affirmation on the receptor-low (TNBC-like) subset",
             fontsize=13, fontweight="bold", color=NAVY, y=1.02)
ns = ", ".join(f"{LBL[c]} n={res[c]['tnbc_n']} ({res[c]['tnbc_pct']}%)" for c in cohorts)
fig.text(0.5, -0.06, "TNBC-like = ESR1, PGR, ERBB2 all low in binarised data (transcriptomic proxy).  "
         + ns + ".  Same capped no-forcing controller, ranked method.",
         ha="center", fontsize=8, color=SLATE)
out = os.path.join(HERE, "preprint", "figures", "figA_tnbc_panel.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=BG)
print("WROTE", out)
