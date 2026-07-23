#!/usr/bin/env python3
"""Regenerate Figure 6 (kernel size + durability, switch vs whole-network) as a
vector PDF for LaTeX. Reads ../setup_a/data/cde_vs_switch_summary.csv.
Run: python3 fig6.py   ->  figures/f06_durability.pdf
"""
import os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
rows = {r["cohort"]: r for r in csv.DictReader(open(os.path.join(HERE, "..", "setup_a", "data", "cde_vs_switch_summary.csv")))}
order = ["TCGA", "METABRIC", "ISPY2"]; labels = ["TCGA", "METABRIC", "I-SPY2"]
S = [float(rows[c]["switch_mean_size"]) for c in order]
W = [float(rows[c]["cde_union_mean_size"]) for c in order]
SD = [float(rows[c]["switch_durable_pct"]) for c in order]
WD = [float(rows[c]["cde_union_durable_pct"]) for c in order]
x = np.arange(3); w = 0.38; SW, WN = "#1b9e77", "#7f7f7f"
fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.5, 3.1))
axA.bar(x-w/2, S, w, label="Switch kernel", color=SW); axA.bar(x+w/2, W, w, label="Whole-network controller", color=WN)
for i in range(3):
    axA.text(x[i]-w/2, S[i]+0.6, f"{S[i]:.2f}", ha="center", fontsize=8)
    axA.text(x[i]+w/2, W[i]+0.6, f"{W[i]:.2f}", ha="center", fontsize=8)
    axA.text(x[i], max(W)+3, f"{W[i]/S[i]:.0f}\u00d7", ha="center", fontsize=9, fontweight="bold", color="#333")
axA.set_xticks(x); axA.set_xticklabels(labels); axA.set_ylim(0, 40); axA.set_ylabel("Mean kernel size (nodes)")
axA.set_title("A  Kernel size", loc="left", fontsize=11, fontweight="bold")
axA.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(0.02, 0.62)); axA.spines[["top","right"]].set_visible(False)
axB.bar(x-w/2, SD, w, color=SW); axB.bar(x+w/2, WD, w, color=WN)
for i in range(3):
    axB.text(x[i]-w/2, SD[i]+0.4, f"{SD[i]:.1f}", ha="center", fontsize=8)
    axB.text(x[i]+w/2, WD[i]+0.4, f"{WD[i]:.1f}", ha="center", fontsize=8)
axB.set_xticks(x); axB.set_xticklabels(labels); axB.set_ylim(80, 100)
axB.set_ylabel("Durability held after\nkernel withdrawal (%)")
axB.set_title("B  Durability at equal reach (death reached 100% either way)", loc="left", fontsize=9.5, fontweight="bold")
axB.spines[["top","right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figures", "f06_durability.pdf"), bbox_inches="tight")
print("wrote figures/f06_durability.pdf")
