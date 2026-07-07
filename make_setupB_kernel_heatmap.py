#!/usr/bin/env python3
"""
Setup B kernel-actuator usage heat map.
Drives the repo's own forward_stab_kernel_design.run(); collects the full node_use
distribution it already computes (the pipeline emits only its argmax as top_node).
No reimplementation of kernel selection. Renders at 300 DPI in the repo style.
"""
import os, sys, io, contextlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
SB = os.path.join(HERE, "setup_b", "code")
SAMP = os.path.join(HERE, "setup_b", "data", "samples")
sys.path.insert(0, SB)
import forward_stab_kernel_design as KD

# repo style
TEAL, GOLD, RED, NAVY, LGREY, BG, SLATE = (
    '#1B6E8B', '#C9992B', '#C0392B', '#0D1F3C', '#ECF0F1', '#FAFAFA', '#4A6170')
CMAP = LinearSegmentedColormap.from_list("teal", ["#FFFFFF", "#CFE2E9", TEAL, NAVY])

COHORTS = [("TCGA", "tcga_brca_switch_inputs.tsv"),
           ("METABRIC", "metabric_switch_inputs.tsv"),
           ("I-SPY2", "ispy2_switch_inputs.tsv")]
METHODS = ["ranked", "stabilize"]


def load(path):
    for kw in (dict(idcols=2, gz=False, header_idcols=None),
               dict(idcols=1, gz=False, header_idcols=None)):
        try:
            s, r = KD.load_matrix(path, **kw)
            if len(s) > 10:
                return s, r
        except Exception:
            continue
    raise RuntimeError("could not load " + path)


# collect node_use and ker_designed for every cohort x method
data = {}   # (method) -> {cohort: (node_use dict, designed)}
for method in METHODS:
    data[method] = {}
for name, fn in COHORTS:
    s, r = load(os.path.join(SAMP, fn))
    for method in METHODS:
        with contextlib.redirect_stdout(io.StringIO()):
            kr = KD.run(name, s, r, None, None, idtrim=lambda x: x, N=None, method=method)
        data[method][name] = (kr["node_use"], kr["ker_designed"])

# node ordering by total usage across all cohorts/methods
totals = {}
for method in METHODS:
    for name, _ in COHORTS:
        nu, _d = data[method][name]
        for nd, c in nu.items():
            totals[nd] = totals.get(nd, 0) + c
nodes = [n for n, _ in sorted(totals.items(), key=lambda kv: -kv[1])]

fig, axes = plt.subplots(1, 2, figsize=(10.2, 0.55 * len(nodes) + 2.2),
                         facecolor=BG, sharey=True)
cohort_names = [c for c, _ in COHORTS]
for ax, method in zip(axes, METHODS):
    M = np.zeros((len(nodes), len(cohort_names)))
    for j, name in enumerate(cohort_names):
        nu, designed = data[method][name]
        for i, nd in enumerate(nodes):
            M[i, j] = 100.0 * nu.get(nd, 0) / designed if designed else 0.0
    im = ax.imshow(M, cmap=CMAP, vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(cohort_names)))
    ax.set_xticklabels(cohort_names, fontsize=10)
    ax.set_yticks(range(len(nodes)))
    ax.set_yticklabels(nodes, fontsize=10)
    ax.set_title(f"{method}", fontsize=12, fontweight="bold", color=NAVY)
    for i in range(len(nodes)):
        for j in range(len(cohort_names)):
            nu, designed = data[method][cohort_names[j]]
            cnt = nu.get(nodes[i], 0)
            val = M[i, j]
            ax.text(j, i, f"{cnt}\n{val:.0f}%", ha="center", va="center",
                    fontsize=8, color=("white" if val > 55 else SLATE))
    ax.set_facecolor(BG)
    for sp in ax.spines.values():
        sp.set_visible(False)

cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
cbar.set_label("% of designed kernels pinning this node", fontsize=9, color=SLATE)
fig.suptitle("Setup B  -  resistance-apoptosis switch: kernel-actuator usage by cohort",
             fontsize=13, fontweight="bold", color=NAVY, y=0.99)
fig.text(0.5, 0.005,
         "Per-box forward-stabilisation kernels (PI3K/AKT axis). Cell = count of designed kernels pinning the node, "
         "and that count as a percentage of all designed kernels in the cohort.",
         ha="center", fontsize=8, color=SLATE)

out = os.path.join(HERE, "preprint", "figures", "figB_kernel_actuator.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=BG)
print("WROTE", out)
print("nodes:", nodes)
for method in METHODS:
    print(method, {c: (data[method][c][0], "designed=%d" % data[method][c][1]) for c in cohort_names})
