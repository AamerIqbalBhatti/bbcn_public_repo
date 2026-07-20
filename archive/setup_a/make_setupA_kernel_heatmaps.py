#!/usr/bin/env python3
"""Render Setup A kernel figures from capture_setupA.json (no recompute)."""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(os.path.dirname(HERE), "preprint", "figures")
sys.path.insert(0, HERE)
from bbcn import harness as H

TEAL, GOLD, RED, NAVY, LGREY, BG, SLATE = (
    '#1B6E8B', '#C9992B', '#C0392B', '#0D1F3C', '#ECF0F1', '#FAFAFA', '#4A6170')
CMAP = LinearSegmentedColormap.from_list("teal", ["#FFFFFF", "#CFE2E9", TEAL, NAVY])

T = json.load(open(os.path.join(HERE, "capture_setupA.json")))
COH = ["TCGA", "METABRIC", "ISPY2"]
LBL = {"TCGA": "TCGA", "METABRIC": "METABRIC", "ISPY2": "I-SPY2"}
METHODS = ["ranked", "stabilize"]
N = T["N"]

# ---------- Figure A1: kernel usage by cohort, two method panels ----------
totals = {}
for key, d in T["node_pat"].items():
    for n, c in d.items():
        totals[n] = totals.get(n, 0) + c
TOPK = 24
nodes = [n for n, _ in sorted(totals.items(), key=lambda kv: -kv[1])[:TOPK]]

fig, axes = plt.subplots(1, 2, figsize=(10.6, 0.42 * len(nodes) + 2.0),
                         facecolor=BG, sharey=True)
for ax, method in zip(axes, METHODS):
    M = np.zeros((len(nodes), len(COH)))
    for j, ch in enumerate(COH):
        d = T["node_pat"].get(f"{ch}|{method}", {})
        for i, nd in enumerate(nodes):
            M[i, j] = 100.0 * d.get(nd, 0) / N[ch]
    im = ax.imshow(M, cmap=CMAP, vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(COH))); ax.set_xticklabels([LBL[c] for c in COH], fontsize=9)
    ax.set_yticks(range(len(nodes))); ax.set_yticklabels(nodes, fontsize=8)
    ax.set_title(method, fontsize=12, fontweight="bold", color=NAVY)
    for i in range(len(nodes)):
        for j in range(len(COH)):
            v = M[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    fontsize=7, color=("white" if v > 55 else SLATE))
    for sp in ax.spines.values():
        sp.set_visible(False)
cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
cbar.set_label("% of patients selecting node as a kernel", fontsize=9, color=SLATE)
fig.suptitle("Setup A  -  135-node controller: kernel usage by cohort (static, isolated/ideal)",
             fontsize=13, fontweight="bold", color=NAVY, y=0.995)
fig.text(0.5, 0.004,
         "Top 24 kernel nodes by total usage. Cell = percent of cohort patients for whom the node is "
         "selected as a stabilising kernel for at least one pathway. ranked reproduces the prior kernels.",
         ha="center", fontsize=8, color=SLATE)
os.makedirs(FIGDIR, exist_ok=True)
p1 = os.path.join(FIGDIR, "figA_kernel_usage.png")
fig.savefig(p1, dpi=300, bbox_inches="tight", facecolor=BG); plt.close(fig)

# ---------- Figure A2: pathway x kernel incidence (canonical = stabilize, pooled) ----------
method = "stabilize"
pooledN = sum(N[c] for c in COH)
pw_node = {}
for ch in COH:
    for kk, c in T["pw_node"].get(f"{ch}|{method}", {}).items():
        pw_node[kk] = pw_node.get(kk, 0) + c

pw_order = [p for p in H.PATHWAYS.keys()]
present_pw = [p for p in pw_order if any(k.split("|", 1)[0] == p for k in pw_node)]
# assign each node to its dominant pathway, order nodes by (pathway, -incidence) for block structure
node_best = {}
for kk, c in pw_node.items():
    pw, nd = kk.split("|", 1)
    if c > node_best.get(nd, (None, -1))[1]:
        node_best[nd] = (pw, c)
all_nodes = sorted(node_best.keys(),
                   key=lambda nd: (present_pw.index(node_best[nd][0]) if node_best[nd][0] in present_pw else 99,
                                   -node_best[nd][1]))

M = np.zeros((len(present_pw), len(all_nodes)))
for i, pw in enumerate(present_pw):
    for j, nd in enumerate(all_nodes):
        M[i, j] = 100.0 * pw_node.get(f"{pw}|{nd}", 0) / (3 * pooledN)

fig2, ax = plt.subplots(figsize=(0.32 * len(all_nodes) + 3, 0.42 * len(present_pw) + 2), facecolor=BG)
im = ax.imshow(M, cmap=CMAP, aspect="auto", vmin=0, vmax=max(1, M.max()))
ax.set_xticks(range(len(all_nodes))); ax.set_xticklabels(all_nodes, rotation=90, fontsize=7)
ax.set_yticks(range(len(present_pw))); ax.set_yticklabels(present_pw, fontsize=9)
ax.set_title("Setup A  -  pathway x kernel incidence (stabilize, three cohorts pooled)",
             fontsize=12, fontweight="bold", color=NAVY)
for sp in ax.spines.values():
    sp.set_visible(False)
cbar = fig2.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label("kernel-selection intensity (% of patient-phenotype problems)", fontsize=8, color=SLATE)
fig2.text(0.5, 0.002,
          "Near block-diagonal structure: each kernel node acts within its own pathway module, "
          "empirical support for the weak-coupling decomposition.",
          ha="center", fontsize=8, color=SLATE)
p2 = os.path.join(FIGDIR, "figA_pathway_kernel.png")
fig2.savefig(p2, dpi=300, bbox_inches="tight", facecolor=BG); plt.close(fig2)

print("WROTE", p1)
print("WROTE", p2)
print("A1 nodes:", nodes)
print("present pathways:", present_pw)
