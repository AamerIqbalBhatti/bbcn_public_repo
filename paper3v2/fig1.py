#!/usr/bin/env python3
"""Figure 1 (the two setups of BBCN), regenerated from code. Column titles are in
separate header bars (not overlapping the content boxes); Setup B's two attractors
are labelled in words as each other's inverse; retired terms replaced.
Run: python3 fig1.py  ->  figures/fig_arch.png
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)

CONT_B, CONT_R = "#dce7f4", "#f7dede"      # content fills
BAR_B,  BAR_R  = "#4472a8", "#b5504f"      # header-bar fills
GREEN, AMBER, GREY = "#dcefe0", "#faf0d7", "#ededed"
EB, ER = "#2f5580", "#8f3a3a"

fig, ax = plt.subplots(figsize=(11.0, 8.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, text, fc, ec="#555", fs=9, bold=False, tc="black", lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35,rounding_size=1.4",
                 fc=fc, ec=ec, lw=lw))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            fontweight=("bold" if bold else "normal"), color=tc, zorder=5)

def arrow(x1, y1, x2, y2, color="#666", lw=1.6):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=13, lw=lw, color=color, shrinkA=1, shrinkB=1, zorder=1))

# ---- shared front end ----
box(33, 93, 34, 5.5, "Patient mRNA  (TCGA / METABRIC / I-SPY2)", GREY, fs=9.5, bold=True)
box(31, 85.5, 38, 5, "robust z-score (median / MAD)  $\\rightarrow$  Boolean state", GREY, fs=9)
arrow(50, 93, 50, 90.5)
arrow(46, 85.5, 26, 83.2)   # to header bar A
arrow(54, 85.5, 74, 83.2)   # to header bar B

# ---- header bars (titles clearly outside the content boxes) ----
box(4, 78.6, 42, 4.4, "SETUP A  \u2014  whole-network controllability (135 nodes)",
    BAR_B, ec=EB, fs=9, bold=True, tc="white")
box(54, 78.6, 42, 4.4, "SETUP B  \u2014  bistable switch (136 nodes: + PHLPP)",
    BAR_R, ec=ER, fs=9, bold=True, tc="white")

# ---- Setup A content ----
ax_x, ax_w = 4, 42
A = [
 (69.5, "BBCN Boolean network on the bus\n135 nodes \u00b7 24 pathways \u00b7 $x(t{+}1)=F(x(t))$"),
 (60.5, "Tier 3  \u2014  phenotype sequence\nResistance-OFF $\\rightarrow$ Apoptosis-ON $\\rightarrow$ Proliferation-OFF"),
 (51.5, "Tier 2 (21-day cycle)  \u2014  pathway selector\nstatic set, or dynamic by failure class"),
 (42.5, "Tier 1 (1-day step)  \u2014  forward-stabilisation kernel design\nminimal kernel ($\\leq$2, escalate 3); pin, hold, relax"),
]
arrow(ax_x+ax_w/2, 78.6, ax_x+ax_w/2, 77)   # bar -> first box
for y, t in A: box(ax_x, y, ax_w, 7.0, t, CONT_B, ec=EB, fs=8.4)
for i in range(len(A)-1): arrow(ax_x+ax_w/2, A[i][0], ax_x+ax_w/2, A[i+1][0]+7.0)
arrow(ax_x+ax_w/2, A[-1][0], ax_x+ax_w/2, 34.5)
box(ax_x, 26, ax_w, 8, "RESULT: phenotypes individually reachable,\njointly near-exclusive ($\\sim$2\u20133% all three)",
    GREEN, ec="#3c7a4e", fs=8.6, bold=True)

# ---- Setup B content ----
bx_x, bx_w = 54, 42
B = [
 (69.5, "Apoptosis core $\\rightarrow$ 9-node AKT1\u2013TP53 switch\nmutual antagonism (bistable)"),
 (60.5, "delayed (autoregressive) multirate engine\nFAST(1) / MID(5) / SLOW(25)  $\\rightarrow$  $X(t{+}1)=CL\\,X(t)$"),
 (51.5, "two attractors of $CL$, each the other's inverse:\nSURVIVAL (AKT$\\uparrow$, p53$\\downarrow$)  \u00b7  APOPTOTIC (p53$\\uparrow$, AKT$\\downarrow$)"),
 (42.5, "forward-stabilisation kernel per box\n(reachability on $CL$ to the apoptotic fixed point)"),
]
arrow(bx_x+bx_w/2, 78.6, bx_x+bx_w/2, 77)
for y, t in B: box(bx_x, y, bx_w, 7.0, t, CONT_R, ec=ER, fs=8.4)
for i in range(len(B)-1): arrow(bx_x+bx_w/2, B[i][0], bx_x+bx_w/2, B[i+1][0]+7.0)
arrow(bx_x+bx_w/2, B[-1][0], bx_x+bx_w/2, 34.5)
box(bx_x, 26, bx_w, 8, "RESULT: 3-cohort routing + hysteresis,\nresistance; PI3K/AKT/mTOR nomination",
    AMBER, ec="#b8860b", fs=8.6, bold=True)

# ---- bottleneck cross-arrow ----
arrow(ax_x+ax_w, 30, bx_x, 71, color="#888", lw=1.7)
ax.text(50, 38, "apoptosis bottleneck\n$\\rightarrow$ re-model the core", ha="center", va="center",
        fontsize=8.0, style="italic", color="#555",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ccc", lw=0.8))

fig.tight_layout(pad=0.4)
out = os.path.join(HERE, "figures", "fig_arch.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", os.path.relpath(out, HERE))
