#!/usr/bin/env python3
"""Figure 1 (the two setups of BBCN), regenerated from code.
Legible fonts, short wrapped lines, no node counts baked in (the count lives in
the caption). Run: python3 fig1.py  ->  figures/fig_arch.png
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)

CONT_B, CONT_R = "#dce7f4", "#f7dede"
BAR_B,  BAR_R  = "#4472a8", "#b5504f"
GREEN, AMBER, GREY = "#dcefe0", "#faf0d7", "#ededed"
EB, ER = "#2f5580", "#8f3a3a"

FIG_W, FIG_H = 9.3, 7.8
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, text, fc, ec="#555", fs=12, bold=False, tc="black", lw=1.3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35,rounding_size=1.4",
                 fc=fc, ec=ec, lw=lw))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            fontweight=("bold" if bold else "normal"), color=tc, zorder=5, linespacing=1.25)

def arrow(x1, y1, x2, y2, color="#666", lw=2.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=16, lw=lw, color=color, shrinkA=1, shrinkB=1, zorder=1))

# ---- shared front end ----
box(30, 90.5, 40, 7.5, "Patient mRNA\n(TCGA / METABRIC / I-SPY2)", GREY, fs=13, bold=True)
box(28, 81.5, 44, 6.5, "robust z-score (median / MAD)\n$\\rightarrow$ Boolean state", GREY, fs=12)
arrow(50, 90.5, 50, 88.0)
arrow(45, 81.5, 25, 78.7)
arrow(55, 81.5, 75, 78.7)

# ---- header bars ----
box(2, 74.0, 46, 4.7, "SETUP A:  whole-network control", BAR_B, ec=EB, fs=13, bold=True, tc="white")
box(52, 74.0, 46, 4.7, "SETUP B:  bistable switch",       BAR_R, ec=ER, fs=13, bold=True, tc="white")

# ---- column layout: stack of (text, fs) in fixed-height boxes ----
def column(x, w, bar_x_center, boxes, content_fc, ec, result, result_fc, result_ec):
    top = 72.0            # top of first content box
    h   = 8.6             # box height
    gap = 1.7             # vertical gap
    y = top
    centers = []
    for text, fs in boxes:
        y -= h
        box(x, y, w, h, text, content_fc, ec=ec, fs=fs)
        centers.append(y + h/2)
        y -= gap
    # arrows: bar -> first, then between boxes
    arrow(bar_x_center, 74.0, bar_x_center, top)
    tops = [top] + [top - i*(h+gap) for i in range(1, len(boxes))]
    for i in range(len(boxes)-1):
        arrow(x+w/2, tops[i]-h, x+w/2, tops[i+1])
    # result box
    ry = 12.0
    box(x, ry, w, 9.0, result, result_fc, ec=result_ec, fs=12.5, bold=True)
    arrow(x+w/2, y+gap, x+w/2, ry+9.0)
    return centers, (x+w/2, ry+9.0/2)

A_boxes = [
 ("BBCN Boolean network on the bus\n24 pathways $\\cdot$ $x(t{+}1)=F(x(t))$", 12),
 ("Tier 3:  phenotype sequence\nResistance-OFF $\\rightarrow$ Apoptosis-ON\n$\\rightarrow$ Proliferation-OFF", 11.5),
 ("Tier 2 (21-day cycle):\npathway selector, static\nor by failure class", 11.5),
 ("Tier 1 (1-day step):\nforward-stabilisation kernel\nminimal kernel ($\\leq$2, escalate 3)", 11.5),
]
B_boxes = [
 ("Apoptosis core $\\rightarrow$ 9-node\nAKT1 / TP53 switch (bistable)", 12),
 ("delayed multirate engine\nFAST(1) / MID(5) / SLOW(25)\n$X(t{+}1)=C\\,L\\,X(t)$", 11.5),
 ("two attractors of $CL$ (inverse):\nSURVIVAL (AKT$\\uparrow$, p53$\\downarrow$)\nAPOPTOTIC (p53$\\uparrow$, AKT$\\downarrow$)", 11.5),
 ("forward-stabilisation kernel\nper box: reach apoptotic\nfixed point on $CL$", 11.5),
]
cA, rA = column(2, 46, 25, A_boxes, CONT_B, EB,
                "RESULT: phenotypes individually\nreachable, jointly near-exclusive\n($\\sim$2\u20133% all three)", GREEN, "#3c7a4e")
cB, rB = column(52, 46, 75, B_boxes, CONT_R, ER,
                "RESULT: 3-cohort routing,\nhysteresis, resistance;\nPI3K/AKT/mTOR nomination", AMBER, "#b8860b")

# ---- bottleneck cross-arrow (A result -> B core) ----
arrow(48, 16.5, 52, cB[0], color="#888", lw=2.0)
ax.text(50, 34, "apoptosis bottleneck\n$\\rightarrow$ re-model the core", ha="center", va="center",
        fontsize=11, style="italic", color="#555",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#ccc", lw=1.0))

fig.tight_layout(pad=0.5)
out = os.path.join(HERE, "figures", "fig_arch.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", os.path.relpath(out, HERE))
