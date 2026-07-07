#!/usr/bin/env python3
"""Appendix B flowcharts: the ranked heuristic and the algebraic stabilize test."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

TEAL, GOLD, RED, NAVY, LGREY, BG, SLATE = (
    '#1B6E8B', '#C9992B', '#C0392B', '#0D1F3C', '#ECF0F1', '#FAFAFA', '#4A6170')
FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "preprint", "figures")


def flowchart(steps, title, fname, accent=TEAL):
    """steps: list of (text, kind) where kind in {io, proc, dec, term}."""
    n = len(steps)
    fig, ax = plt.subplots(figsize=(7.0, 0.95 * n + 1.0), facecolor=BG)
    ax.set_xlim(0, 10); ax.set_ylim(0, n); ax.axis("off")
    ax.set_facecolor(BG)
    y = n - 0.5
    centers = []
    colors = {"io": LGREY, "proc": "#FFFFFF", "dec": "#FBF1D9", "term": accent}
    edges = {"io": SLATE, "proc": accent, "dec": GOLD, "term": NAVY}
    for text, kind in steps:
        w, h = 8.4, 0.66
        x0 = 5 - w / 2
        box = FancyBboxPatch((x0, y - h / 2), w, h,
                             boxstyle="round,pad=0.02,rounding_size=0.12",
                             linewidth=1.6, edgecolor=edges[kind],
                             facecolor=colors[kind])
        ax.add_patch(box)
        tcol = "white" if kind == "term" else NAVY
        ax.text(5, y, text, ha="center", va="center", fontsize=9.2,
                color=tcol, wrap=True)
        centers.append(y)
        y -= 1.0
    for i in range(len(centers) - 1):
        ar = FancyArrowPatch((5, centers[i] - 0.33), (5, centers[i + 1] + 0.33),
                             arrowstyle="-|>", mutation_scale=14,
                             linewidth=1.5, color=SLATE)
        ax.add_patch(ar)
    ax.set_title(title, fontsize=12.5, fontweight="bold", color=NAVY, pad=12)
    out = os.path.join(FIGDIR, fname)
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=BG)
    plt.close(fig)
    print("WROTE", out)


ranked = [
    ("Inputs: pathway rules, target state z, patient state x0", "io"),
    ("Candidate pool = all pathway nodes; size k = 1", "proc"),
    ("Enumerate clamp sets S of size k", "proc"),
    ("For each S: pin S to target, simulate x(t+1)=T.L.x(t) from x0", "proc"),
    ("Does the trajectory reach z?", "dec"),
    ("Rank converging S: impact, delta-match, steps, causal score", "proc"),
    ("Any converging S at this size?  If not, k = k+1 (up to 3)", "dec"),
    ("Select top-ranked S as kernel K", "term"),
]

stabilize = [
    ("Inputs: pathway rules -> transition matrix L, target z", "io"),
    ("Candidate pool = all pathway nodes; size k = 1", "proc"),
    ("Enumerate clamp sets S of size k", "proc"),
    ("Build T_S^z;  L_c = T_S^z . L;  transient power rho", "proc"),
    ("Form Omega = L_c^rho . T_S^z", "proc"),
    ("Every column of Omega = target index?  (global stabilisation)", "dec"),
    ("Keep minimal-size stabilising S; else k = k+1", "dec"),
    ("Tiebreak: least causal score, then lexicographic -> kernel K", "term"),
]

flowchart(ranked, "Appendix B (i)  -  Ranked reachability heuristic", "figC_flow_ranked.png", accent=TEAL)
flowchart(stabilize, "Appendix B (ii)  -  Algebraic global-stabilisation test", "figC_flow_stabilize.png", accent=GOLD)
print("done")
