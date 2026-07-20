#!/usr/bin/env python3
"""Two figures: (1) the Setup B switch-test mechanism as a flow diagram,
(2) the held-vs-durable comparison for Setup A (135-node) and Setup B (9-node core)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- flow diagram
fig, ax = plt.subplots(figsize=(8.2, 10.6)); ax.set_xlim(0, 10); ax.set_ylim(0, 22); ax.axis("off")
def box(y, txt, w=8.6, h=1.5, x=0.7, fc="#eef3fb", ec="#2c4a78", fs=10.5, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.12",
                                fc=fc, ec=ec, lw=1.6))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color="#11243f")
def arrow(y1, y2, x=5.0, label=None):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>", mutation_scale=16,
                                 lw=1.6, color="#2c4a78"))
    if label: ax.text(x + 0.25, (y1 + y2) / 2, label, ha="left", va="center", fontsize=9, color="#555")

ax.text(5, 21.4, "Setup B switch test — one patient", ha="center", fontsize=13, fontweight="bold", color="#11243f")
box(19.4, "Patient transcriptome (full gene matrix)\nthe only place the full data enters", fc="#fde9d9", ec="#b5651d", bold=True)
arrow(19.4, 18.7)
box(17.2, "Per-gene robust z-score (MAD)", h=1.3)
arrow(17.2, 16.5)
box(14.7, "Seed the 9 switch nodes by majority gate\n+ set 8 held inputs (SRC, RHEB, IGF1R, RTK_up,\nCDKN2A, E2F1, ATM, ATR)  + uncoupled flag U", h=2.0, fc="#e8f6ec", ec="#2e7d4f")
arrow(14.7, 14.0)
box(12.2, "Multirate switch:  FAST x1,  MID x5,  SLOW x25\nrun 9-node core to its attractor  (NOT the 135-node net)", h=2.0, fc="#e9eefc", ec="#2c4a78", bold=True)
arrow(12.2, 11.5)
box(10.0, "Route:  APOPTOTIC (TP53=1 & AKT1=0)\n        or  SURVIVAL (AKT1=1)", h=1.9)
arrow(10.0, 9.3, label="if SURVIVAL (resistant)")
box(7.6, "Design apoptosis kernel on FAST/MID boxes\n(forward-stabilization) -> small clamp set\ne.g. AKT1=0, PTEN=1", h=2.0, fc="#fdeef0", ec="#9b2d3a")
arrow(7.6, 6.9)
box(5.1, "PULSE: hold the clamp for P ticks  (drug ON)\nthen RELEASE  (drug OFF)", h=1.9, fc="#fff6da", ec="#9a7b15", bold=True)
arrow(5.1, 4.4)
box(2.6, "Free multirate relaxation to attractor\n(no clamp)", h=1.7)
arrow(2.6, 1.9)
box(0.2, "DURABLE if still APOPTOTIC after release\n= the attractor holds itself, drug off", h=1.6, fc="#e8f6ec", ec="#2e7d4f", bold=True)
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_switch_flow.png"), dpi=140, bbox_inches="tight"); plt.close()

# ---------------------------------------------------------------- durability bars
coh = ["TCGA", "METABRIC", "I-SPY2"]
A_held = [91, 90, 83]; A_dur = [0, 0, 0]          # Setup A, % of all patients
B_held = [24, 18, 17]; B_dur = [24, 18, 17]        # Setup B, % of resistant patients
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6)); x = np.arange(3); w = 0.36
for ax, held, dur, title, sub in [
    (axes[0], A_held, A_dur, "Setup A — full 135-node network", "% of all patients reaching apoptosis"),
    (axes[1], B_held, B_dur, "Setup B — 9-node bistable switch", "% of resistant patients"),
]:
    ax.bar(x - w/2, held, w, label="reached while drug HELD", color="#7aa6d6", ec="#2c4a78")
    ax.bar(x + w/2, dur, w, label="DURABLE after drug withdrawn", color="#2e7d4f", ec="#1c4a2f")
    for i, (h, d) in enumerate(zip(held, dur)):
        ax.text(i - w/2, h + 1.5, f"{h}", ha="center", fontsize=9)
        ax.text(i + w/2, d + 1.5, f"{d}", ha="center", fontsize=9, color="#1c4a2f")
    ax.set_title(title, fontsize=11, fontweight="bold"); ax.set_xticks(x); ax.set_xticklabels(coh)
    ax.set_ylim(0, 100); ax.set_ylabel(sub, fontsize=9); ax.legend(fontsize=8.5, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
axes[0].text(1, 60, "withdraw → collapses to 0\n(no durable attractor)", ha="center", fontsize=9.5, color="#9b2d3a", style="italic")
axes[1].text(1, 60, "withdraw → holds\n(durable = reachable;\nhysteresis latch)", ha="center", fontsize=9.5, color="#1c4a2f", style="italic")
plt.suptitle("Durability of the attractor after the drug is withdrawn", fontsize=13, fontweight="bold")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fig_durability.png"), dpi=140, bbox_inches="tight"); plt.close()
print("wrote", os.path.join(OUT, "fig_switch_flow.png"))
print("wrote", os.path.join(OUT, "fig_durability.png"))
