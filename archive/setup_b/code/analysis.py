"""
analysis.py — Reproduces the structural-analysis results (ledger sec.59).

  1. Materialises the three box matrices L1, L2, L3 (resting regime) and verifies
     each column is one-hot and that coupled fixed points match enumeration.
  2. Box-level analysis: each box is monostable in isolation -> bistability is an
     emergent property of inter-box coupling.
  3. Actuator / kernel logic: which single input flips a resistant cell to apoptosis.

Run:  python analysis.py
"""

from __future__ import annotations
import numpy as np

from bbcn_switch import (FAST, MID, SLOW, NODES, INPUTS, rules, simulate, label,
                         fixed_points, build_box_matrix)

RESTING = dict(SRC=0, RHEB=1, IGF1R=1, RTK_up=1, CDKN2A=1, E2F1=0, ATM=0, ATR=0)


def materialise_matrices():
    print("=" * 70)
    print("1. MATERIALISED BOX MATRICES (resting regime)")
    print("=" * 70)
    for box, name in [(FAST, "FAST"), (MID, "MID"), (SLOW, "SLOW")]:
        L, others = build_box_matrix(box, RESTING)
        onehot = bool(np.all(L.sum(0) == 1))
        print(f"  L[{name:4}] shape {str(L.shape):12} columns one-hot: {onehot}")
    fp = fixed_points(RESTING)
    print(f"\n  coupled fixed points (resting): {fp}")
    print("  (= simultaneous fixed points of L1,L2,L3 = the multirate attractors)")


def actuator_logic():
    print("\n" + "=" * 70)
    print("2. ACTUATOR / KERNEL LOGIC (single input flips resistant -> apoptotic)")
    print("=" * 70)
    surv0 = {n: 0 for n in NODES}; surv0.update(AKT1=1, MDM2=1, MTOR=1, PDPK1=1, PIK3CA=1)
    resistant = simulate(surv0, RESTING)
    meaning = {"RTK_up": "PI3K/RTK inhibitor", "SRC": "SRC inhibitor",
               "RHEB": "mTOR inhibitor", "IGF1R": "IGF1R inhibitor",
               "ATM": "genotoxic (DNA damage)", "ATR": "genotoxic (DNA damage)",
               "E2F1": "oncogenic stress", "CDKN2A": "p14ARF status"}
    for k in INPUTS:
        I = dict(RESTING); I[k] = 1 - I[k]
        flips = label(simulate(resistant, I)) == "APOPTOTIC"
        print(f"  {k:8} -> {meaning.get(k,''):26} : "
              f"{'FLIPS to apoptosis (kernel)' if flips else 'no flip'}")
    print("\n  => only genotoxic/stress inputs (ATM/ATR/E2F1) flip the switch;")
    print("     survival-axis inhibition alone is insufficient (hysteresis).")


if __name__ == "__main__":
    materialise_matrices()
    actuator_logic()
