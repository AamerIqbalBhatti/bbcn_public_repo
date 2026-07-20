#!/usr/bin/env python3
"""Setup B payoff: does a withdrawable pulse latch durably apoptotic via hysteresis?
For each SURVIVAL-routed (resistant) patient: design the apoptosis kernel, clamp it for
P ticks, RELEASE, simulate the switch free, and check it stays APOPTOTIC. Setup A gives 0.
Usage: payoff_pulse.py COHORT [N]"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cohort_pipeline import (load_matrix, NEED, SEED_SIG, AKT_ACT, FOXO_NUC,
                             DAMAGE_SIG, DAMAGE_PCT)
from bbcn_switch import FAST, SLOW, NODES, rules, simulate, label
from forward_stab_kernel_design import design_patient_kernel, APOP_TARGET

TSV = {"TCGA": "tcga_brca_switch_inputs.tsv", "METABRIC": "metabric_switch_inputs.tsv",
       "ISPY2": "ispy2_switch_inputs.tsv"}
cohort = sys.argv[1]; LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 250
HERE = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(HERE, "..", "data", "samples", TSV[cohort])

samples, rows = load_matrix(path, idcols=1)
w = min([len(v) for v in rows.values()] + [len(samples)]); samples = samples[:w]
for g in list(rows): rows[g] = rows[g][:w]
n = min(w, LIMIT)

def z(g):
    if g not in rows: return np.zeros(w)
    v = rows[g]; c = np.nanmedian(v); mad = np.nanmedian(np.abs(v - c)); s = 1.4826 * mad or 1.0
    r = (v - c) / s; r[np.isnan(r)] = 0.0; return r
Z = {g: z(g) for g in NEED}
SEED = {nd: (np.mean([Z[g] for g in SEED_SIG[nd] if g in rows], axis=0) > 0).astype(int)
        if any(g in rows for g in SEED_SIG[nd]) else np.zeros(w, int) for nd in SEED_SIG}
akt = np.mean([Z[g] for g in AKT_ACT if g in rows], axis=0)
fox = np.mean([Z[g] for g in FOXO_NUC if g in rows], axis=0)
U = ((akt > 0) & (fox > 0)).astype(int)
dmgv = np.mean([Z[g] for g in DAMAGE_SIG if g in rows], axis=0); dcut = np.quantile(dmgv, DAMAGE_PCT)

def patient(i):
    rtk = int((Z["GRB2"][i] > 0) or (Z["IRS1"][i] > 0)) if "GRB2" in rows else 0
    I = dict(SRC=int(Z["SRC"][i] > 0) if "SRC" in rows else 0,
             RHEB=int(Z["RHEB"][i] > 0) if "RHEB" in rows else 0,
             IGF1R=int(Z["IGF1R"][i] > 0) if "IGF1R" in rows else 0, RTK_up=rtk,
             CDKN2A=int(Z["CDKN2A"][i] > 0) if "CDKN2A" in rows else 0,
             E2F1=int(Z["E2F1"][i] > 0) if "E2F1" in rows else 0,
             ATM=int(dmgv[i] > dcut), ATR=int(dmgv[i] > dcut))
    x0 = {nd: int(SEED[nd][i]) for nd in NODES}
    return x0, I, int(U[i])

def pulse_then_free(start, I, U, clamp, P, T=300, rate_mid=5, rate_slow=25):
    s = dict(start)
    for t in range(1, T + 1):
        snap = dict(s); nf = rules(snap, I, U)
        for nd in FAST: s[nd] = nf[nd]
        if t % rate_mid == 0: s["MDM2"] = rules(snap, I, U)["MDM2"]
        if t % rate_slow == 0:
            for nd in SLOW: s[nd] = rules(snap, I, U)[nd]
        if t <= P:
            for nd, v in clamp.items(): s[nd] = v     # hold clamp only during the pulse
    return label(s)

Ps = [0, 10, 25, 50, 100, 300]
resistant = held_feasible = 0
durable = {P: 0 for P in Ps}
t0 = time.time()
for i in range(n):
    x0, I, u = patient(i)
    if label(simulate(x0, I, u)) == "APOPTOTIC":
        continue                                       # already apoptotic, not resistant
    resistant += 1
    res, kernels = design_patient_kernel(x0, I, u)
    clamp = {nd: APOP_TARGET[nd] for bx in kernels for nd in kernels[bx]}
    if not clamp:
        continue
    if res == "APOPTOTIC":
        held_feasible += 1                              # reaches apoptosis under continuous hold
    settled = simulate(x0, I, u)
    for P in Ps:
        if pulse_then_free(settled, I, u, clamp, P) == "APOPTOTIC":
            durable[P] += 1

out = {"cohort": cohort, "N": n, "resistant": resistant,
       "held_feasible_pct_of_resistant": (round(100 * held_feasible / resistant) if resistant else 0),
       "durable_after_pulse_pct_of_resistant": {P: (round(100 * durable[P] / resistant) if resistant else 0) for P in Ps}}
json.dump(out, open(os.path.join(HERE, f"payoff_{cohort}.json"), "w"), indent=1)
print(f"[{cohort}] N={n} resistant(SURV-routed)={resistant}  ({time.time()-t0:.0f}s)")
print(f"  feasible under continuous hold: {out['held_feasible_pct_of_resistant']}% of resistant")
print(f"  DURABLE after withdrawable pulse, by pulse length P (ticks):")
print("   P:       " + "  ".join(f"{P:>4}" for P in Ps))
print("   %:       " + "  ".join(f"{out['durable_after_pulse_pct_of_resistant'][P]:>4}" for P in Ps))
print("  (Setup A durable apoptosis = 0% at every duration, for comparison)")
