#!/usr/bin/env python3
"""Per-cohort three-arm comparison: static continuous-hold vs pulsed weekly (1/2/3 days on).
Saves JSON so the cohorts can be run one at a time. Usage: pulsed_compare.py COHORT N"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from bbcn import controller as C, pulsed as P, harness as H

BIN = {"TCGA": "tcga_brca_1082x135.csv", "METABRIC": "metabric_1980x135.csv", "ISPY2": "ispy2_988x135.csv"}
PHEN = ['Apoptosis_ON', 'Proliferation_OFF', 'Resistance_OFF']
cohort = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
HERE = os.path.dirname(os.path.abspath(__file__))

B = pd.read_csv(os.path.join(HERE, "data", "binarized", BIN[cohort]))
cols = [c for c in B.columns if c in set(H.ALL_NODES)]
pts = [{c: int(r[c]) for c in cols} for _, r in B.head(n).iterrows()]
N = len(pts); pct = lambda k: round(100 * k / N)
res = {"cohort": cohort, "N": N}

t0 = time.time()
# Arm 1: static continuous-hold (stabilize). Also measure durable fraction for apoptosis.
iso = [C.run_patient(dict(p), mode="isolated", family="ideal", kernel_method="stabilize") for p in pts]
seq = [C.run_patient(dict(p), mode="sequenced", family="compatible", kernel_method="stabilize") for p in pts]
ach = dur = 0
for p in pts:
    bus0 = {x: 0 for x in H.ALL_NODES}
    for x in C._NODES: bus0[x] = int(p.get(x, 0))
    bus0 = C.evolve(bus0, {})
    a, fb, held = C.control_phenotype(dict(bus0), 'Apoptosis_ON', kernel_method='stabilize')
    if a:
        ach += 1
        if C.is_fixed_point(fb, {}): dur += 1
res["static_hold"] = {ph.split('_')[0].lower(): pct(sum(o[ph] for o in iso)) for ph in PHEN}
res["static_hold"]["all3"] = pct(sum(o['all_three'] for o in seq))
res["static_hold"]["apoptosis_durable_pct_of_achieved"] = (round(100 * dur / ach) if ach else 0)
print(f"[{cohort}] static done {time.time()-t0:.0f}s  apop={res['static_hold']['apoptosis']} durable={res['static_hold']['apoptosis_durable_pct_of_achieved']}")

# Arm 2: pulsed weekly, on_days sweep
res["pulsed"] = {}
for od in (1, 2, 3):
    t1 = time.time()
    pi = [P.run_patient_pulsed(dict(p), mode="isolated", family="ideal", on_days=od, kernel_method="stabilize", druggable=True) for p in pts]
    ps = [P.run_patient_pulsed(dict(p), mode="sequenced", family="compatible", on_days=od, kernel_method="stabilize", druggable=True) for p in pts]
    d = {ph.split('_')[0].lower(): pct(sum(o[ph] for o in pi)) for ph in PHEN}
    d["all3"] = pct(sum(o['all_three'] for o in ps))
    res["pulsed"][f"on{od}"] = d
    print(f"[{cohort}] pulsed on={od} done {time.time()-t1:.0f}s  apop={d['apoptosis']} prolif={d['proliferation']} all3={d['all3']}")

json.dump(res, open(os.path.join(HERE, f"pulsed_cmp_{cohort}.json"), "w"), indent=1)
print(f"[{cohort}] saved.")
