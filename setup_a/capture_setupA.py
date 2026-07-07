#!/usr/bin/env python3
"""
Capture Setup A kernel selections by observing controller._kernel_for (no logic change).
Static isolated / ideal pass, both methods, all three cohorts. Resumable: writes
tallies to capture_setupA.json after each (cohort, method) pair.
"""
import sys, os, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pandas as pd
from bbcn import controller as C, harness as H

A_BIN = {
    "TCGA": os.path.join(HERE, "data/binarized/tcga_brca_1082x135.csv"),
    "METABRIC": os.path.join(HERE, "data/binarized/metabric_1980x135.csv"),
    "ISPY2": os.path.join(HERE, "data/binarized/ispy2_988x135.csv"),
}
OUT = os.path.join(HERE, "capture_setupA.json")
BUDGET = 240.0

orig = C._kernel_for
STATE = {"sel": set()}
def wrapped(bus, pw, ph, method="stabilize"):
    k = orig(bus, pw, ph, method=method)
    for node in k:
        STATE["sel"].add((ph, pw, node))
    return k
C._kernel_for = wrapped

def load(cohort):
    B = pd.read_csv(A_BIN[cohort])
    cols = [c for c in B.columns if c in set(H.ALL_NODES)]
    return [{n: int(r[n]) for n in cols} for _, r in B.iterrows()]

tall = json.load(open(OUT)) if os.path.exists(OUT) else {"node_pat": {}, "pw_node": {}, "N": {}, "done": []}
start = time.time()
for cohort in ("TCGA", "METABRIC", "ISPY2"):
    pts = None
    for method in ("ranked", "stabilize"):
        key = f"{cohort}|{method}"
        if key in tall["done"]:
            continue
        if time.time() - start > BUDGET:
            print("BUDGET reached; re-invoke to continue."); json.dump(tall, open(OUT, "w")); sys.exit(0)
        if pts is None:
            pts = load(cohort)
        tall["N"][cohort] = len(pts)
        np_d = tall["node_pat"].setdefault(key, {})
        pn_d = tall["pw_node"].setdefault(key, {})
        t = time.time()
        for init in pts:
            STATE["sel"] = set()
            C.run_patient(init, mode="isolated", family="ideal", kernel_method=method)
            nodes = set(n for _, _, n in STATE["sel"])
            for n in nodes:
                np_d[n] = np_d.get(n, 0) + 1
            for ph, pw, n in STATE["sel"]:
                kk = f"{pw}|{n}"
                pn_d[kk] = pn_d.get(kk, 0) + 1
        tall["done"].append(key)
        json.dump(tall, open(OUT, "w"))
        print(f"  done {key}: N={len(pts)} in {time.time()-t:.0f}s ; nodes={len(np_d)}")
print("ALL DONE:", tall["done"])
