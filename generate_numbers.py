#!/usr/bin/env python3
r"""
generate_numbers.py — the SINGLE source of every number in the BBCN preprint.

Runs all configurations once and writes three synchronised artefacts:

  results/numbers/numbers.tex   LaTeX \newcommand macros; the preprint \input's this
                                file, so no number is ever hand-typed in the paper.
  results/numbers/numbers.json  the same values, machine-readable.
  results/numbers/all_numbers.csv  one universal table; every row tagged with the
                                exact routine (setup, path, method, cohort, quantity)
                                that produced it.

Configurations
  Setup A : {static controller, dynamic CDE} x {ranked, stabilize}
            on TCGA / METABRIC / I-SPY2, isolated(ideal) + sequenced(compatible).
  Setup B : routing once per cohort (kernel-method-independent) + per-box kernel
            flip rate for ranked vs stabilize, on TCGA / METABRIC / I-SPY2.

Usage
  python generate_numbers.py                # full cohorts
  python generate_numbers.py --fast 60      # reduced-N dry run (wiring check)
  python generate_numbers.py --setup a      # only Setup A   (a|b|all)
"""
from __future__ import annotations
import os, sys, csv, json, time, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
SETUP_A = os.path.join(HERE, "setup_a")
SETUP_B = os.path.join(HERE, "setup_b", "code")
OUTDIR = os.path.join(HERE, "results", "numbers")
UPLOADS = "/mnt/user-data/uploads"

# data locations (raw matrices for Setup B; binarized CSVs for Setup A)
A_BIN = {
    "TCGA":     os.path.join(SETUP_A, "data/binarized/tcga_brca_1082x135.csv"),
    "METABRIC": os.path.join(SETUP_A, "data/binarized/metabric_1980x135.csv"),
    "ISPY2":    os.path.join(SETUP_A, "data/binarized/ispy2_988x135.csv"),
}
A_BIN = {k: v for k, v in A_BIN.items() if os.path.exists(v)}   # auto-detect: run only cohorts present
print("[cohorts present] " + (", ".join(A_BIN) or "NONE"))
B_RAW = {
    "TCGA":     dict(expr=f"{UPLOADS}/data_mrna_seq_v2_rsem.txt",
                     muts=f"{UPLOADS}/data_mutations.txt", gz=False, idtrim=15,
                     idcols=2, header_idcols=None),
    "METABRIC": dict(expr="/tmp/brca_metabric/data_mrna_illumina_microarray.txt",
                     muts="/tmp/brca_metabric/data_mutations.txt", gz=False, idtrim=0,
                     idcols=2, header_idcols=None),
    "ISPY2":    dict(expr=f"{UPLOADS}/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988_txt.gz",
                     muts=None, gz=True, idtrim=0, idcols=1, header_idcols=0),
}

# universal table rows: (setup, path, method, cohort, quantity, value)
ROWS = []
def rec(setup, path, method, cohort, quantity, value):
    ROWS.append(dict(setup=setup, path=path, method=method, cohort=cohort,
                     quantity=quantity, value=value))


# ---------------------------------------------------------------- Setup A
def run_setup_a(fast=None):
    sys.path.insert(0, SETUP_A)
    import pandas as pd
    from bbcn import controller as C, harness as H
    from bbcn.run_cde import run_patient_cde

    def load(cohort):
        B = pd.read_csv(A_BIN[cohort])
        cols = [c for c in B.columns if c in set(H.ALL_NODES)]
        pts = [{n: int(r[n]) for n in cols} for _, r in B.iterrows()]
        return pts[:fast] if fast else pts

    PHEN = ["Resistance_OFF", "Apoptosis_ON", "Proliferation_OFF"]
    for cohort in ("TCGA", "METABRIC", "ISPY2"):
        pts = load(cohort)
        n = len(pts)
        for method in ("ranked", "stabilize"):
            # ---- static controller path: isolated/ideal + sequenced/compatible ----
            iso = {p: 0 for p in PHEN}
            seq = {p: 0 for p in PHEN}; allthree = 0
            t = time.time()
            for init in pts:
                o = C.run_patient(init, mode="isolated", family="ideal", kernel_method=method)
                for p in PHEN: iso[p] += int(o[p])
            for init in pts:
                o = C.run_patient(init, mode="sequenced", family="compatible", kernel_method=method)
                for p in PHEN: seq[p] += int(o[p])
                allthree += int(o.get("all_three", False))
            for p in PHEN:
                rec("A", "static", method, cohort, f"iso_{p}", round(100*iso[p]/n))
                rec("A", "static", method, cohort, f"seq_{p}", round(100*seq[p]/n))
            rec("A", "static", method, cohort, "seq_all_three", round(100*allthree/n))
            rec("A", "static", method, cohort, "N", n)
            print(f"  [A/static/{method}/{cohort}] N={n} "
                  f"iso(R{round(100*iso['Resistance_OFF']/n)}/A{round(100*iso['Apoptosis_ON']/n)}/"
                  f"P{round(100*iso['Proliferation_OFF']/n)}) all3={round(100*allthree/n)}% ({time.time()-t:.0f}s)")

            # ---- dynamic CDE path: sequenced four-stage, held-in-window pass ----
            H._STAB_CACHE.clear()
            cde = {st: 0 for st in ("Resistance_OFF", "Apoptosis_ON", "Proliferation_OFF", "Terminal")}
            t = time.time()
            for init in pts:
                o = run_patient_cde(init, kernel_method=method)
                for st, info in o["stage_summary"].items():
                    key = {"Resistance": "Resistance_OFF", "Apoptosis": "Apoptosis_ON",
                           "Proliferation": "Proliferation_OFF", "Terminal": "Terminal"}.get(st, st)
                    if key in cde and info["passed"]: cde[key] += 1
            for st, v in cde.items():
                rec("A", "cde", method, cohort, f"stagepass_{st}", round(100*v/n))
            print(f"  [A/cde/{method}/{cohort}] N={n} "
                  f"R{round(100*cde['Resistance_OFF']/n)}/A{round(100*cde['Apoptosis_ON']/n)}/"
                  f"P{round(100*cde['Proliferation_OFF']/n)}/T{round(100*cde['Terminal']/n)} ({time.time()-t:.0f}s)")
    sys.path.remove(SETUP_A)


# ---------------------------------------------------------------- Setup B
def run_setup_b(fast=None):
    sys.path.insert(0, SETUP_B)
    import cohort_pipeline as CP
    import forward_stab_kernel_design as KD

    for cohort in ("TCGA", "METABRIC", "ISPY2"):
        cfg = B_RAW[cohort]
        # routing (method-independent), via cohort_pipeline
        s, r = CP.load_matrix(cfg["expr"], idcols=cfg["idcols"], gzipped=cfg["gz"],
                              header_idcols=cfg["header_idcols"])
        tp53 = pik = None
        if cfg["muts"]:
            tp53, pik = CP.load_mutations(cfg["muts"])
        trim = (lambda x: x[:cfg["idtrim"]]) if cfg["idtrim"] else (lambda x: x)
        if fast:
            s = s[:fast]; r = {g: v[:fast] for g, v in r.items()}
        route = CP.run_cohort(cohort, s, r, tp53, pik, idtrim=trim, method="mad")
        for q in ("apop_pct", "surv_pct", "mixed_pct", "surv_uncoupled_pct", "n"):
            rec("B", "routing", "switch", cohort, q, route[q])
        rec("B", "routing", "switch", cohort, "U_pct", round(100*route["U"]))
        if route["resist_surv"] is not None:
            rec("B", "routing", "switch", cohort, "resist_surv_z", round(route["resist_surv"], 2))
            rec("B", "routing", "switch", cohort, "resist_apop_z", round(route["resist_apop"], 2))

        # per-box kernel flip rate: ranked vs stabilize, via forward_stab_kernel_design
        s2, r2 = KD.load_matrix(cfg["expr"], idcols=cfg["idcols"], gz=cfg["gz"],
                                header_idcols=cfg["header_idcols"])
        tp53m = pikm = None
        if cfg["muts"]:
            tp53m, pikm = KD.load_muts(cfg["muts"])
        Ncap = fast if fast else None
        for method in ("ranked", "stabilize"):
            kr = KD.run(cohort, s2, r2, tp53m, pikm, idtrim=trim, N=Ncap, method=method)
            rec("B", "kernel", method, cohort, "flip_pct", kr["flip_pct"])
            rec("B", "kernel", method, cohort, "ker_designed", kr["ker_designed"])
            rec("B", "kernel", method, cohort, "ker_success", kr["ker_success"])
            top = max(kr["node_use"], key=kr["node_use"].get) if kr["node_use"] else "none"
            rec("B", "kernel", method, cohort, "top_node", top)
    sys.path.remove(SETUP_B)


# ---------------------------------------------------------------- emit
def _macro_name(row):
    # \BBCN<Setup><Path><Method><Cohort><Quantity>  (LaTeX-safe: letters only)
    def camel(s):
        return "".join(w.capitalize() for w in str(s).replace("-", " ").replace("_", " ").split())
    return "BBCN" + camel(row["setup"]) + camel(row["path"]) + camel(row["method"]) \
           + camel(row["cohort"]) + camel(row["quantity"])

def emit():
    os.makedirs(OUTDIR, exist_ok=True)
    # numbers.json
    js = {}
    for row in ROWS:
        js[_macro_name(row)] = row["value"]
    with open(os.path.join(OUTDIR, "numbers.json"), "w") as f:
        json.dump(js, f, indent=2, sort_keys=True)
    # all_numbers.csv (universal table, routine-tagged)
    with open(os.path.join(OUTDIR, "all_numbers.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["setup", "path", "method", "cohort", "quantity", "value", "macro"])
        w.writeheader()
        for row in ROWS:
            w.writerow({**row, "macro": _macro_name(row)})
    # numbers.tex (\newcommand macros)
    with open(os.path.join(OUTDIR, "numbers.tex"), "w") as f:
        f.write("% AUTO-GENERATED by generate_numbers.py — DO NOT EDIT BY HAND.\n")
        f.write("% The preprint \\input's this file so paper == code by construction.\n")
        for name in sorted(js):
            v = js[name]
            f.write("\\newcommand{\\%s}{%s}\n" % (name, v))
    print(f"\nwrote {len(ROWS)} numbers -> {OUTDIR}/{{numbers.tex, numbers.json, all_numbers.csv}}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", type=int, default=None, help="reduced N per cohort (wiring dry run)")
    ap.add_argument("--setup", choices=["a", "b", "all"], default="all")
    a = ap.parse_args()
    t0 = time.time()
    if a.setup in ("a", "all"):
        print("=== SETUP A ===")
        run_setup_a(fast=a.fast)
    if a.setup in ("b", "all"):
        print("=== SETUP B ===")
        run_setup_b(fast=a.fast)
    emit()
    print(f"total {time.time()-t0:.0f}s")
