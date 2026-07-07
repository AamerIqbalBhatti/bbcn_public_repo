#!/usr/bin/env python3
"""
run_chunked.py — resumable driver for generate_numbers, sized to the sandbox.

Single core, ~5 min per invocation. This processes the work in patient-chunks,
checkpoints progress to results/numbers/checkpoint/, persists the (x0-independent)
stabilize caches to disk so they stay warm across calls, and exits cleanly after a
wall-time budget. Re-invoke until '--status' shows all tasks done, then '--emit'.

  python run_chunked.py            # do as much as fits in the time budget, then exit
  python run_chunked.py --status   # show per-task progress
  python run_chunked.py --emit     # write numbers.tex/json/csv from completed tasks
"""
from __future__ import annotations
import os, sys, json, time, pickle, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
SETUP_A = os.path.join(HERE, "setup_a")
SETUP_B = os.path.join(HERE, "setup_b", "code")
CKDIR = os.path.join(HERE, "results", "numbers", "checkpoint")
OUTDIR = os.path.join(HERE, "results", "numbers")
UPLOADS = "/mnt/user-data/uploads"
os.makedirs(CKDIR, exist_ok=True)

CKPT = os.path.join(CKDIR, "progress.json")
CACHE = os.path.join(CKDIR, "stab_cache.pkl")
BUDGET = 200          # clean-exit budget; per-chunk commit keeps it safe anyway
CHUNK = 100           # smaller chunks -> tighter deadline adherence

A_BIN = {
    "TCGA":     os.path.join(SETUP_A, "data/binarized/tcga_brca_1082x135.csv"),
    "METABRIC": os.path.join(SETUP_A, "data/binarized/metabric_1980x135.csv"),
    "ISPY2":    os.path.join(SETUP_A, "data/binarized/ispy2_988x135.csv"),
}
B_RAW = {
    "TCGA":     dict(expr=f"{UPLOADS}/data_mrna_seq_v2_rsem.txt",
                     muts=f"{UPLOADS}/data_mutations.txt", gz=False, idtrim=15, idcols=2, hic=None),
    "METABRIC": dict(expr="/tmp/brca_metabric/data_mrna_illumina_microarray.txt",
                     muts="/tmp/brca_metabric/data_mutations.txt", gz=False, idtrim=0, idcols=2, hic=None),
    "ISPY2":    dict(expr=f"{UPLOADS}/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988_txt.gz",
                     muts=None, gz=True, idtrim=0, idcols=1, hic=0),
}
PHEN = ["Resistance_OFF", "Apoptosis_ON", "Proliferation_OFF"]
COHORTS = ("TCGA", "METABRIC", "ISPY2")

# ---- task registry: id -> spec ----
def build_tasks():
    T = {}
    for ch in COHORTS:
        for m in ("ranked", "stabilize"):
            T[f"A|static_iso|{m}|{ch}"] = dict(kind="A_static_iso", cohort=ch, method=m)
            T[f"A|static_seq|{m}|{ch}"] = dict(kind="A_static_seq", cohort=ch, method=m)
            T[f"A|cde|{m}|{ch}"]        = dict(kind="A_cde",        cohort=ch, method=m)
    for ch in COHORTS:
        T[f"B|routing|switch|{ch}"]      = dict(kind="B_routing", cohort=ch, method="switch", atomic=True)
        T[f"B|kernel|ranked|{ch}"]       = dict(kind="B_kernel",  cohort=ch, method="ranked", atomic=True)
        T[f"B|kernel|stabilize|{ch}"]    = dict(kind="B_kernel",  cohort=ch, method="stabilize", atomic=True)
    return T

def load_ckpt():
    if os.path.exists(CKPT):
        return json.load(open(CKPT))
    return {}

def save_ckpt(c):
    json.dump(c, open(CKPT, "w"), indent=1)

# ---- lazy module + data loaders ----
_A = {}
def _amods():
    if not _A:
        sys.path.insert(0, SETUP_A)
        from bbcn import controller as C, harness as H
        from bbcn.run_cde import run_patient_cde
        _A.update(C=C, H=H, cde=run_patient_cde)
        # restore caches from disk
        if os.path.exists(CACHE):
            d = pickle.load(open(CACHE, "rb"))
            H._STAB_CACHE.update(d.get("h", {}))
            C._STAB_KC.update(d.get("c", {}))
    return _A

def _save_caches():
    if _A:
        pickle.dump({"h": _A["H"]._STAB_CACHE, "c": _A["C"]._STAB_KC}, open(CACHE, "wb"))

_PTS = {}
def _patients(cohort):
    if cohort not in _PTS:
        import pandas as pd
        H = _amods()["H"]
        B = pd.read_csv(A_BIN[cohort])
        cols = [c for c in B.columns if c in set(H.ALL_NODES)]
        _PTS[cohort] = [{n: int(r[n]) for n in cols} for _, r in B.iterrows()]
    return _PTS[cohort]

# ---- chunk processors ----
def proc_A_static_iso(st, spec, deadline, commit):
    C = _amods()["C"]; pts = _patients(spec["cohort"]); m = spec["method"]
    a = st["accum"]; i = st["done_n"]
    while i < len(pts) and time.time() < deadline:
        for p in pts[i:i+CHUNK]:
            o = C.run_patient(p, mode="isolated", family="ideal", kernel_method=m)
            for ph in PHEN: a[ph] += int(o[ph])
        i += min(CHUNK, len(pts)-i); st["done_n"] = i; commit()
        if time.time() >= deadline: break
    st["total"] = len(pts)

def proc_A_static_seq(st, spec, deadline, commit):
    C = _amods()["C"]; pts = _patients(spec["cohort"]); m = spec["method"]
    a = st["accum"]; i = st["done_n"]
    while i < len(pts) and time.time() < deadline:
        for p in pts[i:i+CHUNK]:
            o = C.run_patient(p, mode="sequenced", family="compatible", kernel_method=m)
            for ph in PHEN: a[ph] += int(o[ph])
            a["all_three"] += int(o.get("all_three", False))
        i += min(CHUNK, len(pts)-i); st["done_n"] = i; commit()
        if time.time() >= deadline: break
    st["total"] = len(pts)

def proc_A_cde(st, spec, deadline, commit):
    cde = _amods()["cde"]; pts = _patients(spec["cohort"]); m = spec["method"]
    a = st["accum"]; i = st["done_n"]
    keymap = {"Resistance": "Resistance_OFF", "Apoptosis": "Apoptosis_ON",
              "Proliferation": "Proliferation_OFF", "Terminal": "Terminal"}
    while i < len(pts) and time.time() < deadline:
        for p in pts[i:i+CHUNK]:
            o = cde(p, kernel_method=m)
            for stg, info in o["stage_summary"].items():
                k = keymap.get(stg, stg)
                if k in a and info["passed"]: a[k] += 1
        i += min(CHUNK, len(pts)-i); st["done_n"] = i; commit()
        if time.time() >= deadline: break
    st["total"] = len(pts)

def proc_B_routing(st, spec, deadline, commit):
    sys.path.insert(0, SETUP_B)
    import cohort_pipeline as CP
    cfg = B_RAW[spec["cohort"]]
    s, r = CP.load_matrix(cfg["expr"], idcols=cfg["idcols"], gzipped=cfg["gz"], header_idcols=cfg["hic"])
    tp53 = pik = None
    if cfg["muts"]: tp53, pik = CP.load_mutations(cfg["muts"])
    trim = (lambda x: x[:cfg["idtrim"]]) if cfg["idtrim"] else (lambda x: x)
    res = CP.run_cohort(spec["cohort"], s, r, tp53, pik, idtrim=trim, method="mad")
    st["accum"] = {k: res[k] for k in ("apop_pct","surv_pct","mixed_pct","surv_uncoupled_pct","n")}
    st["accum"]["U_pct"] = round(100*res["U"])
    if res["resist_surv"] is not None:
        st["accum"]["resist_surv_z"] = round(res["resist_surv"], 2)
        st["accum"]["resist_apop_z"] = round(res["resist_apop"], 2)
    st["done_n"] = st["total"] = res["n"]

def proc_B_kernel(st, spec, deadline, commit):
    sys.path.insert(0, SETUP_B)
    import forward_stab_kernel_design as KD
    cfg = B_RAW[spec["cohort"]]
    s, r = KD.load_matrix(cfg["expr"], idcols=cfg["idcols"], gz=cfg["gz"], header_idcols=cfg["hic"])
    tp53 = pik = None
    if cfg["muts"]: tp53, pik = KD.load_muts(cfg["muts"])
    trim = (lambda x: x[:cfg["idtrim"]]) if cfg["idtrim"] else (lambda x: x)
    kr = KD.run(spec["cohort"], s, r, tp53, pik, idtrim=trim, N=None, method=spec["method"])
    top = max(kr["node_use"], key=kr["node_use"].get) if kr["node_use"] else "none"
    st["accum"] = dict(flip_pct=kr["flip_pct"], ker_designed=kr["ker_designed"],
                       ker_success=kr["ker_success"], top_node=top)
    st["done_n"] = st["total"] = kr["N"]

PROC = {"A_static_iso": proc_A_static_iso, "A_static_seq": proc_A_static_seq,
        "A_cde": proc_A_cde, "B_routing": proc_B_routing, "B_kernel": proc_B_kernel}

def init_state(spec):
    if spec["kind"] == "A_static_iso":   ac = {p: 0 for p in PHEN}
    elif spec["kind"] == "A_static_seq": ac = {p: 0 for p in PHEN}; ac["all_three"] = 0
    elif spec["kind"] == "A_cde":        ac = {p: 0 for p in PHEN + ["Terminal"]}
    else:                                ac = {}
    return dict(done_n=0, total=None, accum=ac, complete=False)

def run_once():
    tasks = build_tasks(); ck = load_ckpt()
    for tid in tasks:
        ck.setdefault(tid, init_state(tasks[tid]))
    deadline = time.time() + BUDGET
    did = []
    for tid, spec in tasks.items():
        if ck[tid]["complete"]:
            continue
        if time.time() >= deadline:
            break
        commit = lambda: (save_ckpt(ck), _save_caches())
        PROC[spec["kind"]](ck[tid], spec, deadline, commit)
        if ck[tid]["total"] is not None and ck[tid]["done_n"] >= ck[tid]["total"]:
            ck[tid]["complete"] = True
        did.append(f"{tid} {ck[tid]['done_n']}/{ck[tid]['total']}{' DONE' if ck[tid]['complete'] else ''}")
        save_ckpt(ck); _save_caches()
    save_ckpt(ck); _save_caches()
    done = sum(1 for t in ck.values() if t["complete"])
    print(f"this call processed: {len(did)} slice(s)")
    for d in did: print("   ", d)
    print(f"OVERALL: {done}/{len(tasks)} tasks complete")
    return done == len(tasks)

def status():
    tasks = build_tasks(); ck = load_ckpt()
    for tid in tasks:
        s = ck.get(tid, {})
        print(f"  {'OK ' if s.get('complete') else '.. '}{tid:34s} {s.get('done_n',0)}/{s.get('total','?')}")
    done = sum(1 for t in ck.values() if t.get("complete"))
    print(f"OVERALL: {done}/{len(tasks)} complete")

def _macro(setup, path, method, cohort, q):
    def camel(s):
        s = str(s)
        for d, w in zip("0123456789", ["Zero","One","Two","Three","Four","Five","Six","Seven","Eight","Nine"]):
            s = s.replace(d, " " + w + " ")   # digits -> words (LaTeX names are letters-only)
        return "".join(w.capitalize() for w in s.replace("-"," ").replace("_"," ").split())
    return "BBCN"+camel(setup)+camel(path)+camel(method)+camel(cohort)+camel(q)

def emit():
    tasks = build_tasks(); ck = load_ckpt(); ROWS = []
    def rec(setup, path, method, cohort, q, v):
        ROWS.append(dict(setup=setup, path=path, method=method, cohort=cohort, quantity=q, value=v))
    for tid, spec in tasks.items():
        s = ck.get(tid, {})
        if not s.get("complete"):
            print(f"  (skip incomplete: {tid})"); continue
        setup, path, method, cohort = tid.split("|")
        a = s["accum"]; n = s["total"]
        if spec["kind"] == "A_static_iso":
            for ph in PHEN: rec("A","static",method,cohort,f"iso_{ph}", round(100*a[ph]/n))
            rec("A","static",method,cohort,"N", n)
        elif spec["kind"] == "A_static_seq":
            for ph in PHEN: rec("A","static",method,cohort,f"seq_{ph}", round(100*a[ph]/n))
            rec("A","static",method,cohort,"seq_all_three", round(100*a["all_three"]/n))
        elif spec["kind"] == "A_cde":
            for st_ in PHEN+["Terminal"]: rec("A","cde",method,cohort,f"stagepass_{st_}", round(100*a[st_]/n))
        elif spec["kind"] == "B_routing":
            for q,v in a.items(): rec("B","routing","switch",cohort,q, v)
        elif spec["kind"] == "B_kernel":
            for q,v in a.items(): rec("B","kernel",method,cohort,q, v)
    js = {_macro(r["setup"],r["path"],r["method"],r["cohort"],r["quantity"]): r["value"] for r in ROWS}
    os.makedirs(OUTDIR, exist_ok=True)
    json.dump(js, open(os.path.join(OUTDIR,"numbers.json"),"w"), indent=2, sort_keys=True)
    import csv
    with open(os.path.join(OUTDIR,"all_numbers.csv"),"w",newline="") as f:
        w = csv.DictWriter(f, fieldnames=["setup","path","method","cohort","quantity","value","macro"]); w.writeheader()
        for r in ROWS: w.writerow({**r, "macro": _macro(r["setup"],r["path"],r["method"],r["cohort"],r["quantity"])})
    with open(os.path.join(OUTDIR,"numbers.tex"),"w") as f:
        f.write("% AUTO-GENERATED by run_chunked.py --emit — DO NOT EDIT BY HAND.\n")
        for k in sorted(js): f.write("\\newcommand{\\%s}{%s}\n" % (k, js[k]))
    print(f"emitted {len(ROWS)} numbers -> {OUTDIR}/{{numbers.tex,numbers.json,all_numbers.csv}}")

def run_all():
    """Process every task to completion in one process (no per-call budget).
    Use on an unconstrained machine; the sandbox uses repeated run_once() instead."""
    tasks = build_tasks(); ck = load_ckpt()
    for tid in tasks:
        ck.setdefault(tid, init_state(tasks[tid]))
    far = time.time() + 10**9
    for tid, spec in tasks.items():
        if ck[tid]["complete"]:
            continue
        commit = lambda: (save_ckpt(ck), _save_caches())
        PROC[spec["kind"]](ck[tid], spec, far, commit)
        if ck[tid]["total"] is not None and ck[tid]["done_n"] >= ck[tid]["total"]:
            ck[tid]["complete"] = True
        save_ckpt(ck); _save_caches()
        print(f"  done {tid} ({ck[tid]['done_n']}/{ck[tid]['total']})")
    emit()
    return True

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--emit", action="store_true")
    ap.add_argument("--all", action="store_true", help="run every task to completion in one process")
    a = ap.parse_args()
    if a.status: status()
    elif a.emit: emit()
    elif a.all: run_all()
    else: run_once()
