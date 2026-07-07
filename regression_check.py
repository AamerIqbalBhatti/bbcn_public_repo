"""
regression_check.py  -- the Phase-0 safety net.

Runs the two load-bearing experiments (E40 durability, E41 staged controller) on a
FIXED small sample (the first N patients per cohort, deterministic), prints the numbers,
and compares them to a saved baseline. On the SAME code the output is deterministic, so
any difference after an edit means that edit changed behaviour -- which is exactly what
we want to catch, one spot at a time.

EXECUTION PATTERN: functions are DEFINED in the top half; the RUN SECTION at the bottom
(from "t0 = time.time()") actually CALLS run_E40() and run_E41(), prints the tables, and
compares to the baseline. Running `python regression_check.py` fires that bottom section.

Usage (from repo root):
    python regression_check.py                 # run + compare to baseline
    python regression_check.py --sample 40     # bigger sample
    python regression_check.py --update-baseline   # capture current numbers AS the baseline

Baseline file: docs/regression_baseline.json
Full-N acceptance bar (for reference, not re-run here):
    resistance==U 99.98% | commit 92-95% | strict 8-14% | durable FP 15-17%
    | apoptosis 72-78% | proliferation 84-87%
"""
import sys, os, json, time, subprocess, re
sys.path.insert(0, 'setup_a'); sys.path.insert(0, 'setup_b/code')

SAMPLE = 30
if '--sample' in sys.argv: SAMPLE = int(sys.argv[sys.argv.index('--sample') + 1])
UPDATE = '--update-baseline' in sys.argv
BASELINE = 'docs/regression_baseline.json'
COH = {'TCGA': 'setup_a/data/binarized/tcga_brca_1082x135.csv',
       'METABRIC': 'setup_a/data/binarized/metabric_1980x135.csv',
       'ISPY2': 'setup_a/data/binarized/ispy2_988x135.csv'}
COH = {k: v for k, v in COH.items() if os.path.exists(v)}   # auto-detect: run only cohorts present

# ---------- E40 durability (verbatim logic from run_durability_full.py) ----------
import pandas as pd
from bbcn import honest_bbcn as HB
import forward_stab_kernel_design as FSK
from bbcn_switch import FAST, MID, SLOW, rules as sw_rules, label as sw_label
NODES = HB.NODES; lags = HB.make_lags(); APOP = FSK.APOP_TARGET
SW_STATE = ['AKT1','PTEN','MTOR','PDPK1','PIK3CA','MDM2','TP53','FOXO3']

def seed(b):
    U = HB.patient_clamp_off(b); x0 = {nd: int(b.get(nd,0)) for nd in SW_STATE}
    x0['PHLPP'] = int(x0['FOXO3'] and not U)
    I = dict(SRC=int(b.get('SRC',0)), RHEB=int(b.get('RHEB',0)), IGF1R=int(b.get('IGF1R',0)),
             RTK_up=int(b.get('GRB2',0) or b.get('IRS1',0)), CDKN2A=int(b.get('CDKN2A',0)),
             E2F1=int(b.get('E2F1',0)), ATM=int(b.get('ATM',0)), ATR=int(b.get('ATR',0)))
    return x0, I, U

def sw_dyn(x0, I, U, clamp, hold=120, total=200):
    s = dict(x0)
    for nd, v in clamp.items(): s[nd] = v
    tail = []
    for t in range(1, total+1):
        nf = sw_rules(dict(s), I, U)
        for nd in FAST: s[nd] = nf[nd]
        if t % 5 == 0: s['MDM2'] = nf['MDM2']
        if t % 25 == 0:
            for nd in SLOW: s[nd] = nf[nd]
        if t <= hold:
            for nd, v in clamp.items(): s[nd] = v
        tail.append(sw_label(s))
    return sum(l == 'APOPTOTIC' for l in tail[-40:]) / 40 > 0.5

def A_dyn(b, U, clamp, T1=120, T2=180, theta=10, cmax=20):
    bus = {n: int(b.get(n,0)) for n in NODES}
    for k, v in clamp.items(): bus[k] = v
    hist = [bus]
    for t in range(T1+T2):
        nb = HB.step_honest(HB._view(hist, t, lags), U, cascade=True)   # Step 16: real cascade, no brake-zeroing
        if t < T1:
            for k, v in clamp.items(): nb[k] = v
        hist.append(nb)
    tail = hist[-40:]
    strict = sum(s['CASP3'] == 1 and s['AKT1'] == 0 for s in tail) / 40 > 0.5
    commitd = sum(s['CASP3'] == 1 for s in tail) / 40 > 0.5
    return strict, commitd

def run_E40(sample):
    out = {}
    for name, path in COH.items():
        df = pd.read_csv(path, index_col=0); n = min(sample, len(df))
        found = resist = strictc = commitc = swc = 0
        for i in range(n):
            b = {nd: int(df.iloc[i].get(nd, 0)) for nd in NODES}
            x0, I, U = seed(b)
            r, kernels = FSK.design_patient_kernel(x0, I, U, method='stabilize')
            if not kernels and r == 'APOPTOTIC': continue
            resist += 1
            if not kernels: continue
            found += 1
            clamp = {nd: APOP[nd] for bx in kernels for nd in kernels[bx]}
            st, cm = A_dyn(b, U, clamp); strictc += st; commitc += cm
            swc += sw_dyn(x0, I, U, clamp)
        out[name] = dict(n=n, resist=resist, found=found, strict=strictc, commit=commitc, switch=swc)
    return out

# ---------- E41 monolith (run INLINE; no subprocess, no file writing, no hard-coded paths) ----------
from bbcn import harness as _H
from bbcn import repaired_branch as _RB
_RB.apply()                      # bake the repair into harness.PATHWAYS (once)
from bbcn import run_cde as _CDE

def run_E41(sample):
    out = {}
    for name, path in COH.items():
        df = pd.read_csv(path, index_col=0); n = min(sample, len(df))
        apop = prolif = joint = dcommit = cnt = 0
        for i in range(n):
            b = {nd: int(df.iloc[i].get(nd, 0)) for nd in df.columns}
            U = _RB.patient_clamp_off(b)
            init = dict(b); init['CLAMP_OFF'] = U
            init['PHLPP'] = int(b.get('FOXO3', 0) and not U); init['COMMIT'] = 0
            r = _CDE.run_patient_cde(init, kernel_method='stabilize'); cnt += 1
            ss = r['stage_summary']; a = ss['Apoptosis_ON']['passed']
            apop += a; prolif += ss['Proliferation_OFF']['passed']; joint += r['terminal_pass']
            if a:
                fb = r['final_bus']; ffp = _CDE._is_free_fixed_point(fb)
                dcommit += int(ffp and fb.get('CASP3', 0) == 1)
        out[name] = dict(cnt=cnt, apop=apop, prolif=prolif, joint=joint, dcommit=dcommit)
    return out

# ---------- E42 clock durability: idealization (clock-continuous) vs clinical (pulsed, capped) ----------
import durability_clock as _DC
_APOP = FSK.APOP_TARGET

def run_E42_clock(sample):
    out = {}
    for name, path in COH.items():
        df = pd.read_csv(path, index_col=0); n = min(sample, len(df))
        found = cont = puls = capped = 0
        for i in range(n):
            b = {nd: int(df.iloc[i].get(nd, 0)) for nd in NODES}
            x0, I, U = seed(b)
            r, kernels = FSK.design_patient_kernel(x0, I, U, method='stabilize')
            if not kernels:
                continue
            found += 1
            clamp = {nd: _APOP[nd] for bx in kernels for nd in kernels[bx]}
            cont += _DC.A_dyn_clock(b, U, clamp, pulsed=False)               # idealized hold, on the clock
            puls += _DC.A_dyn_clock(b, U, clamp, pulsed=True)                # clinical 4-on/3-off
            capped += _DC.run_capped_on_clock(b, U, kernels, _APOP, pulsed=True)['durable']  # + caps/escalation
        out[name] = dict(found=found, cont=cont, puls=puls, capped=capped)
    return out

# ---------- run, print, compare ----------
def pct(x, d): return round(100*x/d) if d else 0

t0 = time.time()
print(f"Regression check  |  sample = {SAMPLE} patients/cohort  |  {time.strftime('%Y-%m-%d %H:%M')}")
print("=" * 78)
e40 = run_E40(SAMPLE)
print("E40 durability:")
for c, d in e40.items():
    print(f"  {c:9} resistant {pct(d['resist'],d['n']):3}%  kernel {pct(d['found'],d['resist']):3}%  "
          f"switch {pct(d['switch'],d['found']):3}%  strict {pct(d['strict'],d['found']):3}%  "
          f"commit {pct(d['commit'],d['found']):3}%")
e41 = run_E41(SAMPLE)
print("E41 monolith:")
for c, d in e41.items():
    n = d['cnt'] or 1
    print(f"  {c:9} apop {pct(d['apop'],n):3}%  prolif {pct(d['prolif'],n):3}%  "
          f"joint {pct(d['joint'],n):3}%  durable {pct(d['dcommit'],n):3}%")

e42 = run_E42_clock(SAMPLE)
print("E42 clock durability  (idealized clock-continuous | clinical pulsed | clinical capped+escalation):")
for c, d in e42.items():
    f = d['found'] or 1
    print(f"  {c:9} clock-cont {pct(d['cont'],f):3}%  pulsed {pct(d['puls'],f):3}%  capped {pct(d['capped'],f):3}%")

current = {'sample': SAMPLE, 'E40': e40, 'E41': e41, 'E42': e42}
print("=" * 78)

if UPDATE or not os.path.exists(BASELINE):
    json.dump(current, open(BASELINE, 'w'), indent=2)
    print(f"BASELINE CAPTURED -> {BASELINE}  (re-run later to check against this)")
else:
    base = json.load(open(BASELINE))
    drift = []
    if base.get('sample') != SAMPLE:
        print(f"NOTE: baseline sample={base.get('sample')} != current {SAMPLE}; compare like-for-like.")
    for blk in ('E40', 'E41', 'E42'):
        for c in current[blk]:
            for k, v in current[blk][c].items():
                bv = base.get(blk, {}).get(c, {}).get(k)
                if bv is not None and bv != v:
                    drift.append(f"{blk}/{c}/{k}: baseline {bv} -> now {v}")
    if drift:
        print("DRIFT DETECTED (a change moved these raw counts):")
        for d in drift: print("  ! " + d)
    else:
        print("PASS: every raw count matches the baseline exactly. No regression.")
print(f"[{time.time()-t0:.0f}s]")
