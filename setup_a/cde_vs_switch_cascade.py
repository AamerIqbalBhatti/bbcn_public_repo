"""
cde_vs_switch_cascade.py -- Paper 4 closing comparison.

ONE test engine, TWO kernel families. Both kernels are clamped and tested on the
SAME full cascaded delayed BBCN-D engine (honest_bbcn.step_honest(cascade=True)
over the lagged history), with the SAME durability readout. The only thing that
differs between arms is where the kernel came from.

  Switch arm (Setup B / BBCN-S):
      kernel = forward_stab_kernel_design.design_patient_kernel(x0, I, U,
               method='stabilize') on the 9-node switch  (this is exactly E40).
      clamp  = those nodes -> APOP_TARGET.

  CDE arm (Setup A):
      kernel = the nodes the harness/CDE controller pins for the patient, read
               off via the C6 KERNEL_CAPTURE hook while running run_patient_cde.
      clamp  = those nodes -> apoptotic targets (APOP_FULL, see below).

PROVENANCE NOTE (read this):
  The CDE kernel is DESIGNED on the *repaired* harness (BBCN-M), which retains the
  abstract COMMIT latch + COMMIT/CLAMP_OFF nodes (repaired_branch.apply()). The
  switch kernel is designed on the 9-node switch (BBCN-S). BOTH are then TESTED on
  the flag-free cascaded delayed engine (BBCN-D), where commitment is the real
  caspase cascade (PUMA/NOXA + SMAC->XIAP + CASP6->CASP8->CASP3 loop), NO abstract
  flag. So the durability READOUT is flag-free and identical for both arms; only
  the CDE kernel's design model carries the flag. (Decision: proceed + document.)

A_dyn dynamics and the strict/durable readout below are VERBATIM from
regression_check.py A_dyn (lines 66-78). The only addition is `held` -- a second
CASP3 readout taken from the hold-phase tail (for the held-vs-durable figure); it
does not touch the strict/durable computation, so the switch arm still reproduces
E40's locked 'commit' column exactly (Step 2 cross-check).

Outputs:
  setup_a/data/cde_vs_switch_{cohort}.csv   per-patient, both arms
  setup_a/data/cde_vs_switch_summary.csv    cohort x arm summary

Usage (from repo root):
  python setup_a/cde_vs_switch_cascade.py                 # all cohorts, full N
  python setup_a/cde_vs_switch_cascade.py --cohort TCGA   # one cohort (chunking)
  python setup_a/cde_vs_switch_cascade.py --sample 30     # cross-check sample
"""
import os
import sys
import csv
import time
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'setup_b', 'code'))

import pandas as pd
from bbcn import honest_bbcn as HB
from bbcn import harness as H
from bbcn import repaired_branch as RB
from bbcn import run_cde as CDE
import forward_stab_kernel_design as FSK

RB.apply()                                  # CDE arm runs on the repaired harness (as E41 does)

# ----- verbatim E40 scaffolding (regression_check.py) -----
COH = {
    'TCGA':     os.path.join(ROOT, 'setup_a/data/binarized/tcga_brca_1082x135.csv'),
    'METABRIC': os.path.join(ROOT, 'setup_a/data/binarized/metabric_1980x135.csv'),
    'ISPY2':    os.path.join(ROOT, 'setup_a/data/binarized/ispy2_988x135.csv'),
}
COH = {k: v for k, v in COH.items() if os.path.exists(v)}   # auto-detect: run only cohorts present
NODES = HB.NODES
lags = HB.make_lags()
APOP = FSK.APOP_TARGET
SW_STATE = ['AKT1', 'PTEN', 'MTOR', 'PDPK1', 'PIK3CA', 'MDM2', 'TP53', 'FOXO3']

# apoptotic-target map for the CDE arm: harness Apoptosis_ON full-bus target,
# with APOP_TARGET layered on top (it supplies PHLPP=1, absent from the 135-node
# stage target, and agrees with it on all other switch nodes).
H._load_stage_targets('compatible')
APOP_FULL = {**H._STAGE_TARGETS['Apoptosis_ON'], **APOP}

# Cascade-internal executioner/brake nodes that honest_bbcn.step_honest(cascade=True)
# COMPUTES each step (the real caspase machinery). A kernel pins upstream control
# drivers; freezing these downstream nodes to a static target would override the very
# cascade we test (e.g. clamping CASP9/CYCS=0 blocks CASP3 while held, an artifact —
# not biology). The switch kernel only ever pins upstream nodes, so for an
# apples-to-apples comparison the CDE clamp excludes this set and lets the cascade
# drive them, exactly as for the switch arm.
CASCADE_INTERNAL = {'PUMA', 'NOXA', 'MCL1', 'BCL2', 'BCL2L1', 'SMAC', 'XIAP',
                    'CASP3', 'CASP6', 'CASP7', 'CASP8', 'CASP9', 'CYCS', 'BAX', 'BAK1'}


def seed(b):
    U = HB.patient_clamp_off(b); x0 = {nd: int(b.get(nd, 0)) for nd in SW_STATE}
    x0['PHLPP'] = int(x0['FOXO3'] and not U)
    I = dict(SRC=int(b.get('SRC', 0)), RHEB=int(b.get('RHEB', 0)), IGF1R=int(b.get('IGF1R', 0)),
             RTK_up=int(b.get('GRB2', 0) or b.get('IRS1', 0)), CDKN2A=int(b.get('CDKN2A', 0)),
             E2F1=int(b.get('E2F1', 0)), ATM=int(b.get('ATM', 0)), ATR=int(b.get('ATR', 0)))
    return x0, I, U


def A_dyn(b, U, clamp, T1=120, T2=180, theta=10, cmax=20):
    """VERBATIM regression A_dyn dynamics + strict/durable; `held` is an added
    second readout from the hold-phase tail (does not affect strict/durable)."""
    bus = {n: int(b.get(n, 0)) for n in NODES}
    for k, v in clamp.items(): bus[k] = v
    hist = [bus]
    for t in range(T1 + T2):
        nb = HB.step_honest(HB._view(hist, t, lags), U, cascade=True)   # Step 16: real cascade
        if t < T1:
            for k, v in clamp.items(): nb[k] = v
        hist.append(nb)
    tail = hist[-40:]                                   # post-release window (verbatim)
    strict = sum(s['CASP3'] == 1 and s['AKT1'] == 0 for s in tail) / 40 > 0.5
    commitd = sum(s['CASP3'] == 1 for s in tail) / 40 > 0.5
    held_tail = hist[T1 - 40:T1]                        # added: last 40 steps while held
    held = sum(s['CASP3'] == 1 for s in held_tail) / 40 > 0.5
    return held, strict, commitd


def cde_kernel(b, U):
    """Run the harness/CDE control with capture ON; return (union_nodes, apop_nodes,
    pathways) the controller pinned. union = across the whole run; apop = pinned in
    the Apoptosis_ON stage only."""
    init = dict(b)
    init['CLAMP_OFF'] = U
    init['PHLPP'] = int(b.get('FOXO3', 0) and not U)
    init['COMMIT'] = 0
    H.KERNEL_CAPTURE = []
    CDE.run_patient_cde(init, kernel_method='stabilize')
    recs = H.KERNEL_CAPTURE
    H.KERNEL_CAPTURE = None
    union_nodes, apop_nodes, pathways = set(), set(), set()
    for rec in recs:
        pathways.add(rec['pathway'])
        for nd in rec['kernel_nodes']:
            union_nodes.add(nd)
            if rec['stage'] == 'Apoptosis_ON':
                apop_nodes.add(nd)
    return union_nodes, apop_nodes, pathways


FIELDS = ['patient_id', 'cohort', 'resistant', 'switch_found',
          'switch_size', 'switch_nodes', 'switch_held', 'switch_durable',
          'cde_found', 'cde_size_union', 'cde_nodes_union', 'cde_pathways',
          'cde_held_union', 'cde_durable_union',
          'cde_size_apop', 'cde_held_apop', 'cde_durable_apop']


def run_cohort(name, path, sample=None):
    df = pd.read_csv(path, index_col=0)
    n = len(df) if sample is None else min(sample, len(df))
    rows = []
    agg = dict(N=n, resist=0, found=0,
               sw_held=0, sw_dur=0, sw_size=0,
               cde_found=0, cde_held_u=0, cde_dur_u=0, cde_size_u=0,
               cde_held_a=0, cde_dur_a=0, cde_size_a=0)
    t0 = time.time()
    for i in range(n):
        b = {nd: int(df.iloc[i].get(nd, 0)) for nd in NODES}
        pid = str(df.index[i])
        x0, I, U = seed(b)

        # ---- switch arm (exactly E40's design + denominators) ----
        r, kernels = FSK.design_patient_kernel(x0, I, U, method='stabilize')
        if not kernels and r == 'APOPTOTIC':
            continue                                    # not resistant (already apoptotic)
        agg['resist'] += 1
        resistant = True
        if not kernels:
            # resistant but no switch kernel: record and skip the paired comparison
            rows.append({**{f: '' for f in FIELDS}, 'patient_id': pid, 'cohort': name,
                         'resistant': 1, 'switch_found': 0, 'cde_found': 0})
            continue
        agg['found'] += 1
        clamp_S = {nd: APOP[nd] for bx in kernels for nd in kernels[bx]}
        sw_held, sw_strict, sw_dur = A_dyn(b, U, clamp_S)
        agg['sw_held'] += sw_held; agg['sw_dur'] += sw_dur; agg['sw_size'] += len(clamp_S)

        # ---- CDE arm (same patient, kernel from the harness/CDE controller) ----
        union_nodes, apop_nodes, pathways = cde_kernel(b, U)
        clamp_Cu = {nd: APOP_FULL[nd] for nd in union_nodes
                    if nd in APOP_FULL and nd not in CASCADE_INTERNAL}
        clamp_Ca = {nd: APOP_FULL[nd] for nd in apop_nodes
                    if nd in APOP_FULL and nd not in CASCADE_INTERNAL}
        cu_held, cu_strict, cu_dur = A_dyn(b, U, clamp_Cu) if clamp_Cu else (False, False, False)
        ca_held, ca_strict, ca_dur = A_dyn(b, U, clamp_Ca) if clamp_Ca else (False, False, False)
        cde_found = 1 if clamp_Cu else 0
        agg['cde_found'] += cde_found
        agg['cde_held_u'] += cu_held; agg['cde_dur_u'] += cu_dur; agg['cde_size_u'] += len(clamp_Cu)
        agg['cde_held_a'] += ca_held; agg['cde_dur_a'] += ca_dur; agg['cde_size_a'] += len(clamp_Ca)

        rows.append({
            'patient_id': pid, 'cohort': name, 'resistant': 1, 'switch_found': 1,
            'switch_size': len(clamp_S), 'switch_nodes': ';'.join(sorted(clamp_S)),
            'switch_held': int(sw_held), 'switch_durable': int(sw_dur),
            'cde_found': cde_found, 'cde_size_union': len(clamp_Cu),
            'cde_nodes_union': ';'.join(sorted(clamp_Cu)), 'cde_pathways': ';'.join(sorted(pathways)),
            'cde_held_union': int(cu_held), 'cde_durable_union': int(cu_dur),
            'cde_size_apop': len(clamp_Ca), 'cde_held_apop': int(ca_held),
            'cde_durable_apop': int(ca_dur),
        })
    agg['secs'] = round(time.time() - t0)
    print(f"  {name:9} N={n} resist={agg['resist']} found={agg['found']} "
          f"| switch durable={agg['sw_dur']} | cde durable(union)={agg['cde_dur_u']} "
          f"[{agg['secs']}s]")
    return rows, agg


def summarize(name, agg):
    f = agg['found'] or 1
    cf = agg['cde_found'] or 1
    def pct(x, d): return round(100 * x / d, 1) if d else 0.0
    return {
        'cohort': name, 'N': agg['N'], 'resistant': agg['resist'], 'found': agg['found'],
        'cde_found': agg['cde_found'],
        'switch_held_pct': pct(agg['sw_held'], agg['found']),
        'switch_durable_pct': pct(agg['sw_dur'], agg['found']),
        'switch_mean_size': round(agg['sw_size'] / f, 2),
        'cde_union_held_pct': pct(agg['cde_held_u'], agg['cde_found']),
        'cde_union_durable_pct': pct(agg['cde_dur_u'], agg['cde_found']),
        'cde_union_mean_size': round(agg['cde_size_u'] / cf, 2),
        'cde_apop_held_pct': pct(agg['cde_held_a'], agg['cde_found']),
        'cde_apop_durable_pct': pct(agg['cde_dur_a'], agg['cde_found']),
        'cde_apop_mean_size': round(agg['cde_size_a'] / cf, 2),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohort', default=None)
    ap.add_argument('--sample', type=int, default=None)
    args = ap.parse_args()
    cohorts = [args.cohort] if args.cohort else list(COH)
    outdir = os.path.join(ROOT, 'setup_a', 'data')
    os.makedirs(outdir, exist_ok=True)

    print(f"CDE vs switch on cascaded BBCN-D  |  cohorts={cohorts}  |  sample={args.sample or 'full-N'}")
    summary_path = os.path.join(outdir, 'cde_vs_switch_summary.csv')
    summaries = []
    # preserve existing summary rows for cohorts we are not re-running (chunked runs)
    if os.path.exists(summary_path) and args.cohort:
        prev = pd.read_csv(summary_path)
        summaries = [r for _, r in prev.iterrows() if r['cohort'] not in cohorts]
        summaries = [dict(r) for r in summaries]

    for name in cohorts:
        rows, agg = run_cohort(name, COH[name], sample=args.sample)
        out = os.path.join(outdir, f'cde_vs_switch_{name}.csv')
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
        print(f"  wrote {out}")
        summaries.append(summarize(name, agg))

    sdf = pd.DataFrame(summaries)
    order = {c: i for i, c in enumerate(COH)}
    sdf = sdf.sort_values('cohort', key=lambda s: s.map(order)).reset_index(drop=True)
    sdf.to_csv(summary_path, index=False)
    print(f"  wrote {summary_path}")
    print(sdf.to_string(index=False))


if __name__ == '__main__':
    main()
