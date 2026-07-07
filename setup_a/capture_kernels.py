"""
capture_kernels.py -- per-patient kernel-composition capture (read-only).

Runs the CANONICAL staged controller (the exact E41 path used by
regression_check.py: harness + run_cde driving the Delayed-BBCN engine) on each
patient, and records the kernel composition harness ALREADY selects -- which
pathway nodes were pinned, in which stage -- via the read-only
`harness.KERNEL_CAPTURE` hook. It changes no decision: capture is observation
only, and the hook is a no-op whenever KERNEL_CAPTURE is None.

The hook fires inside harness._update_pathway, which is exercised by the
run_cde / E41 path. (E42's clock-durability kernels are designed by
forward_stab_kernel_design and never pass through _update_pathway, so they are
out of scope for this capture by construction.)

Output: one CSV per cohort,
    setup_a/data/kernel_composition_{cohort}.csv
with columns: patient_id, cohort, stage, pathway, kernel_node, kernel_size
(one row per kernel node).

Usage (from repo root):
    python setup_a/capture_kernels.py                  # all cohorts, full N
    python setup_a/capture_kernels.py --cohort TCGA    # one cohort
    python setup_a/capture_kernels.py --sample 50      # first 50 patients/cohort
"""
import os
import sys
import csv
import time
import argparse

# repo root = parent of this file's directory (setup_a/)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)            # so `import bbcn...` resolves like regression_check.py

import pandas as pd
import bbcn.harness as H
from bbcn import repaired_branch as RB
RB.apply()                          # bake the repair into harness.PATHWAYS (once), as E41 does
from bbcn import run_cde as CDE

COH = {
    'TCGA':     os.path.join(ROOT, 'setup_a/data/binarized/tcga_brca_1082x135.csv'),
    'METABRIC': os.path.join(ROOT, 'setup_a/data/binarized/metabric_1980x135.csv'),
    'ISPY2':    os.path.join(ROOT, 'setup_a/data/binarized/ispy2_988x135.csv'),
}

FIELDS = ['patient_id', 'cohort', 'stage', 'pathway', 'kernel_node', 'kernel_size']


def capture_cohort(name, path, sample=None):
    """Run the canonical control on every patient in one cohort, return kernel rows."""
    df = pd.read_csv(path, index_col=0)
    n = len(df) if sample is None else min(sample, len(df))
    rows = []
    t0 = time.time()
    for i in range(n):
        row = df.iloc[i]
        b = {nd: int(row.get(nd, 0)) for nd in df.columns}
        # init EXACTLY as regression_check.run_E41 does (verbatim canonical setup)
        U = RB.patient_clamp_off(b)
        init = dict(b)
        init['CLAMP_OFF'] = U
        init['PHLPP'] = int(b.get('FOXO3', 0) and not U)
        init['COMMIT'] = 0
        pid = str(df.index[i])

        H.KERNEL_CAPTURE = []                       # turn capture ON for this patient
        CDE.run_patient_cde(init, kernel_method='stabilize')
        recs = H.KERNEL_CAPTURE
        H.KERNEL_CAPTURE = None                     # turn capture OFF again

        # _update_pathway records its choice every inner step; the kernel COMPOSITION
        # is the distinct (stage, pathway, node) pinned, so collapse timestep repeats
        # per patient (cross-patient identical compositions stay as separate rows).
        seen = set()
        for rec in recs:
            for node in rec['kernel_nodes']:
                key = (rec['stage'], rec['pathway'], node, rec['kernel_size'])
                if key in seen:
                    continue
                seen.add(key)
                rows.append({
                    'patient_id': pid, 'cohort': name,
                    'stage': rec['stage'], 'pathway': rec['pathway'],
                    'kernel_node': node, 'kernel_size': rec['kernel_size'],
                })
    print(f"  {name:9} {n:5} patients -> {len(rows):6} kernel-node rows "
          f"[{time.time() - t0:.0f}s]")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohort', default=None, help='one of TCGA / METABRIC / ISPY2 (default: all)')
    ap.add_argument('--sample', type=int, default=None, help='first N patients/cohort (default: full N)')
    args = ap.parse_args()

    cohorts = [args.cohort] if args.cohort else list(COH)
    outdir = os.path.join(ROOT, 'setup_a', 'data')
    os.makedirs(outdir, exist_ok=True)

    print(f"Kernel capture  |  cohorts={cohorts}  |  sample={args.sample or 'full-N'}")
    for name in cohorts:
        if name not in COH:
            sys.exit(f"unknown cohort {name!r}; choose from {list(COH)}")
        rows = capture_cohort(name, COH[name], sample=args.sample)
        out = os.path.join(outdir, f'kernel_composition_{name}.csv')
        with open(out, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)
        print(f"  wrote {out}")


if __name__ == '__main__':
    main()
