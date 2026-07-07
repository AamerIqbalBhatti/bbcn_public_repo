#!/usr/bin/env python3
"""
reproduce_combined.py - one entry point for the COMBINED repository (v22).

The combined repo carries TWO networks in one tree:
  * BARE network    (135 nodes, pathways.py untouched)  -> paper 1 numbers  (results/numbers/)
  * REPAIRED network (136 nodes, overlay via repaired_branch.apply) -> paper 3 numbers (results/numbers2/)

The repaired network is NOT a fork. It is the bare network plus three biological
repairs applied as an overlay (PHLPP, COMMIT latch, per-patient U flag). Selecting a
network is a one-line switch:

    from bbcn import repaired_branch as RB; RB.apply()   # turn the controller repaired

Usage:
    python reproduce_combined.py bare        # regenerate paper-1 numbers (numbers.tex)
    python reproduce_combined.py repaired    # regenerate paper-3 numbers (numbers2.tex)
    python reproduce_combined.py all         # both

'repaired' re-runs the four repaired pipelines, then generate_numbers2.py:
    E36 setup_a/run_full_cohorts_repaired.py   (resistance == AKT-FOXO3a uncoupling)
    E38 run_designed_kernel_cohorts.py         (kernel nomination, N=200/cohort)
    E40 run_durability_full.py                 (switch / strict / commit durability)
    E41 setup_a/run_repaired_monolith.py       (staged controller, per-cohort, full N)
The monolith is the slow step (~19 min full N); pass --quick to sample N=200/cohort.
"""
import subprocess, sys, os, shutil

ROOT = os.path.dirname(os.path.abspath(__file__))
SA   = os.path.join(ROOT, 'setup_a')

def run(cmd, cwd):
    print('>>', ' '.join(cmd), '   (cwd=%s)' % os.path.relpath(cwd, ROOT))
    return subprocess.call(cmd, cwd=cwd)

def bare():
    # paper-1 numbers: the existing locked pipeline
    if os.path.exists(os.path.join(ROOT, 'run_chunked.py')):
        run([sys.executable, 'run_chunked.py', '--emit'], ROOT)
    run([sys.executable, 'generate_numbers.py'], ROOT)
    print('bare (paper 1) numbers refreshed in results/numbers/')

def repaired(quick=False):
    outdir = os.path.join(ROOT, 'results', 'numbers2', 'run_outputs')
    os.makedirs(outdir, exist_ok=True)
    # E36 full cohorts (writes to stdout -> capture)
    with open(os.path.join(outdir, 'full_cohorts_E36.log'), 'w') as f:
        subprocess.call([sys.executable, 'run_full_cohorts_repaired.py'], cwd=SA, stdout=f, stderr=subprocess.STDOUT)
    # E38 designed kernel  (writes designed_kernel_agg.json into cwd)
    run([sys.executable, 'run_designed_kernel_cohorts.py'], ROOT)
    if os.path.exists(os.path.join(ROOT, 'designed_kernel_agg.json')):
        shutil.move(os.path.join(ROOT, 'designed_kernel_agg.json'), os.path.join(outdir, 'designed_kernel_E38.json'))
    # E40 durability  (writes durability_full.json into cwd)
    run([sys.executable, 'run_durability_full.py'], ROOT)
    if os.path.exists(os.path.join(ROOT, 'durability_full.json')):
        shutil.move(os.path.join(ROOT, 'durability_full.json'), os.path.join(outdir, 'durability_full_E40.json'))
    # E41 monolith per cohort (full N, or quick sample)
    # NOTE: full N is slow; the locked E41 values are re-verified, not overwritten, by --quick.
    print('E41 monolith: run setup_a/run_repaired_monolith.py per cohort for full N.')
    # finally lock the numbers
    run([sys.executable, 'generate_numbers2.py'], ROOT)
    print('repaired (paper 3) numbers refreshed in results/numbers2/')

if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    quick = '--quick' in sys.argv
    if what in ('bare', 'all'): bare()
    if what in ('repaired', 'all'): repaired(quick=quick)
    if what not in ('bare', 'repaired', 'all'):
        print(__doc__)
