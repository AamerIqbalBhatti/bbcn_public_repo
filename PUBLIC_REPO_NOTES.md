# Public repo: what changed vs the private working repo

This is the public, TCGA-only version. Built from the working repo with:

## Removed (controlled / third-party data, not redistributable)
All METABRIC and I-SPY2 patient-level derivatives (22 files): binarised matrices,
switch inputs, mutation gates, per-cohort cde/kernel-composition results, payoff/
pulsed JSONs. TCGA-BRCA (open access) is kept. A .gitignore prevents re-adding them.

## Removed (redundant / dead code, zero references)
run_repaired_cohorts.py, run_unified_repaired.py, run_multirate_test.py,
run_designed_kernel_AB.py, setup_a/harness.py, setup_a/honest_bbcn.py.
(The referenced repaired runners -- run_full_cohorts_repaired.py,
run_repaired_monolith.py, run_designed_kernel_cohorts.py -- were kept.)

## Removed (bloat)
BBCN_repository_v22.zip (15 MB nested archive).

## Added
- run_all.py : single auto-detecting entry point (+ run_all.sh wrapper).
- Cohort dicts in generate_numbers.py, cde_vs_switch_cascade.py, regression_check.py
  now filter to files that exist, so the code runs on whatever cohorts are present.
- README.md (public), DATA_ACCESS.md, .gitignore.
- paper3v2/fig1.py and paper3v2/fig_pipelines.py (code-based figures).

## Reproducibility contract
Committed numbers*.tex carry published aggregate values for all three cohorts, so
the manuscript always compiles. run_all.py reproduces whatever cohorts are present
(TCGA out of the box). Smoke test: `python run_all.py --quick` -> PASS on TCGA.
