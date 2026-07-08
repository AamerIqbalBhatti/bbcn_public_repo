# BBCN: the Boolean Breast-Cancer Network

Control-theoretic Boolean model of breast-cancer signalling. This repository
reproduces every number, figure, and table of the manuscript from binarised
cohort inputs with a single command.

## Quick start
```bash
git clone https://github.com/AamerIqbalBhatti/bbcn_public_repo
cd bbcn_public_repo
pip install -r requirements.txt
python run_all.py --quick     # fast smoke test (seconds): pipelines wired + reproduce reference
python run_all.py             # full run on every cohort present
python run_all.py --list      # show which cohorts were detected
python run_all.py --pipelines 1,4    # run a subset
```
Plain Python; works on Windows, macOS and Linux. `run_all.sh` is an optional wrapper.

## Cohorts are auto-detected
`run_all.py` runs on whatever binarised inputs are in `setup_a/data/binarized/`.
This repository ships **TCGA-BRCA only** (open access), so a fresh clone
reproduces the TCGA results out of the box. To reproduce all three cohorts, add
the METABRIC and I-SPY2 binarised inputs to the same folder (they are under
controlled / third-party terms and are not redistributed here; see
`DATA_ACCESS.md`). The code is identical either way, behaviour follows the data.

## The four pipelines
| # | Pipeline | Writes | Paper element |
|---|----------|--------|---------------|
| 1 | Baseline whole-network controller (Setup A) | `results/numbers/` (`numbers.tex`) | Table 1; Sec 3.1-3.2; Supplement |
| 2 | Biology-faithful / delayed 136-node network (slow) | `results/numbers2/` (`numbers2.tex`) | Abstract; Table 3; Sec 3.5 |
| 3 | Setup B 9-node bistable switch | routing / hysteresis | Table 2; Fig 4-5; Sec 3.3-3.4 |
| 4 | Switch vs whole-network durability/minimality | `cde_vs_switch_summary.csv` (`numbers_v2.tex`) | Table 4; Fig 7; Sec 3.6 |

All four are required; each writes a different set of macros/figures the manuscript reads.
The committed `numbers*.tex` carry the published (aggregate) values for all three
cohorts, so the manuscript always compiles; `run_all.py` refreshes whichever
cohorts are present.

## Run or test an individual pipeline
```bash
python run_all.py --list          # show detected cohorts (instant)
python run_all.py --quick         # fast smoke test: regression check on the present cohort(s)

python run_all.py --pipelines 1   # baseline whole-network controller  -> results/numbers/
python run_all.py --pipelines 2   # biology-faithful network (SLOW, minutes) -> results/numbers2/
python run_all.py --pipelines 3   # Setup B bistable switch (routing / hysteresis)
python run_all.py --pipelines 4   # switch vs whole-network durability -> cde_vs_switch_summary.csv
python run_all.py --pipelines 3,4 # any comma-separated subset
```
**How to check a run.** `python run_all.py --quick` runs `regression_check.py` and prints
`PASS: every raw count matches the baseline` when the present cohort's numbers match the
committed reference. For an exact match at the reference sample size, run
`python regression_check.py --sample 30`. Pipelines 1 and 2 are the whole-network
computations and take minutes; pipelines 3 and 4 are faster.

## Rebuild the manuscript
```bash
cd paper3v2 && bash assemble.sh      # regenerates figures/numbers, builds BBCN_paper3v2.pdf
```
The `paper3v2/` folder is self-contained and can be synced to Overleaf directly.

## Other run routes
- Direct CLI: `python run_all.py` (above).
- Google Colab / Binder: clone the repo and run the same `python run_all.py`.

## Layout
```
run_all.py            single entry point (auto-detects cohorts)
setup_a/              Setup A: 135/136-node network, controller, kernels (bbcn/ package)
setup_b/              Setup B: 9-node bistable switch, multirate schedule
results/              committed aggregate numbers and run logs
paper3v2/             self-contained manuscript (LaTeX) + figure generators
DATA_ACCESS.md        how to obtain the controlled cohorts
```
