# BBCN: the Boolean Breast-Cancer Network

Control-theoretic Boolean model of breast-cancer signalling. This repository
reproduces every number, figure, and table of the manuscript from binarised
cohort inputs with a single command.

## Quick start
```bash
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
| 1 | Bare Setup A whole-network controller | `results/numbers/` (`numbers.tex`) | Table 1; Sec 3.1-3.2; Supplement |
| 2 | Repaired / delayed 136-node network (slow) | `results/numbers2/` (`numbers2.tex`) | Abstract; Table 3; Sec 3.5 |
| 3 | Setup B 9-node bistable switch | routing / hysteresis | Table 2; Fig 4-5; Sec 3.3-3.4 |
| 4 | Switch vs whole-network durability/minimality | `cde_vs_switch_summary.csv` (`numbers_v2.tex`) | Table 4; Fig 7; Sec 3.6 |

All four are required; each writes a different set of macros/figures the manuscript reads.
The committed `numbers*.tex` carry the published (aggregate) values for all three
cohorts, so the manuscript always compiles; `run_all.py` refreshes whichever
cohorts are present.

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
