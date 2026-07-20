# Paper 3 v3 build

Paper 3 v3 is the merged manuscript (Paper 3 with the Paper 4 additions). Like
`paper3/`, this folder is **self-contained** (local `figures/`, `numbers*.tex`,
`appendix_*.tex`) so it builds in place and syncs to Overleaf as-is.

## Build (self-contained; this is what Overleaf does)
    latexmk -pdf BBCN_paper3v3.tex          # or: bash build.sh   ->  BBCN_paper3v3.pdf

## Re-assemble from the pipeline (regenerate + refresh the committed copies)
    bash assemble.sh                        # run in the repo; needs preprint/, results/, setup_a/

`paper3v3` is a **sibling** of `paper3`, not a copy of it: `assemble.sh` pulls the
shared assets straight from the generators' outputs, the same sources `paper3`
uses, so the two never drift apart.

## Where each committed file comes from
Regenerated here by the two scripts in this folder:
  numbers_v2.tex             <- make_numbers_v2.py  (setup_a/data/cde_vs_switch_summary.csv)
  figures/f06_durability.pdf <- fig6.py             (same CSV)
Copied from the canonical pipeline outputs:
  numbers.tex        <- results/numbers/numbers.tex        (generate_numbers.py)
  numbers2.tex       <- results/numbers2/numbers2.tex      (generate_numbers2.py)
  appendix_rules.tex <- preprint/appendix_rules.tex        (generate_rules_appendix.py)
  appendix_theorem.tex <- preprint/appendix_theorem.tex    (static, prior-work theorem)
  figures/*.png      <- preprint/figures/*                 (make_figures.py, generate_flowcharts.py, ...)

## When to run assemble.sh
Only when upstream results change (rerun of the cde_vs_switch pipeline,
generate_numbers2.py, or a figure generator). Run it locally, commit, and push;
Overleaf then syncs the refreshed files. Overleaf never runs Python; it only
typesets the committed .tex and figures. Plain text edits need no assembly.

## Overleaf
Sync the repo and set the main document to `paper3v3/BBCN_paper3v3.tex`. The .py
and .sh files sit alongside and are ignored by the LaTeX compiler.
