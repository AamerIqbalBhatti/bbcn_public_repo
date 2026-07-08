#!/usr/bin/env bash
# Assemble the self-contained paper3v2 folder from the PIPELINE SOURCES.
# paper3v2 is a SIBLING of paper3: both draw from preprint/ and results/, so
# paper3v2 does not depend on paper3 and never lags behind it.
# Run from inside the repo (needs preprint/, results/, setup_a/ present).
set -e; cd "$(dirname "$0")"
mkdir -p figures

echo "[1/3] regenerate paper3v2-only artifacts (from ../setup_a/data/cde_vs_switch_summary.csv)"
python3 make_numbers_v2.py                                   # -> numbers_v2.tex
python3 fig6.py                                              # -> figures/f06_durability.pdf
python3 fig1.py                                             # -> figures/fig_arch.png (schematic, code-based)

echo "[2/3] pull shared assets from their canonical generators' outputs"
cp ../results/numbers/numbers.tex    numbers.tex             # <- generate_numbers.py
cp ../results/numbers2/numbers2.tex  numbers2.tex            # <- generate_numbers2.py
cp ../preprint/appendix_rules.tex    appendix_rules.tex      # <- generate_rules_appendix.py
cp ../preprint/appendix_theorem.tex  appendix_theorem.tex    # <- static (prior-work theorem)
for f in figA_tnbc_panel fig2_scc figB figB_kernel_actuator \
         figB_tnbc_panel figA_kernel_usage figA_pathway_kernel \
         figC_flow_ranked figC_flow_stabilize; do
  cp "../preprint/figures/$f.png" figures/                   # <- figure generators (fig_arch is now local, from fig1.py)
done

echo "[3/3] build in place"
latexmk -pdf -interaction=nonstopmode -halt-on-error BBCN_paper3v2.tex
latexmk -c
echo "assembled -> BBCN_paper3v2.pdf  (self-contained sibling of paper3)"
