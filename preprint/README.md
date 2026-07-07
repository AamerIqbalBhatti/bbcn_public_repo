# BBCN preprint

LaTeX source of the paper, plus the compiled PDF and a DOCX. Every number, figure, and
appendix is generated from the repository root and `\input` here, so the text cannot drift
from the code. Build with `python reproduce_all.py` (or `--all` for a full regenerate).

## Files
- `BBCN_preprint.tex`   LaTeX source (pdflatex; standard packages only)
- `BBCN_preprint.pdf`   compiled PDF
- `BBCN_preprint.docx`  same content as a Word document (pandoc)
- `numbers.tex`         auto-generated numbers (\input; from `run_chunked.py --emit`)
- `appendix_rules.tex`  auto-generated Appendix A, the digital logic model (\input; from `generate_rules_appendix.py`)
- `appendix_theorem.tex` Appendix C, attributed restatement of the prior-work theorem (static text)
- `figures/`
  - `fig_arch.png`              two-setup architecture
  - `fig2_scc.png`              pathway coupling graph (the strongly-connected core)
  - `figB.png`                  switch routing and biological check
  - `figB_kernel_actuator.png`  Figure 4, Setup B kernel actuators (PI3K/AKT/mTOR)
  - `figA_kernel_usage.png`     Figure 5, Setup A kernel usage by cohort
  - `figA_pathway_kernel.png`   Figure 6, Setup A pathway-by-kernel modularity
  - `figC_flow_ranked.png`      Appendix B, ranked-method flowchart
  - `figC_flow_stabilize.png`   Appendix B, stabilize-method flowchart

## Build
```
cd ..            # repository root
python reproduce_all.py          # build from committed numbers + figures
python reproduce_all.py --all    # regenerate numbers + figures, then build
```
No bibliography tool is needed (references are inline via `thebibliography`).
