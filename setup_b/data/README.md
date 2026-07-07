# Setup B — input data

## Default: small pre-extracted CSVs (shipped, ~0.3-0.7 MB each)
`reproduce_all.py` runs from these automatically. They live in `samples/` and contain only the
59 genes the switch needs, extracted from the public matrices:

| Cohort | File |
|---|---|
| TCGA-BRCA | `samples/tcga_brca_switch_inputs.tsv` |
| METABRIC | `samples/metabric_switch_inputs.tsv` |
| I-SPY2 | `samples/ispy2_switch_inputs.tsv` |

These reproduce the locked numbers exactly. No large downloads needed.

## Optional: full public matrices
To re-run from the originals (e.g. to re-extract or audit), drop these here and the pipeline will
use them if the small CSVs are absent:

| Cohort | File | Source |
|---|---|---|
| TCGA-BRCA | `data_mrna_seq_v2_rsem.txt` | cBioPortal brca_tcga_pan_can_atlas_2018 |
| METABRIC | `data_mrna_illumina_microarray.txt` | cBioPortal brca_metabric |
| I-SPY2 | `GSE194040_..._n988_txt.gz` | GEO GSE194040 |

Phospho front-end uses I-SPY2 RPPA from GEO GSE196093 (SuperSeries GSE196096).
