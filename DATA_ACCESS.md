# Data access

This repository ships binarised inputs for **TCGA-BRCA only**, which is open access.
It does **not** redistribute patient-level data for the other two cohorts:

- **METABRIC** is under controlled access (European Genome-phenome Archive; processed
  data via cBioPortal under its data-use terms). Obtain it from the source, then run
  the binarisation step to regenerate the inputs used here.
- **I-SPY2** expression is available from NCBI GEO, series **GSE194040**. Download it
  and run the same binarisation step to regenerate the inputs.

The binarisation is a deterministic robust z-score (median / MAD) thresholded at zero,
so an approved user reproduces byte-identical binarised matrices and every downstream
number. Only TCGA-BRCA is required to run the pipeline end to end out of the box.
