# Scope — what BBCN does and does not claim

## What it does
- Models a breast tumour as a 135-node Boolean network across 24 signalling pathways.
- Computes minimal drug-target "kernels" that drive the network to target phenotypes
  (resistance-off, apoptosis-on, proliferation-off) as genuine fixed points, no forcing.
- Reproduces results on two independent primary-tumour cohorts (TCGA-BRCA, METABRIC).

## What it does NOT claim
- **No survival benefit is demonstrated.** Controllability did not predict overall survival in
  either cohort (pre-registered test; null under all formulations). These are mechanistic /
  hypothesis-generating results, not validated clinical benefit.
- **No drug-response validation.** The targets the controller recommends (PI3K/AKT/mTOR/SRC axis)
  were essentially not administered to these patients (treated in a chemo/endocrine era), so there
  is no overlap cohort to validate response. This is an era/availability gap, not evidence of effect.
- **Primary tumours only.** Not validated on metastatic disease.
- **Invasion is out of scope as a headline.** An invasion target exists, but the network contains
  only 7 of ~27 canonical EMT drivers (core effectors SNAI/ZEB/CDH1/VIM/MMPs absent), so any
  invasion result is a thin proxy — exploratory at best.
- **Binary, not continuous.** States are binarized (per-gene median). This is a deliberate modelling
  choice for driver-level logic, not a claim that signalling is digital.
- **Sequenced multi-axis control is rare (~3%).** We report this honestly rather than presenting the
  high isolated numbers as if all three phenotypes can be held at once.

## Reproducibility
- Binarized inputs are shipped; `reproduce.py` regenerates all headline numbers and diffs them
  against the committed reference. Raw TCGA/METABRIC matrices are NOT shipped (size + data-use
  terms); `scripts/` documents how to regenerate the binarized inputs from the public sources.
