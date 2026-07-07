# Plain-language captions for the reference outputs

**REFERENCE_TABLE.csv** — For each cohort and regime, the percentage of patients whose tumour
network the controller drove to a *genuine* controlled fixed point (no forcing). "Isolated" treats
each phenotype on its own from the patient's state; "Sequenced + state-preservation" requires
reaching all three in order without undoing earlier gains.

**controllability_TCGA.png / controllability_METABRIC.png** — Each blue bar is the % of patients
for whom that phenotype alone is controllable; each orange bar is the % when phenotypes must be
achieved in sequence while preserving the earlier ones. The tall blue bars (apoptosis ~86%,
proliferation ~94%) versus the short orange bars (single digits; all-three ~3%) show the central
finding: phenotypes are individually controllable but jointly hard. Resistance-off (~4–8%) is the
hardest in every regime.
