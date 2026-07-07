# Consolidation work item: new / overlay nodes to introduce properly

Every node below currently lives as an overlay, an inline assignment, or an abstract flag rather than
as a first-class network node (a `PATHWAYS` entry with its own rule, lag, bus persistence, and a
binarised initial condition from patient mRNA). This list is the formalisation backlog, to be worked
during consolidation (Phase 6), AFTER the current phase closes. Data availability checked in the
files in hand: I-SPY2 (GSE194040 expression) and METABRIC (tarball expression). TCGA marked pending
only because the full expression matrix is not yet in the sandbox; these are standard genes and are
certainly in TCGA-BRCA RNA-seq.

## A. Data-seedable real proteins — promote to proper nodes
Formalise: add a `PATHWAYS` rule, a lag, bus persistence, and binarise the gene per cohort (median
split, matching the existing 135 nodes).

| node | gene symbol | I-SPY2 | METABRIC | TCGA | origin | role |
|---|---|---|---|---|---|---|
| PUMA | BBC3 | yes | yes | pending file | new cascade (this phase) | p53 -| BCL2/BCL-xL |
| NOXA | PMAIP1 | yes | yes | pending file | new cascade (this phase) | p53 -| MCL1 |
| SMAC | DIABLO | yes | yes | pending file | new cascade (this phase) | -| XIAP, unblocks executioner |
| CASP6 | CASP6 | yes | yes | pending file | Phase 5 plan (step 15) | executioner caspase |
| CASP7 | CASP7 | yes | yes | pending file | Phase 5 plan (step 15) | executioner caspase |
| PHLPP | PHLPP1 / PHLPP2 | yes | yes | pending file | repaired-branch overlay (136th node) | closes AKT1-FOXO3 loop |

Note on PHLPP: currently a DERIVED overlay (`PHLPP = FOXO3 and not CLAMP_OFF`), representing
phosphatase activity, not expression. Promoting it to data-seeded (from PHLPP1/PHLPP2 mRNA) is a
modelling choice — functional activity vs measured transcript — and should be decided, not assumed.

## B. Non-measurable scaffolding — retire or keep as derived input
These are not proteins that can be binarised; they are constructs. They do NOT join list A.

| node | what it is | disposition |
|---|---|---|
| COMMIT | abstract theta-accumulator commitment flag | RETIRE once the real PUMA/NOXA/SMAC cascade owns commitment (Phase 5 steps 16-17). Cascade already validated to reproduce it: 98% vs the flag's 89%. |
| CLAMP_OFF (U) | inferred failure of the AKT->14-3-3 nuclear-export clamp on FOXO3 | KEEP as a derived per-patient input. Its cause (14-3-3 phospho-state, JNK/MST1) is not measured in any cohort, so it cannot be data-seeded; it is defined by its baseline effect. Documented as such. |

## Sequencing
1. Finish the current phase (genotoxic closed; cascade built and guarded).
2. Commit-readout switch (Phase 5 step 16): turn the COMMIT flag off, commit through the real cascade
   with an irreversibility (MOMP) latch.
3. Consolidation (Phase 6): formalise every node in list A as a data-seeded `PATHWAYS` node, binarise
   the genes across all three cohorts, re-lock at full N, and update the node count and convention.

Blocking input still needed: the TCGA expression matrix (cBioPortal `brca_tcga` data_mrna, ideally the
same source used for the original 135-node binarisation) so list-A genes binarise consistently across
all three cohorts.
