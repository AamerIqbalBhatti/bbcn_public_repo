# BBCN — Submission kit

Companion to the preprint *The Boolean Breast-Cancer Network (BBCN): structural
controllability of cell-fate signalling and a bistable resistance–apoptosis switch.*
Post the preprint on bioRxiv first; submit the journal version separately. The
cover note names the target journal and is never part of the posted preprint.

---

## 1. bioRxiv portal fields

- Category: Systems Biology (optional cross-list: Cancer Biology)
- License: CC BY 4.0
- Corresponding author: Aamer Iqbal Bhatti, KFUPM (ORCID 0000-0002-3373-3388); enter institutional email in the portal
- Abstract field: paste the paper abstract unchanged (bioRxiv has no length limit)
- Manuscript file: BBCN_preprint.pdf

---

## 2. Abstract trimmed to 200 words (for the npj version)

We model a breast tumour as a Boolean network of signalling pathways and pose
therapy as control: find minimal, druggable interventions that drive the network
to a target cell fate as a genuine fixed point of the free dynamics, without
forcing any node. First, a 135-node network under a three-tier capped controller:
single phenotypes are highly controllable (apoptosis 81–86%, proliferation 94–95%
across three cohorts), yet reaching all three in sequence while preserving earlier
gains is rare (2–3%), because the phenotypes share a dominant strongly-connected
core (16 of 22 pathways). Feasibility does not imply reachability. Second, the
apoptosis core is re-modelled as a bistable resistance–apoptosis switch on
AKT1–TP53 antagonism, a delay-difference (Boolean-ARMA) system over GF(2) on a
multi-timescale schedule. Across TCGA-BRCA, METABRIC and I-SPY2 (N=1082/1980/988)
it routes about half of patients to apoptosis and a quarter to survival;
survival-routed patients are 89–95% in the uncoupled AKT–FOXO3a resistant state
versus 26–37% cohort-wide, a falsifiable alignment. The resistant state is held by
hysteresis and flipped only by genotoxic input. Nominated targets concentrate on
PI3K/AKT/mTOR, robust across three data front-ends. The model is a structural
stratifier and drug-target nominator, not a response predictor.

---

## 3. Cover note (journal submission)

We submit "The Boolean Breast-Cancer Network (BBCN): structural controllability of
cell-fate signalling and a bistable resistance–apoptosis switch" for
consideration. The work treats breast-cancer therapy as a control problem on a
literature-derived Boolean network and makes two linked contributions. Across
three independent patient cohorts it shows that individual cell-fate phenotypes
are controllable while their joint control is obstructed by a dominant
strongly-connected core, that is, feasibility without reachability; and it
resolves the apoptosis bottleneck by re-modelling the core as a bistable AKT1–TP53
switch whose hysteresis explains residual resistance and yields PI3K/AKT/mTOR-
focused druggable targets that are robust across expression, signature, and
phospho-protein readouts. We are explicit about scope: the model is a structural
stratifier and drug-target nominator, validated against curated drug-target
evidence and cohort data, with no survival or pathologic-complete-response claim.
Every number, figure, and table regenerates from a public repository with one
command, and a preprint is posted on bioRxiv. We believe the combination of
semi-tensor-product control theory with patient-resolved oncology fits the
journal's scope.
