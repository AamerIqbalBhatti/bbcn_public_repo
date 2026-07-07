# Paper 3 Repair Plan

Paper 3 is **paper 1, repaired**. It keeps paper 1's IMRAD skeleton, introduction,
cohorts, controller, and switch. It folds in the three biological repairs and replaces
paper 1's negative results with the repaired results. It is **not** built on paper 2;
paper 2's findings enter only as the corrections listed below. Nothing here is built
until you approve this plan.

All repaired numbers below are macro-sourced from `results/numbers2/numbers2.tex`
(regenerated this build; E36/E38/E40 full-N reproduced exactly, E41 re-verified).

---

## Part A. The three repairs (the wiring differences)

These are the only mechanistic changes. Each is one collapsed feedback loop restored to
textbook biology. `pathways.py` is untouched; the repairs apply as an overlay
(`repaired_branch.apply()`), so the bare network of paper 1 is preserved intact.

**R1. p53-MDM2-ATM damage sensor.**
Bare rules: `MDM2 = TP53 or MDM2` (self-locks high), `TP53 = not MDM2` (no damage input).
Result: TP53 is identically 0 whatever the damage. Verified live: under sustained ATM/ATR,
bare TP53 = 0.
Repair: `MDM2 = (TP53 or AKT1) and not(CDKN2A and (E2F1 or ATM or ATR))`,
`TP53 = (not MDM2) or ATM or ATR`. Now damage reaches p53; verified live: repaired TP53 = 1.

**R2. AKT1-FOXO3a-PHLPP loop and the uncoupling flag U.**
Bare: the AKT-FOXO arm has no closing feedback.
Repair: add PHLPP (the 136th node), the phosphatase that closes the loop, plus a per-patient
flag U that marks cells where AKT activity and nuclear FOXO co-occur (the resistant state).
U is computed from each patient's binarised state, not assumed.

**R3. Death-engaged commitment integrator (COMMIT latch).**
Repair: an integrator that accumulates only when death is winning (TP53 high and AKT down
together) and, on crossing threshold, latches a commitment that releases the survival brakes
(MCL1, XIAP, CFLAR, BCL2, BCL2L1). Integrating raw p53 instead breaks bistability, because the
uncoupled resistant state carries high p53 yet survives; the death-engaged form preserves it.

Node count: **135 (bare, paper 1) -> 136 (repaired, paper 3)**; the staged-controller branch
additionally carries COMMIT and U.

---

## Part B. The claim-level corrections (paper 1 -> paper 3)

| # | Paper 1 claim (as posted) | Repaired result | Edit to paper 1 |
|---|---------------------------|-----------------|-----------------|
| C1 | The full network has **no durable cell-fate fixed point**; under synchronous updating it is globally oscillatory, so durability is **relocated to the reduced 9-node switch**. | Repaired staged controller recovers **durable apoptotic fixed points in 15-17%** of patients (TCGA 17, METABRIC 15, I-SPY2 16), where the bare network had **0%**. | Replace the "no durable fixed point / relocate to switch" conclusion with "non-durability was a wiring artifact; repair recovers durable apoptosis." |
| C2 | The resistant attractor is held by hysteresis: **survival-axis inhibition alone cannot flip it; only genotoxic input can.** | **Reversed.** A designed PI3K/AKT/mTOR (survival-axis) kernel **commits 92-95%** of resistant patients; **genotoxic is the weak lever (3-4%)**. | Rewrite this sentence. Survival-axis drug works once the apoptotic execution machinery is repaired; genotoxic is weak. |
| C3 | Survival-routed patients are **89-95% uncoupled vs 26-37% cohort-wide** (a falsifiable prognostic *alignment*, i.e. correlational). | **Deterministic equivalence:** drug resistance **==** AKT-FOXO3a uncoupling for **4049/4050 patients (99.98%)**. | Upgrade the prognostic alignment to a patient-by-patient equivalence. |
| C4 | Individual phenotypes controllable (apoptosis 81-86%, proliferation 94-95%) but **all three jointly only 2-3%**. | Repaired staged controller: apoptosis 72-78%, proliferation 84-87%, **joint 7-14%**, durable 15-17%. | Update the joint-controllability numbers; keep the "feasibility != reachability" framing but note the repair narrows the gap. |
| C5 | Network is **135 nodes**. | Repaired network is **136 nodes** (adds PHLPP). | Update node count and the model schematic (Figure 1). |

Two readouts are kept explicit and honest (this is not a correction, it is preserved from
paper 2's care): **strict** durability (CASP3 on AND survival kinase stays off) reads 8-14%;
**commitment** durability (CASP3 sustained = irreversible caspase activation) reads 92-95%.
The two answer different questions; paper 3 reports both, as paper 2 did.

---

## Part C. Section-by-section edit map for paper 3

Paper 1 file: `preprint/BBCN_preprint.tex`. Paper 3 inherits it and edits:

1. **Title / Abstract.** New title signalling the repair. Abstract: keep setup and cohorts;
   replace the negative half (C1, C2) with the repaired result; fold in C3 as equivalence.
2. **Introduction.** Unchanged through the Boolean-modelling background. Add one paragraph at
   the end of the gap analysis: the bare model's pessimism is traced to three collapsed loops.
3. **Methods.** Keep network, cohorts, binarisation, three-tier controller, switch. **Add 2.x
   "Three biological repairs"** (R1-R3) and the 135->136 node statement.
4. **Results.** Replace the "individually feasible, not jointly reachable / durability only in
   the switch" results with: (i) resistance == uncoupling (C3), (ii) designed kernel commits
   ~90% (C2), (iii) staged controller recovers durable fixed points (C1, C4). Keep the
   strict-vs-commit dual readout.
5. **Discussion.** Reframe: feasibility-vs-reachability gap is real on the bare net but largely
   a wiring artifact; the repair is textbook biology, not tuning. Keep the honest limitation
   (Boolean attractor != pathologic complete response).
6. **Appendix.** Add the repaired update rules (already in paper 2's Appendix A); keep paper 1's
   rule appendix for the bare network so both are documented.
7. **Numbers.** All figures macro-sourced from `numbers2.tex` so paper 3 is drift-proof.

---

## Part D. What stays the same

The model's standing as a structural stratifier and drug-target nominator, the PI3K/AKT/mTOR
nomination, the three cohorts and N (1082/1980/988), the controller design, and the honest
limitation that the model does not predict pathologic complete response. Paper 3 strengthens
paper 1; it does not discard it.
