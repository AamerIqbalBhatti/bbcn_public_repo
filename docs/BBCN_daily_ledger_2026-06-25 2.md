# BBCN — Daily Ledger — 2026-06-25 (Thursday)

Continues `docs/BBCN_daily_ledger_2026-06-22.md` (entries through E74). Phase 6 (Consolidation)
opened today. Timestamping policy (in effect from 2026-06-24): each task close is stamped with
date and UTC time where the clock is available; entry numbers always.

---

## Session execution log (2026-06-25)

### E75. REPO DE-SPACING + SELF-RUNNABLE REGRESSION — closed 02:48 UTC

Problem found at session open: the repo shipped both spaced (latest) and underscore (stale)
copies of the core modules. Python imports the underscore names, so `python regression_check.py`
could not run as committed — it loaded the stale copies and then died on a missing underscore
`durability_clock`.

Root cause was deeper than a rename: two engine lineages were live under the single import name
`honest_bbcn`. The cascade lineage uses `patient_clamp_off` (strict majority: AKT >=3/4 on AND
both FOXO on); the older flag lineage uses `patient_U` (any AKT on AND any FOXO on). These are
**not the same function**, so aliasing one to the other would have silently changed every old
script's numbers. AIB confirmed the direction: `patient_U` ("capital U") was replaced by
`patient_clamp_off` yesterday and is being retired; the majority definition is correct; the
spaced files are the latest.

Fix (plumbing only, no logic changed): promoted the spaced/latest modules to the importable
underscore names and removed the stale shadows, for the regression chain — `honest_bbcn`,
`clinical_clock`, `durability_clock`, `repaired_branch`, `run_cde`, `regression_check`, and
`docs/regression_baseline.json`. Verified: `python regression_check.py` now runs plainly with no
launcher and **PASSES**, every raw count identical to baseline (E40 durability, E41 monolith,
E42 clock). Convention frozen: underscores, no spaces, from here on.

Deferred (gradual, per AIB): the ~8 `patient_U` runner scripts now import the cascade engine and
will error on `patient_U`; retire/migrate them later (tracked under C3). The setup_a-level orphan
spaced files and the spaced docs were left for a later cleanup pass.

### E76. C1a — LIST-A DATA FOUNDATION BINARISED + VERIFIED (all three cohorts)

Blocking input for C1 arrived: raw expression for all three cohorts.
- TCGA: `data_mrna_seq_v2_rsem.txt` (RSEM, 20532 genes x 1082 samples).
- ISPY2: `GSE194040 ... AgilentGeneExp ... geneLevel n988` (Agilent, 988 samples).
- METABRIC: `brca_metabric` tarball -> `data_mrna_illumina_microarray.txt` (microarray, 1980 samples).

Binarisation = per-gene **cohort-median split** (strict `>`), reproducing `mrna_init.py` exactly.

Knife-edge verification (method must reproduce the committed pipeline before trusting it on new
genes): median-split on each uploaded file reproduces the **committed** binarisation —
- TCGA 133/133 mRNA-derived columns (KMT2D, CXCL8 absent, as the known missing-defaults),
- ISPY2 135/135,
- METABRIC 135/135 (after collapsing 218 duplicate gene rows by mean — which the original
  pipeline also did, since the reproduction is exact).
So the uploaded files and the method are the same source as the original 135-node binarisation.

Six list-A genes binarised, no NaNs, aligned to each committed cohort index:
- BBC3 -> PUMA, PMAIP1 -> NOXA, DIABLO -> SMAC, CASP6, CASP7 (each ~50% ON by construction of the
  median split).
- PHLPP node rule chosen = **OR**: `PHLPP = PHLPP1 OR PHLPP2`. A patient is PHLPP-deficient only
  when both isoforms are below median. Real OR-on fractions: TCGA 72.3% / ISPY2 66.3% /
  METABRIC 72.9%. Isoform correlation weak (+0.08 to +0.35), so OR pools independent information
  rather than double-counting one signal.

Files written: `setup_a/data/binarized/{coh}_{N}x7_listA_genes.csv`,
`{coh}_{N}x6_listA_nodes.csv`, and the staged merged `{coh}_{N}x141.csv` (135 originals + 6 new;
originals untouched, verified by collision + index-order checks).

Status: **C1a (data) DONE**. C1b (network formalization + PHLPP seeding) deferred — see Future Work.

### E77. THE QUARTET — why PHLPP, and the AKT / p53 / FOXO3 / PHLPP link to durability

Recorded at AIB's request as the standing rationale for the deferred PHLPP-seeding work.
Reconstructed from the verification record S37–S39 and ledger E73 (numbers verified from code,
not asserted).

**Why PHLPP exists (verification S39).** AIB's question: p53 has one dominant off-switch, MDM2,
so is there an "MDM2-of-AKT" — a single node that turns the AKT hub off? Answer: PHLPP, not PTEN.
PHLPP directly dephosphorylates AKT at Ser473, shutting the hub off (Gao 2005 Mol Cell 18:13;
Brognard 2007 Mol Cell 25:917). PTEN only removes the upstream PIP3 signal, an indirect analog.
PHLPP is the faithful direct off-switch, the true MDM2-of-AKT. Added as the 136th node.

**The quartet / bistable switch (verification S37).** The core is a p53<->AKT mutual antagonism:
a double-negative, hence positive feedback, hence bistable (PLoS One PMC2634840). AKT -> MDM2 -|
p53; p53 -| AKT (via PTEN and via stabilization). Two stable attractors: PRO-SURVIVAL
(low p53 / high AKT = the resistant attractor) and PRO-APOPTOTIC (high p53 / low AKT). FOXO3 and
PHLPP lock the switch: AKT -| FOXO3; FOXO3 -> PHLPP (FoxO3a transcribes PHLPP1; FOXO promoter
sites, PubMed 29775231); PHLPP -| AKT. When AKT drops, FOXO3 rises into the nucleus, FOXO3
induces PHLPP, PHLPP drives AKT lower — self-reinforcing. PHLPP also stabilizes p53 through the
same AKT -> MDM2 edge (less AKT -> less MDM2 -> p53 up; O'Neill 2012 FEBS J 280:572), not as a
direct p53 target.

**Durability = coupling (ledger E73, proven from code).** Model rules:
`PHLPP = FOXO3 and not CLAMP_OFF`; AKT1 is held off only while PHLPP holds it.
- Clamp intact (coupled): the kernel knocks AKT to 0 -> FOXO3 rises -> PHLPP on -> AKT stays off
  even after the drug clears and the PI3K nodes recover (the AKT rule becomes monotone) =
  DURABLE. AKT off -> MDM2 falls -> p53 rises -> the cascade reads out apoptosis.
- Clamp failed (uncoupled = resistance): PHLPP is forced off regardless of FOXO3, so FOXO3 can
  sit nuclear and high but cannot actuate AKT suppression; AKT recovers = NON-DURABLE.

Proof (TCGA, kerneled n=78): PHLPP intact -> durable 100%, AKT-on-at-window-end 0%; PHLPP
uncoupled -> durable 41%, AKT-on 100%. One line: **durability is coupling; PHLPP is the actuator
that turns FOXO3-nuclear into a permanent AKT-off state; resistance is uncoupling, where that
actuator is disabled.**

**Decision — PHLPP seeding (deferred).** Promote PHLPP to *initialized* (AIB terminology; never
"measured") via **Option C**, gating the feedback by expression:
`PHLPP = PHLPP_initialized and FOXO3 and not CLAMP_OFF`, where `PHLPP_initialized` is the
per-patient baseline bit (PHLPP1 OR PHLPP2), held constant like CLAMP_OFF. Option C preserves the
E73 mechanism for PHLPP-expressing patients and adds a **second, independent resistance axis**:
a PHLPP-deficient patient (~30% of the cohort) loses the actuator and becomes non-durable even
when coupled.

Predicted effect, to be tested when C1b is implemented: clock-durable drops from ~91% kerneled
toward roughly the 60s–70s, the size of the drop approximately the PHLPP-deficient fraction among
patients who were durable before.

Why deferred: implementing Option C moves the locked numbers **by design**. The regression
baseline would have to re-lock and the verify cycle would reopen. Per AIB, this surgery is **not**
done today, so Phase 6 can close at frozen numbers for the team handoff. The three engine rules
(`honest_bbcn`, `bbcn_switch`, `repaired_branch`) are unchanged today; the x141 data sits ready
and unused, moving no number.

---

### E78. STEP 18 RESOLVED BY DESIGNATION — canonical engine crowned; C4 landed — 05:33 UTC

Engine reconciliation closed without a merge. Audit finding: the three engines already agree on
the core switch (MDM2, TP53, AKT1, FOXO3, PHLPP identical across honest, repaired-branch, and the
9-node switch). The only divergences were the CLAMP_OFF/U name (C4), the commitment mechanism
(honest=cascade, repaired=COMMIT flag, switch=neither — C3), and the genotoxic channel.

CANONICAL DECISION (AIB). The honest engine is crowned canonical and default, named **Delayed BBCN
(BBCN-D)**. It is the superset: delayed / Boolean-ARMA timescales, the p53-MDM2-ATM repair, the
AKT1-FOXO3-PHLPP switch, CLAMP_OFF, cascade commitment, and the genotoxic input. In the Model axis
of the convention it is the C/D (complete, delayed) model. Everything from here defaults to BBCN-D.

The other two are kept as **frozen, reference-only** baselines — not maintained, not re-locked,
fallback only if ever needed:
- **BBCN-M** — the repaired monolith (still runs the COMMIT flag).
- **BBCN-S** — the reduced 9-node analytical switch.

Why designation, not merge: crowning BBCN-D requires no edit to BBCN-M or BBCN-S, so there is no
flag retirement in repaired (C3 not triggered), no re-lock, and the regression baseline stays
frozen for the team handoff. Step 18 thus collapses from "merge three engines" to "declare the
canonical one and freeze the references." The single-durability-function merge is no longer needed.

FILENAME: kept as `honest_bbcn.py`. NOT renamed to `delayed_bbcn.py`, because that name is already
the generic delay-difference lag engine (a different module), and ~20 files import `honest_bbcn`.
"Delayed BBCN / BBCN-D" is the model designation; the file stays put to avoid churn. A physical
rename, if ever wanted, is its own tracked refactor: pick a non-colliding name, update all
importers, regression-verify.

C4 LANDED (number-safe, verified). The uncoupling bit is now named `CLAMP_OFF` consistently: the
switch's `U` parameter renamed to `CLAMP_OFF` in `bbcn_switch.py` (rules/simulate/fixed_points/
build_box_matrix) and its two keyword callers in `reproduce_all.py` updated; honest and repaired
already shared the definition (the majority rule). Regression PASSES, every raw count identical.

STEP 18 STATUS: core confirmed unified; canonical engine designated (BBCN-D); C4 done. C3 (flag
retirement in BBCN-M) stays deferred — it moves the E41 numbers and is no longer required, since
BBCN-M is now reference-only. Step 18 is effectively closed by designation.

### E79. CONTROL ARCHITECTURE — two independent tiers, defaults vs benchmarks

Clarified the final control architecture as two independent choices made in sequence.

UPPER TIER — pathway selection (which pathway to act on):
- **CDE** (adaptive; clinical decision engine, `run_cde` / `select_active_pathway_set`) = DEFAULT.
- **Static** (max-mismatch; `controller.control_phenotype` via `_pathway_mismatch`, `rank='mismatch'`)
  = BENCHMARK.
Currently these live in two separate drivers, not one switch (the gap C9 fills).

LOWER TIER — kernel selection (which nodes to pin):
- **Stabilization** (algebraic Theorem-1, `forward_stab`) = DEFAULT.
- **Greedy Search** (heuristic add-one-by-rank; the code identifier is `method='ranked'` /
  `kernel_method='ranked'`; "Greedy Search" is its name, NOT "reducer") = BENCHMARK.
This tier ALREADY has a clean switch (`kernel_method='stabilize'` default | `'ranked'`).

The two tiers are independent, giving four configurations: CDE+Stabilization (= canonical
**Delayed BBCN**, the default, frozen numbers), and CDE+Greedy / Static+Stabilization /
Static+Greedy (the three benchmarks, comparison runs only).

NAMING: the ranked kernel method is **Greedy Search**. The code literal stays `'ranked'` for now
(renaming it touches many `kernel_method='ranked'` callers; a minor future cleanup, not done yet).

NEW TASK C9 (number-safe, Phase 6): unify the two pathway drivers into one entry point with
`pathway_mode='cde'` (default) | `'static'`, mirroring the kernel switch. Default path
(CDE + Stabilization) stays byte-identical, so the regression baseline stays frozen; Static
becomes a labelled benchmark. Slotted after C5 and C6 (cleaner once C5 has consolidated the kernel
tier underneath).

## State at end of 2026-06-25

- Phases 0–5 closed (17/19 roadmap steps). Phase 6 (Consolidation) opened.
- Repo is self-runnable (de-spaced); `python regression_check.py` PASSES at the frozen cascade
  baseline.
- C1a data foundation complete and verified across all 4050 patients. No engine touched, no
  number moved this session.

---

## FUTURE WORK (Phase 6 consolidation, and Phase 7 for the number-moving surgery)

Governing rule unchanged: **one core change at a time — freeze, run the E40/E41/E42 regression,
verify against the acceptance bar, commit, then the next.** Each item below either moves the
locked numbers or is the engine surgery, so each is its own gated change.

**C1 — Formalize the six list-A nodes; seed PHLPP via Option C.**
- C1a (data): DONE this session (E76). All three cohorts binarised and verified; x141 staged.
- C1b (network + seeding): PENDING. Promote PUMA, NOXA, SMAC, CASP6, CASP7, PHLPP from inline
  overlay / derived flag to first-class `PATHWAYS` nodes — each with its own rule, lag, bus
  persistence, and an x141-seeded initial condition; grow `ALL_NODES` 135 -> 141; retire COMMIT.
  Seed PHLPP via Option C identically in all three engines (`honest_bbcn`, `bbcn_switch`,
  `repaired_branch`), threading `PHLPP_initialized` from the x141 column. **Moves durability by
  design** (see the E77 prediction); re-lock the baseline only after verifying actual-vs-predicted.
  This is the number-moving surgery -> Phase 7.

**C2 — DRUGGABLE_NODES revision.** Drop PTEN; add MDM2, BCL2, BCL2L1, MCL1, XIAP. Re-capture E40
and E42. Test whether apoptosis controllability rises toward the proliferation ceiling. Moves
numbers (re-capture), so gated.

**C3 — Retire the dead flag scaffolding.** Propagate the one-line cascade switch to the remaining
duplicate / other-engine runners (`fulln_e40`, `fulln_clock`, `fulln_frontload`,
`run_designed_kernel_AB`, `run_designed_kernel_cohorts`, `run_full_cohorts_repaired`,
`run_unified_repaired`); then delete the now-dead scaffolding (`honest_bbcn.commit_signal`,
`BRAKES`, the theta/cmax escalation timer). Touches the regression's own flag-comparison columns,
so it lands cleanly only when the engines merge (Step 18).

**C4 — small-u (CLAMP_OFF) input unification.** Unify how the uncoupling input enters the honest
versus repaired engines — one definition, one entry point — per the engine-reconciliation decision.

**C5 — Merge the kernel-design copies.** Collapse the duplicated kernel-design logic in `harness`
and `controller` into one shared function. This is the deferred roadmap Step 6.

**C6 — Preprint figures.** Instrument per-patient kernel capture (record the full kernel
composition for each patient), then regenerate the TCBB-style kernel heatmaps and the
pathway-burden plots at full N. Needs a full-N run.

**C7 — Deliver the construction ledger to Drive.** Ship the construction ledger (through E77) to
the Drive take-home folder. Pure delivery; the one item with no number risk, doable any time.

**C8 — Repo-wide de-spacing.** Finish the convention sweep: rename the remaining spaced files the
chain fix did not need — the spaced docs (ledgers, roadmap, checkpoints) and the setup_a-level
orphan scripts (`honest bbcn.py`, `repaired branch.py`, `run durability full.py`,
`run repaired monolith.py`, `regression baseline.json`, and the `fulln *.py` set) — to underscore
names, and delete duplicates. No regression impact (these are not on the import chain); pure
hygiene to enforce the no-spaces convention everywhere. Note: the canonical engine file stays
`honest_bbcn.py` (designation BBCN-D); it is NOT renamed to `delayed_bbcn.py` (name already taken
by the lag engine, plus ~20 importers).

**C9 — Finalize the two-tier control opt (number-safe).** Unify the two pathway-selection drivers
(CDE in `run_cde`, Static/max-mismatch in `controller`) into a single entry point with
`pathway_mode='cde'` (default) | `'static'`, the same way the kernel tier already switches
`stabilize` | `'ranked'` (Greedy Search). Default = CDE + Stabilization = the frozen canonical, so
the regression stays locked; Static becomes a labelled benchmark. Slotted after C5 and C6.

**Step 18 — Engine reconciliation.** RESOLVED BY DESIGNATION (E78), number-safe. Honest crowned
canonical and default as **Delayed BBCN (BBCN-D)**; repaired (**BBCN-M**) and the 9-node switch
(**BBCN-S**) frozen as reference-only. C4 (CLAMP_OFF naming) landed. C3 (flag retirement in BBCN-M)
no longer required, since BBCN-M is reference-only — left deferred. The durability-function merge is
moot under designation. Remaining tie-in: C5 (kernel-design merge) is independent and still open.

**Step 19 — Documentation re-lock.** Update the ledger and the convention table with the final
numbers and the final node count, after Step 18.

**Phase 7 — Data-seeded node formalization (new).** The home for changes that intentionally move
the locked numbers: C1b (PHLPP Option C and the six-node formalization) and any further
data-seeding. Kept separate from Phase 6 so Phase 6 can close at frozen numbers for the team
handoff, rather than reopening the verify cycle mid-consolidation.
