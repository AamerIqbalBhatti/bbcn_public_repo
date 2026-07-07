# BBCN — Decisions and Roadmap

**Date: 2026-06-22 (Monday).** Author: Aamer Iqbal Bhatti (KFUPM).
Reference companion to the daily ledger 2026-06-22. This is the clean, human-readable
record of the settled conventions and the implementation plan. No code, no paper.

---

## 1. Naming convention (frozen 2026-06-22)

Every experiment is written **Model / Design / Test**, with a short-form code.

**Axis 1 — Model** (the current network).
- Mechanism: `B` baseline (original literature rules) · `R` repaired (three loops restored) · `C` complete (no abstraction; real inducer and caspase nodes).
- Dynamics: `S` synchronous · `D` delayed (Boolean Delay Equations; each node read at its own lag).

**Axis 2 — Design** (built on the Axis-1 model).
- Pathway handling: `1` single by max-mismatch · `X` CDE (Clinical Decision Engine; coordinated, mode- and phase-aware) · `F` free (no control).
- Kernel method: `r` ranked (greedy forward selection, no L) · `s` stabilising (STP forward-stabilisation, builds L; Rafimanzelat).

**Axis 3 — Test**.
- Test model: `FN` full network · `AS` apoptosis switch · `PS` proliferation switch (may differ from the Axis-1 model; may be older/reduced).
- Data: `T` TCGA-BRCA · `M` METABRIC · `I` I-SPY2 · plus cell-line, clinical chart (future).
- Statistic: `st` structural · `hz` hazard (Cox/KM/log-rank) · `eq` equivalence · `cm` commitment.

Notes: "closed-loop" is not a label (nearly all runs are closed-loop). The biology axis has no
non-biological option by design; `C` is the target.

### Dated experiment table

| Exp | Short form | Date | Where the data lives |
|---|---|---|---|
| Paper 1, Setup A | `BS / Xs / FN·TMI·st+hz` | 2026-06-04 to 06-18 (posted 06-18) | results/numbers/, verification record §1-91, preprint/ |
| Paper 1, Setup B | `BD / Xs / AS·TMI·hz` | 2026-06-04 to 06-18 | results/numbers/, setup_b/, preprint figures |
| E36 | `RD / F / FN·TMI·eq` | 2026-06-18 | ledger 06-18; Drive Take-Home Jun 18; results/numbers2/ |
| E37 | `RD / F / FN·TMI·st` | 2026-06-18 | ledger 06-18; Drive Take-Home Jun 18 |
| E38 | `RD / Xs / FN·TMI·st` | 2026-06-18 | ledger 06-18; results/numbers2/ |
| E40 | `RD / Xs / AS+FN·TMI·cm` | 2026-06-18 | ledger 06-18; results/numbers2/ |
| E41 | `RS / Xs / FN·TMI·st` | 2026-06-18 | ledger 06-18; results/numbers2/ |
| Paper 3 + v22 build | (deliverable, not an experiment) | 2026-06-18/19 | paper3/, results/numbers2/, ledger 06-19 |

---

## 2. Settled decisions (2026-06-22)

**Clinical clock.** Anchored to capivasertib (the AKT inhibitor our kernel nominates),
dosed 4 days on / 3 days off weekly. One binary step = one day. One cycle = one week (4 on, 3 off).
One exposure = one day = one update. Per phenotype: 4 induction cycles + 2 maintenance cycles =
6 treatment weeks, then a separate **4-week fully drug-free observation window** for durability,
sized to the lag-25 slow tier. Total 10 weeks per phenotype. Durability is read at the end of the
washout. The pulsed (withdrawal) controller is the clinically faithful arm.

**Pathway counts per phase.** Induction cap 4, maintenance cap 2; the CDE escalates from 1 on
failure and prefers the fewest pathways that reach the phenotype. (New feature; not yet in code.)

**Apoptosis biology — complete mechanism (point 7).** Retire the COMMIT flag, the brake-zeroing,
and the theta=10 counter. Add 5 real nodes: SMAC/DIABLO (⊣XIAP), NOXA (⊣MCL1), PUMA (with BIM,
BAD ⊣ BCL2, BCL-xL), CASP6 and CASP7 (executioner feedback cascade). The five brakes fall via
named antagonists; commitment emerges from the caspase positive-feedback loop. No MOMP node
(BAX or BAK active = MOMP). U becomes a composite node with a rule, recomputed each step.

**Stimulation input.** Revive the unused MATLAB stimulation channel. Intrinsic quantities are
NODES derived from the cohort (U; baseline ATM/ATR from genomics). Only exogenous treatment
(genotoxic agent, drug) is an input, routed through the stimulation channel.

**Engine reconciliation.** One definition each for CLAMP_OFF (mean-threshold, the renamed U) and commitment (the caspase
cascade), applied identically across the switch, repaired-branch, and honest engines.

**Risk rule (parallel track).** The caspase cascade is biologically correct but not numerically
guaranteed. Keep the working flag version; build the biological version beside it; retire the flag
only after the numbers match.

---

## 3. Acceptance bar (must be preserved within tolerance)

resistance == uncoupling 99.98% (4049/4050) · commitment 92-95% · strict 8-14% ·
durable fixed points 15-17% (bare 0) · apoptosis 72-78% · proliferation 84-87%.
Source of truth: results/numbers2/numbers2.tex.

---

## 4. Implementation roadmap (19 steps, 6 phases)

Governing rule: **one core spot at a time — freeze, run the E40/E41 regression, verify against the
acceptance bar, commit to GitHub, then the next. No bulk changes.**

**Phase 0 — Safety net**
1. Tag current v22 as the frozen baseline; record the acceptance-bar numbers.
2. Add a one-command regression check (E40 + E41 on a fixed small sample). Commit.

**Phase 1 — Convention**
3. Commit CONVENTION.md (the glossary + dated table above). Documentation only.

**Phase 2 — Engine reconciliation**
4. Unify U definition (mean-threshold) across repaired_branch and honest_bbcn.
5. [DONE] U renamed CLAMP_OFF (named 14-3-3 clamp failure), grounded in cited JNK/MST1 biology, held constant (NOT recomputed: dynamic U is circular). Defined by baseline effect.
6. Merge kernel design into one shared function (harness + controller).

**Phase 2.5 — Integrity fixes (found during the redesign audit, E62)**
A. Complete `_pathway_externals` via access-tracking so the kernel cache key captures every
   external that enters L (removes stale-kernel reuse). PREREQUISITE for the clock.
B. AKT1 self-scoring: AKT1 is pinned (real AKT inhibitor) but is also a readout; read phenotype
   verdicts from downstream UNPINNED nodes, never crediting a clamped node. (Measurement decision.)

**Phase 3 — Clinical clock**
7. One binary step = one day.
8. Pulsed schedule: 4 on, 3 off, one-week cycle, as default.
9. Cycles per phenotype = 4 induction + 2 maintenance; add the 4-week drug-free observation; read durability at its end.
10. Add induction/maintenance pathway caps (4 / 2) to the CDE; escalate from 1.

**Phase 4 — Stimulation input**
11. Add the treatment-input vector; route genotoxic through it; init baseline ATM/ATR from genomics.

**Phase 5 — Apoptosis biology (parallel track; flag stays live through step 16)**  [CLOSED 2026-06-24, E72-E74]
12. [DONE, E71] Add SMAC (= BAX or BAK1); XIAP depends on not SMAC.
13. [DONE, E71] Add NOXA (= TP53); MCL1 depends on not NOXA.
14. [DONE, E71] Add PUMA (= TP53); BCL2, BCL-xL depend on not (PUMA or BIM or BAD).
15. [DONE, E72] Add CASP6, CASP7; wire the executioner feedback cascade. (Highest-risk step.)
    Guarded; biologically the point of no return; numerically inert (death state already a stable attractor).
16. [DONE, E74] Turn OFF the COMMIT flag and brake-zeroing; commit only via the real cascade. Decisive
    comparison PASSED (full N: ideal flag 95/cascade 91; clock flag 85/cascade 91, cascade wins the clock).
    Baseline RE-LOCKED to the cascade (flag baseline backed up). Irreversibility PROVEN = the AKT1-FOXO3-PHLPP
    switch, not an apoptosis latch (durable 100% coupled -> 41% uncoupled): see E73.
17. [DONE in the locked path, E74] Flag retired as the commitment mechanism (brake-zeroing removed) in the
    regression, durability_clock, and run_durability_full. RESIDUAL -> CONSOLIDATION (separate task): switch
    the same one line in the remaining duplicate/other-engine runners and DELETE the dead flag scaffolding
    (commit_signal, BRAKES, theta/cmax). Lands cleanly when the engines merge (step 18). See E74 for file list.

**Phase 6 — Consolidation**
18. Reconcile switch / repaired-branch / honest onto the one finalised model; run full-N; re-lock numbers2.
19. Update the ledger and the convention table with the final numbers and node count.

---

## 5. Open items

- All code changes await the go for Phase 0.
- The caspase-cascade commitment risk is resolved only by passing step 16.
- Future data substrates (cell lines, hospital charts) to be added to Axis 3 when available.

---
Phase 6 consolidation backlog -- NEW NODES TO INTRODUCE PROPERLY: see docs/consolidation_new_nodes.md
 List A (data-seedable, formalise as proper nodes): PUMA(BBC3), NOXA(PMAIP1), SMAC(DIABLO),
   CASP6, CASP7, PHLPP(PHLPP1/2). All present in I-SPY2 + METABRIC expression; TCGA expression file pending.
 List B (not measurable, do not formalise): COMMIT (retire when cascade owns commitment),
   CLAMP_OFF/U (keep as derived input; cause not measured).
