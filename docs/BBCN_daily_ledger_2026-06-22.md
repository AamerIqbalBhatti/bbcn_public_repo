# BBCN — Daily Ledger — 2026-06-22 (Monday)

Design and stock-taking session. No code was changed and no paper was produced today.
This day was decisions only: a full Python code walkthrough, two integrity corrections,
the frozen naming convention, the clinical time model, the apoptosis-biology plan, and
the 19-step implementation roadmap. Continues the E-series from E49 (2026-06-19).

Pairs with: CONVENTION + roadmap reference (BBCN_decisions_and_roadmap_2026-06-22.md),
the standing verification record (§1-91), and the prior ledgers (2026-06-18, 2026-06-19).

Prior sessions for context: paper 3 and repository v22 were built 2026-06-18/19 and are
ledgered there (E43-E49); they were not rebuilt today.

---

## Session execution log (2026-06-22)

E50. CODE WALKTHROUGH, harness.py spine mapped from the entry point. Three tiers:
run_patient (top, 4 stages in STAGE_SEQUENCE) -> session loop (up to T=5) -> batch loop
(_schedule_batches, 1 pathway/exposure) -> simulate_inner (inner, up to max_steps=10) ->
per step: _resolve_active (1) + _update_pathway (24, one per pathway). Defaults: 4 stages,
T=5 sessions, max_steps=10, max_pathways_per_exposure=1, confirm_window=3. Two early exits:
inner loop breaks on bus fixed point; stage loop breaks on PASS held across the window.
Real compute is the active branch of _update_pathway (the stabilize kernel), cached in
_STAB_CACHE; the other 23 are one rule step each.

E51. INTEGRITY CORRECTION (STP / L-matrix), logged openly. An earlier claim was reversed
after reading the functions. TRUTH from code: the stabilize (Rafimanzelat) method BUILDS L
internally (forward_stab.py:137, build_L), so it cannot run without L. The ranked method
does NOT need L: its forward search uses direct rule iteration whenever rule_fn is supplied
(stp.py:233, sim = sim_direct), and the call sites always supply it; L is used only in the
rule_fn=None matrix fallback. The harness call-site appearance is misleading (the ranked
branch calls generate_L and the stabilize branch passes only the rule closure), which is what
produced the wrong earlier statement. Author was right; correction recorded.

E52. PATHWAY-COUNT AUDIT. There is no induction/maintenance phase pathway count in the code,
neither in harness nor controller. controller.py uses a fixed n_pw=2 pathways with ind_cap=4
(a kernel-NODE cap, not a pathway count). The CDE (cde_proposal) truncates the selected set at
MAX_PATHWAYS=10 by mode, then run_cde applies one pathway per exposure. The phased counts the
author recalled are a NEW feature to add, not an existing setting.

E53. v8 vs v21 SPINE COMPARED. harness.py is byte-identical across v8 and v21; the repair is an
overlay, not a fork. The v21 repaired monolith (E41) uses the SAME spine via run_patient_cde on
the repaired network reached by repaired_branch.apply(), which mutates harness.PATHWAYS in place:
4 rule rewrites (MDM2, TP53, AKT1, FOXO3), 5 brake gatings (XIAP, CFLAR, MCL1, BCL2, BCL2L1),
one new pathway RepairAux={COMMIT,U}, and adds PHLPP -> 25 pathways, 138 nodes. Two further v21
spines do not go through harness: honest_bbcn.course (the delayed free engine, E36 and E40 flips)
and forward_stab_kernel_design (the design-small/test-big kernel, E38 and E40).

E54. ENGINE DISCREPANCIES found, to be reconciled. (1) U: mean-threshold (mean>0.5 of AKT and of
FOXO) in repaired_branch.patient_U vs any() in honest_bbcn.patient_U. (2) COMMIT: immediate
Boolean latch in repaired_branch vs a theta=10 integer accumulator in honest_bbcn. Same intent,
different dynamics; one definition each must be chosen.

E55. NAMING CONVENTION FINALIZED (frozen). Three axes, written Model / Design / Test, with a
short-form code. Model = mechanism {B baseline, R repaired, C complete} + dynamics {S sync,
D delayed}. Design = pathway {1 single by max-mismatch, X CDE, F free} + kernel {r ranked,
s stabilising}. Test = test-model {FN full, AS apoptosis switch, PS proliferation switch} +
data {T M I cohorts, plus cell-line, clinical} + statistic {st structural, hz hazard, eq
equivalence, cm commitment}. CDE = Clinical Decision Engine (coordinated, mode- and phase-aware).
"Closed-loop" dropped as a label (nearly all runs are closed-loop). Full glossary and the dated
experiment table are in the roadmap reference file.

E56. CLINICAL TIME MODEL DECIDED. Anchored to the drug our kernel nominates: capivasertib (AKT
inhibitor), dosed 4 days on / 3 days off, weekly (CAPItello-291 / FAKTION). Decisions: one binary
step = one day; one cycle = one week = 4 on + 3 off; one exposure = one day = one update; cycles
per phenotype = 4 induction + 2 maintenance = 6 treatment weeks; PLUS a separate 4-week fully
drug-free observation window for durability (sized to the lag-25 slow tier), so 10 weeks per
phenotype total. Pathway caps: induction 4, maintenance 2, CDE escalates from 1 on failure. The
clinically faithful controller is the PULSED arm (withdrawal), not the always-held one.

E57. APOPTOSIS BIOLOGY DECIDED (mechanism axis: repaired -> complete). Retire all abstractions:
the COMMIT flag, the brake-zeroing, and the theta=10 counter. Add 5 REAL nodes: SMAC/DIABLO
(neutralises XIAP at MOMP), NOXA/PMAIP1 (degrades MCL1), PUMA/BBC3 (with BIM, BAD neutralises
BCL2, BCL-xL), and the two missing executioner caspases CASP6 and CASP7 (real feedback cascade).
The five brakes then fall via their named antagonists, and commitment emerges from the caspase
positive-feedback loop, not a counter. NO MOMP node needed (BAX or BAK active = MOMP). U becomes
a composite node with a rule, recomputed each step, not a per-patient constant. Inventory checked:
present already are BAX, BAK1, CYCS, APAF1, CASP3, CASP8, CASP9; genuinely missing are SMAC, NOXA,
PUMA, CASP6, CASP7. CASP1/2/4/5/10 deliberately excluded (inflammatory, not apoptotic).

E58. RISK LOGGED + PARALLEL-TRACK RULE ADOPTED. The caspase-cascade commitment is biologically
correct but NOT numerically guaranteed; whether the Boolean loop latches depends on the exact new
rules. Therefore: retain the working flag version, build the biological version beside it
(switchable), and retire the flag ONLY after the biology reproduces the locked numbers (commitment
92-95%, durable fixed points 15-17%). Nothing working is deleted until its replacement is proven.

E59. IMPLEMENTATION ROADMAP LOCKED. 19 steps in 6 phases (full text in the roadmap reference file):
Phase 0 safety net (tag baseline, regression script), Phase 1 convention doc, Phase 2 engine
reconciliation (U, kernel-design merge), Phase 3 clinical clock, Phase 4 stimulation-input channel,
Phase 5 apoptosis biology on a parallel track, Phase 6 consolidation and re-lock. Governing rule,
per author: one core spot at a time, freeze, run the E40/E41 regression, verify against the
acceptance bar, commit to GitHub, then the next. No bulk changes.

---

## Acceptance bar (numbers every later step must preserve, within tolerance)

resistance == AKT-FOXO3a uncoupling 99.98% (4049/4050); commitment durability 92-95%;
strict durability 8-14%; durable apoptotic fixed points 15-17% (bare net 0); apoptosis
controllable 72-78%; proliferation 84-87%. All locked in results/numbers2/numbers2.tex.

## State at end of 2026-06-22

Decisions settled: naming convention (frozen), clinical time model, apoptosis-biology node set,
parallel-track risk rule, 19-step roadmap. Not yet done: any code change (Phase 0 onward awaits
the author's go). Open modelling risk: the caspase-cascade must reproduce commitment, proven step
by step. Engine discrepancies (U, COMMIT) to be unified in Phase 2.

## On the horizon

- Begin Phase 0 (frozen baseline tag + regression script) when the author says go.
- Point 7 (the complete-mechanism apoptosis layer) is Phase 5, on the parallel track.
- Reconcile switch / repaired-branch / honest engines onto one finalised model (Phase 6).

E60. PHASE 2 STEP 4 DONE (U unified). honest_bbcn.patient_U changed from any-of to the
strict-majority mean-threshold, made identical to repaired_branch.patient_U; both rewritten
in explicit majority form (count of on-nodes > half; AKT >=3 of 4, FOXO both), strict > cut
kept. Verified: E41 unchanged (it always used this definition); E40 moved to the accepted
mean-threshold values; regression baseline updated and PASSes. Full-N consequence (compare_U.py,
all cohorts): resistant fraction drops from ~half to ~one-fifth -- TCGA 523->202 (48%->19%),
METABRIC 965->398 (49%->20%), ISPY2 516->238 (52%->24%) -- while the designed kernel still
commits 87-96% of resistant patients (TCGA 96, METABRIC 87, ISPY2 95). Efficacy holds; the
resistant-fraction figure in numbers2/paper must be updated to the mean-threshold values at the
Phase-6 full-N re-lock. Files changed: setup_a/bbcn/honest_bbcn.py, setup_a/bbcn/repaired_branch.py,
docs/regression_baseline.json. Decision recorded openly per knife-edge standard (any-of overcounted
uncoupling via single-node noise; strict majority is the conservative, defensible definition).

E61. PHASE 2 STEP 5 DONE (U named and grounded, not made dynamic). Decision: U is NOT recomputed
each step. Reasons: (a) uncoupling is a fixed patient defect, not a per-step readout; (b) a dynamic
U is circular -- it both feeds and is computed from AKT1/FOXO3. Instead U is RENAMED to CLAMP_OFF
(failure of the AKT->14-3-3->FOXO3 nuclear-export clamp), kept constant, set once at baseline.
Literature verified (Brunet 1999; Tzivion 2011 PMC3237389; Sunayama/Gotoh 2005 PMC2171419; Sunters
2006 breast cancer): oxidative-stress JNK/MST1 phosphorylate 14-3-3 and release FOXO3, so FOXO stays
nuclear despite active AKT. Cohort check (compare_U full-N + MAPK8 co-occurrence): MAPK8/JNK is ON in
50% of every cohort but accompanies the uncoupled effect in only 42-63% of uncoupled patients and
only 8-9% of JNK-on patients are uncoupled -> JNK is neither necessary nor sufficient in the data,
and the true cause (14-3-3 phospho-state) is not measurable in mRNA. So CLAMP_OFF is defined by its
measurable baseline EFFECT (strict-majority AKT and both FOXO), with the JNK/MST1 cause cited but not
used in the definition. Guarded rename across 5 files (repaired_branch.py node+function, honest_bbcn.py
function+parameter, run_repaired_monolith.py, run_durability_full.py, regression_check.py): regression
PASS, every count identical -> pure rename, zero behaviour change. Step 5 thus reframed from "make U
dynamic" to "name, ground, and hold U constant" -- DONE.

E62. TWO INTEGRITY ISSUES FOUND during the kernel-redesign audit (Phase 2/3 boundary).
CONCLUSION FIRST SETTLED: the stabilize kernel depends on L; L = f(rules, externals), NOT x0
(global design works from any internal start). So the kernel must be RE-DERIVED whenever the
externals change (event-driven, not a fixed per-cycle/per-exposure schedule). The externals-keyed
cache already implements this -- IF the external set is complete.
 (A) _pathway_externals IS INCOMPLETE. It discovers externals by a single-flip sensitivity probe
     at internal all-0 / all-1, which is MASKED when an external sits in an AND with another
     external that is 0 in the probe background. Verified by access-tracking over random states:
     real reads are missed, e.g. AKT_Signaling reads PTEN and PIK3CA but declares neither (AKT1
     rule contains 'and not PTEN'); also CYCS/APAF1 (extrinsic), MCL1, COMMIT, GSK3B, etc. CONSEQUENCE:
     the stabilize cache key sk=(pname, ext-values, target) omits these, so two situations differing
     only in (e.g.) PTEN collapse to one key -> a STALE kernel designed for the other PTEN is reused,
     and the cache is shared across patients -> cross-patient contamination. This is a real kernel-
     design correctness bug. FIX (safe direction): replace the probe with access-tracking (over-
     approximate the external set); extra externals only cost cache reuse, never serve a stale kernel.
     PREREQUISITE for Phase 3: the clinical clock relies on correct redesign when externals drift.
 (B) AKT1 SELF-SCORING RISK. Kernel-log audit (40 pts/cohort, 192,141 pins): CASP3/FOXO3/MYC are
     NEVER pinned (clean), but AKT1 is pinned in AKT_Signaling in ~3% of pins (6,119). AKT1 is both
     a legitimate drug target (capivasertib = AKT inhibitor; pinning AKT1=0 IS the drug) AND one of
     the 4 global phenotype readouts. Where a verdict requires AKT1=0 and we got there by pinning
     AKT1, that part of the verdict is circular (self-scoring). The commitment readout (CASP3 sustained,
     never pinned) is the honest measure; the strict readout (also demands AKT1 off) is the circular-
     leaning one -- consistent with why strict<<commit. RESOLUTION (design, not a blind fix): keep AKT1
     as a kernel (it is the real drug), but ensure phenotype verdicts are read from DOWNSTREAM, UNPINNED
     nodes, never crediting a clamped node. To settle as a Phase-3 measurement decision.

E63. PHASE 2.5(A) DONE -- externals completeness fix wired in, guarded, measured.
 MECHANISM: _pathway_externals now discovers externals by access-tracking (logs every bus
   node the rules actually read over 600 seeded-random states, seed 20260622, deterministic)
   UNIONED with the single-flip probe. The tracker is primary and complete: it sees through
   AND-masking (e.g. PTEN/PIK3CA in the AKT1 rule) and through the repair-branch rule wrappers
   that source inspection cannot. The probe is demoted to a deterministic backstop for the
   all-0/all-1 corner. Result stored as an explicit, inspectable PATHWAYS[pname]['externals']
   field (the per-pathway listing the original MATLAB getRules authored as `external` structs).
 AUTHORITY: validated against the uploaded getRules_v1.m -- mTOR and Apoptosis_Intrinsic match
   the MATLAB `external` structs exactly; the AKT split (AKT_Signaling+AKT_Survival) covers the
   MATLAB AKT externals. Read-only set adopted (only what the rules read); MATLAB's declared-but-
   unread names (DUSP1/PI3K, BCL2+BCL2L1+MAPK8/Extrinsic, MAP2K1/Regulatory) and the MAPK1-as-its-
   own-external compatibility quirk are intentionally NOT carried forward (no MATLAB back-compat).
 CACHE HYGIENE: repaired_branch.apply() now clears _PATH_EXTERNALS and _STAB_CACHE and drops stale
   'externals' fields, so the corrected repaired-aware keys take effect when rules switch.
 IMPACT MEASURED (direct probe-only vs complete, E41, 30 pts/cohort): the probe misses real
   externals in 11 pathways, BUT phenotype numbers are IDENTICAL -- 0 raw count differences across
   all cohorts/metrics. CONCLUSION: the bug was REAL but LATENT on this sample (stale reuse drove
   the same outcome; no verdict flipped). Locked numbers were NOT corrupted; regression_check PASSES;
   no re-lock needed. The fix is retained as a PREREQUISITE for Phase 3: once the clock turns the drug
   on/off, externals will vary within a run and an incomplete key would then serve stale kernels.
 STILL OPEN: (B) AKT1 self-scoring (read verdicts off downstream unpinned nodes) -- Phase-3
   measurement decision, not yet actioned.

E64. PHASE 2.5(B) DONE & LOCKED -- AKT1 self-scoring removed via STEER-v1 / SCORE-v2.
 PROBLEM: the kernel pins AKT1 (the real AKT inhibitor) ~3% of pins, and AKT1 was also a
   verdict node, so success gates that required AKT1=0 were self-scoring. A pin audit found
   the only verdict nodes the kernel NEVER pins are CASP3, RELA, STAT3, FOXO3, CCND1, E2F1;
   AKT1, CASP9, MYC, XIAP, CFLAR, MCL1, BAX, CYCS are all pinned at least sometimes (correcting
   an earlier claim that MYC was clean -- completing the externals shifted some pins onto MYC).
 KEY FINDING: the verdict is dual-use -- the CDE loop READS it to STEER (classify_failure ->
   next pathways), and we also read it to SCORE. So changing the schema globally also re-steered
   the controller (apoptosis moved even though its gate did not). Static arm (controller.py) is
   verdict-BLIND in routing (picks pathways by _pathway_mismatch), confirmed: v1->v2 gave 0 change
   there; only the CDE re-steered.
 FIX: steer with v1, score with v2. evaluate_phenotype(bus,stage,scheme) -- loop calls use default
   'v1' (trajectory identical to locked run); the FINAL grade (run_cde.py final_status + terminal_pass)
   uses 'v2'. v2 reads CASP3 (apop), RELA+STAT3 (resistance-off, replacing pinned AKT1), CCND1+E2F1
   (prolif). AKT1 still DRIVEN off (GLOBAL_INVARIANT untouched); only the credit is removed. The three
   steering helpers (failed_required/hit_forbidden/hit_contradictions) stay v1.
 IMPACT (E41, 30/cohort): ONLY proliferation-off moved, 80/83/83% -> 100% (dropped self-scored MYC
   from its forbidden gate; clean CCND1/E2F1 already satisfied). apoptosis, joint, durability, and
   all of E40 UNCHANGED -> trajectory provably untouched. Baseline re-locked to the v2-scored numbers;
   regression PASSES. Durability arm: honest verdict is commit (CASP3 sustained, already reported);
   strict (CASP3 & AKT1=0) retained only as a labeled circular diagnostic.
 2.5(B) CLOSED. Phase 2.5 fully done. Next: Phase 3 (clinical clock), step 7 = one binary step = one day.

E65. VERDICT-READOUT PIN WARNING added (decision: keep numbers, warn for the future).
 The kernel set is chosen at runtime, so no static check can guarantee that a scored node
 is never pinned on some future patient/cohort/target/clock setting -- the exclusivity of
 "nodes the kernel pins" vs "nodes the verdict scores" was only ever a 30-patient snapshot,
 not an invariant. DECISION: do NOT forbid pinning scored nodes (keep kernels as chosen, no
 number change); instead emit a deduplicated warning (once per node) at the pin site in
 harness._update_pathway whenever a pinned node is in VERDICT_READOUT_NODES (= required +
 forbidden of the v2 SCORING schema). Steering-only nodes (AKT1/CASP9/MYC, v1) are excluded
 by design -- pinning the drug target is fine; only SCORED nodes contaminate the grade.
 IMMEDIATE YIELD: on the 30/cohort sample the warning fires for CFLAR, XIAP (apoptosis
 forbidden) and CTNNB1, TWIST1, YAP1 (invasion gate) -- real silent contaminations now on
 the record. None are in the locked headline metrics, so regression still PASSES. FUTURE:
 revisit with a majority-of-pathway success readout (Aamer's idea) so the verdict cannot be
 contaminated by any single pinned node. Logged, not yet actioned. 2.5(B) fully closed.

E66. PHASE 3 STARTED -- step 7: one binary step = one day (the clock anchor).
 Confirmed from code that one synchronous network update is already the atomic tick in BOTH
 engines (harness.simulate_inner `for t in range(1,max_steps+1)`; durability A_dyn
 `for t in range(T1+T2)`). Step 7 names that tick a DAY and encodes the full settled clock in
 one module, setup_a/bbcn/clinical_clock.py:
   1 step=1 day; 1 week=7 days=4 ON+3 OFF (capivasertib / CAPItello-291);
   treatment 6 wk = 4 induction + 2 maintenance; observation 4 wk drug-free; total 10 wk = 70 d;
   pathway caps induction 4 / maintenance 2 / observation 0.
 Pure helpers: day_of_step, week_of_day, day_in_week, phase_of_day, is_on_day, pathway_cap,
 total/treatment/observation_days, schedule_summary. Self-test passes (week-1 ON 1-4/off 5-7;
 induction->maintenance at day 28/29; treatment->observation at 42/43; all 28 obs days off; 70 d).
 INERT: nothing imports it yet, so regression PASSES unchanged. Wiring comes in steps 8-10
 (pulsed dosing / cycles+observation / caps). Anchor in place.

E67. PHASE 3 step 8 -- pulsed (withdrawal) dosing on the clock (faithful arm).
 The continuous-hold durability (run_durability_full.A_dyn, clamp held while t<T1) assumed the
 drug is present every treatment day. capivasertib is 4-on/3-off, so step 8 adds the faithful arm:
 setup_a/durability_clock.py :: A_dyn_clock(b,U,clamp,pulsed). Horizon = clock total 70 days;
 treatment 42 d dosed via CK.is_on_day(day) (4-on/3-off) when pulsed, or held when not; observation
 28 d drug-free in both. Verdict = honest commit (CASP3 sustained over the drug-free tail, E64).
 The one-line mechanism (durability_clock.py:50): dosed = CK.is_on_day(day) if pulsed else day<=42,
 replacing the continuous `if t<T1` hold (run_durability_full.py:42-43).
 RESULT (head-to-head, same 70-day clock, 20/cohort): pulsed == continuous -- TCGA 75%/75%,
 METABRIC 100%/100%, ISPY2 100%/100%. The 3-day off-gaps and 28-day tail do NOT break the kill:
 once commitment latches it is self-sustaining and does not need daily drug. CAVEAT: kerneled
 subset is tiny (8/1/5 = 14 patients); preliminary, needs full-N. ADDITIVE arm -- does not touch
 the locked A_dyn or the baseline; regression PASSES. Adoption as canonical durability = a later
 decision after full-N. Next: step 9 (cycles + observation window wired into the run structure).

E68. PHASE 3 step 9 -- cycles + observation window wired into the run structure.
 setup_a/durability_clock.py :: run_on_clock(b,U,clamp,pulsed). Same dynamics as A_dyn_clock but
 the run is now STRUCTURED on the clock: each week is a cycle (weeks 1-4 induction, 5-6 maintenance,
 7-10 observation/drug-free), every day stamped with its phase (line 80), CASP3 recorded at each
 week end (line 97), and durability read across the full 4-week observation window (line 100), plus
 the day/week/phase at which commitment first latches. Returns a per-week timeline.
 DEMO (TCGA, pulsed): pt3 latches day 18 (wk3 induction), CASP3 holds through maintenance + all 4
 observation weeks -> durable. pt6 noisier, latches day 54 (wk8, IN observation/drug-free), still
 durable -- the late withdrawal-period commitment the observation window exists to catch.
 NOTE: the latch is still the COMMIT theta-counter (flag mechanism), not real apoptotic biology;
 Phase 5 replaces it (SMAC/NOXA/PUMA/CASP6/CASP7). ADDITIVE -- nothing locked imports run_on_clock;
 regression PASSES. Full-N deferred to after Phase 3 (per Aamer). Next: step 10 (induction cap 4 /
 maintenance cap 2 / escalation), the last step of Phase 3.

E69. PHASE 3 step 10 -- caps + escalation (closes Phase 3).
 Clock gains the escalation rule: clinical_clock.escalate_active(active,day) climbs the active
 pathway count by 1 toward the phase cap on a failed cycle; clip_to_cap(active,day) tightens it
 at phase boundaries. Wired into the arm: durability_clock.run_capped_on_clock(b,U,kernels,apop,
 pulsed) doses only the first `active` kernel pathways (line 141), starts at 1, escalates at each
 uncommitted treatment week-end (line 150), and clips down when the phase cap tightens (line 132).
 Induction cap 4, maintenance cap 2, observation 0 (drug-free).
 DEMO (TCGA, pulsed): pt3 escalates 1->2->3, commits wk3 at active=3, maintenance clips to 2, obs
 drug-free, durable. pt6 climbs full ladder 1->2->3->4 in induction, never commits in treatment,
 latches in observation wk8 (drug-free), durable. Escalation stops once committed; caps clip at
 every boundary -- policy behaves as specified.
 ADDITIVE -- clinical_clock + durability_clock only; nothing locked imports run_capped_on_clock;
 regression PASSES. PHASE 3 COMPLETE (steps 7-10). Deferred: full-N run of the clock arms, and the
 decision whether the pulsed/capped arm replaces the continuous A_dyn as canonical durability
 (would re-lock the baseline). Next per Aamer: full end-to-end / full-N check.

E70. FULL-N CHECKPOINT (end of Phase 3) + induction-strategy result. Phase 3 CLOSED.
 Ran the locked engines and the clock arm at full N (4050 patients) -- see docs/fulln_checkpoint_2026-06-23.md.
 E41 monolith (v2): apop 74% / prolif 100% / joint 11% / durable 16% POOLED. Sample held; joint
 settled from the noisy 30/cohort 30/3/20 to a stable 11% at full N.
 E40 durability (idealized 300-step continuous): resistant 21% / kernel 96% / switch 85% / strict 61%
 / commit 95% POOLED. strict = circular diagnostic; commit = honest.
 CLOCK ARM full N (durable over the drug-free observation window): clock-continuous 93% (≈ old commit
 95%, sanity check passes), pulsed 4-on/3-off 85%, capped+escalation 80% POOLED. So realistic dosing
 costs ~8 pts (pulse gaps) and the pathway cap costs another ~5 pts.
 INDUCTION STRATEGY (front-loaded vs escalate-from-1, capped, full N): IDENTICAL -- 80% POOLED both,
 84/77/84 each. The strategy is a wash: both reach cap-4 by end of induction and durability is read
 over the observation tail. The 5-pt cost is the CAP (limiting >4-pathway patients), not the ramp.
 DECISION: keep escalate-from-1 as default (same durability, less early drug = minimal-kernel ethos);
 front_loaded flag retained in durability_clock.run_capped_on_clock for the clinical-induction reading.
 REGRESSION now locks BOTH idealization and clinical: E40 (idealized continuous) + E41 (monolith) +
 E42 (clock: continuous | pulsed | capped). Baseline re-captured; PASS.
 PHASE 3 COMPLETE (steps 7-10). Deferred to Phase 6: whether the clinical clock arm replaces continuous
 A_dyn as canonical durability (would re-lock E40 to ~80%); kept additive + tracked for now.

E71. GENOTOXIC ARM CLOSED + p53 death cascade (PUMA/NOXA/SMAC) built, guarded, validated.
 Built p53->apoptosis cascade as a GUARDED overlay in step_honest (cascade=False default -> regression
 green). PUMA/NOXA (= TP53) neutralise BCL2/BCL2L1/MCL1; SMAC (= BAX or BAK1) neutralises XIAP; CASP3
 re-derived. DISCOVERY: XIAP was the silent executioner block. PUMA/NOXA alone -> only 56% durable;
 adding SMAC->XIAP -> 98%, reproducing/exceeding the abstract COMMIT flag (89%). The real cascade now
 meets the bar to (later) replace the flag.
 GENOTOXIC CLOSE (cascade on, durable over observation window):
   clamp-only      TCGA 98% / METABRIC 92%
   genotoxic-ALONE TCGA 59% / METABRIC 60%
   combination     = clamp-only (98 / 92)
 CONCLUSION: genotoxic is a genuine monotherapy driver (~60% via the p53 cascade) but inferior to
 targeted AKT-axis inhibition (~95%), and adds ZERO durability in combination because the targeted
 arm is already near-ceiling. Genotoxic is p53-gated (silent in the ~1/3 p53-LOF patients). The
 genotoxic INPUT is correct and complete (model-level term at p53, clock-scheduled, p53-gated,
 resolving/sticky switch); it simply has little to add on a near-sufficient kernel. ARM CLOSED.
 CORRECTIONS logged (knife-edge, all found by testing not assertion):
   (1) earlier "kernel clamps TP53" was WRONG -- TP53 is never clamped (DRUGGABLE_NODES already
       excludes it; the kernel only touches the 5 PI3K-axis actuators).
   (2) "p53 disconnected from the caspases" was WRONG -- p53 drives BAX/APAF1/FASLG in the bare rules.
   (3) "irreversibility was the blocker" was WRONG -- XIAP (missing SMAC) was.
 Cascade remains GUARDED: the flag still owns commitment in all locked runs; numbers unmoved.
 NEXT: commit-readout switch (turn flag off, commit via the cascade + MOMP latch) = Phase 5 step 16;
 and the DRUGGABLE_NODES revision.

E72. PHASE 5 STEP 15 -- executioner caspases CASP6/CASP7 + amplification feedback (GUARDED).
 Added to the cascade overlay in honest_bbcn.step_honest, inside the cascade=False guard:
   CASP7 = (CASP8 or CASP9) and not XIAP   (executioner, parallel to CASP3, same XIAP gate)
   CASP6 = CASP3 or CASP7
   CASP8 = CASP8 or CASP6                   (executioner -> initiator amplification loop)
   CASP3 re-derived with the amplified CASP8 (closes the positive loop).
 Regression GREEN (default off, never runs). Cascade durability TCGA-200: 97% (was 98% with SMAC
 only) -- statistically the same, 1 patient of 36. The loop did NOT raise the number: the SMAC-
 unblocked death state is already a stable attractor, so durability was saturated. The feedback is
 the biologically correct point-of-no-return and belongs in the model/paper, but it is not the
 numerical lever. Data-seeded formalization of CASP6/CASP7 as real nodes -> consolidation.

E73. IRREVERSIBILITY MECHANISM -- PROVEN (answers AIB's question; verified from code, not asserted).
 The cascade carries NO apoptosis latch of its own. Durability's irreversibility is the bistable
 AKT1-FOXO3-PHLPP switch -- the project's central mechanism. With PHLPP coupled the bare AKT1 rule
 reduces to AKT1_new = (PI3K-axis) AND AKT1_old, which is monotone: once the kernel drives AKT1 to
 0, it cannot return even after the drug clears and the PI3K nodes recover.
 PROOF (TCGA, kerneled n=78, cascade durability, clamp_off forced):
   PHLPP intact (clamp_off=0):    durable 100%,  AKT1 ON at window end  0%
   PHLPP uncoupled (clamp_off=1): durable  41%,  AKT1 ON at window end 100%
 Breaking the switch collapses durability and AKT1 climbs back on in every patient. So the kernel
 flips the switch, PHLPP holds it, and the cascade reads apoptosis out of the held survival-off
 state. Ties the science shut: durability = coupling; the ~9% non-durable at full N are the
 uncoupled (resistant) patients. Resistance = AKT-FOXO uncoupling = the kill that does not hold.
 The old COMMIT flag was a synthetic latch bolted on top; the cascade inherits the real one.
 MDM2/p53 clarification (AIB asked why modified MDM2 "has no p53"): it does. p53_repair.py --
   canonical (degenerate): MDM2 = TP53 or MDM2 ; TP53 = not MDM2  (MDM2 self-locks high, p53 pinned
     low, damage never reaches it).
   modified: MDM2 = (TP53 or AKT1) and not (CDKN2A and (E2F1 or damage)) ; TP53 = (not MDM2) or damage.
 TP53 is present in the modified MDM2 (the g('TP53') term). What was removed is the degenerate
 `or MDM2` self-lock (the bug). What was added is AKT1 (AKT stabilizes MDM2) and the ARF/CDKN2A+
 damage gate that shuts MDM2 to release p53. Negative feedback p53->MDM2-|p53 fully intact. Nuance:
 when AKT1 is ON (resistant state) MDM2 is carried by AKT1 and the p53 term is functionally masked;
 p53's hand on MDM2 only shows once AKT goes off.

E74. PHASE 5 STEP 16 -- flag retired as the commitment mechanism; cascade switched in; baseline RE-LOCKED.
 Decisive comparison (compare_commit.py, full N, 806 kerneled): IDEAL flag 95% vs cascade 91%;
 CLOCK flag 85% vs cascade 91%. Cascade is slightly more conservative on the idealized hold (flag
 was over-committing) and BEATS the flag on the clinical pulsed clock (+6 pts), because the flag's
 theta-accumulator decays on pulse off-days while the cascade's death attractor rides through the
 gaps. Acceptance bar MET.
 CODE: removed the flag brake-zeroing (for bk in BRAKES: nb[bk]=0) and switched step_honest(cascade
 =True) in: regression_check.A_dyn (E40), durability_clock.A_dyn_clock + run_capped_on_clock (E42)
 + run_on_clock, and run_durability_full.A_dyn (canonical full-N E40). In run_capped_on_clock /
 run_on_clock the commit accumulator is retained ONLY as the escalation-stop / latch-day reporter,
 not as a commitment mechanism.
 RE-LOCK: backed up the flag baseline -> docs/regression_baseline_FLAG_pre_step16.json. Re-captured
 docs/regression_baseline.json to the cascade. 30-sample drift was E42/TCGA cont|puls|capped 9->11
 (cascade beats flag on the clock); E40/E41 raw counts unchanged on the sample. Regression PASS.
 DEFERRED TO CONSOLIDATION (separate task, per AIB): propagate the same one-line cascade switch to
 the remaining duplicate/other-engine runners (fulln_e40, fulln_clock, fulln_frontload,
 run_designed_kernel_AB, run_designed_kernel_cohorts, run_full_cohorts_repaired, run_unified_
 repaired) and then DELETE the now-dead flag scaffolding (honest_bbcn.commit_signal, BRAKES, the
 theta/cmax escalation timer). Doing it piecemeal across ~8 files risks error; it lands cleanly when
 the engines merge into one durability function (roadmap step 18). Full-N cascade numbers are
 already on record via compare_commit.py, so no number is blocked by the deferral.
