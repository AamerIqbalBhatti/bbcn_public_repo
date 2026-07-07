# BBCN — Daily Ledger — 2026-06-18

Daily working record. Captures every decision, task, and finding from this
session. Pairs with the standing verification record (BBCN_verification_record,
§1-80) and the current manuscript (BBCN_preprint.docx). Successor to the
§81-91 build log that was lost when its container reset.

Session focus: planning and regenerating the kernel heat-map figures for the
manuscript, sourced as designs from the SR preprint and the TCBB paper.

---

## Settled decisions

D1. Model node count is locked at 136. This overrides the "141" that appears in
the §70-80 record and the "118" in the SR and TCBB papers. Per AIB.

D2. Figures will be regenerated on the current 136-node, three-cohort data
(TCGA, METABRIC, I-SPY2). The SR and TCBB figures (118-node, 22-patient,
single basal-like cohort) will not be lifted in. Importing them would mix a
prior model version's empirical results into the three-cohort paper, which is
the version-inconsistency line we do not cross.

D3. Figure designs from the SR preprint and TCBB paper are adopted as
templates. Only the layout and figure type are reused. All data comes from the
current pipeline.

D4. Production figures render at 300 DPI, the npj / Nature Portfolio minimum
for line and combination art at final size. 150 DPI is used only for quick
draft iteration.

D5. Daily-ledger practice adopted. A new dated ledger is started each session,
records every decision, task, and finding, and is delivered as a Markdown file
at end of day.

D6. Cohort scope confirmed (O1 resolved to option a). The main kernel figures run
on all patients across the three cohorts. One basal-like / TNBC subset panel is
added as a re-affirmation that the model recovers sensible TNBC-relevant kernels.

D7. The exact Python pipeline in BBCN repository v2.zip is used as-is. It is not
reimplemented or replaced. Figure generation runs through the existing code. The
pipeline represents significant development effort and is authoritative.

D8. Figures are generated for both experimental setups, setup A and setup B (the
setup column in all_numbers.csv). Confirmed from the repo: Setup A is the
135-node three-tier capped controller (TCGA/METABRIC/I-SPY2); Setup B is the
136-node bistable resistance-apoptosis switch (135 + PHLPP).

---

## Findings

F1. The SR preprint and the TCBB paper are both BBCN118: 118 control nodes (the
TCBB text also references 135 network nodes), a single basal-like TCGA cohort
of N=22 deceased patients, and identical statistics (t=4.95, p=6.7e-5, Cohen's
d=1.06). The TCBB paper is the more developed version of the same study. It
adds a shared Boolean bus, a stage-wise supervisory hierarchy, a scheduling
layer, SAT-based feasibility analysis, and an SCC-based feasibility-versus-
reachability argument.

F2. Source figure inventory.
  - SR Fig 1: per-patient mismatch-reduction bar.
  - SR Fig 2: kernel-usage frequency bar (single cohort).
  - SR Fig 3: pathway x kernel incidence heatmap.
  - SR Fig 4: patient x kernel incidence heatmap.
  - SR Fig 5: pipeline swimlane.
  - TCBB Fig 1: hierarchical control architecture schematic.
  - TCBB Fig 2: pathway-ordered binarization-flip distribution.
  - TCBB Supp S1: binarization-robustness heatmap (node x patient, median vs
    z-score).
  - TCBB Supp S2: phenotype-resolved structural panel (largest SCC, fraction
    cyclic, relevant pathways, condensation DAG depth).
  - TCBB Supp S3: SAT/UNSAT lock-relaxation curves per phenotype.
  - TCBB pipeline outputs (described, not all shown): cohort MPB summaries,
    phenotype evolution curves, representative patient traces, kernel-activity
    heatmaps, pathway-burden profiles.

F3. The current plotting code (plots.py) emits survival and convergence figures
only: Kaplan-Meier convergence, KM kernel-drug alignment, a Cox forest plot,
K_min vs K_max panels, convergence by subtype, V_SCC severity, and a demo-
patient trace. It does not currently emit kernel-incidence heatmaps.

F4. all_numbers.csv holds aggregate quantities only, keyed by setup, path,
method, cohort, and quantity (for example TCGA N=1082). It does not contain
per-patient kernel selections.

F5. Regeneration dependency. The pipeline expects a local data/ directory:
BBCN_cohort_results.csv ("algorithm results per patient"), clinical patient and
sample tables, a timeline file, an mRNA expression matrix, the binarized matrix
bbcn_binarized_135genes.csv, and target JSONs (*_fullbus.json). This data/
directory is not present in the project mount. It must come from
BBCN repository v2.zip on Drive (ID 1wLYShopp2MY8WYa8IMU90CTNcyZ9qWXs) before
any regeneration can run.

F6. The standalone "all numbers.csv" on Drive (9.8 KB) matches the local copy and
holds aggregate quantities only (setup, path, method, cohort, quantity, value).
It carries no per-patient kernel selections. The per-patient kernel data needed
for the heat maps lives inside the data/ directory of BBCN repository v2.zip.

F7. Drive folder inventory (parent 1jP3...). Latest repository is
BBCN repository v2.zip (2026-06-17, 3.20 MB). Also present: the loose
all numbers.csv, BBCN preprint.docx and .pdf, two older BBCN repository.zip
files, BBCN public repo 2026-06-04.zip, the verification record, SR preprint.pdf,
and the TCBB paper and supplementary. v2.zip is the version to use.

---

## Proposed figures to add (regenerated on current data)

P1. Kernel-frequency heatmap, kernels by cohort (TCGA / METABRIC / I-SPY2).
The SR Fig 2 idea upgraded to show all three cohorts in one panel. Candidate
headline kernel heat map.

P2. Pathway x kernel incidence heatmap (SR Fig 3 design). Shows the near
block-diagonal modularity that supports the decomposition claim.

P3. Optional, only if the manuscript keeps a feasibility / reachability thread:
the phenotype structural panel (S2 design) and the lock-relaxation curves (S3
design).

Not carried over: the raw patient x kernel heatmap (SR Fig 4). At roughly 1082
patients per cohort it becomes unreadable. A clustered or per-phenotype-
aggregated version would be used if patient-level detail is wanted.

---

## Open decisions

None open. O1 (cohort scope) is resolved to option (a); see D6.

---

## Session execution log (2026-06-18)

E1. Repository unzipped from the chat upload. Environment: numpy 2.4.4, pandas
3.0.2, matplotlib installed. The exact pipeline runs here.

E2. Reproduction verified. run_chunked.py --status shows 27/27 tasks complete and
checkpointed in the zip. numbers.json matches results/LOCKED_NUMBERS.md (for
example TCGA dynamic stage-pass 14/31/65/4; Setup B routing TCGA 50/36/14).
LOCKED_NUMBERS.md and all_numbers.csv are the authoritative source for caption
numbers.

E3. Heat-map data sources confirmed in code.
  - Setup B: forward_stab_kernel_design.run() returns the full node_use counter;
    the number-generator emits only its argmax (top_node). We collect the full
    distribution. The shipped 59-gene switch-input sample TSVs drive it with no
    downloads.
  - Setup A: per-pathway kernel selections are recorded by the harness kernel_log
    via controller.run_patient. Capturing them needs a per-patient pass.
  - Caveat: the shipped TCGA switch sample loads 1081 rows vs locked N=1082 (one
    row short). This shifts a usage count by at most one and does not move the
    locked flip_pct. Logged, not material to the figure.

E4. Setup B kernel-actuator heat map produced: figB_kernel_usage_setupB.png, 300
DPI, repo style, nodes by cohort, panels ranked and stabilize, driven by KD.run
with no reimplementation. Node usage as percent of designed kernels (stabilize):
    PIK3CA  98 / 97 / 94   (TCGA / METABRIC / I-SPY2)
    AKT1    78 / 76 / 73
    PDPK1   52 / 52 / 54
    MTOR    48 / 38 / 45
    PTEN    27 / 35 / 26
  PI3K/AKT axis dominance, consistent across cohorts and methods, matching the
  paper's PI3K/AKT/mTOR nomination.

E5. Setup A kernel figures produced (static, isolated/ideal; both methods; full N
per cohort; captured by observing controller._kernel_for, no logic change).
  - figA1_kernel_usage_setupA.png: top-24 kernel nodes by cohort, ranked and
    stabilize panels. PTEN, SOS1, NOS3, GRB2 at 100% in every cohort and method.
    Method contrast: MAP2K1 80% ranked vs 0% stabilize; EGF/INS 50% ranked vs
    100% stabilize; ABL1/PIK3CA 100% ranked vs ~79% stabilize.
  - figA2_pathway_kernel_setupA.png: pathway-by-kernel incidence (stabilize,
    pooled). Clean near block-diagonal: RTK_EGFR -> EGF/GRB2/ERBB2/ERBB3;
    RTK_Insulin -> INS/IRS1/INSR; PI3K -> PTEN/NOS3/PDPK1/RAC1; MAPK ->
    SOS1/DUSP1/MAPK8/NF1; JAK_STAT -> JAK2/SOCS3/STAT1; AKT_Signaling ->
    NEDD4L/SGK1/AKT1/GSK3B; mTOR -> TWIST1/RHEB/TSC2. Empirical support for the
    weak-coupling decomposition.

E6. Figures inserted into the manuscript Results with detailed captions, editing
the LaTeX source (preprint/BBCN_preprint.tex; backup .orig kept). New figures
copied to preprint/figures/. Placement: Setup B kernel actuator -> Section 3.4
(Figure 4, the PI3K/AKT/mTOR nomination); Setup A kernel usage -> Section 3.6
(Figure 5); Setup A pathway-by-kernel modularity -> Section 3.6 (Figure 6). Each
caption explains what is shown, the ranked-vs-algebraic contrast, and the link to
the argument. Figure 6 caption is careful to note that kernel modularity is about
where control is applied and does not contradict the strong-coupling obstruction
of Section 3.2.

E7. Methods and Results expanded (Task 2), in the author's style: short
sentences, fully justified, no em-dashes in the new prose. Methods now cover the
cohort rationale and robust z-score, the 135/136 two-setup identity and the bus,
the three-tier controller in six labelled paragraphs (tiers and timescales; caps
and no accumulation; ideal vs compatible targets; algebraic kernel design; the
retained ranked heuristic; success as a genuine fixed point with TRACE
classification), the bistability rationale, and validation with honest scope.
Results 3.1-3.4 and 3.6 each gained interpretation (individual vs joint; what an
SCC is and the BBCN118 boundary; threefold enrichment and falsifiability;
hysteresis and the drug split; kernel sparsity and modularity). Fixed an I-SPY2 N
typo, 986 -> 988, to match the locked N. Every number remains a numbers.tex
macro; nothing is hardcoded.

E8. Rebuilt from source: pdflatex (two passes, clean, no undefined references)
and pandoc. Output is 15 pages; PDF and DOCX regenerated. Table numbers verified
against LOCKED_NUMBERS on the rasterized pages. Two pre-existing cosmetic
overfull-hbox warnings remain (the wide Setup B AKT1 equation; the GitHub URL),
neither introduced here.

E9. Appendix A (the BBCN digital logic model) auto-generated from
setup_a/bbcn/pathways.py by generate_rules_appendix.py, emitting
preprint/appendix_rules.tex which the manuscript inputs. Covers all 24 pathways
and 135 nodes as readable Boolean equations (AND/OR/NOT), with the two auxiliary
input modules shown as pass-throughs. Generating from source means the printed
rules cannot drift from the code.

E10. Appendix B (flowcharts) generated by generate_flowcharts.py: figC_flow_ranked
and figC_flow_stabilize, 300 DPI, repo style, the ranked heuristic and the
algebraic global-stabilisation test as Figures 7 and 8.

E11. Appendix C (decomposition theorem) written as preprint/appendix_theorem.tex:
Theorem 1 with assumptions A1-A3, weighted-Hamming Lyapunov proof sketch,
Corollary 1, and a closing remark naming the assumption BBCN violates (the
dominant 16-of-22 SCC, Section 3.2). Attributed and adapted from the cited prior
preprint (bhatti2026, CC BY 4.0), not presented as new. This attributed form is
the agreed choice and supersedes the older project note about absorbing the
theorem without a citation.

E12. Appendices wired into the manuscript via \appendix before the references;
the float package pins the flowcharts (Figures 7, 8) in Appendix B ahead of the
reference list. Clean compile, 22 pages; PDF and DOCX rebuilt.

E13. Introduction expanded to a full state of the art with verified citations.
Three cited paragraphs under "Boolean models of cancer signalling and their
control": foundations (Kauffman, Thomas, Albert and Othmer); logic-based
signalling (Saez-Rodriguez, Morris, Cell Collective, GINsim); cancer and
breast-cancer Boolean models (Fumia and Martins, von der Heyde, Sgariglia, Taoma,
Zanudo 2017); and control methods (stable motifs, feedback vertex set, STP,
global stabilisation, probabilistic BNs, patient-specific personalisation). 14
new bibitems added (27 total). References verified, most from the author's prior
preprint bibliography; none invented. One Introduction em-dash converted to a
colon. Clean compile, no undefined citations, 24 pages; PDF and DOCX rebuilt.

E14. Repository assembled as a synced, one-button, public-ready tree
(BBCN_repository_v3.zip).
  - Figure/appendix generators made repo-relative and robust (write into
    preprint/ and preprint/figures with the exact names the .tex references):
    generate_rules_appendix.py, generate_flowcharts.py, make_setupB_kernel_heatmap.py,
    setup_a/capture_setupA.py, setup_a/make_setupA_kernel_heatmaps.py.
  - make_figures.py orchestrates every figure/appendix emitter in one call.
  - reproduce_all.py wired with a --figures regen step and an --all one-button
    shortcut; matplotlib auto-installed when needed.
  - Sync map: numbers (run_chunked -> numbers.tex), Appendix A
    (generate_rules_appendix -> appendix_rules.tex), Appendix B flowcharts,
    Figures 4-6 kernel scripts; all are \input or included by the manuscript, so
    text cannot drift from code. Appendix C is static attributed text.
  - Verified the full one-button build end to end (reproduce_all --figures regen):
    figures regenerate, numbers sync, PDF (24pp) and DOCX rebuild, Setup B checks
    pass.
  - README and preprint/README updated; Colab notebook updated (requirements.txt
    + --figures regen); .gitignore extended for LaTeX artifacts; matplotlib added
    to requirements and binder. Backups and build artifacts removed; zip verified
    clean with all 8 figures and all \input targets present.

E15. Introduction subsections completed. Restructured the latter half into 1.2
Gap analysis (the feasibility-versus-reachability gap and the single-clock
timescale gap, sharpened and enumerated), 1.3 Our approach (the two-setup
approach, the four contributions, and the BBCN118 relationship, now pointing to
Appendix C), and 1.4 Paper structure (a roadmap cross-referencing every Results
subsection and the three appendices). Short sentences, justified, no em-dashes.
Clean three-pass compile, no undefined references; 24 pages. Repo rebuilt and
repackaged as BBCN_repository_v4.zip; dated PDF and DOCX refreshed.

E16. Receptor-low (TNBC-like) re-affirmation panel for Setup A added (Figure 2).
Subset defined by low ESR1, PGR, ERBB2 in the binarised data (transcriptomic
proxy for triple-negative status; TCGA 237/22%, METABRIC 363/18%, I-SPY2
279/28%). The same Setup A static controller (ranked) was run on full and subset;
the full-cohort bars reproduce LOCKED_NUMBERS exactly (apoptosis 86/85/81,
all-three 3/2/3), validating the script. In the subset, individual phenotypes stay
highly controllable (apoptosis 75-91%, proliferation 95-98%) while all-three stays
1-3%, re-affirming the individual-versus-joint finding. make_tnbc_panel.py wired
into make_figures.py. Caption is honest that it is a proxy, not IHC-confirmed TNBC.

E17. Em-dash cleanup completed document-wide: the date line, the Setup B
multirate sentence, the two kernel-comparison sentences, and the Table 1 row
labels and empty-cell placeholders (now "n/a"). No prose or table em-dashes remain.

E18. Setup B routing re-affirmation on the TNBC-like subset added (Figure 6). A
backward-compatible per-patient records hook was added to
cohort_pipeline.run_cohort. make_tnbc_setupB.py routes the shipped switch inputs
(no mutation safeguard, so within about two points of the Table 2 routing), joins
receptor status by sample ID (TNBC n=237/363/279), and compares full versus subset
on routing and the survival-routed uncoupled rate. The full-cohort survival-routed
uncoupled rate reproduces the locked headline (96/92/89 vs 95/90/89). In the subset
the uncoupled enrichment holds (81-99%), re-affirming the hysteresis reading across
receptor subtypes. Wired into make_figures.py. Clean three-pass compile, no
undefined references; 25 pages. Repo repackaged as BBCN_repository_v6.zip; PDF and
DOCX refreshed.

---

## Next steps (execution)

Done: kernel figures (both setups); Methods and Results expanded; Appendices A, B,
C; full Introduction (state of the art, Gap analysis, Our approach, Paper
structure); the synced one-button repository; receptor-low (TNBC-like)
re-affirmation panels for both setups (Figures 2 and 6); and document-wide em-dash
cleanup. Current repository is BBCN_repository_v6.zip; PDF and DOCX are 25 pages.

Remaining, optional:
N1. Any reviewer-facing polish before submission to npj Systems Biology
(author-note placeholders in the Discussion still mark prose to be finalised).
N2. If the controlled-access expression and mutation matrices are provided, the
Setup B TNBC routing can be recomputed with the mutation safeguard to match the
Table 2 percentages exactly, rather than the shipped-input approximation.

---

## Session addendum — 2026-06-18 (durability: pulsed arm, scheme-invariance, Setup B payoff)

Successor block to E18. New control arm, durability analysis, an exactness fix to
the existing pipeline, and two manuscript paragraphs. Repackaged as v7.

E19. Pulsed weekly controller added as a third arm (setup_a/bbcn/pulsed.py). One
cycle = one week = 7 Boolean updates, one update = one day. Each cycle acts on ONE
pathway, whose kernel may pin up to ind_cap=3 nodes (a vertical multi-target
cocktail on one axis, distinct from horizontal multi-pathway hold). Exposure is a
drug-on pulse for on_days (default 1, with 2 and 3 compared), then the kernel is
WITHDRAWN and the network relaxes free for the rest of the week and settles.
Success = the phenotype holds as a genuine drug-off (free) fixed point after
relaxation. Cycles up to 5, one pathway each, kernel held within a week and
re-decided at each cycle. Druggable preference added: among equally minimal
stabilising kernels, prefer the most druggable set, weights from drugs.py
DRUG_TABLE (FDA-approved node = 2, trial = 1). drugs.py copied into setup_a/bbcn.

E20. Exactness fix to the kernel cache, applied to the EXISTING pipeline per AIB.
The old external set came from harness._pathway_externals, a functional probe at
only the all-0 and all-1 internal states; verified it misses live externals for
10 of 24 pathways (PTEN for AKT_Signaling and mTOR, MCL1 for Apoptosis_Intrinsic,
TP53 for CellCycle, AKT1 for MAPK, IL6 for JAK_STAT all shown to change the rule
output). Added harness.complete_externals: the complete set of bus nodes the rules
reference, AST-extracted from pathways.py and unioned with the probe. Pointed the
controller stabilize+ranked cache key and the pulsed cache at it. Measured impact:
recomputing each patient from its own full state changes 0 of 120 outcomes;
percentages identical (sample ranked 88/93, stabilize 93/86), so locked numbers
unchanged. The algebraic kernel is x0-independent and reused per external snapshot;
the ranked kernel is x0-dependent and recomputed per cycle per patient (keyed by
externals AND x0).

E21. Static never withdraws — confirmed empirically from control_phenotype on 60
TCGA patients (apoptosis, stabilize). Zero free-mode evolutions during control (86
held evolutions, 0 withdrawn); at success the kernel is non-empty 58/58, holds as
a HELD fixed point 58/58, and as a FREE drug-off fixed point 0/58. So the static
arm is continuous-hold by construction (the old pin-forever paradigm); the
durable-after-withdrawal test is a separate diagnostic, not part of static.

E22. Durability findings (Setup A). Static continuous hold reaches apoptosis
91/90/83% (TCGA/METABRIC/I-SPY2) but 0% durable: of 38/40 achieved on a sample,
all collapse into the apoptosis limit cycle on withdrawal (period>1), none is a
free fixed point. Pulsed weekly (200/cohort): apoptosis, proliferation, resistance,
all-three all 0% durable at on_days 1, 2, 3 in every cohort. Duration sweep
(hold D days then release, D=1..100): apoptosis durable flat 0% at every D;
proliferation saturates at 8% from D>=5; resistance at 2%. The plateau proves the
ceiling is the size of the reachable durable attractor, not an exposure-length
limit. Saved pulsed_cmp_*.json, duration_experiment_TCGA.json.

E23. Scheme-invariance test (answers the synchronous-artifact objection). Fixed
points are invariant to the update scheme, so durability reduces to whether an
apoptotic fixed point exists. Setup A free synchronous dynamics from resting states
settle to a true fixed point in 0/40 patients; all 40 fall into limit cycles
(median period 4). Node-level random-order asynchronous updating does not settle
them either (0/40 from resting; 0/38 apoptosis), so the oscillation is structural,
not a synchronous artifact. Setup B's apoptotic state IS a synchronous fixed point
of its rules for 11/58 resistant patients, and the kernel drives the patient
exactly onto it (11/58). The multirate governs convergence; the durability is
structural. The 19% matches the Setup B durable fraction.

E24. Setup B payoff (setup_b/code/payoff_pulse.py): the latch test. For each
survival-routed (resistant) patient, design the apoptosis kernel, clamp it for P
ticks (drug on), release, simulate the multirate switch free, check it stays
APOPTOTIC. Full N: durable-after-pulse equals reachable-under-hold exactly in every
cohort, TCGA 24% (522/1082 resistant), METABRIC 18% (978/1980), I-SPY2 17%
(480/988), at every pulse length down to 10 ticks, and 0% with no pulse. Setup A
is 0% durable at any duration. The bistable latch keeps everything the kernel can
reach; nothing is lost on withdrawal. Saved payoff_*.json.

E25. Figures (setup_b/code/make_durability_figs.py): fig_switch_flow.png (the
Setup B switch test as a flow diagram, making the 9-node-core scope explicit: 9
simulated nodes + 8 held inputs + U, all patients, NOT the 135-node network) and
fig_durability.png (reached-while-held vs durable-after-withdrawal, per setup per
cohort; Setup A collapses to 0, Setup B bars equal). Scope clarification recorded:
the durability payoff lives in the reduced 9-node bistable core, parameterised per
patient from the full transcriptome; Setup A (135-node) is where durability is
shown impossible. The two setups are complementary, problem and mechanism.

E26. Manuscript additions (preprint/BBCN_preprint.tex), clean 3-pass compile, 26
pages, DOCX rebuilt. (a) Section 2.3, new paragraph "Exactness of the cached
kernel": L is a snapshot of the external inputs, refreshed each cycle; keyed on the
complete external set so the reuse is exact; the ranked heuristic is x0-dependent
and recomputed per cycle. (b) Section 2.4, new paragraph "Why the durable death
state is structural, not a scheduling artifact": fixed-point invariance, the
0-of-40 resting-cycle result with asynchrony not settling either, and the Setup B
synchronous apoptotic fixed point; concludes durability is structural and that
feasibility and durability coincide in Setup B. An author-note flags that the
resting-state attractor census and the switch fixed-point count are to be wired
into numbers.tex and macro-sourced. Clinical reading recorded: withdrawal-loss
matches rebound of a non-curative therapy, but the model's apoptosis is reversible
pro-death signalling, not executed death, so durable benefit means crossing the
Setup B commitment latch, not sustaining Setup A signalling.

---

## Next steps (updated)

Done this block: pulsed weekly arm; exact-externals fix in the live pipeline;
static-no-withdrawal confirmation; Setup A durability and duration sweep; the
scheme-invariance/structural-durability test; the Setup B latch payoff at full N;
the flow and durability figures; and the two Methods paragraphs. Repository is
BBCN_repository_v7.zip; PDF and DOCX are 26 pages.

Remaining, optional:
N3. Wire the resting-state attractor census (Setup A) and the switch
fixed-point/durable-pulse counts (Setup B) into numbers.tex via reproduce_all, so
the new Methods 2.4 paragraph carries exact macro-sourced counts.
N4. If the durability story is to be the headline, add a Results subsection and
place fig_durability.png and fig_switch_flow.png, then regenerate the tables to
report durable-after-withdrawal alongside reached-while-held. Awaiting AIB's call
on prominence.
N5. Run the pulsed weekly arm at full N (currently 200/cohort; saturated at 0, so
the message will not change) via the checkpointed runner.

E27. bioRxiv pre-flight pass. Body confirmed journal-neutral (the only "Albert"
hits are legitimate citations; no target-journal or editor named). Code-and-data
availability section present (GitHub + one-command reproduce + Colab/Binder).
Fixes applied for the posted build: the four internal author-notes are suppressed
(kept in source via a no-op renewcommand), the "working scaffold" date line is
replaced with "June 2026", and the author block now carries a corresponding-author
footnote with ORCID 0000-0002-3373-3388. Clean compile, 26 pages, DOCX rebuilt.
Open: corresponding-author email (entered in the bioRxiv portal, optionally in the
PDF) and abstract trim from 260 to about 200 words for the eventual npj version.
