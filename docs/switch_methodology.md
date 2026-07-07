# BBCN bistable-switch family: design and test methodology

Reusable record for the bistable cell-fate switches built on the Boolean
Breast-Cancer Network (BBCN). Each switch isolates one fate decision as a small
bistable system with timescale separation. This document is the template: to add a
switch, copy Section 4 and fill the same fields. The code modules are the
executable specification and self-verify.

--------------------------------------------------------------------------------
## 1. Why switches exist (purpose)

The full 135-node synchronous network has no durable cell-fate fixed point, and
the three phenotypes share a dominant strongly-connected core (their control sets
union to 67 nodes and all three share the PI3K core PTEN/PDPK1/NOS3/RAC1). So
whole-network control is obstructed: phenotypes are individually reachable but not
jointly. The switch strategy isolates each fate decision as a small bistable
module with two self-consistent attractors, so a transient pulse can flip it and
hysteresis holds it. Design small, test big.

--------------------------------------------------------------------------------
## 2. The eight-step design principle (reusable)

1. Purpose. Isolate one fate decision as a small bistable switch.
2. Decision axis. Find the double-negative (mutual-repression) loop that carries
   the fate choice.
3. Close the loop. Add the minimal arm that completes the antagonism, even if it
   needs one extra node.
4. Collect nodes. Keep the backbone and its immediate regulators; compress the
   rest of the network into held inputs plus one state flag.
5. Assign timescales. Split core nodes by mechanism into fast (phosphorylation),
   medium (turnover), slow (transcription); run on a 1:5:25 multirate schedule.
   The separation breaks the synchronous oscillation and lets the loop settle.
6. Verify bistability. Enumerate the fixed points of the small core; confirm
   exactly two self-consistent attractors.
7. Algebraic form. Lift the delayed rules to delay-free companion form
   (semi-tensor product) so steering is a reachability problem.
8. Design and test kernels. Per box, design the minimal druggable
   forward-stabilization kernel; test by pulse and release (clamp briefly,
   withdraw, simulate free, check the latch holds).

--------------------------------------------------------------------------------
## 3. Timescale boxes (mechanism rule)

Fast = state set by phosphorylation or other post-translational signalling.
Medium = state set by degradation or turnover. Slow = state set by transcription.
Whole network: 81 fast, 10 medium, 44 slow (Setup B adds PHLPP). Each switch
inherits its nodes' classes.

--------------------------------------------------------------------------------
## 4. Built switches

### 4.1 Apoptosis switch (Setup B): survival <-> apoptosis
- Module: setup_b/code/bbcn_switch.py
- Decision axis: AKT1 -| TP53 mutual antagonism; closing arm not-AKT1 -> FOXO3 ->
  PHLPP -> not-AKT1 (PHLPP is the 136th node).
- Core nodes (9): FAST [AKT1, PTEN, MTOR, PDPK1, PIK3CA]; MID [MDM2];
  SLOW [TP53, PHLPP, FOXO3].
- Inputs (8, held, patient-derived): SRC, RHEB, IGF1R, RTK_up, CDKN2A, E2F1, ATM,
  ATR. Plus the uncoupled-state flag U.
- Fixed points: two (SURVIVAL, APOPTOTIC).
- Latch result: survival-axis inhibition alone cannot flip survival to apoptosis;
  only genotoxic (ATM/ATR) input crosses the latch.
- Cohort payoff (durable-after-pulse = reachable-under-hold): TCGA 24% of
  resistant patients (522/1082), METABRIC 18% (978/1980), I-SPY2 17% (480/988);
  zero with no pulse.

### 4.2 Proliferation switch: quiescence <-> proliferation
- Module: setup_b/code/prolif_switch.py
- Decision axis: RB1 -| E2F1; latch E2F1 -> CCNE1 -> CDK2 -| RB1 (and MYC
  autoregulation). The E2F-cyclinE-CDK2-RB positive feedback is the hysteretic
  commitment to division.
- Core nodes (11): FAST [CDK4, CDK6, CDK2, RB1]; MID [CCNE1, CDKN1B];
  SLOW [CCND1, E2F1, MYC, CDKN1A, CDKN2A].
- Inputs (3, held): MITOGEN (growth drive), TP53 (p53 arrest), STRESS
  (senescence/p16).
- Fixed points (resting regime MITOGEN=TP53=STRESS=0): two.
  QUIESCENT: RB1=1, E2F1=0, all CDKs/cyclins=0, p27=1.
  PROLIFERATIVE: RB1=0, E2F1=1, CDK4/6/2=1, cyclin D/E=1, MYC=1, p27=0.
  Mitogen on -> monostable proliferative; p53 on -> bistable with p21 high in the
  quiescent state.
- Latch result (pulse-release from proliferative): restore RB1 latches quiescent;
  knock down E2F1 latches quiescent; CDK4/6 + CDK2 together latches quiescent;
  CDK4/6 inhibition ALONE does not flip (cyclin-E/CDK2 escape re-releases E2F1).
  That reproduces CDK4/6-inhibitor resistance from network structure alone.

--------------------------------------------------------------------------------
## 5. Switches to add (same template)

- Resistance switch: AKT-FOXO survival arm (survival <-> death-primed). Axis around
  AKT1 -| FOXO3 with the survival effectors (RELA, STAT3, XIAP, CFLAR, MCL1).
- Invasion switch: epithelial <-> invasive. Nodes YAP1, CTNNB1, NOTCH1, TWIST1
  (Invasion_OFF schema).

--------------------------------------------------------------------------------
## 6. Delayed full-BBCN test bed (test big)

Goal: run the full 135/136-node network under per-node delays (81/10/44 on the
1:5:25 schedule) as a shared simulation arena, and test whether a kernel designed
on a small switch still drives the full network to the target fate.

Verification (before trusting any result):
1. Rate-1 parity. At rates 1:1:1 the multirate engine must equal one synchronous
   step of the base network. Unit test on the engine.
2. Fixed-point invariance. Any fixed point of the base network must remain fixed
   under any delay schedule (delays change reachability and settling, never the
   existence of fixed points). Enumerate/track and confirm.
3. Switch-embedding consistency. Restrict the full delayed network to a switch
   core, holding the rest as that switch's inputs; it must reproduce the small
   switch's two attractors.

Test protocol: take a patient the small switch routes to the wrong attractor
(survival, or proliferative), apply the small-switch kernel on the full delayed
network from the patient state, simulate free, and check it settles to the target
fate (apoptotic, or quiescent). The pass rate measures transfer of small-switch
design to the full multirate network.

Honest guardrail: a delay does not invent a fixed point the rules do not contain.
Structured slowing of transcription relative to phosphorylation is not uniform
random asynchrony; it can dissolve spurious synchronous limit cycles and let the
network settle, which is why it is worth running, but settling is the most it can
add.

Results (TCGA, reproducible via setup_a/run_multirate_test.py):
- Check 1 rate-1 parity: PASS (10000/10000 comparisons identical to the base step).
- Check 2 fixed-point preservation: holds. A search over 1500 random plus the two
  constant initial states found NO free fixed point: the full network is globally
  oscillatory under both synchronous and 1:5:25 multirate updating, so the engine
  has nothing to move.
- Settle: synchronous 0/120 and multirate 0/120 settle to a fixed point. Structured
  delay does not dissolve the limit cycles; the oscillation is structural.
- Transfer (apoptosis kernel held, then withdrawn): held-reachable 60/60, durable
  after withdrawal 0/60. The delayed net can be driven to apoptosis and held there,
  but has no apoptotic fixed point to fall onto, so nothing survives withdrawal.

Conclusion: the "test big" arena confirms the "design small" necessity. Durability
is not a property of the full network at any timescale; it exists only in the
reduced bistable switch, which is the one place a cell fate is a fixed point.

--------------------------------------------------------------------------------
## 7. Provenance

Executable spec: bbcn_switch.py, prolif_switch.py (each enumerates its own fixed
points). This document is the human-readable companion. Update both together.

--------------------------------------------------------------------------------
## 8. The p53 / DNA-damage axis (delayed full net)

Biology: p53 is a relaxation oscillator on the p53-MDM2 loop driven by ATM.
After double-strand breaks it fires a train of pulses of fixed amplitude and
duration; the NUMBER of pulses scales with damage (frequency/digital code, not
amplitude). The response is excitable (a transient or sustained input triggers a
full pulse). Few pulses -> repair and arrest; sustained/many pulses -> apoptosis
or senescence (Lahav 2004; Batchelor 2008/2011; Loewer 2010; Purvis 2012).

What the full network did: nothing. Its canonical p53 axis is degenerate,
MDM2 = TP53 or MDM2 (self-locks high) and TP53 = not MDM2 with ATM/ATR not wired
in, so MDM2 latches high, TP53 sits low, and damage is ignored (TP53 duty 0.00
with or without damage).

Repair (bbcn/p53_repair.py, an overlay; pathways.py untouched): swap only MDM2
and TP53 for the switch biology,
  MDM2 = (TP53 or AKT1) and not ( CDKN2A and (E2F1 or ATM or ATR) )
  TP53 = (not MDM2) or ATM or ATR.
Effect on the delayed ARMA net (40 TCGA patients): TP53 duty under sustained
damage 0.00 -> 1.00. The sensor is fixed.

But CASP3 (death marker) stays 0.00 even repaired, because TP53 feeds BAX/BAK
which are gated by the survival brakes (MCL1, XIAP, CFLAR, BCL2), and the net
oscillates. And the repaired full net still has zero free fixed points and still
settles 0/80. So:

Would the same repair in Setup A have changed its results? No, not the headlines.
- Durability: unchanged. Still no apoptotic fixed point, still 0% durable.
- Apoptosis execution: p53 high does not reach CASP3 (brakes + oscillation), so
  the damage route adds no durable path to apoptosis; the controller's achieved
  apoptosis comes from direct kernel pinning, not from p53.
- Joint controllability: the repair ADDS coupling (AKT1->MDM2, E2F1/CDKN2A->MDM2,
  p53 into the death arm), tightening the strongly-connected core, so the
  obstruction stands or slightly worsens.
The repair improves biological fidelity of the damage sensor; it does not change
the conclusions. It reinforces why the bistable switch is necessary: even a
fully repaired full net cannot hold apoptosis.

Note (AR vs ARMA): the delayed net is the autoregressive delay-difference form
(lagged states). The biology's pulse-counting is a moving-average term (an
integrator over past damage). Adding such an accumulator node is the natural MA
extension and the Boolean analogue of the p53 pulse counter; not yet built.

--------------------------------------------------------------------------------
## 9. Full repair set + accumulator (revises the Setup A durability answer)

bbcn/honest_bbcn.py carries three repairs as an overlay (pathways.py untouched):
p53-MDM2-ATM; AKT1-FOXO3-PHLPP with the U flag (adds PHLPP, the 136th node);
and the accumulator, an integrator over p53 (the ARMA moving-average term, the
Boolean analogue of pulse-counting) that latches an irreversible commitment past
threshold and releases the survival brakes (MCL1, XIAP, CFLAR, BCL2, BCL2L1).

Results (TCGA, delayed ARMA net, sustained damage unless noted):
- TP53 duty 0.00 -> 1.00 (sensor works).
- CASP3 duty: 0.00 (bare) -> 0.45 (wiring repairs, PHLPP brake lowers AKT1)
  -> 1.00 (accumulator releases brakes).
- Durable CASP3 after damage withdrawal: 19/40 (repairs only) -> 40/40 (accumulator,
  caspase commitment latched).
- Fixed points of the REPAIRED net (no accumulator, no forced damage): a search
  found 448 distinct fixed points and 14/80 patients settle, ALL of them apoptotic.
  The BARE net found 0 and settled 0.

Revised Setup A answer: the wiring repairs DO change the durability result. The
bare full net's global oscillation and 0% durability were substantially an
ARTIFACT of degenerate loops (MDM2 self-lock, missing PHLPP, ATM/ATR unwired), not
pure structure. Restoring the switch's biology gives the full net 448+ apoptotic
fixed points and ~18% of patients settle to them spontaneously; the accumulator's
irreversible commitment then makes apoptosis durable across the board. The switch
remains the tractable DESIGN surface (the algebra), but durable apoptosis is a
property of the repaired AKT-FOXO-PHLPP-p53 feedback structure, not unique to the
reduced switch. The U (uncoupled) arm is in place; sustained genotoxic input
overrides it, matching the switch (survival-axis inhibition alone cannot flip it,
genotoxic input can).

Manuscript implication (for the author to weigh): the Setup A "globally
oscillatory / no durable fixed point" statement should be qualified as a property
of the un-repaired wiring; the biologically-repaired full net does carry durable
apoptotic fixed points.
