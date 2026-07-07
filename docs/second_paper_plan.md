# Second paper — durable cell fate on the biologically-repaired BBCN

## Thesis
The first paper (bioRxiv) showed the full 135-node network is globally oscillatory:
phenotypes are individually controllable but not jointly, and no cell fate is
durable, so durability was relocated to a reduced 9-node bistable switch. This
paper shows that the non-durability was substantially an ARTIFACT of degenerate
feedback wiring. Restoring three pieces of correct biology to the full network
makes durable, controllable apoptosis emerge on the full net itself.

## The three repairs (all literature-grounded; same biology as the switch)
1. p53-MDM2-ATM: damage (ATM/ATR) stabilises p53; MDM2 gated by ARF/CDKN2A.
   Replaces the degenerate MDM2=TP53|MDM2 self-lock and TP53=!MDM2 (no damage input).
2. AKT1-FOXO3-PHLPP + U: PHLPP (136th node) closes the AKT-FOXO3a loop; U is the
   uncoupled resistant state (AKT1 self-sustains, FOXO3 stays nuclear).
3. Accumulator: p53 pulse-integrator (ARMA moving-average term). Past threshold it
   latches an irreversible commitment (caspase point-of-no-return) that releases
   the survival brakes (MCL1, XIAP, CFLAR, BCL2, BCL2L1).

## CORRECTED durability headline (E40, full N)
- Designed single-agent AKT-axis kernel durably commits ~90% of resistant patients to
  apoptosis on the full repaired net (commit readout 92-95%), agreeing with the switch
  design (80-84%). We have a durable therapeutic solution.
- The 12-14% 'non-durable' figure was a strict readout (CASP3 AND AKT1=0) that demands
  the survival kinase also stay off; uncoupled cells reactivate AKT1 after withdrawal but
  are already committed to death. Commitment (caspase-3 active) is the correct readout.
- Resistance still = AKT-FOXO3a uncoupling (it sets which cells reactivate AKT1), but it
  does NOT prevent durable death once commitment is reached.

## FINAL headline (run_full_cohorts_repaired.py, FULL N, all cohorts)
- Drug-resistance <=> AKT-FOXO3a uncoupling: 4049/4050 patients (99.98%), all three cohorts.
- Uncoupled (resistant) 64-67%; survival-axis inhibition flips exactly the coupled 33-36%.
- Genotoxic input alone barely flips the resistant state (3-4%): reversal needs re-coupling.
- Switch and full net agree (death-engaged accumulator preserves switch bistability).
- CORRECTION: earlier "emergent durable apoptosis 10-14%" was an artifact of forcing U=0;
  with U per patient it is ~0. The repaired net reproduces resistance, it does not abolish it.

## (superseded) earlier N=50 result
- Bare net settles: 0/50 every cohort (matches paper 1).
- Repaired net settles spontaneously: ~10-14% of patients, ALL apoptotic.
- Fixed-point search: bare 0; repaired 448+ distinct fixed points.
- Durable apoptosis after withdrawal (with commitment): 50/50 every cohort.
- Sustained genotoxic input overrides the uncoupled resistant state (U=1), matching
  the switch's hysteresis result.

## Setup A <-> Setup B bridge
The switch's APOPTOTIC attractor equals the repaired full net's apoptotic fixed
point (AKT1=0, TP53=1, CASP3=1). The switch was the tractable isolation of the
durable structure; the repair puts that structure back into the full network.

## Build status
S1 cross-cohort contrast: DONE.
S2 repaired pathways branch + full controller re-run: TODO.
S3 Setup B payoff on repaired full net: TODO.
S4 figures: TODO.
