# Relationship to prior work (BBCN118)

This repository builds on the author's earlier preprint:

> Bhatti, A. I. "From Pathway to Patient: Molecular Dysregulation as Basis for Minimal Sequential
> Intervention Strategy in Basal-Like Breast Cancer." Research Square (2026).
> DOI: https://doi.org/10.21203/rs.3.rs-8942485/v1  (CC BY 4.0)

That work introduced **BBCN118**, a 118-node modular Boolean breast-cancer network decomposed into 15
pathway modules, with two contributions:
1. a **Lyapunov-based pathway-decomposition theorem** (Theorem 1 + Corollary 1): under weak
   inter-pathway coupling (max |dFi/dxj| < epsilon) and locally stabilising kernels, sequential
   pathway-wise forward-stabilisation drives the global state to an O(epsilon) neighbourhood of the
   target attractor; and
2. a patient-specific pipeline on 22 deceased basal-like TCGA-BRCA patients, reporting per-patient
   mismatch-fraction reduction (mean 0.508 -> 0.404; paired t, p = 6.7e-5) with kernel hubs
   (JAK2, SOS1, CDKN2A, CDKN1A).

## How the present work differs
- **Scope.** BBCN118 = 118 nodes / 15 pathways on 22 basal-like deceased patients; the present BBCN =
  135 nodes / 24 pathways on full TCGA, METABRIC, and I-SPY2 cohorts (N = 1082/1980/986).
- **Central finding.** BBCN118's decomposition theorem holds **under weak coupling**. The present work
  shows that the full network is in fact dominated by a **strong cyclic core** (16 of 22 pathways), so
  the weak-coupling regime does **not** cover the hard case; the obstruction is structural
  (feasibility != reachability). The decomposition theorem is therefore the *boundary* the present
  work pushes past, not a tool reused here.
- **New mechanism.** The bistable, multi-timescale (Boolean-ARMA) resistance-apoptosis switch is new
  to this work and does not appear in BBCN118.
- **Shared actuator.** Both descend from the same forward-stabilisation kernel design, here
  upgraded from the ranked reachability heuristic (size/impact/steps/causal) to the algebraic
  global-stabilisation test of Rafimanzelat (Theorem 1); both methods are retained and compared,
  so the present controller is continuous with the earlier method at the kernel level.

The decomposition theorem and its Lyapunov proof are preserved in the cited preprint; they are not
reproduced here to keep this paper focused on the strong-coupling regime and its switch-based
resolution. A dedicated methods paper on the decomposition theorem is planned separately.
