# BBCN — Methods and Architecture

This document specifies the BBCN control framework as implemented in this repository: the network,
the controller architecture and its temporal model, the phenotype targets, the efficacy metric, and
the reference results. It is written to be read alongside the code in `bbcn/` and reproduced via
`reproduce.py`.

## 1. The network model
BBCN represents a breast tumour as a **135-node Boolean network** spanning 24 signalling pathways
(`bbcn/pathways.py`, `bbcn/harness.py`). Each node is a gene/protein in {0,1}; pathway update rules
define the synchronous Boolean dynamics. Patient states are initialised by **per-gene median
binarization** of cohort expression (1 if above the gene's cohort median, else 0), restricted to the
135 nodes. Two independent cohorts are used: TCGA-BRCA (1082 patients) and METABRIC (1980), shipped
in `data/binarized/`.

## 2. Control hierarchy
The controller is organised into three tiers, each operating on a distinct unit:

- **Tier 3 — phenotype sequencing.** The target cell-fate phenotypes (status goals), in order:
  Resistance-OFF, Apoptosis-ON, Proliferation-OFF.
- **Tier 2 — cycle (unit: pathways).** Selects a candidate pathway set, ranks pathways
  (upstream-first or by mismatch), and passes 1–2 pathways to Tier 1. No kernel design here.
- **Tier 1 — day / network step (unit: kernels).** Designs the minimal intervention "kernel" for
  one pathway at a time via semi-tensor-product forward selection (`bbcn/stp.py`), schedules the
  1–2 pathways sequentially, and administers.

### 2.1 Temporal model
- Tier-1 interval = 1 day = one network-simulation step.
- Tier-2 interval = 1 cycle = 21 days ≈ 3 weeks. Within a cycle: day 1 = design + administer;
  days 2–21 = network relaxation to a settled state under the held kernel.
- Horizons: primary 8 cycles (~6 months); extended 17 cycles (~1 year).
- Re-decisioning occurs at cycle boundaries; the kernel is held constant within a cycle.

### 2.2 Induction and maintenance; caps
- Induction: up to 2 pathways concurrently → 4–5 target nodes.
- Maintenance: 1 pathway → 2–3 target nodes.
- **No accumulation:** the held kernel is re-decided each cycle and never exceeds the cap; the bus
  state carries forward but prior-cycle pins are not retained.
- **No node forcing:** a phenotype counts as achieved only at a genuine free-dynamics fixed point.

## 3. Phenotype targets (two families)
Targets are full-bus (135-node) reference states in `data/stage_targets/`, in two families:
- **`ideal`** — pure standalone phenotype goals (used for isolated per-phenotype controllability).
- **`compatible`** — protection-aware targets that, by construction, preserve earlier phenotypes
  (Apoptosis compatible with Resistance; Proliferation compatible with both). Used for sequenced
  control with state-preservation.
The loader exposes both via `harness._load_stage_targets(family=...)`.

## 4. Efficacy metric (TRACE)
**TRACE — TRajectory-based Apoptosis Control Efficacy.** For each phenotype, let m(t) be the count of
unsatisfied target nodes (mismatch). Define A = (1/(H·m(0)))·Σ m(t) ∈ [0,1] (0 = instant recovery,
1 = no improvement); TRACE is the weighted sum of A across phenotypes. Outcomes are classified as
RECOVERED (genuine fixed point + phenotype reached), STALLED (mismatch flat for ≥2 cycles ≈ 6 weeks),
or TIMED-OUT (still improving at horizon H).

## 5. Reference results
Controller: capped, no-accumulation, 2 pathways / cap 4 / upstream-first / horizon 8, initialised from
the patient state. Full results in `results/reference/REFERENCE_TABLE.csv`.

| Regime | Resistance | Apoptosis | Proliferation | All three (sequenced) |
|---|---|---|---|---|
| Isolated (ideal targets) | ~4–8% | ~85–87% | ~94–95% | — |
| Sequenced + state-preservation | ~4–8% | single digits | single digits | ~2–4% |

**Interpretation.** Individually, apoptosis and proliferation are highly controllable with a small
(~1-target) maintenance kernel, replicated across both cohorts. Resistance-OFF is the hardest
phenotype in every regime. Achieving all three in sequence while preserving earlier gains is rare:
the phenotypes share a strongly-connected network core and are individually reachable but jointly
near-mutually-exclusive under a realistic targeted-therapy budget. See `docs/SCOPE.md` for what these
results do and do not claim (including null survival and the absence of a drug-response validation
cohort).

## 6. Reproducing
`python reproduce.py` regenerates the table and figures from the shipped data and verifies them
against the committed reference within tolerance. See `README.md`.
