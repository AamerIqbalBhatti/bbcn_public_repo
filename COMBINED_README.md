# BBCN Combined Repository (v22)

This repository is the **superset of v8 and v21**. It carries the complete history of
the project in one tree, with both networks and both result sets reproducible from a
single driver.

## What is merged here

| Source | What it contributed | Status in v22 |
|--------|--------------------|---------------|
| v8  (repository 1) | Paper 1 (bioRxiv preprint), the **bare 135-node** network, paper-1 locked numbers | preserved, unchanged |
| v21 (repository 2) | The **repaired** overlay (PHLPP, COMMIT, U), paper 2, the repaired run-scripts | preserved, unchanged |
| this build | Recovered **repaired numbers** locked into `results/numbers2/`, the unified driver, this note | new |

v21 already contained every file in v8 (102/102 paths). The only shared file that
differed was the daily ledger, where v21 held the newer through-E42 version; the v8
snapshot (through E35) is preserved at `docs/history/`. Nothing from either repo is lost.

## Two networks, one tree

The repaired network is **not a fork**. Every core module is byte-identical to v8.
The repair is an overlay applied at run time:

```python
from bbcn import repaired_branch as RB
RB.apply()        # swaps in: p53-MDM2-ATM sensor, AKT1-FOXO3a-PHLPP loop, COMMIT latch
```

* **Bare network** (paper 1): 135 nodes, `pathways.py` as written. Globally oscillatory,
  no durable apoptotic fixed point, resistance flippable only by genotoxic input.
* **Repaired network** (paper 3): 135 + PHLPP = 136 signalling nodes; the staged-controller
  branch additionally carries the COMMIT latch and the per-patient U flag.

## Result sets (both locked, macro-sourced)

| Set | Network | Paper | Generator | Lock file |
|-----|---------|-------|-----------|-----------|
| `results/numbers/`  | bare     | paper 1 | `run_chunked.py --emit` + `generate_numbers.py` | `preprint/numbers.tex` |
| `results/numbers2/` | repaired | paper 3 | `generate_numbers2.py` (reads `run_outputs/`)   | `results/numbers2/numbers2.tex` |

`results/numbers2/run_outputs/` holds the persisted outputs of the four repaired
pipelines (E36/E38/E40), recovering the results the mid-session container reset had
dropped. The repaired numbers are now drift-proof exactly like the bare numbers.

## One-button reproduce

```bash
python reproduce_combined.py bare       # paper-1 numbers
python reproduce_combined.py repaired   # paper-3 numbers (E41 monolith is the slow step)
python reproduce_combined.py all
```

## Provenance of the locked repaired numbers (regenerated this build)

| Block | Script | Scope | Headline |
|-------|--------|-------|----------|
| E36 | `run_full_cohorts_repaired.py` | full N | resistance == AKT-FOXO3a uncoupling, 4049/4050 (99.98%) |
| E38 | `run_designed_kernel_cohorts.py` | N=200/cohort | kernel found 98-100%, nominates PI3K/AKT/mTOR |
| E40 | `run_durability_full.py` | full N | commitment durability 92-95%; strict 8-14% |
| E41 | `run_repaired_monolith.py` | full N (re-verified N=200/cohort) | durable apoptotic FP 15-17% (bare net 0%) |
