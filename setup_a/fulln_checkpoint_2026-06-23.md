# BBCN Full-N Checkpoint -- end of Phase 3 (2026-06-23)

All three cohorts at full N (TCGA 1082, METABRIC 1980, ISPY2 988 = 4050 patients).

## E41 monolith (v2 scoring)
| cohort | N | apoptosis | prolif-off | joint | durable |
|---|---|---|---|---|---|
| TCGA | 1082 | 76% | 100% | 13% | 17% |
| METABRIC | 1980 | 74% | 100% | 7% | 15% |
| ISPY2 | 988 | 72% | 100% | 14% | 16% |
| POOLED | 4050 | 74% | 100% | 11% | 16% |

## E40 durability (idealized: 300-step continuous hold)
| cohort | N | resistant | kernel found | switch | strict* | commit |
|---|---|---|---|---|---|---|
| TCGA | 1082 | 19% | 99% | 81% | 68% | 97% |
| METABRIC | 1980 | 20% | 94% | 85% | 52% | 93% |
| ISPY2 | 988 | 24% | 97% | 88% | 70% | 98% |
| POOLED | 4050 | 21% | 96% | 85% | 61% | 95% |

\*strict (CASP3 & AKT1=0) is the circular self-scoring diagnostic, kept as a label only; commit (CASP3 sustained) is the honest figure.

## Clock arm (durable = CASP3 across the 4-week drug-free observation window)
| cohort | kerneled | clock-continuous (idealized) | pulsed 4-on/3-off | capped + escalation |
|---|---|---|---|---|
| TCGA | 199 | 95% | 86% | 84% |
| METABRIC | 376 | 90% | 83% | 77% |
| ISPY2 | 231 | 96% | 86% | 84% |
| POOLED | 806 | 93% | 85% | 80% |

Clock-continuous (93%) reproduces the old commit (95%) -> clock engine faithful. Pulsing costs ~8 pts; the cap costs ~5 more.

## Induction strategy (capped, full N) -- front-loaded vs escalate-from-1
| cohort | kerneled | escalate-from-1 | front-loaded |
|---|---|---|---|
| TCGA | 199 | 84% | 84% |
| METABRIC | 376 | 77% | 77% |
| ISPY2 | 231 | 84% | 84% |
| POOLED | 806 | 80% | 80% |

Identical: the induction ramp does not change the durability endpoint. Default kept as escalate-from-1 (parsimony); `front_loaded` flag available.
