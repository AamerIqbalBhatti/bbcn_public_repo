"""
bbcn.cde_proposal — INCREMENT 3b of the CDE port.

The proposal layer that produces `proposed_paths` for the selection switch:
  recommend_control_action  <- recommend_next_control_action_v1.m  (failure_class -> mode + candidate families)
  pathway_priority_score    <- local_pathway_priority_score / local_stage_weight / local_scc_bonus
  select_active_pathway_set <- select_active_pathway_set_v1.m       (rank by score, mode-truncate)

Scoring (verbatim magnitudes):
  score = stage_weight(2.0 strong / 1.0 support / 0, EXACT name match)
        + scc_bonus(0.75 if in stage bigSCC, case-insensitive)
        + current pathway mismatch_fraction
Mode truncation: single_pathway->1 ; coordinated/advance_or_parallel/hold_stage/
protective_recovery->up to MaxPathways(=10) ; advance_stage->0.

Candidate families (EXPECTED/CONTRADICTION_BREAKERS/ENABLING/STAGE_PROTECT) are
reused from harness.py — they are already the faithful ports of
local_expected_pathways_for_phenotype / local_contradiction_breakers /
local_enabling_pathways.
"""
from __future__ import annotations
from typing import List, Dict, Tuple
from bbcn import harness as H

MAX_PATHWAYS = 10

# local_stage_weight: strong (w=2.0) / support (w=1.0), EXACT match (== in source)
_STAGE_STRONG: Dict[str, List[str]] = {
    "Resistance_OFF":    ["RTK", "RTK_EGFR", "PI3K", "AKT", "mTOR", "JAK_STAT", "NFkB"],
    "Apoptosis_ON":      ["Apoptosis_Regulatory", "Apoptosis_Intrinsic", "Apoptosis_Extrinsic"],
    "Proliferation_OFF": ["CellCycle", "CDK_CellCycle"],
    "Invasion_OFF":      ["Wnt", "Hormone", "Others"],
}
_STAGE_SUPPORT: Dict[str, List[str]] = {
    "Resistance_OFF":    ["MAPK", "TF"],
    "Apoptosis_ON":      ["AKT", "mTOR", "NFkB", "JAK_STAT"],
    "Proliferation_OFF": ["MAPK", "RTK", "PI3K", "TF", "Hippo", "Wnt"],
    "Invasion_OFF":      ["NFkB", "Hippo", "TF"],
}
# local_scc_bonus: bigSCC sets (case-insensitive membership -> +0.75)
_SCC_BIG: Dict[str, List[str]] = {
    "Resistance_OFF":    ["RTK", "RTK_EGFR", "PI3K", "AKT", "AKT_Signaling", "AKT_Survival",
                          "mTOR", "JAK_STAT", "NFkB", "TF", "MAPK"],
    "Apoptosis_ON":      ["Apoptosis_Regulatory", "Apoptosis_Intrinsic", "Apoptosis_Extrinsic",
                          "AKT", "AKT_Signaling", "AKT_Survival", "mTOR", "NFkB", "JAK_STAT", "TF"],
    "Proliferation_OFF": ["CellCycle", "CDK_CellCycle", "MAPK", "RTK", "RTK_EGFR", "RTK_Insulin",
                          "PI3K", "TF", "Wnt", "Hippo"],
}


def _stage_weight(pw: str, stage: str) -> float:
    if pw in _STAGE_STRONG.get(stage, []):      # exact == match (source semantics)
        return 2.0
    if pw in _STAGE_SUPPORT.get(stage, []):
        return 1.0
    return 0.0


def _scc_bonus(pw: str, stage: str) -> float:
    big = {b.upper() for b in _SCC_BIG.get(stage, [])}
    return 0.75 if pw.upper() in big else 0.0


def pathway_mismatch_fraction(bus: Dict[str, int], stage: str, pw: str) -> float:
    """Fraction of the pathway's nodes whose bus value != stage target value."""
    nodes = H.PATHWAYS[pw]["nodes"]
    target = H.stage_target_for_pathway(stage, pw)
    if not nodes:
        return 0.0
    return sum(1 for i, nd in enumerate(nodes) if int(bus.get(nd, 0)) != target[i]) / len(nodes)


def pathway_priority_score(bus: Dict[str, int], pw: str, stage: str) -> float:
    return _stage_weight(pw, stage) + _scc_bonus(pw, stage) + pathway_mismatch_fraction(bus, stage, pw)


# recommend_next_control_action_v1: failure_class -> (mode, candidate pathways, protect)
def recommend_control_action(failure_class: str, stage: str) -> Tuple[str, List[str], List[str]]:
    protect = H.STAGE_PROTECT.get(stage, [])
    stage_paths = H.EXPECTED_PATHWAYS.get(stage, [])
    fc = str(failure_class).upper()
    if fc == "STAGE_ACHIEVED":
        return "advance_stage", [], protect
    if fc == "STAGE_PRESERVED":
        return "advance_or_parallel", [], protect
    if fc == "TRANSIENT_STAGE_HIT":
        return "hold_stage", stage_paths, protect
    if fc == "WRONG_PATHWAY_FOCUS":
        return "coordinated_pathway_set", stage_paths, protect
    if fc == "CONTRADICTORY_STATE_PERSISTS":
        return "coordinated_pathway_set", H.CONTRADICTION_BREAKERS.get(stage, []), protect
    if fc == "REQUIRED_MARKERS_NOT_REACHED":
        return "coordinated_pathway_set", H.ENABLING_PATHWAYS.get(stage, []), protect
    if fc == "FEASIBLE_NOT_REACHED":
        return "coordinated_pathway_set", stage_paths, protect
    return "coordinated_pathway_set", stage_paths, protect   # default


# select_active_pathway_set_v1: rank candidates by score (desc, stable), truncate by mode
def select_active_pathway_set(mode: str, candidates: List[str], stage: str,
                              bus: Dict[str, int], max_pathways: int = MAX_PATHWAYS) -> List[str]:
    seen, cand = set(), []
    for c in candidates:
        if c and c.lower() != "auto" and c not in seen:
            seen.add(c); cand.append(c)
    # stable descending sort by score (Python sort is stable; negate score keeps ties in order)
    cand = sorted(cand, key=lambda pw: -pathway_priority_score(bus, pw, stage))
    m = str(mode)
    if m == "single_pathway":
        keepN = min(1, len(cand))
    elif m in ("coordinated_pathway_set", "advance_or_parallel", "hold_stage", "protective_recovery"):
        keepN = min(max_pathways, len(cand))
    elif m == "advance_stage":
        keepN = 0
    else:
        keepN = min(max_pathways, len(cand))
    return cand[:keepN]
