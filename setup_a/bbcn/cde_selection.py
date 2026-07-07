"""
bbcn.cde_selection — INCREMENT 3 of the CDE port.

Faithful port of local_update_control_selected_paths
(single_patient_harness.m:2208-2402): the switch that turns a supervisory mode
plus the session's pathway pools into that session's selected pathway set.

Pools (inputs):
  current_paths    : the set used last session
  proposed_paths   : the CDE mismatch-driven proposal (select_active_pathway_set_v1)
  core_paths       : learned stable pathways for the stage (local_update_core_fringe)
  fringe_paths     : learned exploratory pathways
  protected_core   : union of phenotype_core_map for protected phenotypes
  success_support  : pathways accumulated from sessions that achieved/preserved the stage
  contrad_nodes    : currently-hit contradiction nodes (for PROTECT_THEN_RESOLVE)

maxFringeKeep = 2 (verbatim).

DEPENDENCIES STILL TO PORT (increment 3b / 4 — honest):
  - proposed_paths  <- select_active_pathway_set_v1.m (the mismatch-priority proposal)
  - core/fringe     <- local_update_core_fringe(good_selected_sets)
  - success_support <- accumulated in the session loop on STAGE_PRESERVED/ACHIEVED/
                       STABILIZED_FIXED_KERNEL
  - NODE_TO_PATHWAYS is partial (only the entries read from source so far)
These feeders are required before this switch produces MATLAB-matching selections;
the switch logic itself is complete and faithful.
"""
from __future__ import annotations
from typing import List, Dict, Optional
import math

MAX_FRINGE_KEEP = 2

# local_map_node_to_pathways — partial (entries confirmed from source; extend as read)
NODE_TO_PATHWAYS: Dict[str, List[str]] = {
    "MYC":   ["MAPK", "Wnt", "TF", "RTK_EGFR"],
    "FOXO3": ["AKT_Signaling", "PI3K"],
    "TP53":  ["Apoptosis_Regulatory", "DNA_Repair"],
}


def _ustable(*seqs: List[str]) -> List[str]:
    """unique(..., 'stable'): order-preserving concat + dedup, dropping empties."""
    seen, out = set(), []
    for seq in seqs:
        for x in seq:
            if x and x not in seen:
                seen.add(x); out.append(x)
    return out


def _setdiff_stable(a: List[str], remove: List[str]) -> List[str]:
    rm = set(remove)
    return [x for x in _ustable(a) if x not in rm]


def update_control_selected_paths(super_mode: str,
                                   current_paths: List[str],
                                   proposed_paths: List[str],
                                   core_paths: List[str],
                                   fringe_paths: List[str],
                                   protected_core: List[str],
                                   success_support: List[str],
                                   contrad_nodes: List[str],
                                   max_fringe_keep: int = MAX_FRINGE_KEEP,
                                   stagnant_warn_count: int = 0,
                                   node_to_pathways: Optional[Dict[str, List[str]]] = None
                                   ) -> List[str]:
    cur = _ustable(current_paths)
    prop = _ustable(proposed_paths)
    core = _ustable(core_paths)
    fringe = _ustable(fringe_paths)
    pcore = _ustable(protected_core)
    supp = _ustable(success_support)
    m = str(super_mode).upper()
    n2p = node_to_pathways if node_to_pathways is not None else NODE_TO_PATHWAYS

    if m == "PROMOTE":
        return cur

    if m == "UNDEFINED_STAGE":
        return cur if cur else []

    if m == "EXPLORE":
        return _ustable(prop, supp) if prop else _ustable(cur, supp)

    if m == "STABILIZE":
        if not core:
            return cur if cur else prop
        keep = fringe[:min(max_fringe_keep, len(fringe))]
        return _ustable(core, pcore, keep)

    if m == "PROTECT_RECOVERY":
        base = core if core else cur
        keep = fringe[:min(max_fringe_keep, len(fringe))]
        nxt = _ustable(base, pcore, keep, supp)
        return nxt if nxt else prop

    if m == "STAGNANT_WARN":
        base = core if core else cur
        injectN = min(max_fringe_keep, len(fringe))
        inject: List[str] = []
        if injectN > 0:
            start0 = max(stagnant_warn_count - 1, 0) % injectN   # 1-based idx0 -> 0-based
            inject = fringe[start0:min(start0 + injectN, len(fringe))]
        nxt = _ustable(base, pcore, inject, supp)
        return nxt if nxt else prop

    if m == "ESCAPE_STAGNATION":
        if not core:
            core_keep: List[str] = []
        else:
            keepN = max(1, math.ceil(len(core) / 2))
            core_keep = core[:keepN]
        shell = _setdiff_stable(prop, _ustable(pcore, core_keep))
        if len(shell) > 1:
            shell = shell[1:]                          # drop one repeatedly-used member
        shell = shell[:min(max_fringe_keep + 1, len(shell))]
        nxt = _ustable(pcore, core_keep, shell, supp)
        if not nxt:
            nxt = prop if prop else cur
        return nxt

    if m == "OSCILLATORY_STAGNATION":
        base = prop if prop else cur
        base = _setdiff_stable(base, core)             # drop core dominance
        inject = fringe[:min(3, len(fringe))]          # stronger exploration
        return _ustable(base, inject)

    if m == "PROTECT_THEN_RESOLVE":
        contra_paths: List[str] = []
        for nd in _ustable(contrad_nodes):
            contra_paths = _ustable(contra_paths, n2p.get(str(nd).upper(), []))
        nxt = _ustable(pcore, core, prop, contra_paths)
        if len(nxt) < 10:
            nxt = _ustable(nxt, fringe[:min(3, len(fringe))])
        return nxt if nxt else cur

    return cur  # safety
