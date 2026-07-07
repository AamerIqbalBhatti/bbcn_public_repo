"""
bbcn.cde_supervisor — faithful port of the MATLAB outer-loop supervisor FSM.

Source: +bbcn118/+pipeline/single_patient_harness.m
  - local_supervisory_mode            (line 2068)
  - local_is_good_enough_failure_class
  - local_protected_phenotype_status

This is INCREMENT 1 of the CDE wiring: the supervisory-mode classifier only.
It is a pure function of the per-session memory signals; building those signals
(transient_hit_count, stage_status_hist, warn_signature_hist, selected-set log,
protected-phenotype status) and the per-mode pathway-selection switch are
increments 2 and 3. Nothing here is wired into run_patient yet — that is
deliberate, so this piece can be validated in isolation against the
`supervisory_mode` column of paper_outputs_tcbb/CohortPathwayHistory.csv.

Mode vocabulary (verbatim from source):
  PROMOTE, UNDEFINED_STAGE, OSCILLATORY_STAGNATION, PROTECT_RECOVERY,
  STABILIZE, PROTECT_THEN_RESOLVE, ESCAPE_STAGNATION, STAGNANT_WARN, EXPLORE
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict

# local_is_good_enough_failure_class — verbatim goodSet
_GOOD_ENOUGH = {"UNSTABLE_PASS", "TRANSIENT_STAGE_HIT",
                "STAGE_ACHIEVED", "STAGE_PRESERVED"}


def is_good_enough_failure_class(fc: str) -> bool:
    return str(fc).upper() in _GOOD_ENOUGH


@dataclass
class CDEMid:
    """Mirror of CDE_mid: current-session diagnosis fed to the FSM."""
    recommended_mode: str = ""          # NextControl.recommended_mode
    stage_status_init: str = "UNKNOWN"  # PhenotypeOutcome.stage_status_init
    stage_status_final: str = "UNKNOWN" # PhenotypeOutcome.stage_status_final
    failed_required_markers: List[str] = field(default_factory=list)
    hit_forbidden_markers: List[str] = field(default_factory=list)
    hit_contradictions: List[str] = field(default_factory=list)
    failure_class: str = "UNRESOLVED"
    # protected-phenotype current statuses, e.g. {'Resistance_OFF': 'PASS'}
    protected_status: Dict[str, str] = field(default_factory=dict)


def protected_phenotype_status(cde: CDEMid, protected_phenotypes: List[str]):
    """
    Port of local_protected_phenotype_status. A protected phenotype is 'ok'
    iff its current status is PASS. Returns (all_ok, degraded_list).
    (If no protected phenotypes: all_ok=True, degraded=[].)
    """
    if not protected_phenotypes:
        return True, []
    degraded = [ph for ph in protected_phenotypes
                if str(cde.protected_status.get(ph, "UNKNOWN")).upper() != "PASS"]
    return (len(degraded) == 0), degraded


def supervisory_mode(cde: CDEMid,
                     stage_status_hist: List[str],
                     warn_signature_hist: List[str],
                     selected_pathway_set_log: List[List[str]],
                     k: int,
                     protected_phenotypes: List[str],
                     transient_hit_count: int,
                     stagnant_warn_runlen: int) -> str:
    """
    Faithful port of local_supervisory_mode (single_patient_harness.m:2068).
    `k` is the 1-based session index; selected_pathway_set_log is 0-based here,
    so log[k-2]/log[k-3] correspond to MATLAB's {k-1}/{k-2}.
    """
    # 1) explicit promotion
    if str(cde.recommended_mode).lower() == "advance_stage":
        return "PROMOTE"

    # 2) undefined-stage guard
    si, sf = str(cde.stage_status_init).upper(), str(cde.stage_status_final).upper()
    reqN, forbN, contraN = (len(cde.failed_required_markers),
                            len(cde.hit_forbidden_markers),
                            len(cde.hit_contradictions))
    if si == "UNKNOWN" and sf == "UNKNOWN" and reqN == 0 and forbN == 0 and contraN == 0:
        return "UNDEFINED_STAGE"

    # 3) oscillatory stagnation: repeated transient hits (the apoptosis limit cycle)
    if transient_hit_count >= 3:
        return "OSCILLATORY_STAGNATION"

    # 4) good-enough failure classes
    fc = str(cde.failure_class)
    if is_good_enough_failure_class(fc):
        all_ok, _ = protected_phenotype_status(cde, protected_phenotypes)
        return "STABILIZE" if all_ok else "PROTECT_RECOVERY"

    # 4.5) protected contradiction resolution
    if fc.upper() == "UNRESOLVED":
        all_ok, _ = protected_phenotype_status(cde, protected_phenotypes)
        has_contra = len(cde.hit_contradictions) > 0
        if (not all_ok) and has_contra:
            return "PROTECT_THEN_RESOLVE"

    # 5) stagnant bad regime: last-2 sessions WARN/FAIL, same status+signature+set
    ss = [str(s).upper() for s in stage_status_hist]
    ws = [str(s) for s in warn_signature_hist]
    if len(ss) >= 2 and len(ws) >= 2:
        recent_bad = all(s in ("WARN", "FAIL") for s in ss[-2:])
        same_status = ss[-1] == ss[-2]
        same_sig = ws[-1] == ws[-2]
        same_set = False
        if k >= 3 and len(selected_pathway_set_log) >= 2:
            a = sorted({p for p in selected_pathway_set_log[k - 2] if p})
            b = sorted({p for p in selected_pathway_set_log[k - 3] if p})
            same_set = (a == b)
        if recent_bad and same_status and same_sig and same_set:
            return "ESCAPE_STAGNATION" if stagnant_warn_runlen >= 6 else "STAGNANT_WARN"

    # 6) default
    return "EXPLORE"
