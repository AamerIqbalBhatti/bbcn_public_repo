"""
bbcn.cde_diagnosis — INCREMENT 2 of the CDE port.

Faithful port of the diagnosis + session-memory pieces that feed the
supervisory FSM (bbcn.cde_supervisor):

  - classify_failure          <- local_classify_failure        (analyze_phenotypic_outcome_v1.m:303)
  - stage_stability_class     <- stability block               (analyze_phenotypic_outcome_v1.m:148)
  - outcome_signature         <- local_outcome_signature       (single_patient_harness.m)
  - SessionMemory.update      <- session-loop bookkeeping       (single_patient_harness.m:845-986)

Window sizes (verbatim): warnWindowH = 3 ; confirmH = StageConfirmWindow (default 3).

NOTE on dependencies (honest): classify_failure needs `stage_stability_class`,
which needs kernel-convergence signals from the inner simulator
(kernel_exact_converged, stage_pass_window, kernel_core_ratio vs coreThr). Those
come from the inner free-dynamics run, which run_patient already performs — wiring
them through is part of increment 4. Here they are explicit inputs so this layer
is faithful and testable in isolation.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional

WARN_WINDOW_H = 3
DEFAULT_CONFIRM_H = 3


# ---- stage stability class (analyze_phenotypic_outcome_v1.m:148) ----
def stage_stability_class(stage_status_final: str,
                          stage_pass_window: bool,
                          kernel_exact_converged: bool,
                          kernel_stable_enough: bool) -> str:
    if stage_pass_window and kernel_exact_converged:
        return "STABILIZED_FIXED_KERNEL"
    if stage_pass_window and kernel_stable_enough:
        return "STABILIZED_CORE_KERNEL"
    if str(stage_status_final).upper() == "PASS":   # stage_hit_once
        return "TRANSIENT_STAGE_HIT"
    return "NOT_ACHIEVED"


# ---- failure classification (local_classify_failure, verbatim cascade) ----
def classify_failure(stage_status_init: str,
                     stage_status_final: str,
                     stability_class: str,
                     hit_forbidden_markers: List[str],
                     hit_contradictions: List[str],
                     failed_required_markers: List[str],
                     pathway_focus_match: str) -> str:
    si = str(stage_status_init).upper()
    sf = str(stage_status_final).upper()
    if sf == "PASS":
        if stability_class in ("STABILIZED_FIXED_KERNEL", "STABILIZED_CORE_KERNEL"):
            return "STAGE_ACHIEVED" if si == "FAIL" else "STAGE_PRESERVED"
        if stability_class == "TRANSIENT_STAGE_HIT":
            return "TRANSIENT_STAGE_HIT"
    if sf == "FAIL":
        if hit_forbidden_markers or hit_contradictions:
            return "CONTRADICTORY_STATE_PERSISTS"
        if failed_required_markers:
            return "REQUIRED_MARKERS_NOT_REACHED"
        if str(pathway_focus_match).upper() == "MISMATCH":
            return "WRONG_PATHWAY_FOCUS"
        return "FEASIBLE_NOT_REACHED"
    return "UNRESOLVED"   # WARN / anything else


# ---- outcome signature (local_outcome_signature) ----
# Source builds the signature from (stage final status, #failed_required,
# #forbidden, #contradictions). The literal string format of the assembly line
# was not captured from source; this encoding is behaviourally faithful because
# the FSM only tests signature EQUALITY between consecutive sessions, which this
# preserves exactly (same inputs -> same string). Confirm exact format if the
# logged warn_signature strings are ever compared verbatim.
def outcome_signature(stage_status_final: str,
                      n_failed_required: int,
                      n_forbidden: int,
                      n_contradictions: int) -> str:
    return f"{str(stage_status_final).upper()}|R{n_failed_required}|F{n_forbidden}|C{n_contradictions}"


@dataclass
class SessionMemory:
    """Rolling per-session signals consumed by cde_supervisor.supervisory_mode."""
    confirm_h: int = DEFAULT_CONFIRM_H
    stage_status_hist: List[str] = field(default_factory=list)
    warn_signature_hist: List[str] = field(default_factory=list)
    selected_pathway_set_log: List[List[str]] = field(default_factory=list)
    transient_hit_count: int = 0
    stagnant_warn_runlen: int = 0

    def update_pre_mode(self, *, failure_class: str, curr_stage_status: str,
                        signature: str, selected_paths: List[str]) -> None:
        """
        Updates performed in the session loop BEFORE local_supervisory_mode is
        called (single_patient_harness.m:869-937). transient_hit_count is set
        here; stagnant_warn_runlen is updated AFTER the mode (see update_post_mode),
        matching the one-step lag in the source.
        """
        self.warn_signature_hist.append(signature)
        if len(self.warn_signature_hist) > WARN_WINDOW_H:
            self.warn_signature_hist = self.warn_signature_hist[-WARN_WINDOW_H:]

        self.stage_status_hist.append(str(curr_stage_status))
        if len(self.stage_status_hist) > self.confirm_h:
            self.stage_status_hist = self.stage_status_hist[-self.confirm_h:]

        self.selected_pathway_set_log.append(list(selected_paths))

        if str(failure_class).upper() == "TRANSIENT_STAGE_HIT":
            self.transient_hit_count += 1
        else:
            self.transient_hit_count = 0

    def update_post_mode(self, super_mode: str) -> None:
        """stagnant_warn_runlen update (single_patient_harness.m:976-984)."""
        if super_mode in ("STAGNANT_WARN", "ESCAPE_STAGNATION"):
            self.stagnant_warn_runlen += 1
        else:
            self.stagnant_warn_runlen = 0
