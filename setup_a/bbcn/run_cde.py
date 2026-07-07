"""
bbcn.run_cde — INCREMENT 4b: CDE-driven outer loop, corrected ordering.

Fix vs 4: diagnose each session's OUTCOME (post-control final status) and let it
drive the NEXT session's mode+selection (MATLAB feedback order), and use a real
free-dynamics fixed-point test for kernel_exact_converged (not a one-step proxy).
"""
from __future__ import annotations
from typing import List, Dict
from bbcn import harness as H
from bbcn.cde_supervisor import supervisory_mode, CDEMid, is_good_enough_failure_class
from bbcn.cde_diagnosis import (SessionMemory, classify_failure,
                                stage_stability_class, outcome_signature)
from bbcn.cde_selection import update_control_selected_paths, _ustable
from bbcn.cde_proposal import recommend_control_action, select_active_pathway_set

CONFIRM_WINDOW = 3


def _core_fringe(good_sets):
    if not good_sets:
        return [], []
    sets = [_ustable(s) for s in good_sets]
    core = list(sets[0])
    for s in sets[1:]:
        ss = set(s); core = [p for p in core if p in ss]
    allp = _ustable(*sets); cs = set(core)
    return core, [p for p in allp if p not in cs]


def _is_free_fixed_point(bus: Dict[str, int]) -> bool:
    """Genuine free-dynamics fixed point: one synchronous all-pathway rule update leaves bus unchanged."""
    nb = dict(bus)
    for pname, pdef in H.PATHWAYS.items():
        x = [int(bus.get(nd, 0)) for nd in pdef['nodes']]
        nxt = [int(v) for v in pdef['rules'](x, bus)]
        for i, nd in enumerate(pdef['nodes']):
            nb[nd] = nxt[i]
    return nb == bus


def _find_kernel(bus, pw, stage):
    """Suggested minimal kernel for one pathway toward its stage target (given current bus)."""
    from bbcn import stp
    pdef = H.PATHWAYS[pw]; nd = pdef['nodes']; n = len(nd)
    x0 = [int(bus.get(x, 0)) for x in nd]
    tgt = H.stage_target_for_pathway(stage, pw)
    exts = H._pathway_externals(pw)
    ck = (pw, tuple(int(bus.get(e, 0)) for e in exts))
    L = stp.generate_L(pdef['rules'], n, bus, cache_key=ck)
    sel = (stp.forward_select_kernel(L, tgt, x0, max_k=2, rule_fn=pdef['rules'], ext=bus)
           or stp.forward_select_kernel(L, tgt, x0, max_k=3, rule_fn=pdef['rules'], ext=bus))
    if not sel:
        return {}
    return {nd[s - 1]: tgt[s - 1] for s in sel}


def _evolve_held(bus, stage, held, max_steps=15):
    """Evolve all pathways free, holding the accumulated kernel pins, to convergence."""
    for _ in range(max_steps):
        nb = dict(bus)
        for p, pdef in H.PATHWAYS.items():
            x = [int(bus.get(x_, 0)) for x_ in pdef['nodes']]
            nxt = [int(v) for v in pdef['rules'](x, bus)]
            for i, x_ in enumerate(pdef['nodes']):
                nb[x_] = nxt[i]
        for hn, hv in held.items():
            nb[hn] = hv
        if nb == bus:
            break
        bus = nb
    return bus


def _propose(bus, stage, fc):
    rec_mode, candidates, _ = recommend_control_action(fc, stage)
    return rec_mode, select_active_pathway_set(rec_mode, candidates, stage, bus)


def run_patient_cde(init_bus, T: int = 5, max_steps: int = 10, verbose=False,
                    kernel_method: str = 'stabilize') -> dict:
    H.KERNEL_METHOD = kernel_method
    bus = {nd: int(init_bus.get(nd, 0)) for nd in H.ALL_NODES}
    kernel_log: List[dict] = []
    stage_summary = {}
    pathway_history = []

    for stage in H.STAGE_SEQUENCE:
        mem = SessionMemory(confirm_h=CONFIRM_WINDOW)
        protect = H.STAGE_PROTECT.get(stage, [])
        good_sets: List[List[str]] = []
        success_support: List[str] = []
        final_hist: List[str] = []
        promoted = False
        held: Dict[str, int] = {}   # accumulated kernel pins (batch mode, held)

        stage_init_status = H.evaluate_phenotype(bus, stage)
        # session 1 selection: propose from the stage-initial diagnosis
        fc0 = classify_failure(stage_init_status, stage_init_status,
                               stage_stability_class(stage_init_status, False, False, False),
                               H.hit_forbidden(bus, stage), H.hit_contradictions(bus, stage),
                               H.failed_required(bus, stage), "MATCH")
        _, control_selected = _propose(bus, stage, fc0)

        for k in range(1, T + 1):
            # 1) run THIS session's control — VERBATIM from single_patient_harness:
            # local_schedule_active_pathways splits the selected set into 1-pathway
            # EXPOSURES (MaxActivePathwaysPerExposure=1, sequential), processed IN ORDER
            # with the bus carrying forward. Each exposure runs the inner sim with that
            # ONE pathway selected => that pathway is active and pins its minimal kernel;
            # all 135 nodes update synchronously by rules. MaxSteps=10.
            if control_selected:
                for batch in H._schedule_batches(control_selected, 1):
                    bus = H.simulate_inner(bus, batch, stage, k,
                                           kernel_log, 10, coactive=False)

            # 2) diagnose the OUTCOME
            final_status = H.evaluate_phenotype(bus, stage)
            final_hist.append(final_status)
            fr, hf, hc = (H.failed_required(bus, stage), H.hit_forbidden(bus, stage),
                          H.hit_contradictions(bus, stage))
            pass_window = (len(final_hist) >= CONFIRM_WINDOW and
                           all(s == "PASS" for s in final_hist[-CONFIRM_WINDOW:]))
            converged = _is_free_fixed_point(bus)
            stab = stage_stability_class(final_status, pass_window, converged, converged)
            fc = classify_failure(stage_init_status, final_status, stab, hf, hc, fr, "MATCH")

            # 3) memory + mode (driven by the outcome)
            sig = outcome_signature(final_status, len(fr), len(hf), len(hc))
            mem.update_pre_mode(failure_class=fc, curr_stage_status=final_status,
                                signature=sig, selected_paths=control_selected)
            cde = CDEMid(stage_status_init=stage_init_status, stage_status_final=final_status,
                         failed_required_markers=fr, hit_forbidden_markers=hf,
                         hit_contradictions=hc, failure_class=fc,
                         protected_status={ph: H.evaluate_phenotype(bus, ph) for ph in protect})
            mode = supervisory_mode(cde, mem.stage_status_hist, mem.warn_signature_hist,
                                    mem.selected_pathway_set_log, k, protect,
                                    mem.transient_hit_count, mem.stagnant_warn_runlen)
            pathway_history.append({"stage": stage, "session": k, "mode": mode,
                                    "selected": list(control_selected)})
            mem.update_post_mode(mode)

            if fc in ("STABILIZED_FIXED_KERNEL", "STAGE_PRESERVED", "STAGE_ACHIEVED"):
                success_support = _ustable(success_support, control_selected)
            if is_good_enough_failure_class(fc):
                good_sets.append(list(control_selected))

            # 4) promotion
            if mode == "PROMOTE" or pass_window:
                promoted = True
                break

            # 5) choose NEXT session's selection from this outcome
            rec_mode, proposed = _propose(bus, stage, fc)
            core, fringe = _core_fringe(good_sets)
            control_selected = update_control_selected_paths(
                mode, control_selected, proposed, core, fringe,
                protected_core=[], success_support=success_support,
                contrad_nodes=hc, stagnant_warn_count=mem.stagnant_warn_runlen)
            if not control_selected and proposed:
                control_selected = core if core else proposed

        final_status = H.evaluate_phenotype(bus, stage, scheme='v2')   # SCORE with clean verdict
        stage_summary[stage] = {"final_status": final_status,
                                "passed": (final_status == "PASS") or promoted,
                                "promoted": promoted}

    return {"final_bus": bus, "stage_summary": stage_summary,
            "kernel_history": kernel_log, "pathway_history": pathway_history,
            "terminal_pass": H.evaluate_phenotype(bus, "Terminal", scheme='v2') == "PASS"}
