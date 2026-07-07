"""
bbcn.harness
============
Faithful two-tier port of the MATLAB BBCN control pipeline.

OUTER tier  (single_patient_harness.m): for each of the 4 stages, run up to
    T sessions. Each session: run the clinical decision engine (CDE) to pick the
    active pathway set from the current phenotype diagnosis, schedule the set
    into E_max batches, call the inner simulator, then reassess with the CDE and
    decide stage promotion.

INNER tier  (simulate_bbcn118_sequential.m): MaxSteps free-dynamics steps. Each
    step: resolve the active pathway dynamically by mismatch, build the literal
    STP transition matrix L for it under current externals, search a minimal
    kernel (<=2 escalating to 3), pin the kernel, evolve everything else under
    the true Boolean rules, synchronise the bus. Convergence = bus fixed point.

CDE  (+control package): analyze_phenotypic_outcome_v1 -> recommend_next_control_
    action_v1 -> select_active_pathway_set_v1, plus the phenotype evaluator
    (evaluate_biological_state_bbcn118.m).

Design principles enforced: NO node forcing (only minimal kernels pinned, rest
free), all heuristic knowledge in the model/targets, control technique generic.

Validated against paper_outputs_tcbb/*.csv.
"""

from typing import List, Dict, Tuple, Optional
import json
import os
import numpy as np

from bbcn.pathways import PATHWAYS, ALL_NODES, NODE_TO_PATHWAY
from bbcn.stp import generate_L, forward_select_kernel
from bbcn import forward_stab as FS

# Kernel-selection method for the CDE path: 'ranked' (preserved heuristic) or
# 'stabilize' (algebraic Theorem-1 global stabilization). Set by run_cde before a run.
KERNEL_METHOD = 'stabilize'

# Pathway-tier resolution mode, mirroring KERNEL_METHOD one tier up:
#   'cde'    (default) : adaptive — resolve over the CDE-narrowed set with the
#                        priority weighting (the canonical, frozen behavior).
#   'static'           : pure max-mismatch benchmark — resolve over the full
#                        stage-eligible pool the CDE narrows from, bypassing the
#                        CDE set-narrowing and the priority bonus.
PATHWAY_MODE = 'cde'

# Read-only kernel-composition capture. Default None = no capture (byte-identical
# behavior). Set to a list to record each pinned kernel as _update_pathway selects
# it; this only observes the existing choice, it never changes any decision.
KERNEL_CAPTURE = None

# Persistent cache for the (x0-independent) stabilize kernels, shared across patients.
_STAB_CACHE = {}

State = Dict[str, int]

# ============================================================
# Full-bus stage targets (from +targets/+global/*_fullbus.json)
# Each stage has a complete 135-node target vector; per-pathway targets are
# this vector projected onto the pathway's nodes. mismatch_fraction is then
# computed over ALL pathway nodes vs the full-bus target (matching MATLAB).
# ============================================================
# Two target FAMILIES (both faithful to MATLAB +targets/+global):
#   'compatible' = protection-aware, order-encoded targets. Apoptosis is rewritten
#                  to preserve Resistance; Proliferation to preserve Apoptosis+Resistance.
#                  Use for SEQUENCED control with state-protection (the clinical regime).
#   'ideal'      = pure standalone phenotype goals (no protection baked in).
#                  Use for ISOLATED per-phenotype testing.
_STAGE_TARGET_FILES_COMPATIBLE = {
    'Resistance_OFF':    'resistance_off_ideal_fullbus.json',
    'Apoptosis_ON':      'apoptosis_stage2_resistance_compatible_fullbus.json',
    'Proliferation_OFF': 'proliferation_off_compatible_with_apoptosis_and_resistance_fullbus.json',
    'Terminal':          'four_stage_terminal_ideal_fullbus.json',
}
_STAGE_TARGET_FILES_IDEAL = {
    'Resistance_OFF':    'resistance_off_ideal_fullbus.json',
    'Apoptosis_ON':      'apoptosis_on_ideal_fullbus.json',
    'Proliferation_OFF': 'proliferation_off_ideal_fullbus.json',
    'Invasion_OFF':      'invasion_off_ideal_fullbus.json',
    'Terminal':          'four_stage_terminal_ideal_fullbus.json',
}

# Default family preserves the original behaviour (compatible / sequenced).
TARGET_FAMILY = 'compatible'
_STAGE_TARGETS: Dict[str, Dict[str, int]] = {}


def _load_stage_targets(family: str = None):
    """Load the stage-target vectors. family in {'compatible','ideal'}.
    Re-loads if the requested family differs from what is currently loaded."""
    global TARGET_FAMILY
    fam = family or TARGET_FAMILY
    if _STAGE_TARGETS and fam == TARGET_FAMILY:
        return
    TARGET_FAMILY = fam
    files = _STAGE_TARGET_FILES_IDEAL if fam == 'ideal' else _STAGE_TARGET_FILES_COMPATIBLE
    here = os.path.dirname(__file__)
    tdir = os.path.join(here, '..', 'data', 'stage_targets')
    _STAGE_TARGETS.clear()
    for stage, fname in files.items():
        path = os.path.join(tdir, fname)
        with open(path) as f:
            raw = json.load(f)
        _STAGE_TARGETS[stage] = {k: int(v) for k, v in raw.items()}


def stage_target_for_pathway(stage: str, pname: str) -> List[int]:
    """Project the stage full-bus target onto a pathway's node order."""
    _load_stage_targets()
    fb = _STAGE_TARGETS[stage]
    return [int(fb.get(nd, 0)) for nd in PATHWAYS[pname]['nodes']]


# Externals each pathway reads: probe the rule fn by toggling each non-owned
# node and seeing if any output changes. Computed once per pathway.
_PATH_EXTERNALS: Dict[str, List[str]] = {}


def _pathway_externals(pname: str) -> List[str]:
    """Discover EVERY bus node (outside the pathway) the rules read, completely.

    Access-tracking (logs every external the rules touch over many seeded-random states)
    UNIONED with the single-flip sensitivity probe. Access-tracking catches externals that
    the probe masks via AND-combinations or internal-state gating (ledger E62). The set is a
    safe over-approximation: an extra external only costs cache reuse, never serves a stale
    kernel, whereas a MISSING external silently reuses a wrong kernel.
    """
    if pname in _PATH_EXTERNALS:
        return _PATH_EXTERNALS[pname]
    import random as _random
    pdef = PATHWAYS[pname]
    nodes = pdef['nodes']
    n = len(nodes)
    own = set(nodes)
    alln = list(ALL_NODES)
    others = [nd for nd in alln if nd not in own]
    reads = set()

    # (1) single-flip sensitivity probe at internal all-0 and all-1
    base_bus = {nd: 0 for nd in alln}
    for x in ([0] * n, [1] * n):
        base = pdef['rules'](x, base_bus)
        for nd in others:
            b = dict(base_bus); b[nd] = 1
            if pdef['rules'](x, b) != base:
                reads.add(nd)

    # (2) access-tracking over seeded-random states (deterministic; catches masked externals)
    class _Tracker(dict):
        def __init__(self, *a, **k): super().__init__(*a, **k); self.acc = set()
        def get(self, k, d=0): self.acc.add(k); return super().get(k, d)
        def __getitem__(self, k):
            self.acc.add(k)
            return super().__getitem__(k) if k in self else 0
    rng = _random.Random(20260622)            # fixed seed -> reproducible external set
    for _ in range(600):
        x = [rng.randint(0, 1) for _ in nodes]
        e = _Tracker({nd: rng.randint(0, 1) for nd in alln})
        try:
            pdef['rules'](x, e)
        except Exception:
            pass
        reads |= e.acc

    reads = [nd for nd in alln if nd in reads and nd not in own]   # real, outside-pathway, stable order
    _PATH_EXTERNALS[pname] = reads
    pdef['externals'] = reads        # explicit, inspectable field on the pathway (ledger E63)
    return reads

# ============================================================
# Phenotype evaluator (evaluate_biological_state_bbcn118.m)
# ============================================================
GLOBAL_REQUIRED = {'CASP3': 1, 'AKT1': 0, 'FOXO3': 1, 'MYC': 0}

PHENO_SCHEMA = {
    'Resistance_OFF': {
        'required': {'AKT1': 0},
        'supporting_any': {'FOXO3': 1, 'MTOR': 0, 'RELA': 0, 'STAT3': 0,
                           'XIAP': 0, 'CFLAR': 0, 'MCL1': 0},
        'forbidden_any': {'AKT1': 1, 'RELA': 1, 'STAT3': 1},
        'contradiction_any': {},
    },
    'Apoptosis_ON': {
        'required': {'CASP3': 1},
        'supporting_any': {'CASP8': 1, 'CASP9': 1, 'BAX': 1, 'BAK1': 1,
                           'CYCS': 1, 'APAF1': 1},
        'forbidden_any': {'XIAP': 1, 'CFLAR': 1},
        'contradiction_any': {'CCND1': 1, 'MYC': 1, 'E2F1': 1},
    },
    'Proliferation_OFF': {
        'required': {'CCND1': 0, 'E2F1': 0},
        'supporting_any': {'MYC': 0, 'RB1': 1, 'CDKN1A': 1, 'CDKN2A': 1,
                           'CDK4': 0, 'CDK6': 0, 'CDK2': 0},
        'forbidden_any': {'CCND1': 1, 'E2F1': 1, 'MYC': 1},
        'contradiction_any': {},
    },
    'Invasion_OFF': {
        'required': {},
        'supporting_any': {'YAP1': 0, 'CTNNB1': 0, 'TCF7L2': 0,
                           'NOTCH1': 0, 'TWIST1': 0},
        'forbidden_any': {'YAP1': 1, 'CTNNB1': 1, 'NOTCH1': 1, 'TWIST1': 1},
        'contradiction_any': {},
    },
    'Terminal': {
        'required': {'CASP3': 1, 'CASP9': 1, 'CCND1': 0, 'E2F1': 0,
                     'MYC': 0, 'AKT1': 0, 'STAT3': 0, 'RELA': 0},
        'supporting_any': {'BAX': 1, 'CYCS': 1, 'FOXO3': 1, 'PTEN': 1},
        'forbidden_any': {'AKT1': 1, 'STAT3': 1, 'RELA': 1, 'MYC': 1},
        'contradiction_any': {},
    },
}

# ------------------------------------------------------------------
# Verdict: STEER with v1, SCORE with v2 (ledger E64 / Phase 2.5B).
# The CDE loop STEERS on the original verdict (evaluate_phenotype default scheme='v1',
# which reads AKT1/CASP9/MYC) so the control trajectory is identical to the locked run.
# The FINAL success grade is read once with scheme='v2', which reads only nodes the kernel
# NEVER pins: CASP3 (apoptosis), RELA+STAT3 (resistance-off, replacing the pinned AKT1),
# CCND1+E2F1 (proliferation-off). AKT1 is still DRIVEN off by the drug (GLOBAL_INVARIANT
# untouched); v2 only stops CREDITING a clamped node, removing the self-scoring.
# ------------------------------------------------------------------
PHENO_SCHEMA_V2 = {
    'Resistance_OFF': {
        'required': {'RELA': 0, 'STAT3': 0},          # was {'AKT1': 0}; AKT1 pinned -> read hub
        'supporting_any': {'FOXO3': 1, 'MTOR': 0, 'RELA': 0, 'STAT3': 0,
                           'XIAP': 0, 'CFLAR': 0, 'MCL1': 0},
        'forbidden_any': {'RELA': 1, 'STAT3': 1},     # dropped AKT1: 1
        'contradiction_any': {},
    },
    'Apoptosis_ON': PHENO_SCHEMA['Apoptosis_ON'],     # CASP3 required is already clean -> unchanged
    'Proliferation_OFF': {
        'required': {'CCND1': 0, 'E2F1': 0},          # already clean
        'supporting_any': {'MYC': 0, 'RB1': 1, 'CDKN1A': 1, 'CDKN2A': 1,
                           'CDK4': 0, 'CDK6': 0, 'CDK2': 0},
        'forbidden_any': {'CCND1': 1, 'E2F1': 1},     # dropped MYC: 1
        'contradiction_any': {},
    },
    'Invasion_OFF': PHENO_SCHEMA['Invasion_OFF'],     # unchanged
    'Terminal': {
        'required': {'CASP3': 1, 'CCND1': 0, 'E2F1': 0, 'STAT3': 0, 'RELA': 0},  # dropped AKT1, CASP9, MYC
        'supporting_any': {'BAX': 1, 'CYCS': 1, 'FOXO3': 1, 'PTEN': 1},
        'forbidden_any': {'STAT3': 1, 'RELA': 1},     # dropped AKT1, MYC
        'contradiction_any': {},
    },
}

_SCHEMAS = {'v1': PHENO_SCHEMA, 'v2': PHENO_SCHEMA_V2}

# Nodes the v2 SCORING verdict reads (required + forbidden across stages). If the kernel
# ever PINS one of these, that node's final grade is self-scored for that patient. By
# decision E65 we do NOT forbid this (kernels stay as chosen); instead we WARN once per
# node so any future cohort/target/clock setting surfaces it for us to revisit (e.g. a
# majority-of-pathway success readout). Steering-only nodes (AKT1, CASP9, MYC -- v1 only)
# are intentionally excluded: pinning the drug target is fine, only scored nodes matter.
VERDICT_READOUT_NODES = set()
for _st in PHENO_SCHEMA_V2.values():
    VERDICT_READOUT_NODES |= set(_st.get('required', {}))
    VERDICT_READOUT_NODES |= set(_st.get('forbidden_any', {}))
_VERDICT_PIN_WARNED: set = set()


STAGE_SEQUENCE = ['Resistance_OFF', 'Apoptosis_ON', 'Proliferation_OFF', 'Terminal']


def evaluate_phenotype(bus: State, stage: str, scheme: str = 'v1') -> str:
    """Ternary PASS/WARN/FAIL. scheme='v1' STEERS (locked verdict), 'v2' SCORES (clean readout)."""
    sch = _SCHEMAS[scheme][stage]
    req_fail = any(bus.get(nd, 0) != v for nd, v in sch['required'].items())
    forb_hit = any(bus.get(nd, 0) == v for nd, v in sch['forbidden_any'].items())
    contra_hit = any(bus.get(nd, 0) == v
                     for nd, v in sch['contradiction_any'].items())
    support = sch['supporting_any']
    support_present = len(support) > 0
    support_ok = any(bus.get(nd, 0) == v for nd, v in support.items())
    if req_fail or forb_hit:
        return 'FAIL'
    if contra_hit:
        return 'WARN'
    if support_present and support_ok:
        return 'PASS'
    return 'WARN'


def failed_required(bus: State, stage: str) -> List[str]:
    return [nd for nd, v in PHENO_SCHEMA[stage]['required'].items()
            if bus.get(nd, 0) != v]


def hit_forbidden(bus: State, stage: str) -> List[str]:
    return [nd for nd, v in PHENO_SCHEMA[stage]['forbidden_any'].items()
            if bus.get(nd, 0) == v]


def hit_contradictions(bus: State, stage: str) -> List[str]:
    return [nd for nd, v in PHENO_SCHEMA[stage]['contradiction_any'].items()
            if bus.get(nd, 0) == v]


# ============================================================
# CDE pathway families (verbatim from recommend_next_control_action_v1.m)
# ============================================================
EXPECTED_PATHWAYS = {
    'Resistance_OFF':    ['RTK_EGFR', 'RTK_Insulin', 'PI3K', 'AKT_Signaling',
                          'AKT_Survival', 'mTOR', 'JAK_STAT', 'NFkB', 'TF', 'MAPK'],
    'Apoptosis_ON':      ['Apoptosis_Regulatory', 'Apoptosis_Intrinsic',
                          'Apoptosis_Extrinsic', 'AKT_Signaling', 'AKT_Survival',
                          'mTOR', 'NFkB', 'JAK_STAT'],
    'Proliferation_OFF': ['CellCycle', 'CDK_CellCycle', 'MAPK', 'RTK_EGFR',
                          'RTK_Insulin', 'PI3K', 'TF', 'Hippo', 'Wnt'],
    'Terminal':          ['RTK_EGFR', 'PI3K', 'AKT_Signaling', 'AKT_Survival',
                          'mTOR', 'JAK_STAT', 'NFkB', 'TF', 'Apoptosis_Regulatory',
                          'CellCycle', 'CDK_CellCycle'],
}
ENABLING_PATHWAYS = {
    'Apoptosis_ON':      ['Apoptosis_Regulatory', 'AKT_Signaling', 'AKT_Survival',
                          'mTOR', 'NFkB', 'JAK_STAT', 'Apoptosis_Intrinsic'],
    'Proliferation_OFF': ['CellCycle', 'CDK_CellCycle', 'MAPK', 'RTK_EGFR',
                          'RTK_Insulin', 'TF'],
    'Resistance_OFF':    ['RTK_EGFR', 'RTK_Insulin', 'PI3K', 'AKT_Signaling',
                          'AKT_Survival', 'mTOR', 'JAK_STAT', 'NFkB'],
    'Terminal':          ['RTK_EGFR', 'PI3K', 'AKT_Signaling', 'AKT_Survival',
                          'mTOR', 'JAK_STAT', 'NFkB'],
}
CONTRADICTION_BREAKERS = {
    'Apoptosis_ON':      ['Apoptosis_Regulatory', 'AKT_Signaling', 'AKT_Survival',
                          'mTOR', 'NFkB', 'JAK_STAT'],
    'Resistance_OFF':    ['PI3K', 'AKT_Signaling', 'AKT_Survival', 'mTOR',
                          'JAK_STAT', 'NFkB', 'RTK_EGFR', 'RTK_Insulin'],
    'Proliferation_OFF': ['CellCycle', 'CDK_CellCycle', 'MAPK', 'RTK_EGFR',
                          'RTK_Insulin', 'TF'],
    'Terminal':          ['Apoptosis_Regulatory', 'AKT_Signaling', 'AKT_Survival',
                          'mTOR', 'NFkB', 'JAK_STAT'],
}
STAGE_PROTECT = {
    'Resistance_OFF': [],
    'Apoptosis_ON': ['Resistance_OFF'],
    'Proliferation_OFF': ['Resistance_OFF', 'Apoptosis_ON'],
    'Terminal': ['Resistance_OFF', 'Apoptosis_ON', 'Proliferation_OFF'],
}


# ============================================================
# CDE: diagnosis -> recommended action (recommend_next_control_action_v1.m)
# ============================================================
def cde_recommend(bus: State, stage: str, stage_pass_window: bool) -> dict:
    """
    Faithful port of analyze_phenotypic_outcome_v1 + recommend_next_control_action_v1.
    Returns {mode, pathways, protect, failure_class}.
    """
    status = evaluate_phenotype(bus, stage)
    fr = failed_required(bus, stage)
    hf = hit_forbidden(bus, stage)
    hc = hit_contradictions(bus, stage)

    # failure classification (local_classify_failure)
    if status == 'PASS':
        failure_class = 'STAGE_ACHIEVED' if stage_pass_window else 'TRANSIENT_STAGE_HIT'
    elif status in ('FAIL', 'WARN'):
        if hf or hc:
            failure_class = 'CONTRADICTORY_STATE_PERSISTS'
        elif fr:
            failure_class = 'REQUIRED_MARKERS_NOT_REACHED'
        else:
            failure_class = 'FEASIBLE_NOT_REACHED'
    else:
        failure_class = 'UNRESOLVED'

    stage_paths = EXPECTED_PATHWAYS.get(stage, [])
    protect = STAGE_PROTECT[stage]

    # failure_class -> (mode, pathways)
    if failure_class == 'STAGE_ACHIEVED':
        mode, paths = 'advance_stage', []
    elif failure_class == 'TRANSIENT_STAGE_HIT':
        mode, paths = 'hold_stage', stage_paths
    elif failure_class == 'CONTRADICTORY_STATE_PERSISTS':
        mode, paths = 'coordinated_pathway_set', CONTRADICTION_BREAKERS.get(stage, stage_paths)
    elif failure_class == 'REQUIRED_MARKERS_NOT_REACHED':
        mode, paths = 'coordinated_pathway_set', ENABLING_PATHWAYS.get(stage, stage_paths)
    elif failure_class == 'FEASIBLE_NOT_REACHED':
        mode, paths = 'coordinated_pathway_set', stage_paths
    else:
        mode, paths = 'coordinated_pathway_set', stage_paths

    # XIAP forbidden-hit override (verbatim)
    if 'XIAP' in hf:
        paths = ['Apoptosis_Regulatory', 'AKT_Signaling', 'AKT_Survival', 'mTOR']

    return {'mode': mode, 'pathways': list(paths), 'protect': protect,
            'failure_class': failure_class, 'status': status}


# ============================================================
# INNER tier: simulate_bbcn118_sequential.m
# ============================================================
def _pathway_mismatch(bus: State, pname: str, stage: str) -> float:
    """
    Mismatch fraction over ALL pathway nodes vs the full-bus stage target
    (matching MATLAB: nnz(xor(cur,tgt))/nx).
    """
    nodes = PATHWAYS[pname]['nodes']
    tgt = stage_target_for_pathway(stage, pname)
    mm = sum(1 for i, nd in enumerate(nodes) if int(bus.get(nd, 0)) != tgt[i])
    return mm / max(1, len(nodes))


def _resolve_active(bus: State, selected: List[str], stage: str,
                    priority_weight: float = 2.0) -> Optional[str]:
    """
    Port of local_resolve_active_pathway (auto + useRec mode).
    score(p) = mismatch_fraction(p) + (priority_weight if p in selected else 0),
    restricted to candidates in `selected`. Pick argmax; ties broken by FIRST
    (lowest) global pathway index, exactly like MATLAB max() over the pathway
    array nm. Returns the chosen pathway name, or None if no positive score.

    PATHWAY_MODE switches the resolver:
      'cde'    (default) : the behavior described above, over `selected`.
      'static'           : pure max-mismatch benchmark — resolve over the full
                           stage-eligible pool the CDE narrows from
                           (EXPECTED_PATHWAYS[stage]), ignoring the CDE
                           set-narrowing in `selected` and dropping the priority
                           bonus. argmax of _pathway_mismatch, same first-index
                           tie-break. Reachable only when PATHWAY_MODE='static'.
    """
    if PATHWAY_MODE == 'static':
        pool = set(EXPECTED_PATHWAYS.get(stage, selected))
        best_name, best_score = None, -np.inf
        for pname in PATHWAYS:           # global pathway order = tie-break order
            if pname not in pool:
                continue
            mf = _pathway_mismatch(bus, pname, stage)
            if mf > best_score:          # strict > => first (lowest idx) wins
                best_score, best_name = mf, pname
        return best_name

    sel = set(selected)
    best_name, best_score = None, -np.inf
    for pname in PATHWAYS:               # global pathway order = nm order
        if pname not in sel:
            continue
        mf = _pathway_mismatch(bus, pname, stage)
        score = mf + priority_weight     # bonus uniform within recommended set
        if score > best_score:           # strict > => first max wins (lowest idx)
            best_score, best_name = score, pname
    # MATLAB local_resolve_active_pathway returns the argmax regardless of whether
    # its mismatch is 0 (it does NOT return empty at 0 mismatch). Match that.
    return best_name


def _update_pathway(bus: State, pname: str, is_active: bool,
                    kernel_log: list, stage: str, session: int, sim_step: int,
                    maintain_inv: Dict[str, int] = None
                    ) -> State:
    """
    Port of update_one_pathway.m. If active, search a kernel toward the full-bus
    stage target projected onto this pathway, pin it, evolve the rest free.

    When maintain_inv is non-empty and this pathway owns one of those nodes, the
    target for that node is overridden to the maintained-invariant value, so the
    kernel search drives it there (held by a real kernel, never clamped).
    """
    pdef = PATHWAYS[pname]
    nodes = pdef['nodes']
    n = len(nodes)
    x0 = [int(bus.get(nd, 0)) for nd in nodes]
    target = stage_target_for_pathway(stage, pname)

    if maintain_inv:
        for i, nd in enumerate(nodes):
            if nd in maintain_inv:
                target[i] = maintain_inv[nd]

    kernel_idx = None
    if is_active:
        if KERNEL_METHOD == 'stabilize':
            # Theorem-1 kernel is GLOBAL (x0-independent): depends only on the
            # pathway rules, its externals, and the target -> cache across patients.
            exts = _pathway_externals(pname)
            sk = (pname, tuple(int(bus.get(e, 0)) for e in exts), tuple(target))
            if sk in _STAB_CACHE:
                sel = _STAB_CACHE[sk]
            else:
                rule_closed = lambda x: tuple(int(v) for v in pdef['rules'](list(x), bus))
                sel = FS.stabilize_select_kernel(rule_closed, n, target,
                                                 candidates=list(range(1, n + 1)), max_k=3)
                _STAB_CACHE[sk] = sel
        else:  # ranked (preserved)
            exts = _pathway_externals(pname)
            ckey = (pname, tuple(int(bus.get(e, 0)) for e in exts))
            L = generate_L(pdef['rules'], n, bus, cache_key=ckey)
            sel = forward_select_kernel(L, target, x0, max_k=2,
                                        rule_fn=pdef['rules'], ext=bus)
            if sel is None:
                sel = forward_select_kernel(L, target, x0, max_k=3,
                                            rule_fn=pdef['rules'], ext=bus)
        if sel is not None:
            kernel_idx = [s - 1 for s in sel]   # to 0-based

    # read-only capture of the chosen kernel (no-op unless KERNEL_CAPTURE is a list)
    if KERNEL_CAPTURE is not None and kernel_idx:
        KERNEL_CAPTURE.append({
            'pathway': pname, 'stage': stage,
            'kernel_nodes': [nodes[i] for i in kernel_idx],
            'kernel_size': len(kernel_idx),
        })

    # apply: pin kernel, evolve rest under rules
    x_pinned = list(x0)
    if kernel_idx:
        for i in kernel_idx:
            x_pinned[i] = target[i]
    nxt = [int(v) for v in pdef['rules'](x_pinned, bus)]
    if kernel_idx:
        for i in kernel_idx:
            nxt[i] = target[i]

    new_bus = dict(bus)
    for i, nd in enumerate(nodes):
        new_bus[nd] = nxt[i]

    if kernel_idx:
        knames = [nodes[i] for i in kernel_idx]
        kvals = [target[i] for i in kernel_idx]
        for _nd in knames:
            if _nd in VERDICT_READOUT_NODES and _nd not in _VERDICT_PIN_WARNED:
                _VERDICT_PIN_WARNED.add(_nd)
                import warnings as _warnings
                _warnings.warn(
                    f"[BBCN verdict-integrity] kernel pinned scored node '{_nd}' "
                    f"(pathway {pname}); its v2 grade is self-scored for some patients. "
                    f"Kept by decision E65; revisit later (e.g. majority-of-pathway readout).",
                    stacklevel=2)
        kernel_log.append({
            'stage': stage, 'session': session, 'sim_step': sim_step,
            'pathway': pname, 'kernel_node_names': knames,
            'kernel_values': kvals, 'kernel_size': len(kernel_idx)})
    return new_bus


GLOBAL_INVARIANT = {'AKT1': 0, 'CASP3': 1, 'FOXO3': 1, 'MYC': 0}
_INVARIANT_OWNERS = {'AKT_Signaling', 'Apoptosis_Extrinsic', 'TF'}


def maintained_invariant(stage: str) -> Dict[str, int]:
    """
    Stage-wise maintained invariant M(s), corrected formulation.

        M(s) = GLOBAL_INVARIANT \\ { nodes the current stage's target sets }

    Rationale (the Part-2 finding refined the earlier purely-protective rule):
    the invariant plays TWO roles, not one. (1) PROTECTIVE: holding nodes
    required by already-achieved phenotypes. (2) ENABLING: the structural
    decoupling from holding the invariant (Part 1: freeing JAK_STAT/NFkB,
    suppressing MYC) helps the CURRENT stage even before anything is "achieved".
    The earlier protective-only rule discarded the enabling role and hurt
    Resistance. So we maintain the whole invariant EXCEPT the node(s) the
    current stage is itself trying to set (condition b: release contested nodes
    to the active stage, so maintenance never fights the stage for its own
    target). Nodes are held by re-applied kernels, never clamped.
    """
    stage_nodes = set(PHENO_SCHEMA[stage]['required'].keys())
    return {nd: v for nd, v in GLOBAL_INVARIANT.items() if nd not in stage_nodes}


def _invariant_violated(bus: State, inv: Dict[str, int]) -> bool:
    return any(int(bus.get(nd, 0)) != v for nd, v in inv.items())


def simulate_inner(bus: State, selected_pathways: List[str], stage: str,
                   session: int, kernel_log: list, max_steps: int = 10,
                   maintain_invariant: bool = False,
                   stagewise_invariant: bool = False,
                   coactive: bool = False) -> State:
    """
    Inner simulator. With one pathway per batch that pathway is active each step.
    Synchronous prev->next update; non-active pathways evolve free.

    Invariant maintenance (no clamps, kernels only):
      * maintain_invariant + not stagewise -> hold the full GLOBAL_INVARIANT
      * maintain_invariant + stagewise     -> hold only M(stage) (the formal
        stage-wise maintained set: protected, not currently contested)
    A pathway owning a violated maintained node is made active and searches a
    kernel toward the invariant value. If the search fails, the node is simply
    not held that step (informative, not forced).
    """
    if maintain_invariant:
        inv = maintained_invariant(stage) if stagewise_invariant else GLOBAL_INVARIANT
    else:
        inv = {}

    for t in range(1, max_steps + 1):
        if coactive:
            # MATLAB UseSelectedPathways mode: the ENTIRE selected set is active
            # every step (simulate_bbcn118_sequential.m:230, activeSet = selPaths),
            # each pathway pinning its own kernel toward target simultaneously.
            active_set = set(selected_pathways)
        else:
            active = _resolve_active(bus, selected_pathways, stage)
            active_set = {active} if active else set()

        if inv and _invariant_violated(bus, inv):
            for nd, want in inv.items():
                if int(bus.get(nd, 0)) != want:
                    active_set.add(NODE_TO_PATHWAY[nd])

        new_bus = dict(bus)
        for pname in PATHWAYS:
            is_active = (pname in active_set)
            upd = _update_pathway(bus, pname, is_active,
                                  kernel_log, stage, session, t,
                                  maintain_inv=inv)
            for nd in PATHWAYS[pname]['nodes']:
                new_bus[nd] = upd[nd]
        if new_bus == bus:
            break
        bus = new_bus
    return bus


# ============================================================
# OUTER tier: single_patient_harness.m  (+ stage sequencing)
# ============================================================
def _schedule_batches(paths: List[str], max_per_exposure: int) -> List[List[str]]:
    """local_schedule_active_pathways: split into batches of <=max_per_exposure."""
    seen, uniq = set(), []
    for p in paths:
        if p and p != 'auto' and p not in seen:
            seen.add(p); uniq.append(p)
    if not uniq:
        return []
    m = max(1, int(max_per_exposure))
    return [uniq[i:i + m] for i in range(0, len(uniq), m)]


# Per-stage initial pathway sets, verbatim from run_cohort_harness_from_tcga_cache.m.
# Stage 1 & 3 use UseRecommendedPathways=true (seed lists below); Stage 2 & 4 use
# explicit SelectedPathways. These are the sets fed to the scheduler.
STAGE_SELECTED = {
    'Resistance_OFF':    ['RTK_EGFR', 'PI3K', 'AKT_Signaling', 'AKT_Survival',
                          'mTOR', 'JAK_STAT', 'NFkB'],
    'Apoptosis_ON':      ['AKT_Signaling', 'AKT_Survival', 'Apoptosis_Regulatory',
                          'MAPK', 'JAK_STAT', 'NFkB', 'mTOR', 'PI3K'],
    'Proliferation_OFF': ['CellCycle', 'CDK_CellCycle', 'MAPK', 'RTK_EGFR',
                          'RTK_Insulin', 'PI3K', 'TF'],
    'Terminal':          ['Apoptosis_Regulatory', 'Apoptosis_Intrinsic',
                          'Apoptosis_Extrinsic', 'AKT_Signaling', 'AKT_Survival',
                          'mTOR', 'JAK_STAT', 'NFkB', 'PI3K'],
}


def run_patient(init_bus: State, T: int = 5, max_steps: int = 10,
                max_pathways_per_exposure: int = 1, confirm_window: int = 3,
                maintain_invariant: bool = False,
                stagewise_invariant: bool = False,
                verbose: bool = False) -> dict:
    """
    Full two-tier supervisory run over the 4-stage sequence for one patient.

    Faithful to single_patient_harness.m + run_cohort_harness_from_tcga_cache.m:
      * each stage runs up to T outer sessions
      * the stage's SelectedPathways set is scheduled into batches of
        MaxActivePathwaysPerExposure (=1: one pathway per exposure)
      * each batch runs the inner free-dynamics simulator (MaxSteps)
      * after the session the CDE reassesses; a stage is PROMOTED when its
        phenotype holds PASS across the confirmation window (=3)
      * the bus is carried forward to the next stage (cumulative protection)
    Returns trace with per-stage outcomes, kernel history, and stage summary.
    """
    bus = {nd: int(init_bus.get(nd, 0)) for nd in ALL_NODES}
    kernel_log = []
    stage_summary = {}

    for stage in STAGE_SEQUENCE:
        init_status = evaluate_phenotype(bus, stage)
        selected = list(STAGE_SELECTED[stage])
        pass_history = []
        promoted = False
        promotion_session = None

        for k in range(1, T + 1):
            # schedule the selected set into 1-pathway exposure batches;
            # each batch runs the inner simulator with that single pathway active.
            batches = _schedule_batches(selected, max_pathways_per_exposure)
            for batch in batches:
                bus = simulate_inner(bus, batch, stage, k, kernel_log, max_steps,
                                     maintain_invariant=maintain_invariant,
                                     stagewise_invariant=stagewise_invariant)

            status = evaluate_phenotype(bus, stage)
            pass_history.append(status)

            # stage confirmation: PASS held across the confirmation window -> promote
            if (len(pass_history) >= confirm_window and
                    all(s == 'PASS' for s in pass_history[-confirm_window:])):
                promoted = True
                promotion_session = k
                break

        final_status = evaluate_phenotype(bus, stage)
        # MATLAB final_pass_rate counts a stage as passed if its phenotype is PASS
        # at stage end OR it was confirmed/promoted.
        stage_passed = (final_status == 'PASS') or promoted
        stage_summary[stage] = {
            'init_status': init_status,
            'final_status': final_status,
            'passed': stage_passed,
            'promoted': promoted,
            'promotion_session': promotion_session,
        }

    return {
        'final_bus': bus,
        'stage_summary': stage_summary,
        'kernel_history': kernel_log,
        'terminal_pass': evaluate_phenotype(bus, 'Terminal') == 'PASS',
    }


# ---- complete external dependency set (every e_(e,'X') a pathway's rules reference) ----
import ast as _ast, os as _os
def _build_complete_ext():
    src = open(_os.path.join(_os.path.dirname(__file__), 'pathways.py')).read()
    tree = _ast.parse(src); out = {}
    for fn in tree.body:
        if isinstance(fn, _ast.FunctionDef) and fn.name.startswith('_'):
            nodes = None; refs = set()
            for nd in _ast.walk(fn):
                if isinstance(nd, _ast.Call) and isinstance(nd.func, _ast.Name) and nd.func.id == 'e_':
                    try: refs.add(_ast.literal_eval(nd.args[1]))
                    except Exception: pass
            rets = [s for s in fn.body if isinstance(s, _ast.Return)]
            if rets and isinstance(rets[-1].value, _ast.Dict):
                for k, v in zip(rets[-1].value.keys, rets[-1].value.values):
                    if _ast.literal_eval(k) == 'nodes': nodes = tuple(_ast.literal_eval(v))
            if nodes: out[nodes] = sorted(refs)
    return out
_COMPLETE_EXT = _build_complete_ext()
def complete_externals(pname):
    """Complete external dependency set for a pathway: every bus node its rules read.
    A superset of the functional probe in _pathway_externals; use for EXACT cache keys."""
    own = set(PATHWAYS[pname]['nodes'])
    s = set(_COMPLETE_EXT.get(tuple(PATHWAYS[pname]['nodes']), [])) | set(_pathway_externals(pname))
    return [e for e in ALL_NODES if e in s and e not in own]
