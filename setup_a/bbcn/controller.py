"""
bbcn.controller — the canonical capped, no-accumulation sequential controller.

ONE implementation, with a baseline default and documented options.
Baseline (out of the box): isolated per-phenotype control, ideal targets,
2 pathways / cap 4 / upstream-first / horizon 8 cycles, init from patient, NO accumulation.

Options (documented kwargs):
  family   : 'ideal' (isolated baseline) | 'compatible' (sequenced/state-preserving)
  n_pw     : pathways per induction cycle (1 or 2)
  ind_cap  : induction target cap (3,4,5)  -> ~4-5 == two pathways, ~2-3 == one pathway
  H_cyc    : horizon in cycles (8 ~6mo primary, 17 ~1yr extended)
  rank     : 'upstream' (decisive for apoptosis) | 'mismatch'
  mode     : 'isolated' (each phenotype from patient init) | 'sequenced' (carry bus, preserve prior)

NO FORCING: success requires a genuine free-dynamics fixed point.
NO ACCUMULATION: the held kernel is re-decided each cycle and never exceeds the cap;
the bus state carries forward but prior-cycle pins are NOT retained.
"""
from __future__ import annotations
from bbcn import harness as H
from bbcn import stp
from bbcn import forward_stab as FS

PHENOTYPES = ['Resistance_OFF', 'Apoptosis_ON', 'Proliferation_OFF']
_ORDER = {p: i for i, p in enumerate(H.PATHWAYS)}
_NODES = set(H.ALL_NODES)
_KC = {}
# Persistent cache for the (x0-independent) stabilize kernels: NOT cleared per
# patient, since the Theorem-1 kernel depends only on (family, pathway, externals,
# phenotype), never on the patient's starting state.
_STAB_KC = {}


def _step(bus):
    nb = dict(bus)
    for p, pdef in H.PATHWAYS.items():
        x = [int(bus.get(n, 0)) for n in pdef['nodes']]
        nx = [int(v) for v in pdef['rules'](x, bus)]
        for i, n in enumerate(pdef['nodes']):
            nb[n] = nx[i]
    return nb


def evolve(bus, held, max_steps=25):
    """Tier-1b relaxation: hold the kernel, evolve free dynamics to a fixed point."""
    bus = dict(bus)
    for _ in range(max_steps):
        nb = _step(bus)
        for k, v in held.items():
            nb[k] = v
        if nb == bus:
            break
        bus = nb
    return bus


def is_fixed_point(bus, held):
    """Genuine free-dynamics fixed point under the held kernel (no forcing beyond held)."""
    nb = _step(bus)
    for k, v in held.items():
        nb[k] = v
    return nb == bus


def _kernel_for(bus, pw, ph, method='stabilize'):
    """Minimal kernel for one pathway toward its phenotype target.

    method='ranked'    : preserved heuristic (impact->delta->steps->causal), stp.py.
    method='stabilize' : algebraic Theorem-1 global-stabilization kernel (forward_stab.py).
    Same candidate pool (all pathway nodes) and same size cap (<=3) for an apples-to-apples
    comparison; only the acceptance/selection criterion differs.
    """
    pdef = H.PATHWAYS[pw]; nd = pdef['nodes']; n = len(nd)
    x0 = [int(bus.get(x, 0)) for x in nd]
    tgt = H.stage_target_for_pathway(ph, pw)
    exts = H.complete_externals(pw)   # EXACT: complete external dependency set
    ck = (H.TARGET_FAMILY, pw, tuple(int(bus.get(e, 0)) for e in exts), ph, method)
    if method == 'stabilize':
        if ck in _STAB_KC:
            sel = _STAB_KC[ck]
        else:
            rule_closed = lambda x: tuple(int(v) for v in pdef['rules'](list(x), bus))
            sel = FS.stabilize_select_kernel(rule_closed, n, tgt,
                                             candidates=list(range(1, n + 1)), max_k=3)
            _STAB_KC[ck] = sel
    elif ck in _KC:
        sel = _KC[ck]
    else:  # ranked (preserved)
        L = stp.generate_L(pdef['rules'], n, bus, cache_key=(pw, ck[2]))
        sel = (stp.forward_select_kernel(L, tgt, x0, max_k=2, rule_fn=pdef['rules'], ext=bus)
               or stp.forward_select_kernel(L, tgt, x0, max_k=3, rule_fn=pdef['rules'], ext=bus))
        _KC[ck] = sel
    return {nd[s - 1]: tgt[s - 1] for s in sel} if sel else {}


def control_phenotype(bus, ph, protect=(), n_pw=2, ind_cap=4, H_cyc=8, rank='upstream',
                      kernel_method='stabilize'):
    """Drive one phenotype under the cap, no accumulation. Returns (achieved, bus, held)."""
    held = {}; achieved = False
    def ok(b):
        return (H.evaluate_phenotype(b, ph) == 'PASS' and is_fixed_point(b, held)
                and all(H.evaluate_phenotype(b, pr) == 'PASS' for pr in protect))
    for _ in range(H_cyc):
        if ok(bus):
            achieved = True; break
        avail = [c for c in H.STAGE_SELECTED[ph] if H._pathway_mismatch(bus, c, ph) > 0]
        if not avail:
            break
        if rank == 'upstream':
            avail.sort(key=lambda p: _ORDER[p])
        else:
            avail.sort(key=lambda p: (-H._pathway_mismatch(bus, p, ph), _ORDER[p]))
        new_held = {}
        for pw in avail[:n_pw]:
            for nn, vv in _kernel_for(bus, pw, ph, method=kernel_method).items():
                if len(new_held) >= ind_cap:
                    break
                new_held[nn] = vv
            if len(new_held) >= ind_cap:
                break
        if not new_held:
            break
        held = new_held                 # REPLACE — no accumulation
        bus = evolve(bus, held)
    if not achieved:
        achieved = ok(bus)
    return achieved, bus, held


def run_patient(init_bus, mode='isolated', family='ideal',
                n_pw=2, ind_cap=4, H_cyc=8, rank='upstream', kernel_method='stabilize'):
    """Run the controller for one patient. Returns dict of per-phenotype achieved + all-three."""
    _KC.clear()
    if H.TARGET_FAMILY != family or not H._STAGE_TARGETS:
        H._STAGE_TARGETS.clear(); H._load_stage_targets(family)
    bus0 = {n: 0 for n in H.ALL_NODES}
    for n in _NODES:
        bus0[n] = int(init_bus.get(n, 0))
    bus0 = evolve(bus0, {})
    out = {}
    if mode == 'sequenced':
        bus = dict(bus0); ok_all = True
        for i, ph in enumerate(PHENOTYPES):
            ok, bus, _ = control_phenotype(bus, ph, protect=PHENOTYPES[:i],
                                           n_pw=n_pw, ind_cap=ind_cap, H_cyc=H_cyc, rank=rank,
                                           kernel_method=kernel_method)
            out[ph] = ok; ok_all = ok_all and ok
        out['all_three'] = ok_all
    else:  # isolated
        for ph in PHENOTYPES:
            ok, _, _ = control_phenotype(dict(bus0), ph,
                                         n_pw=n_pw, ind_cap=ind_cap, H_cyc=H_cyc, rank=rank,
                                         kernel_method=kernel_method)
            out[ph] = ok
    return out
