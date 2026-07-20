"""
bbcn.pulsed — the pulsed weekly controller (third control arm).

One cycle = one week = cycle_days Boolean updates (one update = one day).
Each cycle acts on ONE pathway. The pathway's minimal kernel (up to ind_cap
nodes, i.e. a small multi-target cocktail on that one pathway) is applied as a
drug-on pulse for on_days, then removed; the network then evolves FREE for the
remaining relaxation days and settles. Success for a phenotype is that it holds
as a genuine drug-off fixed point after relaxation, not while the drug is held.

This is the realistic, withdrawable schedule. It reuses the exact BBCN dynamics
(controller._step / evolve / is_fixed_point), the exact kernel design
(forward_stab for 'stabilize', stp for 'ranked'), and the exact phenotype
evaluation (harness). The only new logic is the weekly pulse-then-relax loop and
an optional druggable tiebreak over equally minimal stabilizing kernels.
"""
from __future__ import annotations
import os, ast as _ast
from typing import Dict, List, Tuple
from bbcn import harness as H
from bbcn import controller as CTRL
from bbcn import forward_stab as FS
from bbcn import stp
from bbcn import drugs


def _complete_externals() -> Dict[tuple, list]:
    """The COMPLETE set of bus nodes each pathway's rules reference (every e_(e,'X')),
    extracted from source. This is a complete superset of the functional probe used by
    harness._pathway_externals, so keying the kernel cache on it cannot miss a dependency."""
    src = open(os.path.join(os.path.dirname(__file__), 'pathways.py')).read()
    tree = _ast.parse(src); out = {}
    for fn in tree.body:
        if isinstance(fn, _ast.FunctionDef) and fn.name.startswith('_'):
            nodes = None; refs = set()
            for nd in _ast.walk(fn):
                if isinstance(nd, _ast.Call) and isinstance(nd.func, _ast.Name) and nd.func.id == 'e_':
                    try:
                        refs.add(_ast.literal_eval(nd.args[1]))
                    except Exception:
                        pass
            rets = [s for s in fn.body if isinstance(s, _ast.Return)]
            if rets and isinstance(rets[-1].value, _ast.Dict):
                for k, v in zip(rets[-1].value.keys, rets[-1].value.values):
                    if _ast.literal_eval(k) == 'nodes':
                        nodes = tuple(_ast.literal_eval(v))
            if nodes:
                out[nodes] = sorted(refs)
    return out


_COMPLETE_EXT = _complete_externals()


def _ext_for(pw):
    return H.complete_externals(pw)

# druggability weight per node: 2 if a FDA-approved agent targets it, 1 if in trials.
DRUG_WEIGHT: Dict[str, int] = {}
for _d in drugs.DRUG_TABLE:
    _w = 2 if 'FDA' in _d.get('approval', '') else 1
    for _tok in _d['kernel'].split(','):
        _nm = _tok.strip().split('=')[0].strip()
        if _nm:
            DRUG_WEIGHT[_nm] = max(DRUG_WEIGHT.get(_nm, 0), _w)

PHENOTYPES = CTRL.PHENOTYPES
_ORDER = CTRL._ORDER
_DKC: Dict[tuple, tuple] = {}   # druggable-kernel cache (x0-independent for stabilize)
_RKC: Dict[tuple, tuple] = {}   # ranked-kernel cache (x0-dependent: keyed by externals AND x0)


def _druggable_kernel(bus, pw, ph, method='stabilize', ind_cap=3, druggable=True):
    """Minimal kernel for one pathway toward its target; among equally minimal
    stabilizing kernels, prefer the most druggable set. Falls back to the
    controller's deterministic selector for the ranked method."""
    pdef = H.PATHWAYS[pw]; nd = pdef['nodes']; n = len(nd)
    tgt = H.stage_target_for_pathway(ph, pw)
    exts = _ext_for(pw)            # COMPLETE external set (no missed dependency)
    ck = (H.TARGET_FAMILY, pw, tuple(int(bus.get(e, 0)) for e in exts), ph, method, druggable, ind_cap)
    if method == 'stabilize':
        if ck in _DKC:
            sel = _DKC[ck]
        else:
            rule_closed = lambda x: tuple(int(v) for v in pdef['rules'](list(x), bus))
            Lcol = FS.build_L(rule_closed, n)
            hits, _ = FS.stabilizing_kernels(Lcol, list(tgt), n,
                                             candidates=list(range(1, n + 1)),
                                             max_k=min(ind_cap, 3))
            if not hits:
                sel = None
            elif druggable:
                def score(S):
                    drug = sum(DRUG_WEIGHT.get(nd[s - 1], 0) for s in S)
                    return (-drug, sum(S), S)        # most druggable, then the usual tiebreak
                sel = sorted(hits, key=score)[0]
            else:
                sel = sorted(hits, key=lambda S: (sum(S), S))[0]
            _DKC[ck] = sel
        return {nd[s - 1]: tgt[s - 1] for s in sel} if sel else {}
    # ranked: x0-dependent. Recompute each cycle, keyed by (externals, x0); never reuse
    # a kernel designed for a different starting state.
    x0 = tuple(int(bus.get(x, 0)) for x in nd)
    rk = (H.TARGET_FAMILY, pw, tuple(int(bus.get(e, 0)) for e in exts), ph, x0)
    if rk in _RKC:
        sel = _RKC[rk]
    else:
        L = stp.generate_L(pdef['rules'], n, bus, cache_key=(pw, rk[2]))   # L = f(externals) only
        sel = (stp.forward_select_kernel(L, list(tgt), list(x0), max_k=2, rule_fn=pdef['rules'], ext=bus)
               or stp.forward_select_kernel(L, list(tgt), list(x0), max_k=3, rule_fn=pdef['rules'], ext=bus))
        _RKC[rk] = sel
    return {nd[s - 1]: tgt[s - 1] for s in sel} if sel else {}


def _ok_free(bus, ph, protect) -> bool:
    """Phenotype holds as a genuine drug-off (free-dynamics) fixed point."""
    return (H.evaluate_phenotype(bus, ph) == 'PASS'
            and CTRL.is_fixed_point(bus, {})
            and all(H.evaluate_phenotype(bus, pr) == 'PASS' for pr in protect))


def control_phenotype_pulsed(bus, ph, protect=(), on_days=1, cycle_days=7, H_cyc=5,
                             ind_cap=3, rank='upstream', kernel_method='stabilize',
                             druggable=True) -> Tuple[bool, dict, list]:
    """Drive one phenotype by weekly single-pathway pulses; success is drug-off durable."""
    bus = CTRL.evolve(dict(bus), {})                 # start settled (drug off)
    used: List[str] = []
    if _ok_free(bus, ph, protect):
        return True, bus, used
    for _cyc in range(H_cyc):
        avail = [c for c in H.STAGE_SELECTED[ph] if H._pathway_mismatch(bus, c, ph) > 0]
        if not avail:
            break
        if rank == 'upstream':
            avail.sort(key=lambda p: _ORDER[p])
        else:
            avail.sort(key=lambda p: (-H._pathway_mismatch(bus, p, ph), _ORDER[p]))
        kernel, pw = {}, None
        for cand in avail:                            # ONE pathway acts this cycle
            kernel = _druggable_kernel(bus, cand, ph, kernel_method, ind_cap, druggable)
            if kernel:
                pw = cand; break
        if not kernel:
            break
        b = dict(bus)
        for _ in range(on_days):                      # exposure: drug-on pulse
            b = CTRL._step(b)
            for k, v in kernel.items():
                b[k] = v
        for _ in range(max(0, cycle_days - on_days)): # relaxation: drug off
            b = CTRL._step(b)
        bus = CTRL.evolve(b, {})                      # settle buffer (free)
        used.append(pw)
        if _ok_free(bus, ph, protect):
            return True, bus, used
    return _ok_free(bus, ph, protect), bus, used


def run_patient_pulsed(init_bus, mode='isolated', family='ideal', on_days=1, cycle_days=7,
                       H_cyc=5, ind_cap=3, rank='upstream', kernel_method='stabilize',
                       druggable=True) -> dict:
    """Run the pulsed weekly controller for one patient."""
    if H.TARGET_FAMILY != family or not H._STAGE_TARGETS:
        H._STAGE_TARGETS.clear(); H._load_stage_targets(family)
    bus0 = {n: 0 for n in H.ALL_NODES}
    for n in CTRL._NODES:
        bus0[n] = int(init_bus.get(n, 0))
    bus0 = CTRL.evolve(bus0, {})
    out = {}
    if mode == 'sequenced':
        bus = dict(bus0); ok_all = True
        for i, ph in enumerate(PHENOTYPES):
            ok, bus, _ = control_phenotype_pulsed(bus, ph, protect=PHENOTYPES[:i],
                                                  on_days=on_days, cycle_days=cycle_days,
                                                  H_cyc=H_cyc, ind_cap=ind_cap, rank=rank,
                                                  kernel_method=kernel_method, druggable=druggable)
            out[ph] = ok; ok_all = ok_all and ok
        out['all_three'] = ok_all
    else:
        for ph in PHENOTYPES:
            ok, _, _ = control_phenotype_pulsed(dict(bus0), ph, on_days=on_days,
                                                cycle_days=cycle_days, H_cyc=H_cyc, ind_cap=ind_cap,
                                                rank=rank, kernel_method=kernel_method,
                                                druggable=druggable)
            out[ph] = ok
    return out
