"""
bbcn.structural_decoupling
==========================
Part 1 (model-based, NO dynamics, NO forcing): quantify how much the pathway
coupling graph decouples when the global invariant is held constant.

This is a purely structural calculation on the Boolean model:
  * Build the pathway coupling graph G_P: edge p -> q iff pathway q's update
    rules actually read a bus node owned by pathway p (detected robustly by
    probing the rule over many internal states, not just all-zero).
  * Compute the strongly connected components (SCCs) and the dominant SCC.
  * Then FIX the invariant nodes {AKT1=0, CASP3=1, FOXO3=1, MYC=0} as constants
    in the coupling detection (a node that no longer varies cannot transmit a
    perturbation), recompute G_P, and report how the SCC shrinks.

Fixing a node here is an ANALYSIS device — we are asking "if this node were held
constant, would it still carry coupling?" It is NOT a control action and nothing
is clamped during any simulation. It answers: what is the structural ceiling on
decoupling that the gamma_p = 0 invariant could deliver.
"""

from typing import List, Dict, Tuple, Set
import itertools
import numpy as np

from bbcn.pathways import PATHWAYS, ALL_NODES

# The global invariant computed by the gamma_p=0 method (from the graph).
GLOBAL_INVARIANT = {'AKT1': 0, 'CASP3': 1, 'FOXO3': 1, 'MYC': 0}


def _active_pathways() -> List[str]:
    return [p for p in PATHWAYS if PATHWAYS[p].get('layer') != 'Boundary']


def _reads_external(qname: str, ext_node: str, n_samples: int = 24) -> bool:
    """
    Does pathway q's rule output depend on the external bus node `ext_node`?
    Robust detection: sample many random internal states of q and, for each,
    flip ext_node 0->1 with all OTHER externals held at a random fixed value;
    if any output bit changes, q reads ext_node.
    """
    pdef = PATHWAYS[qname]
    nq = len(pdef['nodes'])
    rng = np.random.default_rng(0)
    other_ext = [nd for nd in ALL_NODES if nd not in pdef['nodes']
                 and nd != ext_node]
    for _ in range(n_samples):
        x = [int(b) for b in rng.integers(0, 2, nq)]
        base = {nd: int(b) for nd, b in
                zip(other_ext, rng.integers(0, 2, len(other_ext)))}
        e0 = dict(base); e0[ext_node] = 0
        e1 = dict(base); e1[ext_node] = 1
        if pdef['rules'](x, e0) != pdef['rules'](x, e1):
            return True
    return False


def build_coupling_graph(fixed: Dict[str, int] = None
                         ) -> Tuple[Dict[str, Set[str]], List[str]]:
    """
    Build directed pathway coupling graph. edge p -> q means q depends on a
    bus node owned by p. If `fixed` is given, those nodes are treated as
    constants: a fixed node cannot transmit coupling, so edges carried solely
    by fixed nodes disappear.

    Returns (adjacency, pathway_list) where adjacency[p] = set of q that depend
    on p.
    """
    fixed = fixed or {}
    paths = _active_pathways()
    # node -> owning pathway
    owner = {}
    for p in paths:
        for nd in PATHWAYS[p]['nodes']:
            owner[nd] = p
    adj: Dict[str, Set[str]] = {p: set() for p in paths}
    for p in paths:
        for exported in PATHWAYS[p]['bus_exported']:
            if exported in fixed:
                continue                      # held constant -> carries no coupling
            for q in paths:
                if q == p:
                    continue
                if exported in PATHWAYS[q]['nodes']:
                    continue                  # q owns it, not external to q
                if _reads_external(q, exported):
                    adj[p].add(q)
    return adj, paths


def tarjan_scc(adj: Dict[str, Set[str]], nodes: List[str]) -> List[List[str]]:
    """Standard Tarjan SCC. Returns list of components (each a list of nodes)."""
    index = {}
    low = {}
    onstack = {}
    stack = []
    counter = [0]
    out = []

    import sys
    sys.setrecursionlimit(10000)

    def strong(v):
        index[v] = counter[0]; low[v] = counter[0]; counter[0] += 1
        stack.append(v); onstack[v] = True
        for w in adj[v]:
            if w not in index:
                strong(w); low[v] = min(low[v], low[w])
            elif onstack.get(w):
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = stack.pop(); onstack[w] = False; comp.append(w)
                if w == v:
                    break
            out.append(comp)

    for v in nodes:
        if v not in index:
            strong(v)
    return out


def scc_report(fixed: Dict[str, int] = None) -> dict:
    """Build the graph (optionally with fixed nodes) and summarise its SCCs."""
    adj, paths = build_coupling_graph(fixed=fixed)
    comps = tarjan_scc(adj, paths)
    multi = [c for c in comps if len(c) > 1]
    multi.sort(key=len, reverse=True)
    dominant = multi[0] if multi else []
    n_edges = sum(len(v) for v in adj.values())
    cyclic_pathways = sum(len(c) for c in multi)
    return {
        'n_pathways': len(paths),
        'n_edges': n_edges,
        'n_sccs': len(comps),
        'n_multi_sccs': len(multi),
        'dominant_scc_size': len(dominant),
        'dominant_scc': sorted(dominant),
        'cyclic_fraction': cyclic_pathways / max(1, len(paths)),
        'all_multi_sccs': [sorted(c) for c in multi],
        'adjacency': adj,
    }


def compare_decoupling() -> dict:
    """Baseline SCC vs SCC with the global invariant held constant."""
    base = scc_report(fixed=None)
    inv = scc_report(fixed=GLOBAL_INVARIANT)
    return {'baseline': base, 'invariant_held': inv}


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    r = compare_decoupling()
    b, i = r['baseline'], r['invariant_held']
    print("=" * 60)
    print("STRUCTURAL DECOUPLING ANALYSIS (model-based, no dynamics)")
    print("=" * 60)
    print(f"\nGlobal invariant: {GLOBAL_INVARIANT}\n")
    print(f"{'metric':<28}{'baseline':>14}{'invariant held':>16}")
    print("-" * 58)
    print(f"{'pathways':<28}{b['n_pathways']:>14}{i['n_pathways']:>16}")
    print(f"{'coupling edges':<28}{b['n_edges']:>14}{i['n_edges']:>16}")
    print(f"{'total SCCs':<28}{b['n_sccs']:>14}{i['n_sccs']:>16}")
    print(f"{'multi-pathway SCCs':<28}{b['n_multi_sccs']:>14}{i['n_multi_sccs']:>16}")
    print(f"{'dominant SCC size':<28}{b['dominant_scc_size']:>14}{i['dominant_scc_size']:>16}")
    print(f"{'cyclic fraction':<28}{b['cyclic_fraction']:>14.3f}{i['cyclic_fraction']:>16.3f}")
    print(f"\nBaseline dominant SCC ({b['dominant_scc_size']} pathways):")
    print("  " + ", ".join(b['dominant_scc']))
    print(f"\nInvariant-held dominant SCC ({i['dominant_scc_size']} pathways):")
    print("  " + (", ".join(i['dominant_scc']) if i['dominant_scc'] else "(none / acyclic)"))
    freed = set(b['dominant_scc']) - set(i['dominant_scc'])
    print(f"\nPathways freed from the dominant SCC ({len(freed)}):")
    print("  " + (", ".join(sorted(freed)) if freed else "(none)"))
