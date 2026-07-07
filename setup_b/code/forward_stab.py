"""
forward_stab.py — shared kernel-selection core for BBCN Setup A and Setup B.

Two selectors, chosen by `method`:

  'ranked'    : the PRESERVED heuristic (impact -> delta -> steps -> causal).
                Each setup keeps its own existing ranked path; this module does
                NOT reimplement it, so the previous LOCKED_NUMBERS stay exact.

  'stabilize' : the NEW canonical method = algebraic forward stabilization,
                a verbatim port of the Theorem-1 global-stabilization test of
                  M. R. Rafimanzelat, "Global stabilization of Boolean networks
                  with applications to biomolecular network control,"
                  Scientific Reports 15:15201 (2025).  DOI 10.1038/s41598-025-97684-y
                Applied per pathway (Setup A) or per timescale box (Setup B), on
                the local 2^|nodes| state space (|nodes| <= 9 -> <= 512 states).

The algebraic functions below are copied verbatim from the implementation that
was validated against the paper's worked Examples 1 and 2 (see ledger). The NVR
index convention is identical to bbcn/stp.py, so L/T built here match stp.py.

A clamp set S "stabilizes" the local network to target z iff, after pinning the
nodes in S to their z-values, EVERY initial local state flows to z. This is
strictly stronger than the heuristic (which only needs the patient's own state
to reach z): expect fewer/larger kernels and honest nulls. Kernels remain
patient-specific through the externals captured in rule_fn.
"""

import numpy as np
from itertools import combinations


# ---------- NVR index convention (identical to stp.py) ----------
def state_to_index(x):
    n = len(x); N = 1 << n
    s = 0
    for i, xi in enumerate(x):
        s += xi * (1 << (n - 1 - i))
    return N - s                                  # 1-based

def index_to_state(k, n):
    N = 1 << n; s = N - k
    return tuple((s >> (n - 1 - i)) & 1 for i in range(n))


# ---------- L (Eq. 3-4): column-index vector, Lcol[i-1] = idx(f(state_i)) ----------
def build_L(rule_fn, n):
    N = 1 << n
    Lcol = np.empty(N, dtype=np.int64)
    for i in range(1, N + 1):
        x = index_to_state(i, n)
        Lcol[i - 1] = state_to_index(tuple(rule_fn(x)))
    return Lcol


# ---------- T_S^z (Def. 1, Eq. 5): overwrite S-coords with z ----------
def build_T(S, z, n):
    N = 1 << n; Sset = set(S)
    Tcol = np.empty(N, dtype=np.int64)
    for j in range(1, N + 1):
        x = list(index_to_state(j, n))
        for idx in Sset:
            x[idx - 1] = z[idx - 1]
        Tcol[j - 1] = state_to_index(tuple(x))
    return Tcol


# ---------- logical-matrix algebra by index composition ----------
def mat_mul(Acol, Bcol):
    return Acol[Bcol - 1]                          # (AB)[j] = A[B[j]]

def mat_pow(Mcol, p):
    result = Mcol.copy()
    for _ in range(p - 1):
        result = result[Mcol - 1]
    return result


# ---------- transient period rho (Prop. 1) ----------
def transient_period(Mcol, max_iter=None):
    N = len(Mcol)
    if max_iter is None:
        max_iter = N + 2
    seen = {}; cur = Mcol.copy(); p = 1
    while p <= max_iter:
        key = cur.tobytes()
        if key in seen:
            return seen[key]
        seen[key] = p
        cur = cur[Mcol - 1]
        p += 1
    return p - 1


# ---------- Omega (Def. 3, Eq. 6) and Theorem-1 test (Eq. 7) ----------
def stabilizability_matrix(Lcol, S, z, n):
    Tcol = build_T(S, z, n)
    Lc = mat_mul(Tcol, Lcol)                       # L_c = T_S^z . L
    rho = transient_period(Lc)
    Omega = mat_mul(mat_pow(Lc, rho), Tcol)        # (L_c)^rho . T_S^z
    return Omega, rho

def is_stabilizing(Lcol, S, z, n):
    m = state_to_index(z)
    Omega, rho = stabilizability_matrix(Lcol, S, z, n)
    return bool(np.all(Omega == m)), rho           # COL(Omega) == {e_m}


# ---------- forward selection: all minimal-cardinality stabilizing kernels ----------
def stabilizing_kernels(Lcol, z, n, candidates=None, max_k=None):
    if candidates is None:
        candidates = list(range(1, n + 1))
    if max_k is None:
        max_k = n - 1
    for size in range(1, max_k + 1):
        hits = [tuple(S) for S in combinations(candidates, size)
                if is_stabilizing(Lcol, S, z, n)[0]]
        if hits:
            return hits, size
    return [], None


# ---------- single-kernel selector with deterministic tiebreak ----------
def stabilize_select_kernel(rule_fn, n, z, candidates=None, max_k=3):
    """
    Build L from rule_fn (a function local-state-tuple -> next-local-state-tuple,
    with externals already closed over), then return the chosen stabilizing
    kernel as a 1-based tuple of node indices, or None if none exists up to max_k.

    Among all minimal-cardinality stabilizing kernels the tiebreak is
    (min causal = sum of 1-based indices, then lexicographic) -- a tiebreak only,
    NOT a selection criterion (every candidate already provably stabilizes).
    """
    if list(z) == [int(v) for v in rule_fn(tuple(z))]:
        # z is a fixed point of the free local dynamics; check trivial-clamp first
        pass
    Lcol = build_L(rule_fn, n)
    hits, size = stabilizing_kernels(Lcol, z, n, candidates=candidates, max_k=max_k)
    if not hits:
        return None
    hits.sort(key=lambda S: (sum(S), S))
    return hits[0]


# ---------- unified dispatcher ----------
def design_kernel(method, *, rule_fn, n, z, x0=None, candidates=None,
                  max_k=3, ranked_fn=None):
    """
    method == 'stabilize' -> algebraic Theorem-1 selector (this module).
    method == 'ranked'    -> delegate to ranked_fn() supplied by the caller
                             (each setup keeps its own exact heuristic path).
    Returns a 1-based tuple of node indices (the kernel) or None.
    """
    if method == 'stabilize':
        return stabilize_select_kernel(rule_fn, n, z, candidates=candidates, max_k=max_k)
    elif method == 'ranked':
        if ranked_fn is None:
            raise ValueError("ranked method requires ranked_fn from the caller")
        return ranked_fn()
    else:
        raise ValueError(f"unknown kernel method: {method!r}")
