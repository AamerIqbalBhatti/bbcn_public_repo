"""
bbcn.stp
========
Literal Semi-Tensor-Product (STP) machinery, a faithful port of the MATLAB
utilities:
    generateL_from_rules.m
    normal_vector.m
    targeted_state_mod_matrix.m
    select_best_kernel.m
    new_forward_selection_bn_ranked.m

STP encoding convention (preserved verbatim from MATLAB):
    A length-n Boolean state x (left-MSB bit order) maps to a one-hot
    delta vector of length N = 2^n, with the single 1 at position
        idx = N - bi2de(x, 'left-msb')          (1-based in MATLAB)
    The transition matrix L is N x N with L[row, col] = 1 where
        col = N - state_idx     (state_idx = integer value of x, left-MSB)
        row = N - bi2de(next_x)
    i.e. column j corresponds to the state de2bi(N - j).

We keep MATLAB's 1-based indices internally for exact parity, converting to
0-based only at array access. Node index sets S are 1-based (as in MATLAB),
because select_best_kernel's causal score = sum(S) depends on it.

WARNING: L is dense 2^n x 2^n. BBCN pathways have up to ~9-10 nodes, so the
largest L is 2^10 x 2^10 = 1024x1024, which is fine. The literal construction
is O(2^n) per pathway per externals-configuration.
"""

from typing import List, Dict, Tuple, Callable, Optional
from itertools import combinations
import numpy as np


# ============================================================
# Bit <-> index conversions (left-MSB, matching de2bi/bi2de 'left-msb')
# ============================================================

def bits_to_int_leftmsb(x: List[int]) -> int:
    """bi2de(x, 'left-msb'): x[0] is the most significant bit."""
    v = 0
    for b in x:
        v = (v << 1) | (1 if b else 0)
    return v


def int_to_bits_leftmsb(val: int, n: int) -> List[int]:
    """de2bi(val, n, 'left-msb'): returns length-n list, MSB first."""
    return [(val >> (n - 1 - i)) & 1 for i in range(n)]


def state_to_index(x: List[int]) -> int:
    """
    normal_vector index (1-based, MATLAB):  idx = 2^n - bi2de(x,'left-msb').
    Returned 1-based to match MATLAB; callers subtract 1 for numpy.
    """
    n = len(x)
    N = 1 << n
    return N - bits_to_int_leftmsb(x)


def index_to_state(idx_1based: int, n: int) -> List[int]:
    """Inverse of state_to_index: state at 1-based delta-vector position idx."""
    N = 1 << n
    return int_to_bits_leftmsb(N - idx_1based, n)


# ============================================================
# generateL_from_rules.m  (literal)
# ============================================================

_L_CACHE: Dict[tuple, np.ndarray] = {}


def generate_L(rule_fn: Callable, n: int, ext: Dict[str, int],
               cache_key: Optional[tuple] = None) -> np.ndarray:
    """
    Build the 2^n x 2^n canonical STP transition matrix L for a pathway.

    rule_fn(x, e) -> next-state list (length n), exactly as bbcn.pathways rules.
    ext is the externals dict (bus) read by the rules.

    L[row-1, col-1] = 1 where (1-based):
        col = N - state_idx     ;  state_idx = bi2de(x,'left-msb')
        row = N - bi2de(next_x)

    If cache_key is given, the result is memoized: L depends only on the
    pathway's rules and the externals it reads, so identical keys reuse L.
    """
    if cache_key is not None and cache_key in _L_CACHE:
        return _L_CACHE[cache_key]
    N = 1 << n
    L = np.zeros((N, N), dtype=np.int8)
    for state_idx in range(N):                 # 0..N-1, integer value of x
        x = int_to_bits_leftmsb(state_idx, n)
        nxt = rule_fn(x, ext)
        nxt = [1 if b else 0 for b in nxt]
        row = N - bits_to_int_leftmsb(nxt)     # 1-based
        col = N - state_idx                    # 1-based
        L[row - 1, col - 1] = 1
    if cache_key is not None:
        _L_CACHE[cache_key] = L
    return L


def normal_vector(x: List[int]) -> np.ndarray:
    """2^n x 1 one-hot delta vector for state x (1 at state_to_index(x))."""
    n = len(x)
    N = 1 << n
    vec = np.zeros(N, dtype=np.int8)
    vec[state_to_index(x) - 1] = 1
    return vec


def targeted_state_mod_matrix(S_1based: Tuple[int, ...], z: List[int]) -> np.ndarray:
    """
    T (N x N) that forces nodes in S (1-based indices) to their values in z.
    For each column j (1-based), original state x = de2bi(N-j); set x[S]=z[S];
    new index = N - bi2de(x); T[new-1, j-1] = 1.
    """
    n = len(z)
    N = 1 << n
    T = np.zeros((N, N), dtype=np.int8)
    S0 = [s - 1 for s in S_1based]             # to 0-based for list access
    for j in range(1, N + 1):
        x = int_to_bits_leftmsb(N - j, n)
        for s in S0:
            x[s] = z[s]
        new_idx = N - bits_to_int_leftmsb(x)   # 1-based
        T[new_idx - 1, j - 1] = 1
    return T


# ============================================================
# Forward kernel search  (new_forward_selection_bn_ranked.m, literal STP)
# ============================================================

def _delta_improvement(Lc: np.ndarray, x0_vec: np.ndarray,
                       z: List[int], n: int) -> int:
    """
    #incorrect bits of x0 corrected after one step (literal de2bi math).
    x1_idx = find(Lc * x0_vec); x1 = de2bi(N - x1_idx); delta vs z.
    """
    N = 1 << n
    x1_vec = Lc @ x0_vec
    x1_idx = int(np.argmax(x1_vec)) + 1            # 1-based find
    x1 = int_to_bits_leftmsb(N - x1_idx, n)
    x0_idx = int(np.argmax(x0_vec)) + 1
    x0 = int_to_bits_leftmsb(N - x0_idx, n)
    return sum(1 for i in range(n)
               if x0[i] != z[i] and x1[i] == z[i])


def select_best_kernel(results: List[dict]) -> int:
    """
    Faithful port of select_best_kernel.m.
    Ranking (size NOT filtered — idx1 keeps all, per the MATLAB):
        max impact -> max delta -> min steps -> min causal_score
    Returns index (into results) of the chosen kernel (first survivor).
    """
    if not results:
        return -1
    idx1 = list(range(len(results)))
    max_impact = max(results[i]['impact'] for i in idx1)
    idx2 = [i for i in idx1 if results[i]['impact'] == max_impact]
    max_delta = max(results[i]['delta'] for i in idx2)
    idx3 = [i for i in idx2 if results[i]['delta'] == max_delta]
    min_steps = min(results[i]['steps'] for i in idx3)
    idx4 = [i for i in idx3 if results[i]['steps'] == min_steps]
    min_causal = min(results[i]['causal'] for i in idx4)
    idx5 = [i for i in idx4 if results[i]['causal'] == min_causal]
    return idx5[0]


def forward_select_kernel(L: np.ndarray, z: List[int], x0: List[int],
                          max_k: int = 2,
                          rule_fn=None, ext=None) -> Optional[Tuple[int, ...]]:
    """
    Literal port of new_forward_selection_bn_ranked.m.

    Enumerate all node subsets of size 1..max_k. For each S, pin nodes S to
    their target values and simulate forward (up to min(2^(n-|S|),200) steps,
    cycle-detected); a subset matches if the trajectory reaches z. Collect all
    matches, rank with select_best_kernel, return chosen S (1-based) or None.

    The forward simulation is done by direct rule iteration when `rule_fn` is
    supplied (O(n) per step), which is behaviourally identical to multiplying by
    Lc = T@L but avoids constructing the 2^n x 2^n matrix. If rule_fn is None it
    falls back to the literal matrix simulation using L.
    """
    n = len(z)

    def sim_matrix(S):
        from numpy import argmax
        T = targeted_state_mod_matrix(S, z)
        Lc = T @ L
        z_idx = int(argmax(normal_vector(z)))
        x_curr = normal_vector(x0).astype(np.int64)
        r = min(2 ** (n - len(S)), 200)
        visited = set()
        for t in range(1, r + 1):
            ic = int(argmax(x_curr))
            if ic in visited:
                break
            visited.add(ic)
            x_next = Lc @ x_curr
            if int(argmax(x_next)) == z_idx:
                return True, t
            x_curr = x_next
        return (int(argmax(x_curr)) == z_idx), r

    def sim_direct(S):
        # Lc = T @ L applies L (rules) FIRST, then T (pin). Each step: step the
        # rules, THEN force nodes in S to target. The initial x0 is NOT pre-pinned
        # (the matrix starts from normal_vector(x0) unmodified), matching MATLAB.
        S0 = [s - 1 for s in S]
        x = list(x0)
        r = min(2 ** (n - len(S)), 200)
        visited = set()
        for t in range(1, r + 1):
            key = tuple(x)
            if key in visited:
                break
            visited.add(key)
            nxt = [int(v) for v in rule_fn(x, ext)]
            for i in S0:
                nxt[i] = z[i]
            if nxt == z:
                return True, t
            x = nxt
        return (x == z), r

    sim = sim_direct if rule_fn is not None else sim_matrix

    def delta_direct(S):
        # one step of Lc from x0: rules first, then pin S to target
        S0 = [s - 1 for s in S]
        x1 = [int(v) for v in rule_fn(x0, ext)]
        for i in S0:
            x1[i] = z[i]
        return sum(1 for i in range(n) if x0[i] != z[i] and x1[i] == z[i])

    results = []
    for ksize in range(1, max_k + 1):
        for S in combinations(range(1, n + 1), ksize):   # 1-based, ascending
            matched, steps = sim(S)
            if matched:
                impact = sum(1 for s in S if x0[s - 1] != z[s - 1])
                if rule_fn is not None:
                    delta = delta_direct(S)
                else:
                    delta = _delta_improvement(targeted_state_mod_matrix(S, z) @ L,
                                               normal_vector(x0), z, n)
                causal = sum(S)
                results.append({'nodes': S, 'size': ksize, 'impact': impact,
                                'delta': delta, 'steps': steps, 'causal': causal})

    if not results:
        return None
    best = select_best_kernel(results)
    return results[best]['nodes']
