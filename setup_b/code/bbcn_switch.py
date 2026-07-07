"""
bbcn_switch.py — The BBCN bistable resistance<->apoptosis switch.

Single source of truth for:
  * the 9-node core switch rules (with the biology-grounded uncoupled-state amendment, ledger sec.70)
  * the three-box asynchronous multirate engine (rates 1:5:25; FAST/MID/SLOW)
  * fixed-point enumeration
  * the materialised STP structure matrices L1, L2, L3

All update rules are exact Boolean maps. Rate separation lives entirely in the
update SCHEDULE (a slower box fires less often); there are no rate constants,
no Hill functions, and no per-node persistence counters.

Reference: see paper/bbcn_minipaper.tex and docs/verification_record.md.
"""

from __future__ import annotations
import itertools
import numpy as np

# ---------------------------------------------------------------------------
# Node partition into the three rate boxes
# ---------------------------------------------------------------------------
FAST = ["AKT1", "PTEN", "MTOR", "PDPK1", "PIK3CA"]   # rate 1  (phosphorylation)
MID  = ["MDM2"]                                       # rate 5  (degradation/turnover)
SLOW = ["TP53", "PHLPP", "FOXO3"]                     # rate 25 (transcription)
NODES = FAST + MID + SLOW

# Held, patient-derived inputs (plus CLAMP_OFF, the uncoupled Akt-FOXO3a flag).
INPUTS = ["SRC", "RHEB", "IGF1R", "RTK_up", "CDKN2A", "E2F1", "ATM", "ATR"]

RATES = {"FAST": 1, "MID": 5, "SLOW": 25}


def rules(s: dict, I: dict, CLAMP_OFF: int = 0) -> dict:
    """One synchronous Boolean update of all 9 nodes.

    s : current state dict over NODES
    I : input dict over INPUTS
    CLAMP_OFF : uncoupled Akt-FOXO3a flag (0 = normal/coupled, 1 = uncoupled resistant state)

    The three CLAMP_OFF-gated terms (ledger sec.70) are the only difference from the base switch:
      AKT1 gains a self-sustain term  | (CLAMP_OFF & AKT1)
      PHLPP death-brake is released   = FOXO3 & !CLAMP_OFF
      FOXO3 stays nuclear             = !AKT1 | CLAMP_OFF
    When CLAMP_OFF == 0 these collapse to the base switch exactly.
    """
    n = {}
    self_sustain = 1 if (CLAMP_OFF and s["AKT1"] == 1) else 0
    n["AKT1"] = int(((s["MTOR"] and s["PDPK1"] and s["PIK3CA"])
                     and (not s["PTEN"]) and (not s["PHLPP"])) or self_sustain)
    n["PTEN"]   = int((not I["SRC"]) or s["TP53"])
    n["MTOR"]   = int(I["RHEB"] or (s["PIK3CA"] and not s["PTEN"]))
    n["PDPK1"]  = int((s["PIK3CA"] and not s["PTEN"]) or I["IGF1R"])
    n["PIK3CA"] = int(I["RTK_up"])
    n["MDM2"]   = int((s["TP53"] or s["AKT1"])
                      and not (I["CDKN2A"] and (I["E2F1"] or I["ATM"] or I["ATR"])))
    n["TP53"]   = int((not s["MDM2"]) or I["ATM"] or I["ATR"])
    n["PHLPP"]  = int(s["FOXO3"] and not CLAMP_OFF)
    n["FOXO3"]  = int((not s["AKT1"]) or CLAMP_OFF)
    return n


def simulate(x0: dict, I: dict, CLAMP_OFF: int = 0, T: int = 200,
             rate_mid: int = 5, rate_slow: int = 25) -> dict:
    """Exact three-box multirate simulation (no approximation).

    Every tick the FAST box updates; the MID box updates every `rate_mid` ticks;
    the SLOW box every `rate_slow` ticks. All firing boxes read the SAME
    pre-update snapshot and commit together.
    """
    s = dict(x0)
    for t in range(1, T + 1):
        snap = dict(s)
        nxt = rules(snap, I, CLAMP_OFF)
        for nd in FAST:
            s[nd] = nxt[nd]
        if t % rate_mid == 0:
            s["MDM2"] = rules(snap, I, CLAMP_OFF)["MDM2"]
        if t % rate_slow == 0:
            for nd in SLOW:
                s[nd] = rules(snap, I, CLAMP_OFF)[nd]
    return s


def label(s: dict) -> str:
    """Phenotype label of a settled state."""
    if s["TP53"] == 1 and s["AKT1"] == 0:
        return "APOPTOTIC"
    if s["AKT1"] == 1:
        return "SURVIVAL"
    return "mixed"


def fixed_points(I: dict, CLAMP_OFF: int = 0) -> list:
    """Enumerate the synchronous fixed points (x = f(x)) over all 2^9 states."""
    fps = []
    for bits in itertools.product([0, 1], repeat=len(NODES)):
        s = dict(zip(NODES, bits))
        if rules(s, I, CLAMP_OFF) == s:
            fps.append(tuple(bits))
    return fps


# ---------------------------------------------------------------------------
# Materialised STP structure matrices for the three boxes (fixed input regime)
# ---------------------------------------------------------------------------
def _box_inputs(box):
    """External nodes each box reads (other-box nodes + true inputs)."""
    deps = {
        "AKT1": ["MTOR", "PDPK1", "PIK3CA", "PTEN", "PHLPP"],
        "PTEN": ["SRC", "TP53"], "MTOR": ["RHEB", "PIK3CA", "PTEN"],
        "PDPK1": ["PIK3CA", "PTEN", "IGF1R"], "PIK3CA": ["RTK_up"],
        "MDM2": ["TP53", "AKT1", "CDKN2A", "E2F1", "ATM", "ATR"],
        "TP53": ["MDM2", "ATM", "ATR"], "PHLPP": ["FOXO3"], "FOXO3": ["AKT1"],
    }
    ext = []
    for nd in box:
        for d in deps[nd]:
            if d not in box and d not in ext:
                ext.append(d)
    return ext


def build_box_matrix(box, I, CLAMP_OFF=0):
    """Materialise a box's STP structure matrix L (2^k x 2^(k+m)), column one-hot.

    Returns (L, others) where `others` are the external nodes read by the box.
    """
    others = [o for o in _box_inputs(box) if o not in INPUTS]
    k, m = len(box), len([o for o in _box_inputs(box)])
    full_ext = _box_inputs(box)
    L = np.zeros((2 ** k, 2 ** (k + m)), dtype=int)
    for bs in range(2 ** k):
        bbits = [(bs >> (k - 1 - i)) & 1 for i in range(k)]
        for os in range(2 ** m):
            obits = [(os >> (m - 1 - i)) & 1 for i in range(m)]
            s = {nd: b for nd, b in zip(box, bbits)}
            for nd, b in zip(full_ext, obits):
                if nd in NODES:
                    s[nd] = b
                else:
                    I = dict(I); I[nd] = b
            for nd in NODES:
                s.setdefault(nd, 0)
            nb = rules(s, I, CLAMP_OFF)
            ncol = 0
            for nd in box:
                ncol = (ncol << 1) | nb[nd]
            L[ncol, bs * (2 ** m) + os] = 1
    return L, others


if __name__ == "__main__":
    # quick self-test: resting regime fixed point + bistability demo
    rest = dict(SRC=0, RHEB=1, IGF1R=1, RTK_up=1, CDKN2A=1, E2F1=0, ATM=0, ATR=0)
    print("resting fixed points:", fixed_points(rest))
    surv0 = {n: 0 for n in NODES}; surv0.update(AKT1=1, MDM2=1, MTOR=1, PDPK1=1, PIK3CA=1)
    apop0 = {n: 0 for n in NODES}; apop0.update(TP53=1, PHLPP=1, FOXO3=1, PTEN=1)
    print("surv-seed ->", label(simulate(surv0, rest)))
    print("apop-seed ->", label(simulate(apop0, rest)))
