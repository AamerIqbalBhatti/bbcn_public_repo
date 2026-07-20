"""
prolif_switch.py — The BBCN bistable proliferation (restriction-point) switch.

Companion to bbcn_switch.py (the apoptosis switch), built on the SAME template:
a small bistable core with timescale separation, run on a 1:5:25 multirate
schedule. Here the decision axis is the RB1-E2F1 restriction point. The
double-negative loop is RB1 -| E2F1 and E2F1 -> CCNE1 -> CDK2 -| RB1; the E2F1 ->
CCNE1 -> CDK2 -> (RB phosphorylation) -> E2F1 positive feedback is the hysteretic
commitment to division, the proliferation analogue of the apoptosis latch.

Node convention: 1 = active/high. RB1=1 means active (hypophosphorylated) RB that
represses E2F1; RB1=0 means phosphorylated/inactive RB.
"""
from __future__ import annotations
import itertools

# three rate boxes (same mechanism split as the global classification)
FAST = ["CDK4", "CDK6", "CDK2", "RB1"]                 # rate 1  (phosphorylation)
MID  = ["CCNE1", "CDKN1B"]                             # rate 5  (turnover: cyclin E, p27)
SLOW = ["CCND1", "E2F1", "MYC", "CDKN1A", "CDKN2A"]    # rate 25 (transcription)
NODES = FAST + MID + SLOW

# held, patient-derived inputs
INPUTS = ["MITOGEN", "TP53", "STRESS"]                 # growth drive; p53 arrest; senescence/p16
RATES = {"FAST": 1, "MID": 5, "SLOW": 25}


def rules(s: dict, I: dict) -> dict:
    """One synchronous Boolean update of all 11 nodes."""
    n = {}
    # SLOW transcriptional layer
    n["CCND1"]  = int(I["MITOGEN"] or s["MYC"])                  # cyclin D from mitogen / MYC
    n["MYC"]    = int(I["MITOGEN"] or s["E2F1"])                 # MYC from mitogen / E2F autoreg
    n["E2F1"]   = int(not s["RB1"])                             # E2F free iff RB inactive
    n["CDKN1A"] = int(I["TP53"] and not s["MYC"])               # p21 from p53, repressed by MYC
    n["CDKN2A"] = int(I["STRESS"])                              # p16 from oncogenic/replicative stress
    # MID turnover layer
    n["CCNE1"]  = int(s["E2F1"])                                # cyclin E is an E2F target
    n["CDKN1B"] = int(not (s["CDK2"] or s["MYC"]))              # p27 degraded by CDK2 / MYC-Skp2
    # FAST phospho layer
    n["CDK4"]   = int(s["CCND1"] and not s["CDKN2A"] and not s["CDKN1A"])
    n["CDK6"]   = int(s["CCND1"] and not s["CDKN2A"] and not s["CDKN1A"])
    n["CDK2"]   = int(s["CCNE1"] and not (s["CDKN1A"] or s["CDKN1B"]))
    n["RB1"]    = int(not (s["CDK4"] or s["CDK6"] or s["CDK2"]))  # RB inactivated by any active CDK
    return n


def label(s: dict) -> str:
    if s["RB1"] == 1 and s["E2F1"] == 0: return "QUIESCENT"      # proliferation OFF
    if s["E2F1"] == 1 and s["RB1"] == 0: return "PROLIFERATIVE"  # proliferation ON
    return "mixed"


def simulate(x0: dict, I: dict, T: int = 300, rate_mid: int = 5, rate_slow: int = 25) -> dict:
    s = dict(x0)
    for t in range(1, T + 1):
        snap = dict(s); nxt = rules(snap, I)
        for nd in FAST: s[nd] = nxt[nd]
        if t % rate_mid == 0:
            for nd in MID: s[nd] = nxt[nd]
        if t % rate_slow == 0:
            for nd in SLOW: s[nd] = nxt[nd]
    return s


def fixed_points(I: dict) -> list:
    fps = []
    for bits in itertools.product([0, 1], repeat=len(NODES)):
        s = dict(zip(NODES, bits))
        if rules(s, I) == s:
            fps.append(tuple(bits))
    return fps
