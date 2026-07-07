"""
delayed_bbcn.py — genuine delay-difference (lag) engine for the full BBCN.

Unlike multirate.py (which commits slow nodes on a CADENCE), here EVERY node
updates every tick, but each regulator is SEEN AT ITS OWN LAG: a slow regulator
contributes its value from d_slow steps ago, a medium one from d_mid ago, a fast
one from the current step. This is the literal Boolean-ARMA / delay-difference
form x(t+1) = f( x(t), x(t-d_mid), x(t-d_slow) ); the lag lives on the source node.

Reuses the 81/10/44 timescale classification from multirate.py.
"""
from __future__ import annotations
from bbcn import harness as H
from bbcn.controller import _step
from bbcn.multirate import FAST, MID, SLOW

def make_lags(d_fast=0, d_mid=1, d_slow=2):
    return {nd: (d_fast if nd in FAST else d_mid if nd in MID else d_slow)
            for nd in H.ALL_NODES}

def _view(history, t, lags):
    """Each node shows its value from `lag` steps ago (delayed visibility)."""
    return {nd: history[max(0, t - lags[nd])][nd] for nd in H.ALL_NODES}

def step_lagged(history, t, lags, held=None, step_fn=_step):
    nb = step_fn(_view(history, t, lags))    # all rules read the lagged view
    if held:
        for k, v in held.items(): nb[k] = v
    return nb

def run(bus, lags=None, T=400, held=None, step_fn=_step):
    lags = lags or make_lags()
    bus = {nd: int(bus.get(nd, 0)) for nd in H.ALL_NODES}
    if held:
        for k, v in held.items(): bus[k] = v
    history = [bus]; win = max(lags.values()) + 2
    for t in range(0, T):
        nb = step_lagged(history, t, lags, held=held, step_fn=step_fn)
        history.append(nb)
        if t >= win and all(history[-1] == history[-1-j] for j in range(1, win+1)):
            return nb, True, t, history
    return history[-1], False, T, history
