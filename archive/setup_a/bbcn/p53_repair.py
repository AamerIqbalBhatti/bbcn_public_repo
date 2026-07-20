"""
p53_repair.py — lift the switch's repaired p53-MDM2-ATM wiring into the full net.

The canonical full network (pathways._akt_survival) has a degenerate p53 axis:
  MDM2 = TP53 or MDM2     (self-locks high)
  TP53 = not MDM2         (no damage input; ATM/ATR never reach it)
so MDM2 latches high, TP53 stays low, and DNA damage does nothing.

This overlay replaces ONLY the MDM2 and TP53 rules with the biology-grounded
form used in the Setup B switch (bbcn_switch.rules):
  MDM2 = (TP53 or AKT1) and not ( CDKN2A and (E2F1 or ATM or ATR) )
  TP53 = (not MDM2) or ATM or ATR
i.e. ATM/ATR (DNA damage) stabilise p53, and the ARF/CDKN2A arm releases MDM2
under oncogenic/damage stress. Everything else in the network is unchanged.

It is an overlay (a swapped step function), so the locked Setup A pipeline and
pathways.py are untouched.
"""
from bbcn.controller import _step

def step_repaired(bus):
    nb = _step(bus)
    g = lambda k: int(bus.get(k, 0))
    nb['MDM2'] = int((g('TP53') or g('AKT1'))
                     and not (g('CDKN2A') and (g('E2F1') or g('ATM') or g('ATR'))))
    nb['TP53'] = int((not g('MDM2')) or g('ATM') or g('ATR'))
    return nb
