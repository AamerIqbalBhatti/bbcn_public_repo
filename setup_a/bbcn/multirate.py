"""
multirate.py — delayed/asynchronous full-BBCN engine (the 'test big' arena).

Wraps the exact base-network synchronous step (controller._step: every pathway's
rules(x, bus) applied with the full bus as externals) and commits each node's
update only on the ticks of its timescale class. Fast every tick, medium every 5,
slow every 25 (1:5:25), using the mechanism classification 81/10/44.

A genuine fixed point is unchanged under the base step, so it is invariant under
ANY schedule; the engine never moves it. That is verification check 2 by
construction; check 1 (rate-1 parity) and the settle behaviour are tested in the
companion runner.
"""
from __future__ import annotations
from bbcn import harness as H
from bbcn.controller import _step

FAST = set("EGFR GRB2 GRB7 ERBB2 ERBB3 INSR IRS1 PIK3CA ABL1 PTEN PDPK1 RAC1 NOS3 HRAS MAP2K1 RAF1 KRAS SOS1 NF1 MAPK8 MAPK1 DUSP1 DVL1 FZD3 RND1 ABL2 EIF4EBP1 JAK2 LATS1 NF2 AKT1 SGK1 GSK3B RHEB MTOR TSC2 EIF4E IKBKB TLR2 TLR4 XIAP CFLAR PAK1 PRKACA SRC FAS FADD TRADD CASP8 BID CASP3 CASP9 BAD BAK1 BAX BCL2 BCL2L1 BCL2L11 CYCS APAF1 RB1 CDK4 CDK6 CDK2 ATM ATR BRCA1 BRCA2 PARP1 CHEK1 CHEK2 RAD51 HSP90AA1 IGF1R IL6R LCK PIM1 PIM2 PIM3 PRKAA2 RPS6KB1".split())
MID  = set("APC AXIN1 CTNNB1 SOCS3 NEDD4L NFKBIA MDM2 MCL1 CCNE1 CDKN1B".split())
SLOW = set("EGF INS AR PGR ESR1 KMT2D TCF7L2 WNT1 DLL1 HES1 HEY1 MYOD1 MYOG NOTCH1 STAT3 STAT1 YAP1 FOXO1 FOXO3 TWIST1 CEBPB CREB1 GATA3 CXCL8 MYC RELA IL1A IL1B TP53 JUN FASLG CDKN2A CCND1 E2F1 CDKN1A TEAD1 PDCD1 CD274 IL6 NRG1 PAX7 IFNG PPARG TNF".split())

def rate_of(nd, r_fast=1, r_mid=5, r_slow=25):
    if nd in FAST: return r_fast
    if nd in MID:  return r_mid
    return r_slow  # SLOW

def multirate_step(bus, t, held=None, rates=(1,5,25)):
    """One multirate tick. Fast/mid/slow commit on their cadence; held pins override."""
    nxt = _step(bus)                      # exact base synchronous image of all nodes
    nb = dict(bus)
    for nd in H.ALL_NODES:
        if t % rate_of(nd, *rates) == 0:
            nb[nd] = nxt[nd]
    if held:
        for k, v in held.items(): nb[k] = v
    return nb

def is_fixed(bus):
    return _step(bus) == bus

def run(bus, T=400, held=None, rates=(1,5,25), settle_window=25):
    """Run T ticks. Detect settling = state constant over a full schedule period."""
    bus = {nd: int(bus.get(nd,0)) for nd in H.ALL_NODES}
    history = [tuple(bus[n] for n in H.ALL_NODES)]
    for t in range(1, T+1):
        bus = multirate_step(bus, t, held=held, rates=rates)
        history.append(tuple(bus[n] for n in H.ALL_NODES))
        if t >= settle_window and len(set(history[-settle_window-1:])) == 1:
            return bus, True, t   # settled to a fixed point
    return bus, is_fixed(bus), T
