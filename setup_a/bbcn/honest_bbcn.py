"""
honest_bbcn.py — biologically-honest delayed full BBCN (ARMA form + repairs + sensor).

Carries, as an overlay on the base step (pathways.py untouched):
  1. p53-MDM2-ATM repair  : MDM2/TP53 rewired to the switch biology (damage drives p53).
  2. AKT1-FOXO3-PHLPP + CLAMP_OFF : the uncoupled-resistance repair, adding PHLPP (the 136th
     node) and the CLAMP_OFF flag (0 = clamp intact / coupled, 1 = clamp failed / uncoupled).
  3. Accumulator (sensor)  : an integrator over p53 (the moving-average / ARMA-MA term,
     the Boolean analogue of p53 pulse-counting). When the integral crosses theta it
     LATCHES a commitment (caspase point-of-no-return) that releases the survival
     brakes (MCL1, XIAP, CFLAR, BCL2, BCL2L1), letting BAX/BAK->CASP fire.

Genuine lag (ARMA) timescales: fast 0, mid d_mid, slow d_slow on the source node.

EXECUTION PATTERN: library. Nothing runs on import. Defines patient_clamp_off (read once at
patient init), step_honest (one delayed Boolean update), and run_honest (a full trajectory).
Other scripts import and call these, e.g. run_durability_full.py and regression_check.py call
step_honest inside their own simulation loops.
"""
from __future__ import annotations
from bbcn import harness as H
from bbcn.controller import _step
from bbcn.multirate import FAST, MID, SLOW

NODES  = list(H.ALL_NODES) + ['PHLPP']          # 136 nodes
BRAKES = ['MCL1', 'XIAP', 'CFLAR', 'BCL2', 'BCL2L1']

def make_lags(d_fast=0, d_mid=5, d_slow=25):
    def lag(nd):
        if nd == 'PHLPP': return d_slow
        if nd in FAST: return d_fast
        if nd in MID:  return d_mid
        return d_slow
    return {nd: lag(nd) for nd in NODES}

# CLAMP_OFF (AKT->14-3-3->FOXO3 nuclear-export clamp has failed): FOXO stays nuclear despite
# active AKT. Named cause = oxidative-stress JNK/MST1 phosphorylating 14-3-3 (cohort cannot measure
# it); defined here by its measurable baseline EFFECT. A passed value, not a bus node, in this engine.
AKT_ACT = ['AKT1', 'MTOR', 'RPS6KB1', 'EIF4EBP1']
FOXO_NUC = ['FOXO3', 'FOXO1']
def patient_clamp_off(bus):
    # Identical definition to repaired_branch.patient_clamp_off.
    # Strict majority of AKT-activity nodes AND strict majority of FOXO-nuclear nodes on:
    #   AKT (4 nodes): at least 3 of 4 on.    FOXO (2 nodes): both on.   (count > half = mean > 0.5)
    akt_on = sum(int(bus.get(k, 0)) for k in AKT_ACT)
    fox_on = sum(int(bus.get(k, 0)) for k in FOXO_NUC)
    return 1 if (akt_on > len(AKT_ACT) / 2 and fox_on > len(FOXO_NUC) / 2) else 0

def commit_signal(s):
    # death-engaged: commit only when p53 is high AND the survival hub is actually lost,
    # NOT on raw p53 (the uncoupled resistant state has high p53 but survives).
    return 1 if (s['TP53'] and s['AKT1'] == 0) else 0

def step_honest(bus, clamp_off=0, genotoxic=0, cascade=False):
    nb = _step(bus)
    g = lambda k: int(bus.get(k, 0))
    # DNA-damage signal = endogenous sensors (ATM/ATR) UNION the genotoxic treatment input.
    # genotoxic defaults to 0, so dmg == (ATM or ATR) and every legacy caller is unchanged.
    # A genotoxic-aware caller passes a per-day, p53-gated value; it is OR-ed with the natural
    # signal, never overriding it. This is the first instance of the input-union (designed +
    # environmental) that small-u will formalise.
    dmg = int(g('ATM') or g('ATR') or genotoxic)
    # 1. p53-MDM2-ATM
    nb['MDM2'] = int((g('TP53') or g('AKT1'))
                     and not (g('CDKN2A') and (g('E2F1') or dmg)))
    nb['TP53'] = int((not g('MDM2')) or dmg)
    # 2. AKT1-FOXO3-PHLPP + CLAMP_OFF  (FOXO nuclear-export clamp failure)
    phlpp = int(g('FOXO3') and not clamp_off)
    nb['PHLPP'] = phlpp
    nb['AKT1']  = int(((g('MTOR') and g('PDPK1') and g('PIK3CA') and not g('PTEN'))
                       and not phlpp) or (clamp_off and g('AKT1')))
    nb['FOXO3'] = int((not g('AKT1')) or clamp_off)
    # 3. p53 -> PUMA/NOXA -> guardian neutralisation  (Phase 5 cascade, GUARDED: default off).
    # The p53-induced BH3-only proteins PUMA (BBC3) and NOXA (PMAIP1) neutralise the anti-apoptotic
    # guardians: NOXA -| MCL1, PUMA -| BCL2/BCL-xL. With the guardians relieved by p53 itself, the
    # existing bare chain BAX/BAK1 -> CYCS -> CASP9 -> CASP3 fires. This is the BIOLOGICAL route the
    # abstract COMMIT flag currently fakes by zeroing brakes. cascade=False -> none of this runs, so
    # every locked number holds; the flag still owns commitment until we deliberately switch it off.
    if cascade:
        puma = int(nb['TP53'])
        noxa = int(nb['TP53'])
        nb['PUMA'] = puma
        nb['NOXA'] = noxa
        nb['MCL1']   = int(nb.get('MCL1', 0)   and not noxa)
        nb['BCL2']   = int(nb.get('BCL2', 0)   and not puma)
        nb['BCL2L1'] = int(nb.get('BCL2L1', 0) and not puma)
        # SMAC/DIABLO, released from the mitochondria with cytochrome c when BAX/BAK fire,
        # neutralises XIAP and unblocks the executioner. Without this arm XIAP silently gates
        # CASP3 off and the cascade stalls at ~56%; with it the real biology reaches ~98%,
        # reproducing (and slightly exceeding) the abstract COMMIT flag.
        smac = int(nb.get('BAX', 0) or nb.get('BAK1', 0))
        nb['SMAC'] = smac
        nb['XIAP'] = int(nb.get('XIAP', 0) and not smac)
        nb['CASP3'] = int((nb.get('CASP8', 0) or nb.get('CASP9', 0)) and not nb['XIAP'])
        # Executioner caspases + amplification feedback (Phase 5 step 15). CASP7 runs parallel to
        # CASP3 under the same XIAP gate; CASP6 is cleaved by the executioners; CASP6 then back-
        # cleaves the initiator CASP8. The resulting CASP6 -> CASP8 -> CASP3 -> CASP6 positive loop
        # is the molecular point of no return: once an executioner fires it self-sustains across
        # steps, so a transient drug pulse that lights the cascade stays committed through the
        # drug-free observation window without any abstract flag.
        nb['CASP7'] = int((nb.get('CASP8', 0) or nb.get('CASP9', 0)) and not nb['XIAP'])
        nb['CASP6'] = int(nb['CASP3'] or nb['CASP7'])
        nb['CASP8'] = int(nb.get('CASP8', 0) or nb['CASP6'])
        nb['CASP3'] = int((nb.get('CASP8', 0) or nb.get('CASP9', 0)) and not nb['XIAP'])
    return nb

def _view(history, t, lags):
    return {nd: history[max(0, t - lags[nd])][nd] for nd in NODES}

def run_honest(bus, clamp_off=0, lags=None, T=400, held=None, theta=10, cmax=20, accumulator=True):
    lags = lags or make_lags()
    bus = {nd: int(bus.get(nd, 0)) for nd in NODES}
    if held:
        for k, v in held.items(): bus[k] = v
    c = 0; commit = 0; history = [bus]
    for t in range(T):
        nb = step_honest(_view(history, t, lags), clamp_off)
        if accumulator:
            c = min(c + 1, cmax) if commit_signal(nb) else max(c - 1, 0)
            commit = commit or (1 if c >= theta else 0)
            if commit:
                for b in BRAKES: nb[b] = 0
        if held:
            for k, v in held.items(): nb[k] = v
        history.append(nb)
    return history[-1], commit, history
