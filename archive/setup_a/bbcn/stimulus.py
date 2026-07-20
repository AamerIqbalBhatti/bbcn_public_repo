"""
stimulus.py -- treatment-input vector for the BBCN clinical model (Phase 4, step 11).

WHY. Until now the only "drug" in the durability / clock arms is the kernel CLAMP: it forces the
designed apoptosis-target nodes (the targeted AKT-inhibitor / capivasertib reading). Real therapy
also includes GENOTOXIC modalities (chemo, radiation) that act UPSTREAM by causing DNA damage,
which activates the ATM/ATR sensors, which drive p53, which drives apoptosis. This module adds a
small, explicit TREATMENT-INPUT VECTOR so genotoxic stress can be routed through ATM/ATR on the
clinical clock, alongside the targeted clamp.

VERIFIED MECHANICS (DNA_Repair pathway + honest_bbcn.step_honest):
  - ATM and ATR are NOT pure read-only inputs. They are nodes of the DNA_Repair pathway whose
    update rule is the IDENTITY (ATM <- ATM, ATR <- ATR): a self-sustaining latch. So whatever
    value they are SEEDED to from the patient's binarised genomics PERSISTS through the run on its
    own -- no `held` dict needed. (This is why the genomic baseline, part (c), already holds.)
  - Nothing upstream currently sets ATM/ATR. A patient seeded ATM=0 keeps ATM=0 forever.
    Genotoxic therapy is precisely the missing upstream: a stimulus that ENGAGES the latch (ATM<-1).
  - step_honest reads the sensors in the p53 arm:  TP53 = (not MDM2) or ATM or ATR.
    So engaging the damage latch pushes p53 up, feeding the commitment accumulator.

ONE-CHANGE DISCIPLINE. Step 11 ships in two inert-by-default sub-steps:
  (i)  THIS module -- the input-vector structure, the genomic-baseline reader, and the genotoxic
       routing logic, unit-tested standalone but wired into NO runner. Locked numbers cannot move.
  (ii) WIRING -- route the stimulus into the clock runners behind an opt-in arg that defaults OFF,
       then measure the genotoxic arm at full N.

EXECUTION PATTERN: library. Nothing runs on import.
"""
from __future__ import annotations

# The DNA-damage sensors a genotoxic dose can engage. Both are p53 drivers in step_honest.
DAMAGE_NODES = ('ATM', 'ATR')


def genomic_baseline(b):
    """Untreated DNA-damage sensor state, read from the patient's binarised genomics.
    This is the value ATM/ATR latch to in the absence of any genotoxic stimulus (step 11, part c)."""
    return {nd: int(b.get(nd, 0)) for nd in DAMAGE_NODES}


def make_stimulus(genotoxic=False, hits=DAMAGE_NODES):
    """Build a treatment-input vector.

    genotoxic=False reproduces current behaviour exactly (no overrides -> locked numbers hold).
    genotoxic=True  means a genotoxic agent is part of the regimen; on a dosing on-day it engages
                    the damage sensors named in `hits` (default both ATM and ATR).
    `hits` lets a caller model an agent that drives only double-strand-break sensing (ATM) or only
    replication-stress sensing (ATR); the default engages both.
    """
    hits = tuple(h for h in hits if h in DAMAGE_NODES)
    return {'genotoxic': bool(genotoxic), 'hits': hits}


DEFAULT_STIMULUS = make_stimulus()   # genotoxic OFF: fully inert


def stimulus_overrides(stim, dosing_day, p53_lof=False):
    """Node overrides contributed by the treatment-input vector on this step.

    Returns a dict {node: 1, ...} of damage sensors to force, or {} when nothing applies.
    Genotoxic only acts on a DOSING on-day (the same gate the clamp uses); off-days contribute
    nothing, so the latch then simply holds whatever it reached.

    NARROW p53 GATE (Phase 4). Genotoxic's only death route in this network is ATM/ATR -> p53 ->
    commitment. In a p53 loss-of-function tumour that route is dead, so genotoxic produces no death
    signal: we return {} for a p53_lof patient. This is the narrow scope -- it touches ONLY the
    genotoxic arm; the targeted-clamp readout and the global TP53 rule are left untouched. The gate
    status per patient comes from load_p53_gate (mutation-derived for TCGA/METABRIC, proxy for I-SPY2).
    """
    if not stim or not stim.get('genotoxic') or not dosing_day:
        return {}
    if p53_lof:
        return {}                      # genotoxic cannot rescue a p53-null tumour
    return {nd: 1 for nd in stim.get('hits', ())}


def load_p53_gate(path):
    """Load a per-patient p53 loss-of-function gate into {sample_id: 0/1}.

    Accepts either a mutation-derived gate (column TP53_LOF) or the I-SPY2 proxy (column
    TP53_LOF_PROXY). 1 = p53 non-functional (genotoxic death route disabled for that patient).
    """
    import csv
    gate = {}
    with open(path, newline='') as fh:
        r = csv.DictReader(fh)
        col = 'TP53_LOF' if 'TP53_LOF' in r.fieldnames else 'TP53_LOF_PROXY'
        for row in r:
            gate[row['sample_id']] = int(row[col])
    return gate


if __name__ == '__main__':
    # standalone unit test: prove the latch, the genomic baseline, and genotoxic engagement,
    # all WITHOUT touching any runner. Uses zero lags to isolate the rule logic from ARMA delays.
    import sys
    sys.path.insert(0, '.'); sys.path.insert(0, '..'); sys.path.insert(0, '../setup_b/code')
    from bbcn import honest_bbcn as HB

    lags = {nd: 0 for nd in HB.NODES}

    def short_run(seed_atm, stim, steps=20):
        b = {n: 0 for n in HB.NODES}
        b['ATM'] = seed_atm
        hist = [dict(b)]
        for t in range(steps):
            nb = HB.step_honest(HB._view(hist, t, lags), 0)
            for k, v in stimulus_overrides(stim, dosing_day=True).items():
                nb[k] = v
            hist.append(nb)
        return hist[-1]

    # (c) genomic baseline reader
    assert genomic_baseline({'ATM': 1, 'ATR': 0}) == {'ATM': 1, 'ATR': 0}

    # latch: untreated ATM stays where it was seeded
    off = short_run(0, DEFAULT_STIMULUS)
    on  = short_run(1, DEFAULT_STIMULUS)
    print(f"latch  seed ATM=0 -> ATM={off['ATM']} (expect 0) | seed ATM=1 -> ATM={on['ATM']} (expect 1)")
    assert off['ATM'] == 0 and on['ATM'] == 1

    # genotoxic engages the latch from cold and lifts p53 (p53-functional patient)
    g = short_run(0, make_stimulus(genotoxic=True))
    print(f"genotoxic seed ATM=0 -> ATM={g['ATM']} ATR={g['ATR']} TP53={g['TP53']} (expect 1,1,1)")
    assert g['ATM'] == 1 and g['ATR'] == 1 and g['TP53'] == 1

    # NARROW GATE: a p53-LOF patient gets no genotoxic override -> stays cold
    assert stimulus_overrides(make_stimulus(genotoxic=True), dosing_day=True, p53_lof=True) == {}
    assert stimulus_overrides(make_stimulus(genotoxic=True), dosing_day=True, p53_lof=False) == {'ATM': 1, 'ATR': 1}
    print("narrow p53 gate: genotoxic silent for p53-LOF, active for p53-functional -> OK")

    # gate loader round-trips a real gate file (mutation-derived)
    import os
    gp = os.path.join(os.path.dirname(__file__), '..', 'data', 'mutation_gate', 'tp53_lof_tcga.csv')
    if os.path.exists(gp):
        gate = load_p53_gate(gp)
        print(f"load_p53_gate: {len(gate)} patients, {sum(gate.values())} p53-LOF")
        assert set(gate.values()) <= {0, 1}

    # default vector is inert
    assert stimulus_overrides(DEFAULT_STIMULUS, dosing_day=True) == {}
    print("stimulus.py self-test PASSED")
