"""
durability_clock.py -- pulsed (withdrawal) durability on the clinical clock (Phase 3, step 8).

This is the FAITHFUL dosing arm. The continuous-hold durability (run_durability_full.A_dyn,
clamp held while t<T1 then released) assumed the drug is on every day of treatment. Real
capivasertib is 4-on / 3-off. Step 8 replaces the continuous hold with the clock's pulse:

    horizon      = clinical_clock.total_days()        = 70 days
    treatment    = clinical_clock.treatment_days()    = 42 days (first 6 weeks)
    observation  = clinical_clock.observation_days()  = 28 days (last 4 weeks, drug-free)
    dosing       = clinical_clock.is_on_day(day)      (ON first 4 of each treatment week)

The durability verdict is the honest one (E64): CASP3 sustained over the drug-free observation
tail. The single switch `pulsed` flips between the faithful pulse and the old continuous hold,
on the SAME 70-day clock, so a head-to-head isolates exactly the cost of the off-day gaps.

This module does not touch run_durability_full.A_dyn or the regression baseline; it is the new
arm, run for comparison. See ledger E67.
"""
import sys
sys.path.insert(0, 'setup_a'); sys.path.insert(0, 'setup_b/code')
from bbcn import honest_bbcn as HB
from bbcn import clinical_clock as CK

NODES = HB.NODES
lags = HB.make_lags()


def A_dyn_clock(b, U, clamp, pulsed=True, genotoxic=False, p53_lof=0, sticky=False,
                theta=10, cmax=20):
    """Full-network durability over the 70-day clock.

    pulsed=True  -> drug dosed on clinical_clock.is_on_day(day)  (4-on/3-off, faithful)
    pulsed=False -> drug held every treatment day                (continuous, old assumption)

    genotoxic=True adds a genotoxic agent on the same dosing schedule. It enters the model as the
    `genotoxic` damage input to step_honest (unioned with endogenous ATM/ATR at the p53 node).
    p53_lof=1 (from the mutation/proxy gate) disables it for that patient -- genotoxic cannot rescue
    a p53-null tumour. sticky=False (default) is resolving damage: the genotoxic signal recovers on
    off-days. sticky=True keeps it engaged once triggered (persistent unrepaired damage); the two
    bracket the durability, since commitment, once latched, carries most kills either way.

    Observation (days 43-70) is drug-free for the clamp in both. Returns True if CASP3 stays
    committed through the drug-free tail (durable kill).
    """
    bus = {n: int(b.get(n, 0)) for n in NODES}
    for k, v in clamp.items():
        bus[k] = v
    c = 0; commit = 0; hist = [bus]; triggered = 0
    treat = CK.treatment_days()                      # 42
    for t in range(CK.total_days()):                 # 70 days
        day = t + 1
        dosed = CK.is_on_day(day) if pulsed else (day <= treat)   # <-- the pulse vs the hold
        gtx = int((dosed or (sticky and triggered)) and genotoxic and not p53_lof)
        if gtx:
            triggered = 1
        nb = HB.step_honest(HB._view(hist, t, lags), U, genotoxic=gtx, cascade=True)
        # Step 16: commitment is owned by the real p53 -> PUMA/NOXA -> SMAC -> executioner cascade.
        # Durability is carried by the bistable AKT1-FOXO3-PHLPP switch the kernel latches off
        # (proven: breaking PHLPP collapses durability 100% -> 41%). The old COMMIT-flag accumulator
        # and brake-zeroing are gone; the cascade reads apoptosis out of the held survival-off state.
        if dosed:
            for k, v in clamp.items():
                nb[k] = v
        hist.append(nb)
    tail = hist[-CK.observation_days():]             # the 28 drug-free observation days
    return sum(s['CASP3'] == 1 for s in tail) / len(tail) > 0.5


def run_on_clock(b, U, clamp, pulsed=True, theta=10, cmax=20):
    """Step 9 -- the run lives on the clock, week by week, with explicit phases.

    Same dynamics as A_dyn_clock, but structured into the clinical timeline and reported:
      * each WEEK is a cycle; weeks 1-4 induction, 5-6 maintenance, 7-10 observation (drug-free);
      * we record, per week, the phase, how many days were dosed, and whether CASP3 is committed
        at week end, plus the day commitment first latched;
      * durability = CASP3 sustained across the whole 4-week observation window (not just a tail).
    Returns a dict with the per-week timeline and the durability verdict.
    """
    bus = {n: int(b.get(n, 0)) for n in NODES}
    for k, v in clamp.items():
        bus[k] = v
    c = 0; commit = 0; hist = [bus]
    treat = CK.treatment_days()
    latch_day = None
    weeks = {w: {'phase': None, 'dosed_days': 0, 'casp3_end': 0} for w in range(1, CK.TOTAL_WEEKS + 1)}

    for t in range(CK.total_days()):
        day = t + 1
        w = CK.week_of_day(day)
        weeks[w]['phase'] = CK.phase_of_day(day)
        nb = HB.step_honest(HB._view(hist, t, lags), U, cascade=True)   # Step 16: real cascade
        c = min(c + 1, cmax) if HB.commit_signal(nb) else max(c - 1, 0)
        if c >= theta and not commit:
            commit = 1
            latch_day = day
        elif c >= theta:
            commit = 1
        # commit/latch_day retained for the per-week latch reporting only; no brake-zeroing.
        dosed = CK.is_on_day(day) if pulsed else (day <= treat)
        if dosed:
            for k, v in clamp.items():
                nb[k] = v
            weeks[w]['dosed_days'] += 1
        hist.append(nb)
        if CK.day_in_week(day) == CK.DAYS_PER_WEEK or day == CK.total_days():
            weeks[w]['casp3_end'] = int(nb.get('CASP3', 0))

    obs_days = hist[-CK.observation_days():]
    durable = sum(s['CASP3'] == 1 for s in obs_days) / len(obs_days) > 0.5
    return {
        'weeks': weeks,
        'latch_day': latch_day,
        'latch_week': (CK.week_of_day(latch_day) if latch_day else None),
        'latch_phase': (CK.phase_of_day(latch_day) if latch_day else None),
        'durable': durable,
    }


def run_capped_on_clock(b, U, kernels, apop, pulsed=True, front_loaded=False,
                        genotoxic=False, p53_lof=0, sticky=False, theta=10, cmax=20):
    """Step 10 -- phase caps + escalation govern how many pathways are dosed.

    Two induction philosophies, selected by `front_loaded`:
      front_loaded=False (escalate-from-1): start at 1 active pathway and climb by 1 at the end of
          any uncommitted treatment week (minimal-effective-dose / parsimonious).
      front_loaded=True (induction blast): start at the induction cap (4) from day one and let the
          phase boundaries de-escalate it (maintenance clips to 2, observation to 0) -- the classic
          clinical meaning of induction-then-maintenance.
    Observation is always drug-free (cap 0). Returns the per-week timeline and the durability verdict.
    """
    boxes = list(kernels.keys())
    bus = {n: int(b.get(n, 0)) for n in NODES}
    c = 0; commit = 0; hist = [bus]; triggered = 0
    treat = CK.treatment_days()
    active = CK.INDUCTION_PATHWAY_CAP if front_loaded else 1   # blast vs titrate
    weeks = {}

    for t in range(CK.total_days()):
        day = t + 1
        w = CK.week_of_day(day)
        cap = CK.pathway_cap(day)
        active = CK.clip_to_cap(active, day)             # tighten when phase tightens (de-escalation)
        dosing_day = (CK.is_on_day(day) if pulsed else (day <= treat))
        gtx = int((dosing_day or (sticky and triggered)) and genotoxic and not p53_lof)
        if gtx:
            triggered = 1
        nb = HB.step_honest(HB._view(hist, t, lags), U, genotoxic=gtx, cascade=True)
        c = min(c + 1, cmax) if HB.commit_signal(nb) else max(c - 1, 0)
        commit = commit or (c >= theta)   # Step 16: retained ONLY as the escalation-stop signal
        # (has the AKT1-FOXO3 switch flipped). The flag's brake-zeroing is gone; commitment and
        # durability are now read from the real cascade via step_honest(cascade=True).
        if dosing_day and cap > 0 and active > 0:
            for bx in boxes[:active]:                     # only the active pathways are dosed
                for nd in kernels[bx]:
                    nb[nd] = apop[nd]
        hist.append(nb)
        # end of a week = end of a cycle: record, then (escalate-from-1 mode only) climb if uncommitted
        if CK.day_in_week(day) == CK.DAYS_PER_WEEK or day == CK.total_days():
            weeks[w] = {'phase': CK.phase_of_day(day), 'cap': cap,
                        'active': (active if cap > 0 else 0), 'committed': int(commit)}
            if not front_loaded and CK.phase_of_day(day) != 'observation' and not commit:
                active = CK.escalate_active(active, day)   # FAILURE -> escalate from 1 upward

    obs_days = hist[-CK.observation_days():]
    durable = sum(s['CASP3'] == 1 for s in obs_days) / len(obs_days) > 0.5
    return {'weeks': weeks, 'n_boxes': len(boxes), 'durable': durable}


if __name__ == '__main__':
    # Head-to-head: pulsed (faithful) vs continuous hold, same 70-day clock, same kernels.
    import pandas as pd
    import forward_stab_kernel_design as FSK
    APOP = FSK.APOP_TARGET
    SW_STATE = ['AKT1', 'PTEN', 'MTOR', 'PDPK1', 'PIK3CA', 'MDM2', 'TP53', 'FOXO3']
    COH = {'TCGA': 'setup_a/data/binarized/tcga_brca_1082x135.csv',
           'METABRIC': 'setup_a/data/binarized/metabric_1980x135.csv',
           'ISPY2': 'setup_a/data/binarized/ispy2_988x135.csv'}
    SAMPLE = int(sys.argv[1]) if len(sys.argv) > 1 else 25

    def seed(b):
        U = HB.patient_clamp_off(b)
        x0 = {nd: int(b.get(nd, 0)) for nd in SW_STATE}
        x0['PHLPP'] = int(x0['FOXO3'] and not U)
        I = dict(SRC=int(b.get('SRC', 0)), RHEB=int(b.get('RHEB', 0)), IGF1R=int(b.get('IGF1R', 0)),
                 RTK_up=int(b.get('GRB2', 0) or b.get('IRS1', 0)), CDKN2A=int(b.get('CDKN2A', 0)),
                 E2F1=int(b.get('E2F1', 0)), ATM=int(b.get('ATM', 0)), ATR=int(b.get('ATR', 0)))
        return x0, I, U

    print(f"clock durability  |  {SAMPLE}/cohort  |  continuous-hold vs pulsed 4-on/3-off (70-day clock)")
    for name, path in COH.items():
        df = pd.read_csv(path, index_col=0)
        n = min(SAMPLE, len(df))
        found = cont = puls = 0
        for i in range(n):
            b = {nd: int(df.iloc[i].get(nd, 0)) for nd in NODES}
            x0, I, U = seed(b)
            r, kernels = FSK.design_patient_kernel(x0, I, U, method='stabilize')
            if not kernels:
                continue
            found += 1
            clamp = {nd: APOP[nd] for bx in kernels for nd in kernels[bx]}
            cont += A_dyn_clock(b, U, clamp, pulsed=False)
            puls += A_dyn_clock(b, U, clamp, pulsed=True)
        pc = lambda x: f"{round(100*x/found)}%" if found else "-"
        print(f"  {name:9} kerneled={found:3} | continuous {pc(cont):>4} | pulsed {pc(puls):>4}")
