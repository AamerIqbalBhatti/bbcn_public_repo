"""
clinical_clock.py -- canonical clinical time for BBCN (Phase 3, step 7).

ANCHOR (the one decision this module encodes):
    one binary network update  =  one clinical DAY.

Everything downstream is built on that tick. Dosing follows capivasertib (an AKT
inhibitor; CAPItello-291), which is given on a 4-days-on / 3-days-off weekly cycle.

    1 step       = 1 day
    1 week       = 7 days = 4 ON + 3 OFF
    treatment    = 6 weeks = 4 induction + 2 maintenance
    observation  = 4 weeks, drug-free (the durability tail)
    total        = 10 weeks = 70 days

This module DEFINES the clock only. Wiring it into the controllers -- pulsed dosing
(step 8), cycles + observation window (step 9), and the induction/maintenance pathway
caps (step 10) -- comes later. Nothing imports it yet, so it is behaviorally inert and
the regression bar is untouched. See ledger E66.
"""

# --- the anchor ---
DAYS_PER_STEP = 1
DAYS_PER_WEEK = 7
ON_DAYS = 4          # capivasertib: dosed the first 4 days of each treatment week
OFF_DAYS = 3         # then 3 days off
assert ON_DAYS + OFF_DAYS == DAYS_PER_WEEK

# --- windows (weeks) ---
INDUCTION_WEEKS = 4
MAINTENANCE_WEEKS = 2
TREATMENT_WEEKS = INDUCTION_WEEKS + MAINTENANCE_WEEKS    # 6
OBSERVATION_WEEKS = 4
TOTAL_WEEKS = TREATMENT_WEEKS + OBSERVATION_WEEKS         # 10

# --- pathway caps per phase (used at step 10) ---
INDUCTION_PATHWAY_CAP = 4
MAINTENANCE_PATHWAY_CAP = 2
OBSERVATION_PATHWAY_CAP = 0      # drug-free: no new control during observation


def day_of_step(step: int) -> int:
    """1-based step index -> 1-based day index (1 step = 1 day)."""
    return step * DAYS_PER_STEP


def week_of_day(day: int) -> int:
    """1-based day -> 1-based week."""
    return (day - 1) // DAYS_PER_WEEK + 1


def day_in_week(day: int) -> int:
    """1-based day -> position within its week, 1..7."""
    return (day - 1) % DAYS_PER_WEEK + 1


def phase_of_day(day: int) -> str:
    """'induction' | 'maintenance' | 'observation' for a 1-based day."""
    w = week_of_day(day)
    if w <= INDUCTION_WEEKS:
        return 'induction'
    if w <= TREATMENT_WEEKS:
        return 'maintenance'
    return 'observation'


def is_on_day(day: int) -> bool:
    """True if the drug is dosed on this day.

    Dosed only during the treatment window, and within a treatment week only on the
    first ON_DAYS days. The observation window is entirely drug-free.
    """
    if phase_of_day(day) == 'observation':
        return False
    return day_in_week(day) <= ON_DAYS


def pathway_cap(day: int) -> int:
    """Induction cap 4, maintenance cap 2, observation 0 (drug-free)."""
    return {
        'induction': INDUCTION_PATHWAY_CAP,
        'maintenance': MAINTENANCE_PATHWAY_CAP,
        'observation': OBSERVATION_PATHWAY_CAP,
    }[phase_of_day(day)]


def escalate_active(active: int, day: int) -> int:
    """Escalation rule: the controller starts at 1 active pathway and, on a cycle that fails to
    commit, climbs by 1 toward the day's phase cap (induction up to 4, maintenance up to 2).
    """
    return min(active + 1, pathway_cap(day))


def clip_to_cap(active: int, day: int) -> int:
    """Clip the active count to the current phase cap (used at phase boundaries)."""
    return min(active, pathway_cap(day))


def total_days() -> int:
    """Full horizon in days (treatment + observation)."""
    return TOTAL_WEEKS * DAYS_PER_WEEK   # 70


def treatment_days() -> int:
    return TREATMENT_WEEKS * DAYS_PER_WEEK   # 42


def observation_days() -> int:
    return OBSERVATION_WEEKS * DAYS_PER_WEEK  # 28


def schedule_summary() -> str:
    """Human-readable day-by-day schedule, for inspection/tests."""
    rows = []
    for d in range(1, total_days() + 1):
        rows.append((d, week_of_day(d), day_in_week(d), phase_of_day(d),
                     'ON' if is_on_day(d) else 'off', pathway_cap(d)))
    return rows
