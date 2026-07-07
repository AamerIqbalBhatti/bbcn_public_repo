"""
bbcn.survival
================
Cox PH, Kaplan-Meier, and log-rank survival analysis.

All models correspond to the clinical validation in BBCN Paper III.

Functions
---------
cox_os(df)              : Cox PH — Overall Survival, fully adjusted
cox_dfs(df)             : Cox PH — Disease-Free Survival
cox_vscc(df)            : Cox PH — V_SCC as continuous predictor
cox_drug_overlap(df)    : Cox PH — kernel drug alignment
km_convergence(df)      : KM curves — converged vs non-converged
km_drug_overlap(df)     : KM curves — high vs low drug alignment
run_all_models(df)      : run all models, return dict of results
"""

from typing import Dict, Any, Optional, Tuple, List
import pandas as pd
import numpy as np

try:
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False
    print("WARNING: lifelines not installed. Run: pip install lifelines")


# ============================================================
# COX PH MODELS
# ============================================================

def cox_os(
    df: pd.DataFrame,
    penalizer: float = 0.1,
    verbose: bool = True,
) -> Optional[CoxPHFitter]:
    """
    Cox PH — Overall Survival.
    Covariates: converged, init_vscc, AGE, stage_num,
                subtype dummies, TMB_NONSYNONYMOUS, ANEUPLOIDY_SCORE

    Required columns: OS_MONTHS, OS_event, converged, AGE, stage_num,
                      subtype_Basal, subtype_Her2, subtype_LumB,
                      TMB_NONSYNONYMOUS, ANEUPLOIDY_SCORE
    """
    _check_lifelines()

    cols = [
        'OS_MONTHS', 'OS_event', 'converged', 'init_vscc',
        'AGE', 'stage_num', 'subtype_Basal', 'subtype_Her2',
        'subtype_LumB', 'subtype_Normal', 'TMB_NONSYNONYMOUS', 'ANEUPLOIDY_SCORE',
    ]
    sub = _prep(df, cols)
    if verbose:
        print(f"\nCox OS: N={len(sub)}, Events={int(sub['OS_event'].sum())}")

    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(sub, duration_col='OS_MONTHS', event_col='OS_event')
    if verbose:
        cph.print_summary(decimals=3)
        print(f"C-index: {cph.concordance_index_:.3f}")
    return cph


def cox_dfs(
    df: pd.DataFrame,
    penalizer: float = 0.1,
    verbose: bool = True,
) -> Optional[CoxPHFitter]:
    """
    Cox PH — Disease-Free Survival.
    Required columns: DFS_MONTHS, DFS_event + same covariates as cox_os.
    """
    _check_lifelines()

    cols = [
        'DFS_MONTHS', 'DFS_event', 'converged', 'init_vscc',
        'AGE', 'stage_num', 'subtype_Basal', 'subtype_Her2',
        'subtype_LumB', 'subtype_Normal', 'TMB_NONSYNONYMOUS', 'ANEUPLOIDY_SCORE',
    ]
    sub = _prep(df, cols)
    if verbose:
        print(f"\nCox DFS: N={len(sub)}, Events={int(sub['DFS_event'].sum())}")

    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(sub, duration_col='DFS_MONTHS', event_col='DFS_event')
    if verbose:
        cph.print_summary(decimals=3)
    return cph


def cox_vscc(
    df: pd.DataFrame,
    penalizer: float = 0.1,
    verbose: bool = True,
) -> Optional[CoxPHFitter]:
    """Cox PH — V_SCC as continuous predictor of OS."""
    _check_lifelines()
    cols = [
        'OS_MONTHS', 'OS_event', 'init_vscc', 'AGE', 'stage_num',
        'subtype_Basal', 'subtype_Her2', 'subtype_LumB', 'subtype_Normal',
    ]
    sub = _prep(df, cols)
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(sub, duration_col='OS_MONTHS', event_col='OS_event')
    if verbose:
        cph.print_summary(decimals=3)
    return cph


def cox_drug_overlap(
    df: pd.DataFrame,
    penalizer: float = 0.1,
    verbose: bool = True,
) -> Optional[CoxPHFitter]:
    """
    Cox PH — kernel drug alignment score predicting OS.
    Requires drug overlap columns added by cohort.add_drug_overlap().
    """
    _check_lifelines()
    cols = [
        'OS_MONTHS', 'OS_event', 'kernel_drug_overlap',
        'AGE', 'stage_num', 'subtype_Basal', 'subtype_Her2',
        'subtype_LumB', 'subtype_Normal',
    ]
    sub = _prep(df, cols)
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(sub, duration_col='OS_MONTHS', event_col='OS_event')
    if verbose:
        cph.print_summary(decimals=3)
    return cph


# ============================================================
# KAPLAN-MEIER
# ============================================================

def km_convergence(df: pd.DataFrame) -> Dict[str, Any]:
    """
    KM curves — converged vs non-converged patients.

    Returns dict with KMF objects, log-rank p-value,
    and median OS for each group.
    """
    _check_lifelines()
    km_data = df[['OS_MONTHS', 'OS_event', 'converged']].dropna()

    results = {}
    for c, label in [(1, 'Converged'), (0, 'Non-converged')]:
        sub = km_data[km_data['converged'] == c]
        kmf = KaplanMeierFitter()
        kmf.fit(sub['OS_MONTHS'], sub['OS_event'], label=label)
        results[label] = {
            'kmf': kmf, 'n': len(sub),
            'events': int(sub['OS_event'].sum()),
            'median_os': kmf.median_survival_time_,
        }

    ga = km_data[km_data['converged'] == 1]
    gb = km_data[km_data['converged'] == 0]
    lr = logrank_test(
        ga['OS_MONTHS'], gb['OS_MONTHS'],
        event_observed_A=ga['OS_event'], event_observed_B=gb['OS_event'],
    )
    results['logrank_p']    = lr.p_value
    results['logrank_stat'] = lr.test_statistic

    return results


def km_drug_overlap(df: pd.DataFrame) -> Dict[str, Any]:
    """
    KM curves — high (≥2 kernel drugs) vs low (0–1) alignment.
    Requires drug overlap columns.
    """
    _check_lifelines()
    km_data = df[['OS_MONTHS', 'OS_event', 'high_overlap']].dropna()

    results = {}
    labels  = {1: '≥2 kernel drugs', 0: '0–1 kernel drugs'}
    for v, label in labels.items():
        sub = km_data[km_data['high_overlap'] == v]
        kmf = KaplanMeierFitter()
        kmf.fit(sub['OS_MONTHS'], sub['OS_event'], label=label)
        results[label] = {
            'kmf': kmf, 'n': len(sub),
            'events': int(sub['OS_event'].sum()),
            'median_os': kmf.median_survival_time_,
        }

    ga = km_data[km_data['high_overlap'] == 1]
    gb = km_data[km_data['high_overlap'] == 0]
    lr = logrank_test(
        ga['OS_MONTHS'], gb['OS_MONTHS'],
        event_observed_A=ga['OS_event'], event_observed_B=gb['OS_event'],
    )
    results['logrank_p']    = lr.p_value
    results['logrank_stat'] = lr.test_statistic

    return results


# ============================================================
# RUN ALL MODELS
# ============================================================

def run_all_models(df: pd.DataFrame, verbose: bool = True) -> Dict[str, Any]:
    """
    Run all survival models. Returns dict of results.

    Keys: cox_os, cox_dfs, cox_vscc, cox_drug_overlap,
          km_convergence, km_drug_overlap
    """
    results = {}

    print("\n" + "=" * 60)
    print("RUNNING ALL SURVIVAL MODELS")
    print("=" * 60)

    for name, fn in [
        ('cox_os',           cox_os),
        ('cox_dfs',          cox_dfs),
        ('cox_vscc',         cox_vscc),
        ('cox_drug_overlap', cox_drug_overlap),
    ]:
        try:
            results[name] = fn(df, verbose=verbose)
        except Exception as ex:
            print(f"  {name}: SKIPPED — {ex}")
            results[name] = None

    for name, fn in [
        ('km_convergence',  km_convergence),
        ('km_drug_overlap', km_drug_overlap),
    ]:
        try:
            results[name] = fn(df)
            r = results[name]
            print(f"\n{name}: log-rank p={r['logrank_p']:.4f}")
        except Exception as ex:
            print(f"  {name}: SKIPPED — {ex}")
            results[name] = None

    return results


# ============================================================
# HELPERS
# ============================================================

def _check_lifelines():
    if not HAS_LIFELINES:
        raise ImportError(
            "lifelines is required for survival analysis.\n"
            "Install it with:  pip install lifelines"
        )


def _prep(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    """Select available columns and drop rows with missing values."""
    available = [c for c in cols if c in df.columns]
    missing   = [c for c in cols if c not in df.columns]
    if missing:
        print(f"  WARNING: Missing columns (skipped): {missing}")
    return df[available].dropna()
