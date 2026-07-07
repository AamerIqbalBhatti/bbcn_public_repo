"""
plots_cde_vs_switch.py -- Paper 4 figures: CDE vs switch kernels on the cascaded
delayed BBCN-D engine.

Reads the comparison written by cde_vs_switch_cascade.py
(setup_a/data/cde_vs_switch_summary.csv and cde_vs_switch_{cohort}.csv) and saves
to setup_a/figures/ (PNG + SVG):

  1. cde_vs_switch_durability.(png|svg)  -- HEADLINE: durable% per cohort, CDE vs
                                            switch, grouped bars.
  2. cde_vs_switch_held_vs_durable.(png) -- held% vs durable% per arm per cohort
                                            (transient vs lasting).
  3. cde_vs_switch_kernels.(png)         -- mean kernel size per arm + CDE pathway
                                            burden (C6 style).

Both arms are tested on the SAME flag-free cascaded delayed engine
(honest_bbcn.step_honest(cascade=True)); only the kernel design differs. The CDE
kernel is designed on the repaired harness (BBCN-M, COMMIT-latched) and excludes
the cascade-internal executioner nodes (which the engine computes); the switch
kernel is the 9-node BBCN-S design.

Usage (from repo root):
    python setup_a/plots_cde_vs_switch.py
"""
import os
import sys
from collections import Counter

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATADIR = os.path.join(ROOT, 'setup_a', 'data')
FIGDIR = os.path.join(ROOT, 'setup_a', 'figures')

CDE_C = '#4C72B0'    # CDE arm color
SW_C = '#C44E52'     # switch arm color


def _save(fig, stem):
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ('png', 'svg'):
        out = os.path.join(FIGDIR, f'{stem}.{ext}')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"  wrote {out}")
    plt.close(fig)


def load_summary():
    p = os.path.join(DATADIR, 'cde_vs_switch_summary.csv')
    if not os.path.exists(p):
        sys.exit(f"missing {p}; run cde_vs_switch_cascade.py first")
    return pd.read_csv(p)


def load_patients():
    import glob
    paths = sorted(glob.glob(os.path.join(DATADIR, 'cde_vs_switch_*.csv')))
    paths = [p for p in paths if not p.endswith('summary.csv')]
    return pd.concat((pd.read_csv(p) for p in paths), ignore_index=True)


def _grouped(ax, cohorts, cde_vals, sw_vals, ylabel, title):
    import numpy as np
    x = np.arange(len(cohorts)); w = 0.38
    b1 = ax.bar(x - w / 2, cde_vals, w, label='CDE (Setup A)', color=CDE_C)
    b2 = ax.bar(x + w / 2, sw_vals, w, label='Switch (Setup B)', color=SW_C)
    ax.set_xticks(x); ax.set_xticklabels(cohorts)
    ax.set_ylabel(ylabel); ax.set_title(title)
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h, f'{h:.0f}',
                    ha='center', va='bottom', fontsize=8)
    ax.legend()


def fig_durability(s):
    fig, ax = plt.subplots(figsize=(8, 5))
    _grouped(ax, s['cohort'], s['cde_union_durable_pct'], s['switch_durable_pct'],
             'durable %  (CASP3 sustained after release)',
             'Durability on the cascaded delayed BBCN-D: CDE vs switch kernel')
    ax.set_ylim(0, 105)
    fig.text(0.5, -0.02,
             'Same flag-free cascade engine + readout for both arms; only the kernel design differs.',
             ha='center', fontsize=8, style='italic')
    _save(fig, 'cde_vs_switch_durability')


def fig_held_vs_durable(s):
    import numpy as np
    cohorts = list(s['cohort']); x = np.arange(len(cohorts)); w = 0.2
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x - 1.5 * w, s['cde_union_held_pct'], w, label='CDE held', color=CDE_C, alpha=0.55)
    ax.bar(x - 0.5 * w, s['cde_union_durable_pct'], w, label='CDE durable', color=CDE_C)
    ax.bar(x + 0.5 * w, s['switch_held_pct'], w, label='Switch held', color=SW_C, alpha=0.55)
    ax.bar(x + 1.5 * w, s['switch_durable_pct'], w, label='Switch durable', color=SW_C)
    ax.set_xticks(x); ax.set_xticklabels(cohorts)
    ax.set_ylabel('% of kerneled patients'); ax.set_ylim(0, 105)
    ax.set_title('Apoptosis reached (held) vs durable (after release), per arm')
    ax.legend(ncol=2)
    _save(fig, 'cde_vs_switch_held_vs_durable')


def fig_kernels(s, pat):
    import numpy as np
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # (a) mean kernel size per arm per cohort
    cohorts = list(s['cohort']); x = np.arange(len(cohorts)); w = 0.38
    ax1.bar(x - w / 2, s['cde_union_mean_size'], w, label='CDE (union)', color=CDE_C)
    ax1.bar(x + w / 2, s['switch_mean_size'], w, label='Switch', color=SW_C)
    ax1.set_xticks(x); ax1.set_xticklabels(cohorts)
    ax1.set_ylabel('mean kernel size (nodes pinned)')
    ax1.set_title('Kernel size: CDE vs switch')
    for arr, off in ((s['cde_union_mean_size'], -w / 2), (s['switch_mean_size'], w / 2)):
        for xi, v in zip(x, arr):
            ax1.text(xi + off, v, f'{v:.1f}', ha='center', va='bottom', fontsize=8)
    ax1.legend()

    # (b) CDE pathway burden (pooled over found patients)
    cnt = Counter()
    for s_ in pat.loc[pat['switch_found'] == 1, 'cde_pathways'].dropna():
        for pw in str(s_).split(';'):
            if pw:
                cnt[pw] += 1
    if cnt:
        items = sorted(cnt.items(), key=lambda kv: kv[1])
        names = [k for k, _ in items]; vals = [v for _, v in items]
        ax2.barh(names, vals, color=CDE_C)
        for yi, v in enumerate(vals):
            ax2.text(v, yi, f' {v}', va='center', fontsize=8)
    ax2.set_xlabel('# patients whose CDE kernel touches the pathway')
    ax2.set_title('CDE pathway burden  (switch arm: fixed 9-node BBCN-S)')

    fig.suptitle('Kernel composition: CDE (Setup A) vs switch (Setup B)', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, 'cde_vs_switch_kernels')


def main():
    s = load_summary()
    pat = load_patients()
    print(f"summary cohorts: {list(s['cohort'])}; patients rows: {len(pat)}")
    fig_durability(s)
    fig_held_vs_durable(s)
    fig_kernels(s, pat)
    print("done.")


if __name__ == '__main__':
    main()
