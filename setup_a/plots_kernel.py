"""
plots_kernel.py -- figures from the per-patient kernel-composition capture.

Reads the CSVs written by capture_kernels.py
(setup_a/data/kernel_composition_{cohort}.csv) and saves to setup_a/figures/:

  1. kernel_heatmap.(png|svg)   -- pathway x kernel-node frequency: how often each
                                    node is pinned, grouped by the pathway that pins it
                                    (fraction of patients, pooled across cohorts).
  2. pathway_burden.(png|svg)   -- (a) distribution of kernel size, and
                                    (b) which pathways carry the kernel load
                                    (distinct patient-kernels contributed per pathway).

Usage (from repo root):
    python setup_a/plots_kernel.py
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')                       # headless
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATADIR = os.path.join(ROOT, 'setup_a', 'data')
FIGDIR = os.path.join(ROOT, 'setup_a', 'figures')


def load() -> pd.DataFrame:
    paths = sorted(glob.glob(os.path.join(DATADIR, 'kernel_composition_*.csv')))
    if not paths:
        sys.exit(f"no kernel_composition_*.csv in {DATADIR}; run capture_kernels.py first")
    df = pd.concat((pd.read_csv(p) for p in paths), ignore_index=True)
    print(f"loaded {len(df)} rows from {len(paths)} cohort file(s): "
          f"{df['patient_id'].nunique()} patients, {df['cohort'].nunique()} cohorts")
    return df


def _save(fig, stem):
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ('png', 'svg'):
        out = os.path.join(FIGDIR, f'{stem}.{ext}')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"  wrote {out}")
    plt.close(fig)


def kernel_heatmap(df: pd.DataFrame):
    """pathway x node: fraction of patients whose kernel pins that node in that pathway."""
    n_pat = df['patient_id'].nunique()
    # one count per (patient, pathway, node) -> patient coverage, not step/stage repeats
    uniq = df.drop_duplicates(['patient_id', 'pathway', 'kernel_node'])
    mat = (uniq.groupby(['pathway', 'kernel_node']).size()
                .unstack(fill_value=0).astype(float) / n_pat)
    # order rows by total load, columns by total frequency, for readability
    mat = mat.loc[mat.sum(1).sort_values(ascending=False).index,
                  mat.sum(0).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(max(8, 0.32 * mat.shape[1]),
                                    max(4, 0.42 * mat.shape[0])))
    im = ax.imshow(mat.values, aspect='auto', cmap='viridis', vmin=0, vmax=1)
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_xlabel('kernel node'); ax.set_ylabel('pathway')
    ax.set_title(f'Kernel composition: P(node pinned by pathway)  (n={n_pat} patients)')
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label='fraction of patients')
    _save(fig, 'kernel_heatmap')


def pathway_burden(df: pd.DataFrame):
    """(a) kernel-size distribution; (b) per-pathway kernel load."""
    # one record per (patient, stage, pathway) = one selected kernel
    kernels = df.drop_duplicates(['patient_id', 'stage', 'pathway', 'kernel_size'])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # (a) distribution of kernel size
    sizes = kernels['kernel_size'].value_counts().sort_index()
    ax1.bar(sizes.index.astype(str), sizes.values, color='#4C72B0')
    ax1.set_xlabel('kernel size (nodes pinned)'); ax1.set_ylabel('# selected kernels')
    ax1.set_title('Kernel-size distribution')
    for x, v in zip(sizes.index.astype(str), sizes.values):
        ax1.text(x, v, f'{v}', ha='center', va='bottom', fontsize=8)

    # (b) which pathways carry the load (distinct patient-kernels per pathway)
    load = kernels['pathway'].value_counts().sort_values(ascending=True)
    ax2.barh(load.index, load.values, color='#C44E52')
    ax2.set_xlabel('# selected kernels (patient x stage)'); ax2.set_ylabel('pathway')
    ax2.set_title('Pathway burden: who carries the kernel load')
    for y, v in zip(range(len(load)), load.values):
        ax2.text(v, y, f' {v}', va='center', fontsize=8)

    fig.suptitle('Pathway burden of the harness kernel selection', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, 'pathway_burden')


def main():
    df = load()
    kernel_heatmap(df)
    pathway_burden(df)
    print("done.")


if __name__ == '__main__':
    main()
