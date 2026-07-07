#!/usr/bin/env python3
"""
reproduce.py — ONE command to regenerate BBCN's headline numbers and figures.

    python reproduce.py            # regenerate, write results/generated/, diff vs reference
    python reproduce.py --write-reference   # (maintainers) refresh the committed reference

What it does, in plain terms:
  1. Loads the two shipped binarized cohorts (TCGA 1082, METABRIC 1980).
  2. Runs the canonical controller in three regimes and counts, per phenotype,
     the % of patients driven to a GENUINE fixed point (no forcing, no accumulation).
  3. Writes a results table (CSV) and bar-chart figures with plain captions.
  4. Compares the freshly-generated numbers against the committed reference and
     prints, in words, whether they match.

A non-expert can simply read results/reference/REFERENCE_TABLE.csv and the figures.
"""
import sys, os, csv, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from bbcn import controller as C
from bbcn import harness as H

HERE = os.path.dirname(os.path.abspath(__file__))
COHORTS = {
    'TCGA': os.path.join(HERE, 'data/binarized/tcga_brca_1082x135.csv'),
    'METABRIC': os.path.join(HERE, 'data/binarized/metabric_1980x135.csv'),
}
PHEN = ['Resistance_OFF', 'Apoptosis_ON', 'Proliferation_OFF']
TOL = 2.0  # percentage-point tolerance for the verification diff


def _rows(path):
    B = pd.read_csv(path); cols = [c for c in B.columns if c in set(H.ALL_NODES)]
    for _, r in B.iterrows():
        yield {n: int(r[n]) for n in cols}


def run_regime(path, mode, family):
    agg = {p: 0 for p in PHEN}; allc = 0; N = 0
    for init in _rows(path):
        o = C.run_patient(init, mode=mode, family=family)
        for p in PHEN:
            agg[p] += int(o[p])
        if mode == 'sequenced':
            allc += int(o.get('all_three', False))
        N += 1
    res = {p: round(100 * agg[p] / N) for p in PHEN}
    if mode == 'sequenced':
        res['all_three'] = round(100 * allc / N)
    res['N'] = N
    return res


def generate():
    regimes = [
        ('Isolated / ideal targets',        'isolated',  'ideal'),
        ('Isolated / compatible targets',   'isolated',  'compatible'),
        ('Sequenced + state-preservation',  'sequenced', 'compatible'),
    ]
    table = []
    for cohort, path in COHORTS.items():
        for label, mode, family in regimes:
            r = run_regime(path, mode, family)
            row = {'cohort': cohort, 'regime': label, 'N': r['N'],
                   'Resistance_%': r['Resistance_OFF'], 'Apoptosis_%': r['Apoptosis_ON'],
                   'Proliferation_%': r['Proliferation_OFF'],
                   'AllThree_%': r.get('all_three', '')}
            table.append(row); print("  ", row)
    return table


def write_csv(table, path):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)


def make_figures(table, outdir):
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for cohort in COHORTS:
        iso = [r for r in table if r['cohort'] == cohort and r['regime'].startswith('Isolated / ideal')][0]
        seq = [r for r in table if r['cohort'] == cohort and r['regime'].startswith('Sequenced')][0]
        fig, ax = plt.subplots(figsize=(7, 4.2))
        labels = ['Resistance\nOFF', 'Apoptosis\nON', 'Proliferation\nOFF', 'ALL THREE\n(sequenced)']
        iso_v = [iso['Resistance_%'], iso['Apoptosis_%'], iso['Proliferation_%'], None]
        seq_v = [seq['Resistance_%'], seq['Apoptosis_%'], seq['Proliferation_%'], seq['AllThree_%']]
        x = range(len(labels)); w = 0.38
        ax.bar([i - w/2 for i in x], [v if v is not None else 0 for v in iso_v], w,
               label='Isolated (each phenotype alone)', color='#2c7fb8')
        ax.bar([i + w/2 for i in x], seq_v, w,
               label='Sequenced + state-preservation', color='#d95f0e')
        ax.set_ylabel('% of patients reaching a genuine\ncontrolled fixed point'); ax.set_ylim(0, 100)
        ax.set_xticks(list(x)); ax.set_xticklabels(labels, fontsize=9)
        ax.set_title(f'BBCN controllability — {cohort} (N={iso["N"]})', fontsize=11)
        ax.legend(fontsize=8, loc='upper center')
        for i, v in enumerate(seq_v):
            ax.text(i + w/2, v + 1.5, f'{v}%', ha='center', fontsize=8)
        for i, v in enumerate(iso_v):
            if v is not None: ax.text(i - w/2, v + 1.5, f'{v}%', ha='center', fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(outdir, f'controllability_{cohort}.png'), dpi=130)
        plt.close(fig)


def verify(gen, ref):
    def key(r): return (r['cohort'], r['regime'])
    refmap = {key(r): r for r in ref}; ok = True; msgs = []
    for g in gen:
        rr = refmap.get(key(g))
        if not rr:
            msgs.append(f"  NEW row not in reference: {key(g)}"); ok = False; continue
        for col in ['Resistance_%', 'Apoptosis_%', 'Proliferation_%', 'AllThree_%']:
            gv, rv = g[col], rr[col]
            if gv == '' and (rv == '' or rv is None): continue
            if abs(float(gv) - float(rv)) > TOL:
                msgs.append(f"  MISMATCH {key(g)} {col}: generated {gv} vs reference {rv}"); ok = False
    return ok, msgs


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--write-reference', action='store_true')
    args = ap.parse_args()
    os.makedirs(os.path.join(HERE, 'results/generated'), exist_ok=True)
    os.makedirs(os.path.join(HERE, 'results/figures'), exist_ok=True)
    print("Regenerating BBCN headline numbers (this runs the full cohorts)...")
    table = generate()
    gen_csv = os.path.join(HERE, 'results/generated/RESULTS_TABLE.csv'); write_csv(table, gen_csv)
    make_figures(table, os.path.join(HERE, 'results/figures'))
    ref_csv = os.path.join(HERE, 'results/reference/REFERENCE_TABLE.csv')
    if args.write_reference:
        write_csv(table, ref_csv); print(f"\nReference written to {ref_csv}")
    else:
        if not os.path.exists(ref_csv):
            print("\nNo committed reference yet; run with --write-reference once to create it."); sys.exit(0)
        ref = list(csv.DictReader(open(ref_csv)))
        ok, msgs = verify(table, ref)
        print("\n=== VERIFICATION ===")
        if ok:
            print("  PASS — regenerated numbers match the committed reference (within +/-2 pts).")
        else:
            print("  DIFFERENCES vs reference:"); [print(m) for m in msgs]
            sys.exit(1)
